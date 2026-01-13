//! Metrics collection and reporting.
//!
//! Provides a simple recorder that collects histogram statistics and
//! prints a summary to stderr at the end of execution.

use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::sync::{Arc, Mutex};

use metrics::{
    Counter, Gauge, Histogram, HistogramFn, Key, KeyName, Metadata, Recorder, SetRecorderError,
    SharedString, Unit,
};

mod sketch;
use sketch::DDSketch;

/// A histogram that tracks count, sum, min, max and quantiles via DDSketch.
struct SketchHistogram {
    count: u64,
    sum: f64,
    sum_squared: f64,
    min: f64,
    max: f64,
    sketch: DDSketch,
}

impl Default for SketchHistogram {
    fn default() -> Self {
        Self {
            count: 0,
            sum: 0.0,
            sum_squared: 0.0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
            sketch: DDSketch::new(0.01), // 1% relative accuracy
        }
    }
}

impl SketchHistogram {
    fn record(&mut self, value: f64) {
        if self.count == 0 {
            self.min = value;
            self.max = value;
        } else {
            self.min = self.min.min(value);
            self.max = self.max.max(value);
        }
        self.count += 1;
        self.sum += value;
        self.sum_squared += value * value;
        self.sketch.add(value);
    }

    fn mean(&self) -> f64 {
        if self.count == 0 {
            0.0
        } else {
            self.sum / self.count as f64
        }
    }

    fn stddev(&self) -> f64 {
        if self.count == 0 {
            0.0
        } else {
            let mean = self.mean();
            let variance = (self.sum_squared / self.count as f64) - (mean * mean);
            if variance < 0.0 { 0.0 } else { variance.sqrt() }
        }
    }

    fn quantile(&self, q: f64) -> f64 {
        self.sketch.quantile(q).unwrap_or(0.0)
    }
}

/// Snapshot of histogram data for reporting
struct HistogramSnapshot {
    count: u64,
    mean: f64,
    stddev: f64,
    min: f64,
    max: f64,
    quantiles: [f64; 10],
}

/// Thread-safe wrapper for SketchHistogram
struct AtomicHistogram {
    inner: Mutex<SketchHistogram>,
}

impl AtomicHistogram {
    fn new() -> Self {
        Self {
            inner: Mutex::new(SketchHistogram::default()),
        }
    }

    fn snapshot(&self) -> HistogramSnapshot {
        let guard = self.inner.lock().unwrap();
        HistogramSnapshot {
            count: guard.count,
            mean: guard.mean(),
            stddev: guard.stddev(),
            min: if guard.count > 0 { guard.min } else { 0.0 },
            max: if guard.count > 0 { guard.max } else { 0.0 },
            quantiles: [
                guard.quantile(0.0),
                guard.quantile(0.10),
                guard.quantile(0.20),
                guard.quantile(0.30),
                guard.quantile(0.40),
                guard.quantile(0.50),
                guard.quantile(0.60),
                guard.quantile(0.70),
                guard.quantile(0.80),
                guard.quantile(0.90),
            ],
        }
    }
}

impl HistogramFn for AtomicHistogram {
    fn record(&self, value: f64) {
        self.inner.lock().unwrap().record(value);
    }
}

/// Storage for our metrics
struct SummaryStorage {
    histograms: Mutex<HashMap<String, Arc<AtomicHistogram>>>,
}

impl SummaryStorage {
    fn new() -> Self {
        Self {
            histograms: Mutex::new(HashMap::new()),
        }
    }

    fn get_or_create_histogram(&self, name: &str) -> Arc<AtomicHistogram> {
        let mut histograms = self.histograms.lock().unwrap();
        histograms
            .entry(name.to_string())
            .or_insert_with(|| Arc::new(AtomicHistogram::new()))
            .clone()
    }

    fn snapshot(&self) -> HashMap<String, HistogramSnapshot> {
        let histograms = self.histograms.lock().unwrap();
        histograms
            .iter()
            .map(|(k, v)| (k.clone(), v.snapshot()))
            .collect()
    }
}

/// A recorder that collects summary statistics and prints them on drop.
pub struct SummaryRecorder {
    storage: Arc<SummaryStorage>,
}

impl SummaryRecorder {
    /// Create and install a new SummaryRecorder as the global recorder.
    pub fn install() -> Result<SummaryHandle, SetRecorderError<Self>> {
        let storage = Arc::new(SummaryStorage::new());
        let recorder = Self {
            storage: storage.clone(),
        };
        metrics::set_global_recorder(recorder)?;
        Ok(SummaryHandle { storage })
    }
}

impl Recorder for SummaryRecorder {
    fn describe_counter(&self, _key: KeyName, _unit: Option<Unit>, _description: SharedString) {}
    fn describe_gauge(&self, _key: KeyName, _unit: Option<Unit>, _description: SharedString) {}
    fn describe_histogram(&self, _key: KeyName, _unit: Option<Unit>, _description: SharedString) {}

    fn register_counter(&self, _key: &Key, _metadata: &Metadata<'_>) -> Counter {
        // We don't track counters in this simple implementation
        Counter::noop()
    }

    fn register_gauge(&self, _key: &Key, _metadata: &Metadata<'_>) -> Gauge {
        // We don't track gauges in this simple implementation
        Gauge::noop()
    }

    fn register_histogram(&self, key: &Key, _metadata: &Metadata<'_>) -> Histogram {
        let name = key.name().to_string();
        let histogram = self.storage.get_or_create_histogram(&name);
        Histogram::from_arc(histogram)
    }
}

/// Handle to access the collected metrics summary.
pub struct SummaryHandle {
    storage: Arc<SummaryStorage>,
}

impl SummaryHandle {
    /// Print a summary of all collected metrics to stderr as a tab-separated table.
    pub fn print_summary(&self) {
        let histograms = self.storage.snapshot();

        if histograms.is_empty() {
            return;
        }

        let out = File::create("parallax-stats.tsv").unwrap();
        let mut writer = std::io::BufWriter::new(out);

        let quantile_names = [
            "q10", "q20", "q30", "q40", "q50", "q60", "q70", "q80", "q90",
        ];

        // Header row
        writeln!(
            writer,
            "metric\tcount\tmean\tstddev\tmin\tmax\t{}",
            quantile_names.join("\t")
        )
        .unwrap();

        // Sort by name for consistent output
        let mut names: Vec<_> = histograms.keys().collect();
        names.sort();

        for name in names {
            let h = &histograms[name];
            if h.count > 0 {
                let quantiles: String = h
                    .quantiles
                    .iter()
                    .skip(1)
                    .map(|q| format!("\t{:.2}", q))
                    .collect();
                writeln!(
                    writer,
                    "{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}{}",
                    name, h.count, h.mean, h.stddev, h.min, h.max, quantiles
                )
                .unwrap();
            }
        }

        writer.flush().unwrap();
    }
}
