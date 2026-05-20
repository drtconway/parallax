use std::sync::{Arc, Mutex};

use super::*;

pub struct SimpleSummarizer {
    n: u64,
    s: f64,
    s2: f64,
    min: f64,
    max: f64,
}

impl SimpleSummarizer {
    pub fn new() -> Self {
        SimpleSummarizer {
            n: 0,
            s: 0.0,
            s2: 0.0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
        }
    }

    pub fn record(&mut self, value: f64) {
        self.n += 1;
        self.s += value;
        self.s2 += value * value;
        self.min = self.min.min(value);
        self.max = self.max.max(value);
    }

    pub fn report(&self) -> SummaryReport {
        let n = self.n;
        let mean = self.s / (n as f64);
        let std_dev = self.s2 / (n as f64) - mean * mean;
        let min = self.min;
        let max = self.max;

        SummaryReport {
            n,
            mean,
            std_dev,
            min,
            max,
        }
    }
}

pub struct SummaryReport {
    n: u64,
    mean: f64,
    std_dev: f64,
    min: f64,
    max: f64,
}

impl Reporter for SummaryReport {
    fn columns(&self) -> &[&str] {
        static COLUMS: [&'static str; 5] = ["n", "mean", "std_dev", "min", "max"];
        &COLUMS
    }

    fn rows(&self) -> Vec<Vec<Value>> {
        vec![vec![
            Value::from(self.n),
            Value::from(self.mean),
            Value::from(self.std_dev),
            Value::from(self.min),
            Value::from(self.max),
        ]]
    }
}

pub struct SimpleSummaryRecorder {
    inner: Arc<Mutex<SimpleSummarizer>>,
}

impl SimpleSummaryRecorder {
    pub fn new() -> Self {
        SimpleSummaryRecorder {
            inner: Arc::new(Mutex::new(SimpleSummarizer::new())),
        }
    }
}

impl Recorder for SimpleSummaryRecorder {
    fn record(&self, value: Value) {
        if let Some(value) = value.as_float() {
            let mut inner = self.inner.lock().unwrap();
            inner.record(value);
        }
    }

    fn report(&self) -> Box<dyn Reporter> {
        let inner = self.inner.lock().unwrap();
        Box::new(inner.report())
    }
}
