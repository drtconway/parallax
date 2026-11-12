use std::sync::{Arc, Mutex};

use super::*;

pub struct Histogram {
    histogram: HashMap<i64, usize>,
}

impl Histogram {
    pub fn new() -> Self {
        Histogram {
            histogram: HashMap::new(),
        }
    }

    pub fn record(&mut self, value: Value) {
        if let Some(int_value) = value.as_int() {
            *self.histogram.entry(int_value).or_insert(0) += 1;
        }
    }

    pub fn report(&self) -> Histo {
        let mut histogram: Vec<(i64, usize)> = self.histogram.iter().map(|(k, v)| (*k, *v)).collect();
        histogram.sort_unstable();
        Histo {
            histogram,
        }
    }
}

pub struct Histo {
    histogram: Vec<(i64, usize)>,
}

impl Reporter for Histo {
    fn columns(&self) -> &[&str] {
        static COLUMS: [&str; 2] = ["value", "count"];
        &COLUMS
    }

    fn rows(&self) -> Vec<Vec<Value>> {
        let res: Vec<Vec<Value>> = self
            .histogram
            .iter()
            .map(|(k, v)| vec![Value::from(*k), Value::from(*v)])
            .collect();
        res
    }
}

pub struct HistogramRecorder {
    inner: Arc<Mutex<Histogram>>,
}

impl HistogramRecorder {
    pub fn new() -> Self {
        HistogramRecorder {
            inner: Arc::new(Mutex::new(Histogram::new())),
        }
    }

    pub fn new_registered(key: &str) -> &'static Self {
        let recorder = Box::leak(Box::new(Self::new()));
        super::registry().register(key, recorder);
        recorder
    }
}

impl Recorder for HistogramRecorder {
    fn record_value(&self, value: Value) {
        let mut inner = self.inner.lock().unwrap();
        inner.record(value);
    }

    fn report(&self) -> Box<dyn Reporter> {
        let inner = self.inner.lock().unwrap();
        Box::new(inner.report())
    }
}
