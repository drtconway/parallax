

use crate::utils::telemetry::sampler::{Sample, Sampler};

use super::*;

pub struct Quantiles {
    n: usize,
    sampler: Sampler,
}

impl Quantiles {
    pub fn new(n: usize, f: usize, seed: u64) -> Self {
        Quantiles {
            n,
            sampler: Sampler::new(n * f, seed),
        }
    }

    fn record(&mut self, value: Value) {
        self.sampler.record(value);
    }
}


pub struct QuantileReport {
    pub quantiles: Vec<(f64, f64)>
}

impl QuantileReport {
    pub fn new(n: usize, sample: Sample) -> Self {
        let mut values: Vec<f64> = sample.values.iter().filter_map(|v| v.as_float()).collect();
        values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let quantiles = (0..n).map(|i| {
            let q = (i as f64 + 1.0) / (n as f64);
            let idx = (q * (values.len() as f64 - 1.0)).floor() as usize;
            (q, values[idx])
        }).collect();
        QuantileReport { quantiles }
    }
}

impl Reporter for QuantileReport {
    fn columns(&self) -> &[&str] {
        static COLUMS: [&str; 2] = ["quantile", "value"];
        &COLUMS
    }

    fn rows(&self) -> Vec<Vec<Value>> {
        self.quantiles
            .iter()
            .map(|(q, v)| vec![Value::Flt(*q), Value::Flt(*v)])
            .collect()
    }
}

pub struct QuantileRecorder {
    quantiles: Mutex<Quantiles>,
}

impl QuantileRecorder {
    pub fn new(n: usize, f: usize, seed: u64) -> Self {
        QuantileRecorder {
            quantiles: Mutex::new(Quantiles::new(n, f, seed)),
        }
    }
    
    pub fn new_registered(name: &str, n: usize, f: usize, seed: u64) -> &'static Self {
        let recorder = Box::leak(Box::new(Self::new(n, f, seed)));
        super::registry().register(name, recorder);
        recorder
    }
}

impl Recorder for QuantileRecorder {
    fn record_value(&self, value: Value) {
        if let Some(_) = value.as_float() {
            self.quantiles.lock().unwrap().record(value);
        }
    }
    
    fn report(&self) -> Box<dyn Reporter> {
        let quantiles = self.quantiles.lock().unwrap();
        Box::new(QuantileReport::new(10, quantiles.sampler.report()))
    }
}