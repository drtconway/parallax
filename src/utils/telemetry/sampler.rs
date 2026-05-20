use std::sync::{Arc, Mutex};

use rand::{RngExt, SeedableRng, rngs::StdRng};

use super::*;

pub struct Sampler {
    n: usize,
    count: usize,
    values: Vec<Value>,
    rng: StdRng,
}

impl Sampler {
    pub fn new(n: usize, seed: u64) -> Self {
        let rng = StdRng::seed_from_u64(seed);
        Sampler {
            n,
            count: 0,
            values: vec![],
            rng,
        }
    }

    pub fn record(&mut self, value: Value) {
        self.count += 1;
        if self.values.len() < self.n {
            self.values.push(value);
        } else {
            let u = self.rng.random_range(0..self.count);
            if u < self.n {
                self.values[u] = value;
            }
        }
    }

    pub fn report(&self) -> Sample {
        Sample {
            values: self.values.clone(),
        }
    }
}

pub struct Sample {
    values: Vec<Value>,
}

impl Reporter for Sample {
    fn columns(&self) -> &[&str] {
        static COLUMS: [&str; 2] = ["n", "value"];
        &COLUMS
    }

    fn rows(&self) -> Vec<Vec<Value>> {
        let res: Vec<Vec<Value>> = self
            .values
            .iter()
            .enumerate()
            .map(|(n, v)| vec![Value::from(n), v.clone()])
            .collect();
        res
    }
}

pub struct SamplingRecorder {
    inner: Arc<Mutex<Sampler>>,
}

impl SamplingRecorder {
    pub fn new(n: usize, seed: u64) -> Self {
        SamplingRecorder {
            inner: Arc::new(Mutex::new(Sampler::new(n, seed))),
        }
    }
}

impl Recorder for SamplingRecorder {
    fn record(&self, value: Value) {
        let mut inner = self.inner.lock().unwrap();
        inner.record(value);
    }

    fn report(&self) -> Box<dyn Reporter> {
        let inner = self.inner.lock().unwrap();
        Box::new(inner.report())
    }
}
