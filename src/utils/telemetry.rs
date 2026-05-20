use std::collections::HashMap;

#[derive(Debug, Clone)]
pub enum Value {
    Int(i64),
    Flt(f64),
    Str(String),
}

impl Value {
    pub fn as_int(&self) -> Option<i64> {
        if let Self::Int(value) = self {
            Some(*value)
        } else {
            None
        }
    }

    pub fn as_float(&self) -> Option<f64> {
        match self {
            Value::Int(value) => Some(*value as f64),
            Value::Flt(value) => Some(*value),
            Value::Str(_value) => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        if let Self::Str(value) = self {
            Some(value)
        } else {
            None
        }
    }
}

impl From<i64> for Value {
    fn from(value: i64) -> Self {
        Value::Int(value)
    }
}

impl From<u64> for Value {
    fn from(value: u64) -> Self {
        Value::Int(value as i64)
    }
}

impl From<usize> for Value {
    fn from(value: usize) -> Self {
        Value::Int(value as i64)
    }
}

impl From<f64> for Value {
    fn from(value: f64) -> Self {
        Value::Flt(value)
    }
}

impl From<&str> for Value {
    fn from(value: &str) -> Self {
        Value::Str(value.into())
    }
}

impl From<String> for Value {
    fn from(value: String) -> Self {
        Value::Str(value)
    }
}

pub trait Recorder {
    fn record(&self, value: Value);

    fn report(&self) -> Box<dyn Reporter>;
}

pub trait Reporter {
    fn columns(&self) -> &[&str];

    fn rows(&self) -> Vec<Vec<Value>>;
}

pub struct Registry {
    recorders: HashMap<String, &'static dyn Recorder>,
}

impl Registry {
    pub fn register(&mut self, label: impl Into<String>, recorder: &'static impl Recorder) {
        self.recorders.insert(label.into(), recorder);
    }

    pub fn iter(&self) -> impl Iterator<Item = (&str, &'static dyn Recorder)> {
        let mut keys: Vec<&str> = self.recorders.keys().map(|s| s as &str).collect();
        keys.sort();
        keys.into_iter()
            .map(|k| (k, *self.recorders.get(k).unwrap()))
    }
}

pub mod sampler;
pub mod summary;
