use std::{
    collections::HashMap,
    io::Write,
    sync::{Arc, LazyLock, Mutex},
};

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

pub trait Recorder: Sync {
    fn record_usize(&self, value: usize) {
        self.record(Value::from(value));
    }

    fn record_f64(&self, value: f64) {
        self.record(Value::from(value));
    }
    
    fn record(&self, value: Value);

    fn report(&self) -> Box<dyn Reporter>;
}

pub trait Reporter {
    fn columns(&self) -> &[&str];

    fn rows(&self) -> Vec<Vec<Value>>;
}

pub struct Registry {
    recorders: Arc<Mutex<HashMap<String, &'static dyn Recorder>>>,
}

impl Registry {
    pub fn new() -> Self {
        Registry {
            recorders: Arc::new(Mutex::new(HashMap::new())),
        }
    }

    pub fn register(&self, label: impl Into<String>, recorder: &'static impl Recorder) {
        let mut recorders = self.recorders.lock().unwrap();
        recorders.insert(label.into(), recorder);
    }

    /// Write out the reports in a multi-section TSV
    pub fn report(&self, out: &mut impl Write) -> std::io::Result<()> {
        let recorders = self.recorders.lock().unwrap();
        let mut keys: Vec<&str> = recorders.keys().map(|s| s as &str).collect();
        keys.sort();
        for label in keys.into_iter() {
            let recorder = recorders.get(label).unwrap();
            let report = recorder.report();
            writeln!(out, "[{}]", label)?;
            for (i, col) in report.columns().iter().enumerate() {
                if i > 0 {
                    write!(out, "\t")?;
                }
                write!(out, "{}", *col)?;
            }
            writeln!(out, "")?;
            for row in report.rows().into_iter() {
                for (i, value) in row.into_iter().enumerate() {
                    if i > 0 {
                        write!(out, "\t")?;
                    }
                    match value {
                        Value::Int(value) => write!(out, "{}", value),
                        Value::Flt(value) => write!(out, "{}", value),
                        Value::Str(value) => write!(out, "{}", value),
                    }?;
                }
                writeln!(out, "")?;
            }
        }
        Ok(())
    }
}

pub fn registry() -> &'static Registry {
    static REGISTRY: LazyLock<Registry> = LazyLock::new(|| Registry::new());
    &*REGISTRY
}

pub mod sampler;
pub mod summary;
