use std::time::Instant;

pub struct RateProgressLabels {
    pub item: String,
    pub unit: String,
}

impl Default for RateProgressLabels {
    fn default() -> Self {
        RateProgressLabels {
            item: "items".to_string(),
            unit: "units".to_string(),
        }
    }
}

pub struct RateProgressConfig {
    pub labels: RateProgressLabels,
    pub interval: f64,
    pub formatter: Box<dyn Fn(&RateProgressLabels, &RateProgressView) + Send + Sync>,
}

impl RateProgressConfig {
    pub fn with_item<S: Into<String>>(mut self, item: S) -> Self {
        self.labels.item = item.into();
        self
    }

    pub fn with_unit<S: Into<String>>(mut self, unit: S) -> Self {
        self.labels.unit = unit.into();
        self
    }

    pub fn with_interval(mut self, interval: f64) -> Self {
        self.interval = interval;
        self
    }

    pub fn with_formatter<F>(mut self, formatter: F) -> Self
    where
        F: Fn(&RateProgressLabels, &RateProgressView) + Send + Sync + 'static,
    {
        self.formatter = Box::new(formatter);
        self
    }
}

impl Default for RateProgressConfig {
    fn default() -> Self {
        RateProgressConfig {
            labels: RateProgressLabels::default(),
            interval: 10.0,
            formatter: Box::new(default_formatter),
        }
    }
}

pub struct RateProgress {
    config: RateProgressConfig,
    start: Instant,
    since: Instant,
    view: RateProgressView,
}

impl RateProgress {
    pub fn new() -> Self {
        RateProgress {
            config: RateProgressConfig::default(),
            start: Instant::now(),
            since: Instant::now(),
            view: RateProgressView::default(),
        }
    }

    pub fn with_config(config: RateProgressConfig) -> Self {
        RateProgress {
            config,
            start: Instant::now(),
            since: Instant::now(),
            view: RateProgressView::default(),
        }
    }

    pub fn record(&mut self, amount: u64) {
        let now = Instant::now();
        let elapsed = (now - self.since).as_secs_f64();

        self.view.update(amount, elapsed);
        if elapsed >= self.config.interval {
            self.view.update_rates();
            (self.config.formatter)(&self.config.labels, &self.view);
            self.view.reset_recent();
            self.since = now;
        }
    }

    pub fn finish(&mut self) {
        self.view.update_rates();
        (self.config.formatter)(&self.config.labels, &self.view);
    }
}

pub struct RateProgressView {
    total_count: usize,
    total_amount: u64,
    total_time: f64,
    recent_count: usize,
    recent_amount: u64,
    recent_time: f64,
    recent_count_rate: f64,
    recent_amount_rate: f64,
}

impl Default for RateProgressView {
    fn default() -> Self {
        RateProgressView {
            total_count: 0,
            total_amount: 0,
            total_time: 0.0,
            recent_count: 0,
            recent_amount: 0,
            recent_time: 0.0,
            recent_count_rate: 0.0,
            recent_amount_rate: 0.0,
        }
    }
}

impl RateProgressView {
    pub fn update(&mut self, amount: u64, elapsed: f64) {
        self.total_count += 1;
        self.total_amount += amount;

        self.recent_count += 1;
        self.recent_amount += amount;
        self.recent_time = elapsed;
    }

    pub fn update_rates(&mut self) {
        self.total_time += self.recent_time;
        
        let count_rate = if self.recent_time > 0.0 {
            self.recent_count as f64 / self.recent_time
        } else {
            0.0
        };
        let amount_rate = if self.recent_time > 0.0 {
            self.recent_amount as f64 / self.recent_time
        } else {
            0.0
        };
        self.recent_count_rate = 0.2 * self.recent_count_rate + 0.8 * count_rate;
        self.recent_amount_rate = 0.2 * self.recent_amount_rate + 0.8 * amount_rate;
    }

    pub fn reset_recent(&mut self) {
        
        self.recent_count = 0;
        self.recent_amount = 0;
        self.recent_time = 0.0;
    }
}

pub fn humanize(x: f64) -> (f64, &'static str) {
    let factors = [
        (1e9, "G"),
        (1e6, "M"),
        (1e3, "K"),
        (1.0, ""),
        //(1e-3, "m"),
        //(1e-6, "µ"),
        //(1e-9, "n"),
    ];
    let abs = if x < 0.0 { -x } else { x };
    for (factor, suffix) in factors {
        if abs >= factor {
            return (x / factor, suffix);
        }
    }
    (x, "")
}

pub fn default_formatter(labels: &RateProgressLabels, view: &RateProgressView) {
    let (count_rate, count_suffix) = humanize(view.recent_count_rate);
    let (amount_rate, amount_suffix) = humanize(view.recent_amount_rate);
    let item = &labels.item;
    let unit = &labels.unit;
    log::info!(
        "Processed {} {item} in {:.0}s [{:.2} {}{item}/s, {:.2} {}{unit}/s]",
        view.total_count,
        view.total_time,
        count_rate,
        count_suffix,
        amount_rate,
        amount_suffix
    );
}
