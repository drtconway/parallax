use std::time::Instant;

pub struct RateProgressConfig {
    pub item: String,
    pub unit: String,
    pub interval: f64,
}

impl RateProgressConfig {
    pub fn with_item<S: Into<String>>(mut self, item: S) -> Self {
        self.item = item.into();
        self
    }

    pub fn with_unit<S: Into<String>>(mut self, unit: S) -> Self {
        self.unit = unit.into();
        self
    }

    pub fn with_interval(mut self, interval: f64) -> Self {
        self.interval = interval;
        self
    }
}

impl Default for RateProgressConfig {
    fn default() -> Self {
        RateProgressConfig {
            item: "items".to_string(),
            unit: "units".to_string(),
            interval: 10.0,
        }
    }
}

pub struct RateProgress {
    config: RateProgressConfig,
    counter: usize,
    amount: u64,
    start: Instant,
    recent_counter: usize,
    recent_amount: u64,
    since: Instant,
}

impl RateProgress {
    pub fn new() -> Self {
        RateProgress {
            config: RateProgressConfig::default(),
            counter: 0,
            amount: 0,
            start: Instant::now(),
            recent_counter: 0,
            recent_amount: 0,
            since: Instant::now(),
        }
    }

    pub fn with_config(config: RateProgressConfig) -> Self {
        RateProgress {
            config,
            counter: 0,
            amount: 0,
            start: Instant::now(),
            recent_counter: 0,
            recent_amount: 0,
            since: Instant::now(),
        }
    }

    pub fn record(&mut self, amount: u64) {
        self.counter += 1;
        self.recent_counter += 1;
        let total_count = self.counter;
        let recent_count = self.recent_counter;

        self.amount += amount;
        self.recent_amount += amount;
        let recent_amount = self.recent_amount;

        let now = Instant::now();
        let total_elapsed = (now - self.start).as_secs_f64();
        let recent_elapsed = (now - self.since).as_secs_f64();
        if recent_elapsed >= self.config.interval {
            let items = &self.config.item;
            let (count_rate, count_suffix) = humanize(recent_count as f64 / recent_elapsed);
            let (amount_rate, amount_suffix) = humanize(recent_amount as f64 / recent_elapsed);
            let unit = &self.config.unit;
            log::info!(
                "Processed {} {items} in {:.0}s [{:.3} {}{items}/s, {:.3} {}{unit}/s]",
                total_count,
                total_elapsed,
                count_rate,
                count_suffix,
                amount_rate,
                amount_suffix
            );
            self.recent_counter = 0;
            self.recent_amount = 0;
            self.since = now;
        }
    }

    pub fn finish(&self) {
        let count = self.counter;
        let total = self.amount;
        let elapsed = (Instant::now() - self.start).as_secs_f64();
        if count > 0 {
            let items = &self.config.item;
            let (count_rate, count_suffix) = humanize(count as f64 / elapsed);
            let (total_rate, total_suffix) = humanize(total as f64 / elapsed);
            let unit = &self.config.unit;
            log::info!(
                "Finished processing {} {items} in {:.0}s [{:.3} {}{items}/s, {:.3} {}{unit}/s]",
                count,
                elapsed,
                count_rate,
                count_suffix,
                total_rate,
                total_suffix
            );
        }
    }
}

fn humanize(x: f64) -> (f64, &'static str) {
    let factors = [(1e9, "G"), (1e6, "M"), (1e3, "K"), (1.0, ""), (1e-3, "m"), (1e-6, "µ"), (1e-9, "n")];
    let abs = if x < 0.0 { -x } else { x };
    for (factor, suffix) in factors {
        if abs >= factor {
            return (x / factor, suffix);
        }
    }
    (x, "")
}