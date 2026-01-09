
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Selection<T, U> {
    Left(T),
    Right(U),
    Both(T, U),
}

pub mod table;

pub struct GroupByKey<'a, F: Fn(&'a T) -> K, T, K: PartialEq> {
    items: &'a [T],
    key_fn: F,
    begin: usize,
    end: usize,
}

impl<'a, F: Fn(&'a T) -> K, T, K: PartialEq> GroupByKey<'a, F, T, K> {
    pub fn new(iter: &'a [T], key_fn: F) -> Self {
        GroupByKey {
            items: iter,
            key_fn,
            begin: 0,
            end: 0,
        }
    }
}

impl<'a, F: Fn(&'a T) -> K, T, K: PartialEq> Iterator for GroupByKey<'a, F, T, K> {
    type Item = (K, &'a [T]);

    fn next(&mut self) -> Option<Self::Item> {
        if self.end >= self.items.len() {
            return None;
        }

        self.begin = self.end;
        let key = (self.key_fn)(&self.items[self.begin]);
        self.end += 1;

        while self.end < self.items.len() {
            let next_key = (self.key_fn)(&self.items[self.end]);
            if &next_key != &key {
                break;
            }
            self.end += 1;
        }

        Some((key, &self.items[self.begin..self.end]))
    }
}

pub trait GroupByTrait<F: Fn(&Self::Item) -> K, K: PartialEq> {
    type Item;

    fn group_by(&'_ self, key_fn: F) -> GroupByKey<'_, F, Self::Item, K>;
}

impl<T, F: Fn(&T) -> K, K: PartialEq> GroupByTrait<F, K> for [T] {
    type Item = T;

    fn group_by(&'_ self, key_fn: F) -> GroupByKey<'_, F, T, K> {
        GroupByKey::new(self, key_fn)
    }
}