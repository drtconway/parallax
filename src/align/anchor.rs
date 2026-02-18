#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Anchor {
    pub query_pos: usize,
    pub ref_pos: usize,
    pub length: usize,
}

impl Anchor {
    pub fn new(query_pos: usize, ref_pos: usize, length: usize) -> Self {
        Self {
            query_pos,
            ref_pos,
            length,
        }
    }

    pub fn diagonal(&self) -> isize {
        self.ref_pos as isize - self.query_pos as isize
    }

    #[allow(dead_code)]
    pub fn order_by_length(a: &Anchor, b: &Anchor) -> std::cmp::Ordering {
        let res = b.length.cmp(&a.length);
        if res == std::cmp::Ordering::Equal {
            a.query_pos.cmp(&b.query_pos)
        } else {
            res
        }
    }

    pub fn order_by_query_pos(a: &Anchor, b: &Anchor) -> std::cmp::Ordering {
        a.query_pos
            .cmp(&b.query_pos)
            .then(b.ref_pos.cmp(&a.ref_pos))
    }
}
