use std::sync::Arc;

#[derive(Clone)]
pub enum Rope<'a> {
    Atom(&'a str),
    Concat(Arc<Rope<'a>>, Arc<Rope<'a>>),
    Substr(Arc<Rope<'a>>, usize, usize),
    Invert(Arc<Rope<'a>>),
}

impl<'a> Rope<'a> {
    pub fn len(&self) -> usize {
        match self {
            Rope::Atom(s) => s.len(),
            Rope::Concat(left, right) => left.len() + right.len(),
            Rope::Substr(rope, start, len) => *len,
            Rope::Invert(rope) => rope.len(),
        }
    }

    pub fn slice(&self, range: impl RopeSlice) -> Rope<'a> {
        let (start, len) = range.to_start_len(self.len());
        Rope::Substr(Arc::new(self.clone()), start, len)
    }

    pub fn invert(&self) -> Rope<'a> {
        Rope::Invert(Arc::new(self.clone()))
    }

    fn materialize(&self) -> String {
        let mut pieces: Vec<Piece<'a>> = Vec::new();
        let n = self.len();
        Rope::materialize_inner(self, &mut pieces, 0, n, true);
        let mut result = String::with_capacity(n);
        for piece in &pieces {
            match piece {
                Piece::Forward(s) => result.push_str(s),
                Piece::Reverse(s) => result.extend(s.chars().rev()),
            }
        }
        result
    }

    fn materialize_inner(
        rope: &Rope<'a>,
        pieces: &mut Vec<Piece<'a>>,
        start: usize,
        len: usize,
        forward: bool,
    ) {
        assert!(
            start + len <= rope.len(),
            "Range must be within the bounds of the rope"
        );
        match rope {
            Rope::Atom(s) => {
                let end = start + len;
                if forward {
                    pieces.push(Piece::Forward(&s[start..end]));
                } else {
                    pieces.push(Piece::Reverse(&s[start..end]));
                }
            }
            Rope::Concat(left, right) => {
                let left_len = left.len();
                let end = start + len;
                if forward {
                    if start < left_len {
                        let left_part_len = usize::min(len, left_len - start);
                        Rope::materialize_inner(left, pieces, start, left_part_len, forward);
                        if left_part_len < len {
                            Rope::materialize_inner(right, pieces, 0, len - left_part_len, forward);
                        }
                    } else {
                        Rope::materialize_inner(right, pieces, start - left_len, len, forward);
                    }
                } else {
                    if end > left_len {
                        let right_start = start.saturating_sub(left_len);
                        let right_end = end - left_len;
                        Rope::materialize_inner(right, pieces, right_start, right_end - right_start, forward);
                    }
                    if start < left_len {
                        let left_end = usize::min(end, left_len);
                        Rope::materialize_inner(left, pieces, start, left_end - start, forward);
                    }
                }
            }
            Rope::Substr(inner_rope, inner_start, inner_len) => {
                let substr_start = inner_start + start;
                let substr_end = substr_start + len;
                assert!(
                    substr_end <= inner_start + inner_len,
                    "Substr range must be within the bounds of the inner rope"
                );
                Rope::materialize_inner(inner_rope, pieces, substr_start, len, forward);
            }
            Rope::Invert(inner_rope) => {
                let new_start = inner_rope.len() - start - len;
                Rope::materialize_inner(inner_rope, pieces, new_start, len, !forward);
            }
        }
    }
}

pub trait RopeSlice {
    fn to_start_len(&self, rope_len: usize) -> (usize, usize);
}

impl RopeSlice for std::ops::Range<usize> {
    fn to_start_len(&self, rope_len: usize) -> (usize, usize) {
        assert!(self.end <= rope_len && self.start <= self.end);
        (self.start, self.end - self.start)
    }
}

impl RopeSlice for std::ops::RangeFrom<usize> {
    fn to_start_len(&self, rope_len: usize) -> (usize, usize) {
        assert!(self.start <= rope_len);
        (self.start, rope_len - self.start)
    }
}

impl RopeSlice for std::ops::RangeTo<usize> {
    fn to_start_len(&self, rope_len: usize) -> (usize, usize) {
        assert!(self.end <= rope_len);
        (0, self.end)
    }
}

impl RopeSlice for std::ops::RangeFull {
    fn to_start_len(&self, rope_len: usize) -> (usize, usize) {
        (0, rope_len)
    }
}

impl<'a> std::ops::Add for Rope<'a> {
    type Output = Rope<'a>;
    fn add(self, rhs: Rope<'a>) -> Rope<'a> {
        Rope::Concat(Arc::new(self), Arc::new(rhs))
    }
}

impl<'a> From<&'a str> for Rope<'a> {
    fn from(s: &'a str) -> Self {
        Rope::Atom(s)
    }
}

impl<'a> From<Rope<'a>> for String {
    fn from(rope: Rope<'a>) -> Self {
        rope.materialize()
    }
}

enum Piece<'a> {
    Forward(&'a str),
    Reverse(&'a str),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atom_len() {
        let r = Rope::from("hello");
        assert_eq!(r.len(), 5);
    }

    #[test]
    fn atom_materialize() {
        let r = Rope::from("hello");
        assert_eq!(String::from(r), "hello");
    }

    #[test]
    fn concat_len() {
        let r = Rope::from("hello") + Rope::from(" world");
        assert_eq!(r.len(), 11);
    }

    #[test]
    fn concat_materialize() {
        let r = Rope::from("hello") + Rope::from(" world");
        assert_eq!(String::from(r), "hello world");
    }

    #[test]
    fn nested_concat() {
        let ab = Rope::from("aaa") + Rope::from("bbb");
        let abc = ab + Rope::from("ccc");
        assert_eq!(abc.len(), 9);
        assert_eq!(String::from(abc), "aaabbbccc");
    }

    #[test]
    fn slice_range() {
        let r = Rope::from("hello world");
        let s = r.slice(0..5);
        assert_eq!(s.len(), 5);
        assert_eq!(String::from(s), "hello");
    }

    #[test]
    fn slice_range_middle() {
        let r = Rope::from("hello world");
        let s = r.slice(6..11);
        assert_eq!(String::from(s), "world");
    }

    #[test]
    fn slice_range_from() {
        let r = Rope::from("hello world");
        let s = r.slice(6..);
        assert_eq!(String::from(s), "world");
    }

    #[test]
    fn slice_range_to() {
        let r = Rope::from("hello world");
        let s = r.slice(..5);
        assert_eq!(String::from(s), "hello");
    }

    #[test]
    fn slice_range_full() {
        let r = Rope::from("hello world");
        let s = r.slice(..);
        assert_eq!(String::from(s), "hello world");
    }

    #[test]
    fn slice_empty() {
        let r = Rope::from("hello");
        let s = r.slice(2..2);
        assert_eq!(s.len(), 0);
        assert_eq!(String::from(s), "");
    }

    #[test]
    fn slice_of_concat() {
        let r = Rope::from("hello") + Rope::from(" world");
        let s = r.slice(3..8);
        assert_eq!(String::from(s), "lo wo");
    }

    #[test]
    fn slice_of_slice() {
        let r = Rope::from("hello world");
        let s1 = r.slice(2..9);
        assert_eq!(String::from(s1.clone()), "llo wor");
        let s2 = s1.slice(1..5);
        assert_eq!(String::from(s2), "lo w");
    }

    #[test]
    fn invert_atom() {
        let r = Rope::from("hello");
        let inv = r.invert();
        assert_eq!(inv.len(), 5);
        assert_eq!(String::from(inv), "olleh");
    }

    #[test]
    fn invert_concat() {
        let r = Rope::from("abc") + Rope::from("def");
        let inv = r.invert();
        assert_eq!(String::from(inv), "fedcba");
    }

    #[test]
    fn double_invert() {
        let r = Rope::from("hello");
        let inv2 = r.invert().invert();
        assert_eq!(String::from(inv2), "hello");
    }

    #[test]
    fn invert_of_slice() {
        let r = Rope::from("abcdef");
        let s = r.slice(1..5);
        assert_eq!(String::from(s.clone()), "bcde");
        let inv = s.invert();
        assert_eq!(String::from(inv), "edcb");
    }

    #[test]
    fn slice_of_invert() {
        let r = Rope::from("abcdef");
        let inv = r.invert();
        assert_eq!(String::from(inv.clone()), "fedcba");
        let s = inv.slice(1..4);
        assert_eq!(String::from(s), "edc");
    }

    #[test]
    fn concat_with_invert() {
        let r = Rope::from("abc") + Rope::from("def").invert();
        assert_eq!(String::from(r), "abcfed");
    }

    #[test]
    fn empty_atom() {
        let r = Rope::from("");
        assert_eq!(r.len(), 0);
        assert_eq!(String::from(r), "");
    }

    #[test]
    fn single_char_atom() {
        let r = Rope::from("x");
        assert_eq!(r.len(), 1);
        assert_eq!(String::from(r.clone()), "x");
        assert_eq!(String::from(r.invert()), "x");
    }

    #[test]
    #[should_panic]
    fn slice_out_of_bounds() {
        let r = Rope::from("hello");
        r.slice(0..6);
    }

    #[test]
    #[should_panic]
    fn slice_start_after_end() {
        let r = Rope::from("hello");
        r.slice(3..2);
    }

    #[test]
    fn complex_tree() {
        // Build: ((ab + cd) + (ef + gh)), then slice and invert
        let ab = Rope::from("ab") + Rope::from("cd");
        let ef = Rope::from("ef") + Rope::from("gh");
        let full = ab + ef;
        assert_eq!(full.len(), 8);
        assert_eq!(String::from(full.clone()), "abcdefgh");

        let mid = full.slice(2..6);
        assert_eq!(String::from(mid.clone()), "cdef");

        let inv_mid = mid.invert();
        assert_eq!(String::from(inv_mid), "fedc");
    }
}