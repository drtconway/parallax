//! Persistent rope data structure for efficient slicing, concatenation, and
//! reversal of symbol sequences.
//!
//! A [`Rope`] represents a logical sequence built by composing borrowed slices
//! without copying. Supported operations:
//!
//! - **Concatenation** (`+`): joins two ropes into a single logical sequence.
//! - **Slicing** ([`Rope::slice`]): extracts a sub-range lazily.
//! - **Inversion** ([`Rope::invert`]): reverses the element order lazily.
//!
//! All operations are O(1) and allocation-free (apart from the `Arc` node).
//! The underlying data is only copied when the rope is materialised into an
//! owned value via `String::from(rope)` or `Vec::<u8>::from(rope)`.
//!
//! The rope is generic over any type implementing [`Symbols`], with built-in
//! support for `str` (UTF-8 text) and `[u8]` (raw bytes).

use std::sync::Arc;

// ── Symbols trait ────────────────────────────────────────────────────────────

/// Abstraction over contiguous, sliceable symbol types.
///
/// This trait lets [`Rope`] work generically over different element types
/// without code duplication. It is implemented for [`str`] (character-level
/// reversal) and `[u8]` (byte-level reversal).
pub trait Symbols {
    /// The owned, growable counterpart (e.g. `String` for `str`, `Vec<u8>` for `[u8]`).
    type Owned;

    /// Number of atomic units (bytes for both `str` and `[u8]`).
    fn seq_len(&self) -> usize;

    /// Return the sub-slice `[start..end)`. Must be valid for the type
    /// (e.g. on UTF-8 boundaries for `str`).
    fn subseq(&self, start: usize, end: usize) -> &Self;

    /// Create a new empty owned value with the given byte capacity.
    fn new_owned(capacity: usize) -> Self::Owned;

    /// Append `slice` to `owned` in forward (natural) order.
    fn append_forward(owned: &mut Self::Owned, slice: &Self);

    /// Append `slice` to `owned` in reverse order.
    fn append_reversed(owned: &mut Self::Owned, slice: &Self);
}

impl Symbols for str {
    type Owned = String;
    fn seq_len(&self) -> usize {
        self.len()
    }
    fn subseq(&self, start: usize, end: usize) -> &Self {
        &self[start..end]
    }
    fn new_owned(capacity: usize) -> String {
        String::with_capacity(capacity)
    }
    fn append_forward(owned: &mut String, slice: &str) {
        owned.push_str(slice);
    }
    fn append_reversed(owned: &mut String, slice: &str) {
        owned.extend(slice.chars().rev());
    }
}

impl Symbols for [u8] {
    type Owned = Vec<u8>;
    fn seq_len(&self) -> usize {
        self.len()
    }
    fn subseq(&self, start: usize, end: usize) -> &Self {
        &self[start..end]
    }
    fn new_owned(capacity: usize) -> Vec<u8> {
        Vec::with_capacity(capacity)
    }
    fn append_forward(owned: &mut Vec<u8>, slice: &[u8]) {
        owned.extend_from_slice(slice);
    }
    fn append_reversed(owned: &mut Vec<u8>, slice: &[u8]) {
        owned.extend(slice.iter().rev());
    }
}

// ── Rope ─────────────────────────────────────────────────────────────────────

/// A persistent, lazily-evaluated rope over borrowed symbol slices.
///
/// `Rope` records concatenations, sub-slicing, and inversions as lightweight
/// tree nodes.  No data is copied until [`materialize`](Rope::materialize) is
/// called (implicitly via `String::from` / `Vec::<u8>::from`).
///
/// # Type parameter
///
/// `S` is the underlying slice type, constrained by [`Symbols`].  In
/// practice this is either `str` or `[u8]`.
///
/// # Examples
///
/// ```
/// use parallax::utils::rope::Rope;
///
/// let greeting = Rope::from("hello") + Rope::from(" world");
/// assert_eq!(String::from(greeting.slice(0..5)), "hello");
///
/// let bytes = Rope::from(b"ACGT".as_slice());
/// assert_eq!(Vec::from(bytes.invert()), b"TGCA");
/// ```
pub enum Rope<'a, S: Symbols + ?Sized> {
    /// A leaf node borrowing a contiguous slice.
    Atom(&'a S),
    /// The concatenation of two child ropes.
    Concat(Arc<Rope<'a, S>>, Arc<Rope<'a, S>>),
    /// A sub-range `[start .. start+len)` of a child rope.
    Substr(Arc<Rope<'a, S>>, usize, usize),
    /// The reversal of a child rope.
    Invert(Arc<Rope<'a, S>>),
}

impl<'a, S: Symbols + ?Sized> Clone for Rope<'a, S> {
    fn clone(&self) -> Self {
        match self {
            Rope::Atom(s) => Rope::Atom(s),
            Rope::Concat(l, r) => Rope::Concat(Arc::clone(l), Arc::clone(r)),
            Rope::Substr(r, s, l) => Rope::Substr(Arc::clone(r), *s, *l),
            Rope::Invert(r) => Rope::Invert(Arc::clone(r)),
        }
    }
}

impl<'a, S: Symbols + ?Sized> Rope<'a, S> {
    /// Returns the length (in bytes) of the logical sequence.
    pub fn len(&self) -> usize {
        match self {
            Rope::Atom(s) => s.seq_len(),
            Rope::Concat(left, right) => left.len() + right.len(),
            Rope::Substr(_rope, _start, len) => *len,
            Rope::Invert(rope) => rope.len(),
        }
    }

    /// Lazily extract a sub-range of this rope.
    ///
    /// Accepts any standard range syntax (`a..b`, `a..`, `..b`, `..`).
    ///
    /// # Panics
    ///
    /// Panics if the range is out of bounds.
    pub fn slice(&self, range: impl RopeSlice) -> Rope<'a, S> {
        let (start, len) = range.to_start_len(self.len());
        Rope::Substr(Arc::new(self.clone()), start, len)
    }

    /// Lazily reverse the element order of this rope.
    ///
    /// For `str` ropes this reverses at the `char` level; for `[u8]` ropes
    /// it reverses individual bytes.
    pub fn invert(&self) -> Rope<'a, S> {
        Rope::Invert(Arc::new(self.clone()))
    }

    /// Flatten the rope tree into an owned value (`String` or `Vec<u8>`).
    ///
    /// This is the only operation that copies data. It is called
    /// automatically by the `From<Rope<..>>` conversions.
    fn materialize(&self) -> S::Owned {
        let mut pieces: Vec<Piece<'a, S>> = Vec::new();
        let n = self.len();
        Rope::materialize_inner(self, &mut pieces, 0, n, true);
        let mut result = S::new_owned(n);
        for piece in &pieces {
            match piece {
                Piece::Forward(s) => S::append_forward(&mut result, s),
                Piece::Reverse(s) => S::append_reversed(&mut result, s),
            }
        }
        result
    }

    /// Recursively walk the rope tree, collecting leaf slices into `pieces`.
    ///
    /// `start` and `len` describe the active window within `rope`.
    /// `forward` tracks the current orientation; an `Invert` node flips it.
    fn materialize_inner(
        rope: &Rope<'a, S>,
        pieces: &mut Vec<Piece<'a, S>>,
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
                let slice = s.subseq(start, end);
                if forward {
                    pieces.push(Piece::Forward(slice));
                } else {
                    pieces.push(Piece::Reverse(slice));
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
                        Rope::materialize_inner(
                            right,
                            pieces,
                            right_start,
                            right_end - right_start,
                            forward,
                        );
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

/// Helper trait that converts Rust range types into a `(start, len)` pair
/// for use by [`Rope::slice`].
pub trait RopeSlice {
    /// Convert into a `(start, length)` pair, given the total rope length.
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

impl<'a, S: Symbols + ?Sized> std::ops::Add for Rope<'a, S> {
    type Output = Rope<'a, S>;
    fn add(self, rhs: Rope<'a, S>) -> Rope<'a, S> {
        Rope::Concat(Arc::new(self), Arc::new(rhs))
    }
}

impl<'a> From<&'a str> for Rope<'a, str> {
    fn from(s: &'a str) -> Self {
        Rope::Atom(s)
    }
}

impl<'a> From<&'a [u8]> for Rope<'a, [u8]> {
    fn from(s: &'a [u8]) -> Self {
        Rope::Atom(s)
    }
}

impl<'a> From<Rope<'a, str>> for String {
    fn from(rope: Rope<'a, str>) -> Self {
        rope.materialize()
    }
}

impl<'a> From<Rope<'a, [u8]>> for Vec<u8> {
    fn from(rope: Rope<'a, [u8]>) -> Self {
        rope.materialize()
    }
}

/// A leaf fragment collected during materialisation, tagged with its
/// orientation so that reversal can be applied when building the output.
enum Piece<'a, S: Symbols + ?Sized> {
    Forward(&'a S),
    Reverse(&'a S),
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

    // ── [u8] tests ───────────────────────────────────────────────────────

    #[test]
    fn bytes_atom() {
        let r = Rope::from(b"hello".as_slice());
        assert_eq!(r.len(), 5);
        assert_eq!(Vec::from(r), b"hello");
    }

    #[test]
    fn bytes_concat() {
        let r = Rope::from(b"hello".as_slice()) + Rope::from(b" world".as_slice());
        assert_eq!(Vec::from(r), b"hello world");
    }

    #[test]
    fn bytes_slice() {
        let r = Rope::from(b"hello world".as_slice());
        let s = r.slice(6..11);
        assert_eq!(Vec::from(s), b"world");
    }

    #[test]
    fn bytes_invert() {
        let r = Rope::from(b"abcdef".as_slice());
        let inv = r.invert();
        assert_eq!(Vec::from(inv), b"fedcba");
    }

    #[test]
    fn bytes_complex() {
        let ab = Rope::from(b"ab".as_slice()) + Rope::from(b"cd".as_slice());
        let ef = Rope::from(b"ef".as_slice()) + Rope::from(b"gh".as_slice());
        let full = ab + ef;
        assert_eq!(full.len(), 8);

        let mid = full.slice(2..6);
        assert_eq!(Vec::from(mid.clone()), b"cdef");
        assert_eq!(Vec::from(mid.invert()), b"fedc");
    }
}
