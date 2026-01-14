
/// Complement a single nucleotide
#[inline]
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        _ => b'N',
    }
}

/// Reverse complement a sequence into the provided buffer
pub fn reverse_complement_into(seq: &[u8], buf: &mut Vec<u8>) {
    buf.clear();
    buf.reserve(seq.len());
    for &base in seq.iter().rev() {
        buf.push(complement(base));
    }
}