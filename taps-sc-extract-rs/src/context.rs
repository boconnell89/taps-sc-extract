//! Reference cytosine context (CpG / CHG / CHH and CG / CH).

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TriContext {
    Cpg,
    Chg,
    Chh,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Context {
    Cg,
    Ch,
}

fn byte_at(seq: &[u8], pos: isize) -> u8 {
    if pos < 0 {
        return b'N';
    }
    let i = pos as usize;
    if i >= seq.len() {
        b'N'
    } else {
        seq[i].to_ascii_uppercase()
    }
}

/// Trinucleotide context at 0-based `pos`. None if out of bounds or not C/G.
pub fn classify_trinucleotide(seq: &[u8], pos: usize) -> Option<TriContext> {
    if pos >= seq.len() {
        return None;
    }
    let base = seq[pos].to_ascii_uppercase();
    match base {
        b'C' => {
            let n1 = byte_at(seq, pos as isize + 1);
            if n1 == b'G' {
                Some(TriContext::Cpg)
            } else {
                let n2 = byte_at(seq, pos as isize + 2);
                if n2 == b'G' && matches!(n1, b'A' | b'C' | b'T') {
                    Some(TriContext::Chg)
                } else {
                    Some(TriContext::Chh)
                }
            }
        }
        b'G' => {
            let p1 = byte_at(seq, pos as isize - 1);
            if p1 == b'C' {
                Some(TriContext::Cpg)
            } else {
                let p2 = byte_at(seq, pos as isize - 2);
                if p2 == b'C' && matches!(p1, b'A' | b'G' | b'T') {
                    Some(TriContext::Chg)
                } else {
                    Some(TriContext::Chh)
                }
            }
        }
        _ => None,
    }
}

/// CG vs CH at 0-based `pos`.
pub fn classify_context(seq: &[u8], pos: usize) -> Option<Context> {
    match classify_trinucleotide(seq, pos)? {
        TriContext::Cpg => Some(Context::Cg),
        TriContext::Chg | TriContext::Chh => Some(Context::Ch),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn top_and_bottom_strand() {
        assert_eq!(classify_trinucleotide(b"CGA", 0), Some(TriContext::Cpg));
        assert_eq!(classify_context(b"CGA", 0), Some(Context::Cg));
        assert_eq!(classify_trinucleotide(b"CAG", 0), Some(TriContext::Chg));
        assert_eq!(classify_context(b"CAG", 0), Some(Context::Ch));
        assert_eq!(classify_trinucleotide(b"CTA", 0), Some(TriContext::Chh));
        assert_eq!(classify_trinucleotide(b"CCC", 0), Some(TriContext::Chh));

        assert_eq!(classify_trinucleotide(b"TCG", 2), Some(TriContext::Cpg));
        assert_eq!(classify_context(b"TCG", 2), Some(Context::Cg));
        assert_eq!(classify_trinucleotide(b"CTG", 2), Some(TriContext::Chg));
        assert_eq!(classify_trinucleotide(b"TAG", 2), Some(TriContext::Chh));

        assert_eq!(classify_context(b"C", 0), Some(Context::Ch));
        assert_eq!(classify_context(b"G", 0), Some(Context::Ch));
        assert_eq!(classify_context(b"A", 0), None);
        assert_eq!(classify_context(b"CGA", 5), None);
    }
}
