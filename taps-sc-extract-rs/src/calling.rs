//! Astair-compatible SAM flag classification and TAPS mC-to-T calling.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Ot,
    Ob,
}

/// Classify a SAM flag into OT / OB. Anything else is excluded (unmapped,
/// orphan, secondary, supplementary, improper pair).
pub fn classify_strand(flag: u16) -> Option<Strand> {
    match flag {
        99 | 147 => Some(Strand::Ot),
        83 | 163 => Some(Strand::Ob),
        _ => None,
    }
}

/// TAPS mC-to-T: 1 = methylated, 0 = unmethylated, None = non-informative.
///
/// OT + ref C: read T meth, read C unmeth.
/// OB + ref G: read A meth, read G unmeth.
pub fn call_mctot(strand: Strand, ref_base: u8, read_base: u8) -> Option<u8> {
    let r = ref_base.to_ascii_uppercase();
    let b = read_base.to_ascii_uppercase();
    match (strand, r, b) {
        (Strand::Ot, b'C', b'T') => Some(1),
        (Strand::Ot, b'C', b'C') => Some(0),
        (Strand::Ob, b'G', b'A') => Some(1),
        (Strand::Ob, b'G', b'G') => Some(0),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flag_ot_ob() {
        assert_eq!(classify_strand(99), Some(Strand::Ot));
        assert_eq!(classify_strand(147), Some(Strand::Ot));
        assert_eq!(classify_strand(83), Some(Strand::Ob));
        assert_eq!(classify_strand(163), Some(Strand::Ob));
        for f in [0u16, 16, 4, 355, 2048, 77, 141, 65, 129, 113, 177] {
            assert_eq!(classify_strand(f), None);
        }
    }

    #[test]
    fn mctot_ot_ob() {
        assert_eq!(call_mctot(Strand::Ot, b'C', b'T'), Some(1));
        assert_eq!(call_mctot(Strand::Ot, b'C', b'C'), Some(0));
        assert_eq!(call_mctot(Strand::Ot, b'C', b't'), Some(1));
        assert_eq!(call_mctot(Strand::Ot, b'C', b'A'), None);
        assert_eq!(call_mctot(Strand::Ot, b'G', b'A'), None);
        assert_eq!(call_mctot(Strand::Ob, b'G', b'A'), Some(1));
        assert_eq!(call_mctot(Strand::Ob, b'G', b'G'), Some(0));
        assert_eq!(call_mctot(Strand::Ob, b'G', b'T'), None);
        assert_eq!(call_mctot(Strand::Ob, b'C', b'T'), None);
    }
}
