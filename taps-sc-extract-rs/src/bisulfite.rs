//! Bisulfite / EM-seq calling (unmethylated C reads as T). Kept apart from the
//! TAPS rules in `calling.rs`, which are the inverse; see `chemistry.rs`.

use crate::calling::Strand;

/// Single-end directional libraries (BSBolt `W_C2T` / `C_C2T`): primary
/// forward reads are the original top strand, primary reverse reads the
/// original bottom strand. Paired, secondary and supplementary flags are not
/// classified (supplementary reads are filtered before the pileup).
#[inline]
pub fn classify_strand(flag: u16) -> Option<Strand> {
    match flag {
        0 => Some(Strand::Ot),
        16 => Some(Strand::Ob),
        _ => None,
    }
}

/// OT at reference C: read C = methylated (1), read T = unmethylated (0).
/// OB at reference G: read G = methylated, read A = unmethylated.
#[inline]
pub fn call(strand: Strand, ref_base: u8, read_base: u8) -> Option<u8> {
    let r = ref_base.to_ascii_uppercase();
    let b = read_base.to_ascii_uppercase();
    match (strand, r, b) {
        (Strand::Ot, b'C', b'C') => Some(1),
        (Strand::Ot, b'C', b'T') => Some(0),
        (Strand::Ob, b'G', b'G') => Some(1),
        (Strand::Ob, b'G', b'A') => Some(0),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bisulfite_is_the_inverse_of_taps() {
        for (s, r, b) in [(Strand::Ot, b'C', b'C'), (Strand::Ot, b'C', b'T'), (Strand::Ob, b'G', b'G'), (Strand::Ob, b'G', b'A')] {
            let taps = crate::calling::call_mctot(s, r, b).unwrap();
            assert_eq!(call(s, r, b), Some(1 - taps));
        }
        assert_eq!(call(Strand::Ot, b'C', b'c'), Some(1));
        assert_eq!(call(Strand::Ot, b'C', b'N'), None);
        assert_eq!(call(Strand::Ot, b'G', b'G'), None);
        assert_eq!(call(Strand::Ob, b'G', b'T'), None);
    }

    #[test]
    fn single_end_strands_only() {
        assert_eq!(classify_strand(0), Some(Strand::Ot));
        assert_eq!(classify_strand(16), Some(Strand::Ob));
        for f in [99u16, 147, 83, 163, 2048, 2064, 256, 4, 1] {
            assert_eq!(classify_strand(f), None, "{f}");
        }
        // TAPS is unchanged: SE flags are never TAPS strands.
        assert_eq!(crate::calling::classify_strand(0), None);
        assert_eq!(crate::calling::classify_strand(16), None);
    }
}
