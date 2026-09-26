//! Astair-compatible SAM flag classification and TAPS mC-to-T calling.

/// Default 5' sequenced-base trims (R2 M-bias; R1 typically clean).
pub const DEFAULT_TRIM_R1: u32 = 0;
pub const DEFAULT_TRIM_R2: u32 = 10;

const BAM_FREVERSE: u16 = 16;
const BAM_FREAD2: u16 = 128;

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

/// Keep this pileup base if it is past the 5' sequenced-base trim.
///
/// `qpos` is 0-based into BAM SEQ (reference-oriented). Reverse-strand records
/// store SEQ reverse-complemented, so the original sequenced 5' end is at the
/// high end of SEQ. R2 uses `trim_r2` (default 10); everything else uses
/// `trim_r1` (default 0).
#[inline]
pub fn keep_after_5prime_trim(
    qpos: usize,
    query_len: usize,
    flags: u16,
    trim_r1: u32,
    trim_r2: u32,
) -> bool {
    let trim = if flags & BAM_FREAD2 != 0 {
        trim_r2
    } else {
        trim_r1
    };
    if trim == 0 {
        return true;
    }
    let n = trim as usize;
    if n >= query_len {
        return false;
    }
    if flags & BAM_FREVERSE != 0 {
        qpos < query_len - n
    } else {
        qpos >= n
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

    #[test]
    fn five_prime_trim_r2_forward() {
        // flag 163 = R2 forward. First 10 SEQ bases are the sequenced 5'.
        let flag = 163u16;
        let len = 49usize;
        for q in 0..10 {
            assert!(!keep_after_5prime_trim(q, len, flag, 0, 10), "qpos={q}");
        }
        assert!(keep_after_5prime_trim(10, len, flag, 0, 10));
        assert!(keep_after_5prime_trim(48, len, flag, 0, 10));
    }

    #[test]
    fn five_prime_trim_r2_reverse() {
        // flag 147 = R2 reverse. Sequenced 5' is at the high end of SEQ.
        let flag = 147u16;
        let len = 49usize;
        assert!(keep_after_5prime_trim(0, len, flag, 0, 10));
        assert!(keep_after_5prime_trim(38, len, flag, 0, 10));
        for q in 39..49 {
            assert!(!keep_after_5prime_trim(q, len, flag, 0, 10), "qpos={q}");
        }
    }

    #[test]
    fn five_prime_trim_r1_default_keeps_all() {
        for flag in [99u16, 83] {
            for q in 0..39 {
                assert!(keep_after_5prime_trim(q, 39, flag, 0, 10));
            }
        }
    }

    #[test]
    fn five_prime_trim_r1_reverse() {
        // flag 83 = R1 reverse, trim 5: skip last 5 SEQ bases.
        let flag = 83u16;
        assert!(keep_after_5prime_trim(0, 39, flag, 5, 10));
        assert!(keep_after_5prime_trim(33, 39, flag, 5, 10));
        assert!(!keep_after_5prime_trim(34, 39, flag, 5, 10));
        assert!(!keep_after_5prime_trim(38, 39, flag, 5, 10));
    }

    #[test]
    fn five_prime_trim_zero_and_full() {
        assert!(keep_after_5prime_trim(0, 49, 163, 0, 0));
        assert!(!keep_after_5prime_trim(0, 10, 163, 0, 10));
        assert!(!keep_after_5prime_trim(9, 10, 163, 0, 10));
    }
}
