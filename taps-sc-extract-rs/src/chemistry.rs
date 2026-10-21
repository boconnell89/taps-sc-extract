//! Which conversion chemistry the reads carry, and the strand and call rules
//! that go with it. TAPS (`calling`) and bisulfite (`bisulfite`) keep separate
//! tables and strand classifiers; nothing outside this dispatch mixes them.

use crate::calling::Strand;

#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum Chemistry {
    /// TAPS (mC -> T), paired-end directional (astair flags 99/147/83/163).
    Taps,
    /// Bisulfite / EM-seq (unmethylated C -> T), single-end directional
    /// (BSBolt: primary flag 0 = OT, 16 = OB).
    Bisulfite,
}

impl Chemistry {
    /// Strand of a record, or `None` if this chemistry never calls it.
    #[inline]
    pub fn strand(self, flag: u16) -> Option<Strand> {
        match self {
            Chemistry::Taps => crate::calling::classify_strand(flag),
            Chemistry::Bisulfite => crate::bisulfite::classify_strand(flag),
        }
    }

    /// 1 = methylated, 0 = unmethylated, `None` = no call.
    #[inline]
    pub fn call(self, strand: Strand, ref_base: u8, read_base: u8) -> Option<u8> {
        match self {
            Chemistry::Taps => crate::calling::call_mctot(strand, ref_base, read_base),
            Chemistry::Bisulfite => crate::bisulfite::call(strand, ref_base, read_base),
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Chemistry::Taps => "taps",
            Chemistry::Bisulfite => "bisulfite",
        }
    }
}
