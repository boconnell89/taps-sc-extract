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

    /// Reads one molecule can put at a position: both mates for paired-end
    /// TAPS, one read for single-end bisulfite (`--max-depth auto`).
    pub fn reads_per_molecule(self) -> u32 {
        match self {
            Chemistry::Taps => 2,
            Chemistry::Bisulfite => 1,
        }
    }
}

/// `--chemistry`: what the user asked for (`auto` resolves from the BAM).
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum ChemistryArg {
    Auto,
    Taps,
    Bisulfite,
}

/// `--strand-source` for bisulfite: `auto` takes the strand from the flag and
/// requires every read's BSBolt `YS` tag to agree (the run stops at the first
/// disagreement); `flag` skips the YS check (aligners without YS tags).
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum StrandSource {
    Auto,
    Flag,
}

/// What a sample of the BAM says about its chemistry.
#[derive(Debug, Default)]
pub struct Detected {
    /// Primary mapped reads sampled.
    pub sampled: usize,
    /// ... of which paired (`BAM_FPAIRED`).
    pub paired: usize,
    /// ... of which carry a BSBolt `YS` tag.
    pub with_ys: usize,
    /// BSBolt version from the header (`@PG VN:BSB-x.y.z-...` or a bsbolt command line).
    pub bsbolt: Option<String>,
}

/// Reads sampled for detection (primary, mapped).
pub const DETECT_SAMPLE: usize = 10_000;

/// Sample the start of the BAM. With `verify_ys`, errors on any `*_G2A` YS
/// (non-directional / PBAT libraries are not supported) or a YS that
/// disagrees with the flag.
pub fn detect(bam: &std::path::Path, verify_ys: bool) -> anyhow::Result<Detected> {
    use rust_htslib::bam::Read;
    let mut r = crate::input::open_indexed(bam)?;
    let header = String::from_utf8_lossy(r.header().as_bytes()).into_owned();
    let mut d = Detected {
        bsbolt: header.lines().filter(|l| l.starts_with("@PG")).find_map(|l| {
            l.split('\t').find_map(|f| f.strip_prefix("VN:BSB-").map(|v| v.split('-').next().unwrap_or(v).to_string()))
                .or_else(|| l.to_ascii_lowercase().contains("bsbolt").then(|| "unknown".to_string()))
        }),
        ..Default::default()
    };
    r.fetch(rust_htslib::bam::FetchDefinition::All).map_err(|e| anyhow::anyhow!("read {}: {e}", bam.display()))?;
    let mut rec = rust_htslib::bam::Record::new();
    while d.sampled < DETECT_SAMPLE {
        match r.read(&mut rec) {
            None => break,
            Some(Err(e)) => anyhow::bail!("read {}: {e}", bam.display()),
            Some(Ok(())) => {}
        }
        let f = rec.flags();
        if f & (0x4 | 0x100 | 0x800) != 0 {
            continue;
        }
        d.sampled += 1;
        if f & 0x1 != 0 {
            d.paired += 1;
        }
        if let Some(ys) = ys_tag(&rec) {
            d.with_ys += 1;
            if verify_ys {
                check_ys(&rec, &ys)?;
            }
        }
    }
    Ok(d)
}

/// The BSBolt `YS` tag, if present.
pub fn ys_tag(rec: &rust_htslib::bam::Record) -> Option<String> {
    match rec.aux(b"YS") {
        Ok(rust_htslib::bam::record::Aux::String(s)) => Some(s.to_string()),
        _ => None,
    }
}

/// A directional YS (`W_C2T` / `C_C2T`) that matches the flag's strand.
pub fn check_ys(rec: &rust_htslib::bam::Record, ys: &str) -> anyhow::Result<()> {
    let name = || String::from_utf8_lossy(rec.qname()).into_owned();
    if ys.contains("G2A") {
        anyhow::bail!(
            "read {} has YS:Z:{ys}: non-directional or PBAT bisulfite libraries are not supported",
            name()
        );
    }
    let reverse = rec.flags() & 0x10 != 0;
    let ok = matches!((ys, reverse), ("W_C2T", false) | ("C_C2T", true));
    if !ok {
        anyhow::bail!(
            "read {} has YS:Z:{ys} but flag {} ({}): the strand tag and the alignment disagree \
             (use --strand-source flag to trust the flag)",
            name(),
            rec.flags(),
            if reverse { "reverse" } else { "forward" }
        );
    }
    Ok(())
}

/// Resolve `--chemistry` against what the BAM looks like.
pub fn resolve(arg: ChemistryArg, d: &Detected, strand_source: StrandSource) -> anyhow::Result<Chemistry> {
    let bs_tags = d.with_ys > 0;
    let all_se = d.sampled > 0 && d.paired == 0;
    match arg {
        ChemistryArg::Auto => {
            if bs_tags {
                Ok(Chemistry::Bisulfite)
            } else if d.paired > 0 {
                Ok(Chemistry::Taps)
            } else {
                anyhow::bail!(
                    "cannot tell the chemistry from the BAM ({} single-end reads sampled, no BSBolt YS tags); \
                     pass --chemistry taps or --chemistry bisulfite",
                    d.sampled
                )
            }
        }
        ChemistryArg::Taps => {
            if bs_tags {
                anyhow::bail!("--chemistry taps, but the BAM carries BSBolt bisulfite tags (YS); use --chemistry bisulfite");
            }
            if all_se {
                anyhow::bail!("--chemistry taps needs paired-end reads, but all {} sampled reads are single-end", d.sampled);
            }
            Ok(Chemistry::Taps)
        }
        ChemistryArg::Bisulfite => {
            if d.paired > 0 {
                anyhow::bail!("--chemistry bisulfite supports single-end reads only, but the BAM has paired reads");
            }
            if !bs_tags && d.sampled > 0 && strand_source == StrandSource::Auto {
                anyhow::bail!(
                    "--chemistry bisulfite, but no read carries a BSBolt YS strand tag; \
                     pass --strand-source flag to take the strand from the alignment flag"
                );
            }
            Ok(Chemistry::Bisulfite)
        }
    }
}

/// Per-chemistry defaults, used for every option the user did not set.
pub struct Defaults {
    pub min_mapq: u8,
    pub min_baseq: u8,
    pub max_depth: u32,
    pub trim_r1: u32,
    pub trim_r2: u32,
    pub trim_3p: u32,
    pub baq: bool,
}

impl Chemistry {
    pub fn defaults(self) -> Defaults {
        match self {
            // TAPS: unchanged historical defaults.
            Chemistry::Taps => Defaults {
                min_mapq: 0,
                min_baseq: 20,
                max_depth: crate::extract::DEFAULT_MAX_DEPTH,
                trim_r1: crate::calling::DEFAULT_TRIM_R1,
                trim_r2: crate::calling::DEFAULT_TRIM_R2,
                trim_3p: 0,
                baq: true,
            },
            // Bisulfite SE (BSBolt, Ultima): drop MAPQ 0 (not bisulfite-unique)
            // and 1-9; trim 10 bp off the sequenced 3' end (Tn5 9 bp gap,
            // filled with unmethylated C at end repair); no BAQ (BSBolt,
            // Bismark and MethylDackel do not BAQ). Candidates until the filter
            // report on deduplicated data (code/bisulfite-extraction-plan.md).
            Chemistry::Bisulfite => Defaults {
                min_mapq: 10,
                min_baseq: 20,
                max_depth: crate::extract::DEFAULT_MAX_DEPTH,
                trim_r1: 0,
                trim_r2: 0,
                trim_3p: 10,
                baq: false,
            },
        }
    }
}

/// The FASTA must be the unconverted genome (never a C->T / G->A converted
/// index): every contig used must exist with the BAM's length, and about
/// 100 kb of A/C/G/T from the middle of the first must have both C and G
/// at >= 5% (unconverted mammalian ~20%, converted ~0%).
pub fn check_reference(
    fai: &crate::fasta::FastFaiReader,
    bam: &std::path::Path,
    contigs: &[&str],
) -> anyhow::Result<()> {
    use rust_htslib::bam::Read;
    let r = crate::input::open_indexed(bam)?;
    let h = r.header();
    for &c in contigs {
        let Some(tid) = h.tid(c.as_bytes()) else { continue };
        let bam_len = h.target_len(tid).unwrap_or(0);
        match fai.reference_length(c) {
            None => anyhow::bail!("contig {c} is in the BAM but not in the FASTA"),
            Some(l) if l != bam_len => anyhow::bail!("contig {c}: FASTA length {l} differs from the BAM header ({bam_len})"),
            _ => {}
        }
    }
    let Some(&first) = contigs.first() else { return Ok(()) };
    let Some(len) = fai.reference_length(first) else { return Ok(()) };
    let (mut acgt, mut c, mut g) = (0u64, 0u64, 0u64);
    let mut pos = len / 2;
    while acgt < 100_000 && pos < len {
        let end = (pos + 100_000).min(len);
        for b in fai.fetch(first, pos, end)? {
            match b {
                b'C' | b'c' => { c += 1; acgt += 1 }
                b'G' | b'g' => { g += 1; acgt += 1 }
                b'A' | b'a' | b'T' | b't' => acgt += 1,
                _ => {}
            }
        }
        pos = end;
    }
    if acgt >= 10_000 && (c.min(g) as f64) < 0.05 * acgt as f64 {
        anyhow::bail!(
            "the FASTA looks bisulfite-converted on {first} (C {:.1}%, G {:.1}% of {acgt} bases): \
             pass the unconverted genome the aligner index was built from",
            100.0 * c as f64 / acgt as f64,
            100.0 * g as f64 / acgt as f64
        );
    }
    Ok(())
}
