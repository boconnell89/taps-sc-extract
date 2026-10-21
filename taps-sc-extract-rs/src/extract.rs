//! Column-wise window extraction matching Python `_process_chunk`.
//!
//! rust-htslib's `IndexedReader::pileup()` is a single-file `bam_plp` iterator:
//! it does not run BAQ and does not call `bam_plp_init_overlaps`. This module
//! pushes records into `bam_plp` after samtools-style filters, matching pysam
//! `stepper="samtools"`:
//!
//! - flag filter: UNMAP | SECONDARY | QCFAIL | DUP
//! - `ignore_orphans`: skip paired-but-not-proper
//! - `compute_baq`: `sam_prob_realn(b, full_contig, len, BAQ_APPLY|BAQ_EXTEND)`
//!   (flag 3). The reference must start at contig coordinate 0. With
//!   the TAPS modes the reference has the strand's convertible base masked to
//!   N, and the default `BaqMode::Taps` skips reads that already match it
//!   exactly; see `BaqMode`.
//! - `ignore_overlaps`: `bam_plp_init_overlaps` (zeros the lower-quality
//!   overlapping mate base; pysam then drops it via `min_base_quality`)
//!
//! Python's current extractor sets `compute_baq=True` but never passes
//! `fastafile=` to `AlignmentFile.pileup()`, so BAQ is a no-op there. Rust
//! actually applies BAQ when `compute_baq` is true because we always have FASTA.

use crate::accumulate::{ensure_cell, merge_window_cells, BarcodeIntern, CellMaps};
use crate::baq_simd::{self, BaqJob};
use crate::barcode::barcode_from_qname_bytes;
use crate::calling::{
    call_mctot, classify_strand, keep_after_5prime_trim, DEFAULT_TRIM_R1, DEFAULT_TRIM_R2,
};
use crate::context::{classify_trinucleotide, Context, TriContext};
use crate::fasta::FastFaiReader;
use crate::window::Window;
use anyhow::{Context as AnyhowContext, Result};
use rust_htslib::bam::pileup::Alignment;
use rust_htslib::bam::{IndexedReader, Read, Record};
use rust_htslib::htslib;
use std::cell::Cell;
use std::collections::HashSet;
use std::os::raw::{c_char, c_int};
use std::path::Path;
use std::time::Instant;

#[derive(Clone, Debug, Default)]
pub struct CallStats {
    pub cpg_c: u64,
    pub cpg_t: u64,
    pub chg_c: u64,
    pub chg_t: u64,
    pub chh_c: u64,
    pub chh_t: u64,
    pub cg_c: u64,
    pub cg_t: u64,
    pub ch_c: u64,
    pub ch_t: u64,
}

impl CallStats {
    pub fn add(&mut self, tri: TriContext, ctx: Context, meth: u8) {
        let (c, t) = if meth == 1 { (1, 0) } else { (0, 1) };
        match tri {
            TriContext::Cpg => {
                self.cpg_c += c;
                self.cpg_t += t;
            }
            TriContext::Chg => {
                self.chg_c += c;
                self.chg_t += t;
            }
            TriContext::Chh => {
                self.chh_c += c;
                self.chh_t += t;
            }
        }
        match ctx {
            Context::Cg => {
                self.cg_c += c;
                self.cg_t += t;
            }
            Context::Ch => {
                self.ch_c += c;
                self.ch_t += t;
            }
        }
    }

    pub fn merge(&mut self, other: &CallStats) {
        self.cpg_c += other.cpg_c;
        self.cpg_t += other.cpg_t;
        self.chg_c += other.chg_c;
        self.chg_t += other.chg_t;
        self.chh_c += other.chh_c;
        self.chh_t += other.chh_t;
        self.cg_c += other.cg_c;
        self.cg_t += other.cg_t;
        self.ch_c += other.ch_c;
        self.ch_t += other.ch_t;
    }

    pub fn cg_pct(&self) -> f64 {
        let tot = self.cg_c + self.cg_t;
        if tot == 0 {
            0.0
        } else {
            100.0 * self.cg_c as f64 / tot as f64
        }
    }

    pub fn ch_pct(&self) -> f64 {
        let tot = self.ch_c + self.ch_t;
        if tot == 0 {
            0.0
        } else {
            100.0 * self.ch_c as f64 / tot as f64
        }
    }
}

pub struct WindowResult {
    pub stats: CallStats,
    pub intern: BarcodeIntern,
    /// Parallel to `intern.names()`. Empty when `ExtractParams.accumulate` is false.
    pub cells: Vec<CellMaps>,
    /// Filled when `ExtractParams.profile` is true (record counts always).
    pub timing: WindowTiming,
    /// Per-window `(chunk_id, timing)` collected by the parallel drivers when profiling.
    pub window_timings: Vec<(u32, WindowTiming)>,
}

/// Where one window's time went (seconds on the worker thread).
#[derive(Clone, Debug, Default)]
pub struct WindowTiming {
    /// BAM records read in the window (after `fetch`).
    pub records: u64,
    /// Records that passed the filters and entered the pileup.
    pub pushed: u64,
    /// Records that ran the BAQ HMM.
    pub baq_reads: u64,
    /// Records the `--max-depth` cap dropped (read-centric calling only).
    pub depth_dropped: u64,
    /// Read batches whose BAQ was sent to the GPU (the rest ran on the CPU).
    pub baq_gpu_batches: u64,
    /// Read batches that had BAQ work.
    pub baq_batches: u64,
    pub total_s: f64,
    /// FASTA window fetch plus BAM index seek.
    pub setup_s: f64,
    /// `bam.read` (BGZF inflate, record decode) and the flag/MAPQ/orphan filters.
    pub read_s: f64,
    /// BAQ (`baq_simd::Workspace::run`: SIMD kernel or `sam_prob_realn`).
    pub baq_s: f64,
    /// `bam_plp_push` / `bam_plp_next` and overlap clipping, excluding calling.
    pub pileup_s: f64,
    /// Per-column calling callback: context, strand, trim, barcode, cell maps.
    pub call_s: f64,
}

impl WindowTiming {
    pub fn add(&mut self, o: &WindowTiming) {
        self.records += o.records;
        self.pushed += o.pushed;
        self.baq_reads += o.baq_reads;
        self.depth_dropped += o.depth_dropped;
        self.baq_gpu_batches += o.baq_gpu_batches;
        self.baq_batches += o.baq_batches;
        self.total_s += o.total_s;
        self.setup_s += o.setup_s;
        self.read_s += o.read_s;
        self.baq_s += o.baq_s;
        self.pileup_s += o.pileup_s;
        self.call_s += o.call_s;
    }
}

impl WindowResult {
    pub fn empty() -> Self {
        Self {
            stats: CallStats::default(),
            intern: BarcodeIntern::default(),
            cells: Vec::new(),
            timing: WindowTiming::default(),
            window_timings: Vec::new(),
        }
    }

    pub fn n_cells(&self) -> usize {
        self.intern.len()
    }

    /// Merge counts and interned barcode names (not per-pos maps).
    pub fn absorb_counts(&mut self, other: &WindowResult) {
        self.stats.merge(&other.stats);
        self.timing.add(&other.timing);
        for name in other.intern.names() {
            self.intern.intern_str(name);
        }
    }

    /// Merge counts and per-pos cell maps (barcode-string intern remap).
    pub fn absorb_maps(&mut self, other: &WindowResult) {
        self.stats.merge(&other.stats);
        self.timing.add(&other.timing);
        merge_window_cells(&mut self.intern, &mut self.cells, &other.intern, &other.cells);
    }
}

#[derive(Clone, Debug)]
pub struct ExtractParams {
    pub min_baseq: u8,
    pub min_mapq: u8,
    pub max_depth: u32,
    pub ignore_overlaps: bool,
    pub ignore_orphans: bool,
    pub compute_baq: bool,
    /// Reference BAQ realigns against (when `compute_baq`).
    pub baq_mode: BaqMode,
    /// BAQ kernel: SIMD width or htslib.
    pub baq_width: baq_simd::Width,
    /// Where the BAQ HMM runs (`--baq-device`); the GPU service must be started.
    pub baq_device: crate::baq_gpu::BaqDevice,
    /// Contexts to extract; columns of the others are skipped.
    pub context: ContextSel,
    pub decomp_threads: usize,
    /// Fill per-barcode (pos, t, c) maps. Off for stats-only (saves RAM).
    pub accumulate: bool,
    /// Drop this many 5' sequenced bases of read 1 (default 0).
    pub trim_r1: u32,
    /// Drop this many 5' sequenced bases of read 2 (default 10).
    pub trim_r2: u32,
    /// Time each stage of every window (`WindowResult.timing`).
    pub profile: bool,
    /// How calls are made (`--caller`).
    pub caller: Caller,
}

/// `--caller`: column pileup (htslib `bam_plp`) or read-centric calling
/// (`readcall`), which replays the pileup's depth cap and mate-overlap rules
/// and gives the same calls.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum Caller {
    /// One read at a time (`readcall::ReadPileup`).
    Reads,
    /// htslib `bam_plp` columns.
    Pileup,
}

impl Default for ExtractParams {
    fn default() -> Self {
        Self {
            min_baseq: 20,
            min_mapq: 0,
            max_depth: DEFAULT_MAX_DEPTH,
            ignore_overlaps: true,
            ignore_orphans: true,
            compute_baq: true,
            baq_mode: BaqMode::Taps,
            baq_width: baq_simd::Width::cpu_best(),
            baq_device: crate::baq_gpu::BaqDevice::Cpu,
            context: ContextSel::Both,
            decomp_threads: 1,
            accumulate: false,
            trim_r1: DEFAULT_TRIM_R1,
            trim_r2: DEFAULT_TRIM_R2,
            profile: false,
            caller: Caller::Reads,
        }
    }
}

/// `--context`: which methylation contexts to extract.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum ContextSel {
    /// CG and CH (default).
    Both,
    /// CpG only: CH columns are skipped before any per-read work, and no /CH
    /// group is written.
    #[value(name = "CG", alias = "cg")]
    Cg,
    /// CH only (CHG + CHH); no /CG group is written.
    #[value(name = "CH", alias = "ch")]
    Ch,
}

impl ContextSel {
    #[inline]
    pub fn wants(self, c: Context) -> bool {
        match self {
            ContextSel::Both => true,
            ContextSel::Cg => c == Context::Cg,
            ContextSel::Ch => c == Context::Ch,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            ContextSel::Both => "CG and CH",
            ContextSel::Cg => "CG only",
            ContextSel::Ch => "CH only",
        }
    }

    /// Share of methylation sites in the selected contexts (mm10 test data:
    /// CG ~4% of C/G cell-sites, CH ~96%), for memory estimates.
    pub fn site_share(self) -> f64 {
        match self {
            ContextSel::Both => 1.0,
            ContextSel::Cg => 0.04,
            ContextSel::Ch => 0.96,
        }
    }
}

/// Reference that BAQ (`sam_prob_realn`) realigns each read against.
///
/// In TAPS a methylated C is read as T (A on the bottom strand), so against the
/// plain reference every methylated C is a mismatch. htslib's HMM reads that as
/// possible misalignment and lowers base qualities at and around methylated
/// sites, which drops methylated calls more often than unmethylated ones
/// (chr19: -3.0% methylated vs -0.4% unmethylated CpG calls).
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum BaqMode {
    /// TAPS-masked BAQ (see `TapsFull`), skipped for reads that are all-M and
    /// match the masked reference at every base (~72% of reads). Full-genome
    /// test BAM: 172 of 25.9M CG cell-sites differ from `TapsFull`, processing
    /// ~2x faster.
    Taps,
    /// Mask the convertible base: OT reads see every reference C as N, OB reads
    /// every G (all contexts, CpG and CH). N scores as neither match nor
    /// mismatch in the HMM (`probaln.c`), so methylation state does not affect
    /// BAQ, while indels and real mismatches still do. Every read runs the HMM.
    TapsFull,
    /// TAPS-masked BAQ only for reads with an indel or clip in the CIGAR
    /// (~8%). ~4x faster than `TapsFull`; 0.06% of CG cell-sites differ.
    TapsPartial,
    /// htslib's plain reference (samtools/asTair behaviour; biased for TAPS).
    Htslib,
}

impl BaqMode {
    pub fn label(self) -> &'static str {
        match self {
            BaqMode::Taps => "taps",
            BaqMode::TapsFull => "taps-full",
            BaqMode::TapsPartial => "taps-partial",
            BaqMode::Htslib => "htslib",
        }
    }

    /// Run the HMM on this read? `r` covers the read; `mask` is the base BAQ
    /// sees as N for it (`baq_mask`).
    #[inline]
    fn wants(self, rec: &Record, r: &RefSlice, mask: u8) -> bool {
        match self {
            BaqMode::TapsFull | BaqMode::Htslib => true,
            BaqMode::TapsPartial => !cigar_all_match(rec),
            BaqMode::Taps => !(cigar_all_match(rec) && matches_masked_ref(rec, r, mask)),
        }
    }
}

/// The reference base BAQ sees as N for this read: C on OT reads, G on OB
/// reads in the TAPS modes; none for `htslib`. The SIMD kernel and the htslib
/// fallback apply it while reading the reference, so no masked copy is stored.
#[inline]
pub(crate) fn baq_mask(mode: BaqMode, flags: u16) -> u8 {
    match mode {
        BaqMode::Htslib => 0,
        _ if baq_is_ob(flags) => b'G',
        _ => b'C',
    }
}

/// Unmasked, uppercase reference bases `[lo, hi)` of one contig; `hi` is the
/// contig length when the slice reaches the contig end.
#[derive(Clone, Copy)]
pub(crate) struct RefSlice<'a> {
    pub seq: &'a [u8],
    pub lo: u64,
    pub hi: u64,
}

impl RefSlice<'_> {
    #[inline]
    fn covers(&self, lo: u64, hi: u64) -> bool {
        lo >= self.lo && hi <= self.hi
    }
}

/// Reference span `sam_prob_realn` can read for this read, clipped to the
/// contig: `[pos - slack, end + slack)` with slack = 2 * read length + ref span
/// + 16. Must match `taps_baq_need` in csrc/baq_simd.c.
pub(crate) fn baq_need(rec: &Record, contig_len: u64) -> (u64, u64) {
    let span: u64 = rec
        .raw_cigar()
        .iter()
        .filter(|&&c| matches!(c & 0xf, 0 | 2 | 3 | 7 | 8)) // M D N = X
        .map(|&c| u64::from(c >> 4))
        .sum();
    let slack = 2 * rec.seq_len() as u64 + span + 16;
    let pos = rec.pos().max(0) as u64;
    (pos.saturating_sub(slack), (pos + span + slack).min(contig_len))
}

/// CIGAR is only M operations (no clip, indel or skip).
#[inline]
fn cigar_all_match(rec: &Record) -> bool {
    rec.raw_cigar().iter().all(|&c| c & 0xf == htslib::BAM_CMATCH)
}

/// BAM 4-bit code (A=1 C=2 G=4 T=8) of each reference byte, 0 where the base
/// matches anything: non-ACGT bytes and the masked base.
const fn ref_codes(mask: u8) -> [u8; 256] {
    let mut t = [0u8; 256];
    t[b'A' as usize] = 1;
    t[b'C' as usize] = 2;
    t[b'G' as usize] = 4;
    t[b'T' as usize] = 8;
    if mask != 0 {
        t[mask as usize] = 0;
    }
    t
}
static REF_CODES_NONE: [u8; 256] = ref_codes(0);
static REF_CODES_C: [u8; 256] = ref_codes(b'C');
static REF_CODES_G: [u8; 256] = ref_codes(b'G');

/// Every read base equals the reference base where the reference is A/C/G/T
/// and not the masked base (masked positions match anything; a read N is a
/// mismatch). Assumes an all-M CIGAR, so read base i sits at `pos + i`.
/// Branch-free over the read: this runs on every read when BAQ is on.
#[inline]
fn matches_masked_ref(rec: &Record, r: &RefSlice, mask: u8) -> bool {
    let n = rec.seq_len();
    let Ok(pos) = u64::try_from(rec.pos()) else {
        return false;
    };
    if !r.covers(pos, pos + n as u64) {
        return false;
    }
    let codes = match mask {
        b'C' => &REF_CODES_C,
        b'G' => &REF_CODES_G,
        _ => &REF_CODES_NONE,
    };
    let bases = &r.seq[(pos - r.lo) as usize..(pos - r.lo) as usize + n];
    let raw = rec.seq().encoded;
    let mut bad = 0u8;
    // Two read bases per SEQ byte, high nibble first.
    for (pair, &byte) in bases.chunks_exact(2).zip(raw) {
        let (c0, c1) = (codes[pair[0] as usize], codes[pair[1] as usize]);
        bad |= u8::from(c0 != 0 && c0 != byte >> 4) | u8::from(c1 != 0 && c1 != byte & 0xf);
    }
    if n % 2 == 1 {
        let c = codes[bases[n - 1] as usize];
        bad |= u8::from(c != 0 && c != raw[n / 2] >> 4);
    }
    bad == 0
}

/// Bottom-strand (OB) read for TAPS-aware BAQ: read 1 reverse or read 2
/// forward (flags 83/163; unpaired reverse reads too). Everything else,
/// including 99/147, is OT.
#[inline]
pub(crate) fn baq_is_ob(flags: u16) -> bool {
    let reverse = flags & htslib::BAM_FREVERSE as u16 != 0;
    let read2 = flags & htslib::BAM_FPAIRED as u16 != 0 && flags & htslib::BAM_FREAD2 as u16 != 0;
    reverse != read2
}

/// pysam / samtools default: skip unmapped, secondary, QC fail, duplicate.
const FLAG_FILTER: u16 = (htslib::BAM_FUNMAP
    | htslib::BAM_FSECONDARY
    | htslib::BAM_FQCFAIL
    | htslib::BAM_FDUP) as u16;

/// Default `--max-depth`: htslib's per-position read cap, counted across all
/// cells. Single-cell depth at a site is roughly the number of cells covering
/// it, so 250 dropped whole cells at well-covered loci (805k reads on the
/// test BAM); 100,000 only guards against pathological repeats.
pub const DEFAULT_MAX_DEPTH: u32 = 100_000;

/// The cap passed to htslib / `ReadPileup`: 0 means unlimited. htslib takes a
/// signed int and treats a negative cap as "drop", so never pass `u32::MAX`.
pub fn effective_max_depth(max_depth: u32) -> i32 {
    if max_depth == 0 || max_depth > i32::MAX as u32 {
        i32::MAX
    } else {
        max_depth as i32
    }
}

/// Records read, filtered and BAQ'd together before they are pushed in order.
/// Large enough that most reads find same-shape partners for the SIMD lanes
/// (chr19 at the default skip: ~92% of 4-lane slots filled at 8,192).
const READ_BATCH: usize = 8192;

thread_local! {
    /// This worker's push-and-call seconds per `READ_BATCH` batch (GPU BAQ
    /// sharing: how long it is busy before it needs a batch back).
    static CALL_BATCH_S: Cell<f64> = const { Cell::new(0.0) };
}

struct PlpGuard(htslib::bam_plp_t);

impl Drop for PlpGuard {
    fn drop(&mut self) {
        if !self.0.is_null() {
            unsafe { htslib::bam_plp_destroy(self.0) }
        }
    }
}

// Not in the public htslib header (only bam_mplp_init_overlaps is), but
// present in libhts.a.
extern "C" {
    fn bam_plp_init_overlaps(iter: htslib::bam_plp_t) -> std::os::raw::c_int;
}

/// `--profile-json` times one record (and one pileup column) in this many.
const PROFILE_STRIDE: u64 = 16;

/// Filters. `None`: the read is not pushed; `Some(true)`: pushed and a BAQ
/// candidate (BAQ on and the cell kept); `Some(false)`: pushed without BAQ.
fn accept_record(
    rec: &Record,
    params: &ExtractParams,
    whitelist: Option<&HashSet<String>>,
) -> Option<bool> {
    let flags = rec.flags();
    if flags & FLAG_FILTER != 0 {
        return None;
    }
    // MAPQ and orphan tests only reject. A rejected read is never pushed, so
    // testing before BAQ leaves calls unchanged and skips the HMM for it.
    if rec.mapq() < params.min_mapq {
        return None;
    }
    if params.ignore_orphans
        && (flags & htslib::BAM_FPAIRED as u16) != 0
        && (flags & htslib::BAM_FPROPER_PAIR as u16) == 0
    {
        return None;
    }
    // A read from a filtered-out cell is still pushed (depth caps and overlap
    // pairing stay the same), but its BAQ is skipped: BAQ only changes this
    // read's own qualities, and overlap tweaks pair it only with its mate,
    // which has the same barcode. Kept cells' calls are unchanged.
    let excluded = whitelist.is_some_and(|wl| {
        !wl.contains(std::str::from_utf8(barcode_from_qname_bytes(rec.qname())).unwrap_or(""))
    });
    Some(params.compute_baq && !excluded)
}

/// Padding around a window's reference slice when BAQ is on. A short read's
/// BAQ span (`baq_need`) reaches a few hundred bases past the window; reads
/// that reach further get their own fetch.
const BAQ_PAD: u64 = 16_384;

/// Process one genomic window `[start, end)` with column-wise `bam_plp` pileup.
///
/// The window's reference (padded by `BAQ_PAD` when BAQ is on) serves both
/// context calls and BAQ; `sam_prob_realn`'s contig-coordinate contract is
/// kept by passing slice bounds (`BaqJob::ref_lo`/`ref_hi`).
pub fn process_window(
    bam: &mut IndexedReader,
    fasta: &FastFaiReader,
    window: &Window,
    whitelist: Option<&HashSet<String>>,
    params: &ExtractParams,
) -> Result<WindowResult> {
    let profile = params.profile;
    let t_window = Instant::now();
    let mut timing = WindowTiming::default();
    let pad = if params.compute_baq { BAQ_PAD } else { 2 };
    let ref_len = fasta
        .reference_length(&window.contig)
        .ok_or_else(|| anyhow::anyhow!("contig {} missing in FASTA", window.contig))?;
    let ref_start = window.start.saturating_sub(pad);
    let ref_end = (window.end + pad).min(ref_len);
    let ref_seq = fasta
        .fetch(&window.contig, ref_start, ref_end)
        .with_context(|| format!("FASTA fetch {}:{ref_start}-{ref_end}", window.contig))?;

    bam.fetch((window.contig.as_str(), window.start, window.end))
        .map_err(|e| anyhow::anyhow!("BAM fetch: {e}"))?;

    if profile {
        timing.setup_s = t_window.elapsed().as_secs_f64();
    }

    // func=NULL: we push records ourselves so overlap_push sees the same
    // bam1_t that rust-htslib just filled (no extra bam_copy1 into iter->b).
    let plp = unsafe { htslib::bam_plp_init(None, std::ptr::null_mut()) };
    if plp.is_null() {
        anyhow::bail!("bam_plp_init failed");
    }
    let _guard = PlpGuard(plp);
    unsafe {
        htslib::bam_plp_set_maxcnt(plp, effective_max_depth(params.max_depth));
        if params.ignore_overlaps && bam_plp_init_overlaps(plp) < 0 {
            anyhow::bail!("bam_plp_init_overlaps failed");
        }
    }

    let mut stats = CallStats::default();
    let mut intern = BarcodeIntern::default();
    let mut cells: Vec<CellMaps> = Vec::new();
    // Read-centric calling (`--caller reads`): same calls, no columns.
    let read_mode = params.caller == Caller::Reads;
    let mut rp = crate::readcall::ReadPileup::new(effective_max_depth(params.max_depth), params.ignore_overlaps);
    let cw = crate::readcall::CallWindow {
        ref_seq: &ref_seq,
        ref_start,
        win_start: window.start as u32,
        win_end: window.end as u32,
        params,
        whitelist,
    };
    let mut sink = crate::readcall::CallSink::default();
    let win_start = window.start as u32;
    let win_end = window.end as u32;

    // Cell id per pileup read, indexed by `bam1_t.id`: htslib numbers each read
    // it stores in this iterator 0, 1, 2, ... and never reuses a number. (The
    // record's address is reused once htslib recycles its node, so it cannot
    // identify a read.) Resolved on the read's first call.
    const CELL_UNRESOLVED: u32 = u32::MAX;
    const CELL_NONE: u32 = u32::MAX - 1;
    let mut read_cell: Vec<u32> = Vec::new();

    let mut consume_column = |pos: u32, n_plp: i32, col: *const htslib::bam_pileup1_t| {
        if pos < win_start || pos >= win_end || n_plp <= 0 || col.is_null() {
            return;
        }
        let rel = (pos as u64).saturating_sub(ref_start) as usize;
        if rel >= ref_seq.len() {
            return;
        }
        let ref_b = ref_seq[rel];
        let Some(tri) = classify_trinucleotide(&ref_seq, rel) else {
            return;
        };
        let ctx = match tri {
            TriContext::Cpg => Context::Cg,
            _ => Context::Ch,
        };
        // --context: skip the other context's columns before any per-read work.
        if !params.context.wants(ctx) {
            return;
        }
        let pile = unsafe { std::slice::from_raw_parts(col, n_plp as usize) };
        for p in pile {
            let aln = Alignment::new(p);
            if aln.is_del() || aln.is_refskip() {
                continue;
            }
            let Some(qpos) = aln.qpos() else {
                continue;
            };
            let rec = aln.record();
            let seq = rec.seq();
            if qpos >= seq.len() {
                continue;
            }
            if !keep_after_5prime_trim(
                qpos,
                seq.len(),
                rec.flags(),
                params.trim_r1,
                params.trim_r2,
            ) {
                continue;
            }
            let qual = rec.qual().get(qpos).copied().unwrap_or(0);
            if qual < params.min_baseq {
                continue;
            }
            let Some(strand) = classify_strand(rec.flags()) else {
                continue;
            };
            let Some(call) = call_mctot(strand, ref_b, seq[qpos]) else {
                continue;
            };

            let read_id = unsafe { (*p.b).id } as usize;
            if read_id >= read_cell.len() {
                read_cell.resize(read_id + 1, CELL_UNRESOLVED);
            }
            if read_cell[read_id] == CELL_UNRESOLVED {
                let bc_bytes = barcode_from_qname_bytes(rec.qname());
                let valid = match whitelist {
                    Some(wl) => wl.contains(std::str::from_utf8(bc_bytes).unwrap_or("")),
                    None => true,
                };
                read_cell[read_id] = if valid {
                    intern.intern_bytes(bc_bytes)
                } else {
                    CELL_NONE
                };
            }
            let id_opt = match read_cell[read_id] {
                CELL_NONE => None,
                id => Some(id),
            };

            let Some(id) = id_opt else {
                continue;
            };

            if params.accumulate {
                ensure_cell(&mut cells, id).add(ctx, pos + 1, call);
            }
            stats.add(tri, ctx, call);
        }
    };

    // Calling time when profiling: every PROFILE_STRIDE-th column is timed and
    // scaled up, so the clock is not read twice per column (that alone cost
    // ~70% extra wall time on a 64-core server).
    let call_ns = Cell::new(0u64);
    let n_col = Cell::new(0u64);
    let mut consume_column = |pos: u32, n_plp: i32, col: *const htslib::bam_pileup1_t| {
        if profile && n_col.replace(n_col.get() + 1) % PROFILE_STRIDE == 0 {
            let t = Instant::now();
            consume_column(pos, n_plp, col);
            call_ns.set(call_ns.get() + t.elapsed().as_nanos() as u64 * PROFILE_STRIDE);
        } else {
            consume_column(pos, n_plp, col);
        }
    };

    let n_read = Cell::new(0u64);
    let mut call_one = |r: &Record| {
        if profile && n_read.replace(n_read.get() + 1) % PROFILE_STRIDE == 0 {
            let t = Instant::now();
            crate::readcall::call_read(&cw, r, &mut sink);
            call_ns.set(call_ns.get() + t.elapsed().as_nanos() as u64 * PROFILE_STRIDE);
        } else {
            crate::readcall::call_read(&cw, r, &mut sink);
        }
    };

    let drain_ready = |plp: htslib::bam_plp_t, consume: &mut dyn FnMut(u32, i32, *const htslib::bam_pileup1_t)| -> Result<()> {
        loop {
            let mut tid = 0i32;
            let mut pos = 0i32;
            let mut n_plp = 0i32;
            let col = unsafe { htslib::bam_plp_next(plp, &mut tid, &mut pos, &mut n_plp) };
            if col.is_null() {
                if n_plp < 0 {
                    anyhow::bail!("bam_plp_next error");
                }
                break;
            }
            if pos >= 0 {
                consume(pos as u32, n_plp, col);
            }
        }
        Ok(())
    };

    // Records are read in batches: filter, BAQ the whole batch (SIMD kernel
    // across reads), then push and drain each record in file order exactly as
    // one-at-a-time processing would. BAQ changes only the read's own
    // qualities, so running it ahead of the pushes leaves the pileup unchanged.
    // Batches in flight. On the CPU one: read and filter a batch, BAQ it, push
    // it. With the GPU several: while the GPU runs BAQ on the newest batches,
    // the oldest finished batch is pushed and called, so calling overlaps the
    // device. Batches are pushed in file order either way.
    let gpu = params.compute_baq && params.baq_device == crate::baq_gpu::BaqDevice::Gpu;
    let depth = if gpu { crate::baq_gpu::depth() } else { 1 };
    struct Slot {
        pool: Vec<Record>,
        keep: Vec<bool>,
        n: usize,
        jobs: Vec<BaqJob>,
        // Reads whose BAQ span reaches past the padded window (rare): own
        // fetches, alive until their batch's BAQ has run.
        extras: Vec<Vec<u8>>,
        ws: baq_simd::Workspace,
        pending: Option<crate::baq_gpu::Pending>,
    }
    let mut slots: Vec<Slot> = (0..depth)
        .map(|_| Slot {
            pool: (0..READ_BATCH).map(|_| Record::new()).collect(),
            keep: vec![false; READ_BATCH],
            n: 0,
            jobs: Vec::with_capacity(READ_BATCH),
            extras: Vec::new(),
            ws: baq_simd::Workspace::new(),
            pending: None,
        })
        .collect();
    let mut free: Vec<usize> = (0..depth).rev().collect();
    let mut inflight: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
    let win_ref = RefSlice { seq: &ref_seq, lo: ref_start, hi: ref_end };
    // Stage clock: `lap` returns seconds since the previous lap when profiling.
    let mut mark = Instant::now();
    let mut lap = |on: bool| -> f64 {
        if !on {
            return 0.0;
        }
        let now = Instant::now();
        let d = now.duration_since(mark).as_secs_f64();
        mark = now;
        d
    };
    let mut plp_s = 0.0f64;
    let mut eof = false;
    while !(eof && inflight.is_empty()) {
        if !eof && !free.is_empty() {
            // Read, filter and start BAQ on a new batch.
            lap(profile);
            let si = free.pop().expect("checked");
            let slot = &mut slots[si];
            slot.jobs.clear();
            slot.extras.clear();
            let mut n = 0;
            while n < READ_BATCH {
                match bam.read(&mut slot.pool[n]) {
                    None => {
                        eof = true;
                        break;
                    }
                    Some(Err(e)) => anyhow::bail!("BAM read: {e}"),
                    Some(Ok(())) => {
                        timing.records += 1;
                        let acc = accept_record(&slot.pool[n], params, whitelist);
                        slot.keep[n] = acc.is_some();
                        if acc == Some(true) {
                            let rec = &mut slot.pool[n];
                            let (lo, hi) = baq_need(rec, ref_len);
                            let r = if win_ref.covers(lo, hi) {
                                win_ref
                            } else {
                                slot.extras.push(fasta.fetch(&window.contig, lo, hi).with_context(|| {
                                    format!("FASTA BAQ fetch {}:{lo}-{hi}", window.contig)
                                })?);
                                RefSlice { seq: slot.extras.last().expect("pushed above"), lo, hi }
                            };
                            let mask = baq_mask(params.baq_mode, rec.flags());
                            if params.baq_mode.wants(rec, &r, mask) {
                                slot.jobs.push(BaqJob {
                                    b: rec.inner_mut(),
                                    reference: r.seq.as_ptr() as *const c_char,
                                    ref_lo: r.lo as i64,
                                    ref_hi: r.hi as i64,
                                    mask: c_int::from(mask),
                                });
                            }
                        }
                        n += 1;
                    }
                }
            }
            slot.n = n;
            timing.read_s += lap(profile);
            timing.baq_reads += slot.jobs.len() as u64;
            timing.baq_batches += u64::from(!slot.jobs.is_empty());
            // With the GPU, a batch goes to it while it keeps up; otherwise
            // this worker runs the batch's BAQ itself (`baq_gpu::take_batch`).
            // Slack: this worker calls its other in-flight batches before it
            // needs this one back.
            let slack_s = CALL_BATCH_S.get() * inflight.len() as f64;
            if gpu && crate::baq_gpu::take_batch(slack_s) {
                slot.pending = Some(crate::baq_gpu::submit(&mut slot.ws, &slot.jobs));
                timing.baq_gpu_batches += 1;
            } else {
                slot.ws.run(&slot.jobs, params.baq_width);
            }
            timing.baq_s += lap(profile);
            if n > 0 {
                inflight.push_back(si);
            } else {
                free.push(si);
            }
            continue;
        }
        // Finish the oldest batch and push it into the pileup in file order.
        let si = inflight.pop_front().expect("loop condition");
        let slot = &mut slots[si];
        lap(profile);
        if let Some(p) = slot.pending.take() {
            crate::baq_gpu::complete(p, &mut slot.ws, &slot.jobs, params.baq_width);
        }
        timing.baq_s += lap(profile);
        let t_call = Instant::now();
        for (rec, _) in slot.pool[..slot.n].iter().zip(&slot.keep[..slot.n]).filter(|(_, k)| **k) {
            timing.pushed += 1;
            if read_mode {
                rp.push(rec, &mut call_one)?;
                continue;
            }
            if unsafe { htslib::bam_plp_push(plp, rec.inner()) } < 0 {
                anyhow::bail!("bam_plp_push failed");
            }
            drain_ready(plp, &mut consume_column)?;
        }
        if gpu {
            // Moving average of push-and-call time per full batch (for `take_batch`).
            let t = t_call.elapsed().as_secs_f64() * READ_BATCH as f64 / slot.n.max(1) as f64;
            CALL_BATCH_S.set(if CALL_BATCH_S.get() == 0.0 { t } else { 0.9 * CALL_BATCH_S.get() + 0.1 * t });
        }
        plp_s += lap(profile);
        free.push(si);
    }
    if profile {
        lap(true);
    }
    if read_mode {
        rp.finish(&mut call_one);
        timing.depth_dropped = rp.dropped;
        stats = sink.stats;
        intern = sink.intern;
        cells = sink.cells;
    } else {
        if unsafe { htslib::bam_plp_push(plp, std::ptr::null()) } < 0 {
            anyhow::bail!("bam_plp_push(EOF) failed");
        }
        drain_ready(plp, &mut consume_column)?;
    }
    plp_s += lap(profile);
    for cell in cells.iter_mut() {
        cell.normalize();
    }
    if profile {
        timing.call_s = call_ns.get() as f64 * 1e-9;
        timing.pileup_s = (plp_s - timing.call_s).max(0.0);
        timing.total_s = t_window.elapsed().as_secs_f64();
    }

    Ok(WindowResult {
        stats,
        intern,
        cells,
        timing,
        window_timings: Vec::new(),
    })
}

/// Mapped reads on `contigs`, and an estimated share of the BAM file those reads occupy.
pub struct BamTraffic {
    pub mapped_reads: u64,
    pub estimated_bytes: u64,
}

pub fn estimate_bam_traffic(bam_path: &Path, contigs: &[&str]) -> Result<BamTraffic> {
    let mut bam = crate::input::open_indexed(bam_path)
        .map_err(|e| anyhow::anyhow!("open BAM for idxstats {}: {e}", bam_path.display()))?;
    let stats = bam
        .index_stats()
        .map_err(|e| anyhow::anyhow!("BAM index stats: {e}"))?;
    let names: Vec<String> = bam
        .header()
        .target_names()
        .into_iter()
        .map(|n| String::from_utf8_lossy(n).into_owned())
        .collect();
    let want: std::collections::HashSet<&str> = contigs.iter().copied().collect();
    let mut mapped_reads = 0u64;
    let mut records = 0u64;
    for (tid, _len, mapped, unmapped) in stats {
        records += mapped + unmapped;
        if tid < 0 {
            continue;
        }
        if let Some(name) = names.get(tid as usize) {
            if want.contains(name.as_str()) {
                mapped_reads += mapped;
            }
        }
    }
    let file_len = std::fs::metadata(bam_path)
        .with_context(|| format!("stat {}", bam_path.display()))?
        .len();
    let estimated_bytes = if records == 0 {
        file_len
    } else {
        ((file_len as u128) * u128::from(mapped_reads) / u128::from(records)) as u64
    };
    Ok(BamTraffic {
        mapped_reads,
        estimated_bytes,
    })
}

pub fn open_bam(path: &Path, decomp_threads: usize) -> Result<IndexedReader> {
    let mut bam = crate::input::open_indexed(path)
        .map_err(|e| anyhow::anyhow!("open BAM {}: {e}", path.display()))?;
    if decomp_threads > 0 {
        bam.set_threads(decomp_threads)
            .map_err(|e| anyhow::anyhow!("BAM set_threads: {e}"))?;
    }
    Ok(bam)
}

pub fn extract_regions_stats(
    bam_path: &Path,
    fasta_path: &Path,
    windows: &[Window],
    whitelist: Option<&HashSet<String>>,
    params: &ExtractParams,
) -> Result<WindowResult> {
    let fasta = FastFaiReader::open(fasta_path)?;
    let mut bam = open_bam(bam_path, params.decomp_threads)?;
    let mut total = WindowResult::empty();
    for (i, w) in windows.iter().enumerate() {
        let r = process_window(&mut bam, &fasta, w, whitelist, params)?;
        total.stats.merge(&r.stats);
        if params.accumulate {
            merge_window_cells(&mut total.intern, &mut total.cells, &r.intern, &r.cells);
        } else {
            for name in r.intern.names() {
                total.intern.intern_str(name);
            }
        }
        if (i + 1) % 5 == 0 || i + 1 == windows.len() {
            eprintln!(
                "Rust extract: [{}/{}] cells={} CG={:.2}% CH={:.3}%",
                i + 1,
                windows.len(),
                total.n_cells(),
                total.stats.cg_pct(),
                total.stats.ch_pct()
            );
        }
    }
    Ok(total)
}

/// Single-thread stream extract: accumulate one window, write shard temp files, drop maps.
pub fn extract_windows_to_temp(
    bam_path: &Path,
    fasta_path: &Path,
    windows: &[Window],
    whitelist: Option<&HashSet<String>>,
    params: &ExtractParams,
    temp_dir: &Path,
    n_shards: usize,
) -> Result<WindowResult> {
    let fasta = FastFaiReader::open(fasta_path)?;
    let mut bam = open_bam(bam_path, params.decomp_threads)?;
    let mut total = WindowResult::empty();
    let mut n_files = 0usize;
    for (i, w) in windows.iter().enumerate() {
        let mut p = params.clone();
        p.accumulate = true;
        let r = process_window(&mut bam, &fasta, w, whitelist, &p)?;
        n_files += crate::shard_io::write_window_shards(temp_dir, w, &r.intern, &r.cells, n_shards)?;
        total.stats.merge(&r.stats);
        for name in r.intern.names() {
            total.intern.intern_str(name);
        }
        if (i + 1) % 5 == 0 || i + 1 == windows.len() {
            eprintln!(
                "Rust extract: [{}/{}] cells={} files={n_files} CG={:.2}% CH={:.3}%",
                i + 1,
                windows.len(),
                total.n_cells(),
                total.stats.cg_pct(),
                total.stats.ch_pct()
            );
        }
    }
    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn max_depth_zero_is_unlimited_and_never_negative() {
        assert_eq!(effective_max_depth(0), i32::MAX);
        assert_eq!(effective_max_depth(250), 250);
        assert_eq!(effective_max_depth(u32::MAX), i32::MAX);
        assert_eq!(DEFAULT_MAX_DEPTH, 100_000);
    }

    #[test]
    fn baq_strand_matches_calling_strand() {
        use crate::calling::Strand;
        for flag in [99u16, 147, 83, 163] {
            let ob = matches!(classify_strand(flag), Some(Strand::Ob));
            assert_eq!(baq_is_ob(flag), ob, "flag {flag}");
        }
        // Supplementary copies follow their primary's strand; unpaired reads use orientation.
        assert!(!baq_is_ob(99 | 2048));
        assert!(baq_is_ob(83 | 2048));
        assert!(!baq_is_ob(0));
        assert!(baq_is_ob(16));
    }

    #[test]
    fn baq_mask_by_strand_and_mode() {
        assert_eq!(baq_mask(BaqMode::Taps, 99), b'C');
        assert_eq!(baq_mask(BaqMode::TapsFull, 147), b'C');
        assert_eq!(baq_mask(BaqMode::TapsPartial, 83), b'G');
        assert_eq!(baq_mask(BaqMode::Taps, 163), b'G');
        assert_eq!(baq_mask(BaqMode::Htslib, 99), 0);
    }

    #[test]
    fn baq_need_covers_band_and_clips_to_contig() {
        let r = record("2S96M1D2M", &[b'A'; 100], 1_000);
        let (lo, hi) = baq_need(&r, 1_000_000);
        // span 99, slack 2*100 + 99 + 16 = 315
        assert_eq!((lo, hi), (685, 1_000 + 99 + 315));
        assert_eq!(baq_need(&r, 1_200), (685, 1_200));
        assert_eq!(baq_need(&record("4M", b"ACGT", 3), 100).0, 0);
    }

    fn record(cigar: &str, seq: &[u8], pos: i64) -> Record {
        let mut r = Record::new();
        let c = rust_htslib::bam::record::CigarString::try_from(cigar).unwrap();
        r.set(b"q", Some(&c), seq, &vec![30u8; seq.len()]);
        r.set_pos(pos);
        r
    }

    #[test]
    fn taps_baq_skips_only_perfect_all_match_reads() {
        // OT read: reference C is masked.
        let seq = b"AACGTTCGAA";
        let full = RefSlice { seq, lo: 0, hi: 10 };
        let m = b'C';
        let perfect = record("4M", b"ATGT", 1); // ref ANGT: T vs masked C is fine
        let mismatch = record("4M", b"AAGA", 1);
        let read_n = record("4M", b"ATNT", 1); // N against reference G
        let clipped = record("1S3M", b"ATGT", 2);
        for (rec, taps, partial) in [(&perfect, false, false), (&mismatch, true, false), (&read_n, true, false), (&clipped, true, true)] {
            assert_eq!(BaqMode::Taps.wants(rec, &full, m), taps, "{:?}", rec.seq().as_bytes());
            assert_eq!(BaqMode::TapsPartial.wants(rec, &full, m), partial, "{:?}", rec.seq().as_bytes());
            assert!(BaqMode::TapsFull.wants(rec, &full, m));
        }
        // Without the mask (htslib) the T against C is a mismatch.
        assert!(BaqMode::Taps.wants(&perfect, &full, 0));
        // A slice that starts later: same answer, indexed by contig position.
        let part = RefSlice { seq: &seq[1..], lo: 1, hi: 10 };
        assert!(!BaqMode::Taps.wants(&perfect, &part, m));
        // Running off the reference end: run BAQ.
        assert!(BaqMode::Taps.wants(&record("4M", b"AAAA", 8), &full, m));
    }
}
