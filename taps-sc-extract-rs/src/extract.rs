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
use rustc_hash::FxHashMap;
use std::cell::Cell;
use std::collections::HashSet;
use std::os::raw::c_char;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
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
    pub decomp_threads: usize,
    /// Fill per-barcode (pos, t, c) maps. Off for stats-only (saves RAM).
    pub accumulate: bool,
    /// Drop this many 5' sequenced bases of read 1 (default 0).
    pub trim_r1: u32,
    /// Drop this many 5' sequenced bases of read 2 (default 10).
    pub trim_r2: u32,
    /// Time each stage of every window (`WindowResult.timing`).
    pub profile: bool,
}

impl Default for ExtractParams {
    fn default() -> Self {
        Self {
            min_baseq: 20,
            min_mapq: 0,
            max_depth: 250,
            ignore_overlaps: true,
            ignore_orphans: true,
            compute_baq: true,
            baq_mode: BaqMode::Taps,
            baq_width: baq_simd::Width::cpu_best(),
            decomp_threads: 1,
            accumulate: false,
            trim_r1: DEFAULT_TRIM_R1,
            trim_r2: DEFAULT_TRIM_R2,
            profile: false,
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

    /// Run the HMM on this read? (`seq` is the reference BAQ would use.)
    #[inline]
    fn wants(self, rec: &Record, seq: &[u8]) -> bool {
        match self {
            BaqMode::TapsFull | BaqMode::Htslib => true,
            BaqMode::TapsPartial => !cigar_all_match(rec),
            BaqMode::Taps => !(cigar_all_match(rec) && matches_masked_ref(rec, seq)),
        }
    }
}

/// CIGAR is only M operations (no clip, indel or skip).
#[inline]
fn cigar_all_match(rec: &Record) -> bool {
    rec.raw_cigar().iter().all(|&c| c & 0xf == htslib::BAM_CMATCH)
}

/// Every read base equals the reference base where the reference is A/C/G/T
/// (masked N positions match anything; a read N is a mismatch). Assumes an
/// all-M CIGAR, so read base i sits at `pos + i`.
#[inline]
fn matches_masked_ref(rec: &Record, refseq: &[u8]) -> bool {
    let n = rec.seq_len();
    let Ok(pos) = usize::try_from(rec.pos()) else {
        return false;
    };
    let Some(r) = refseq.get(pos..pos + n) else {
        return false;
    };
    let raw = rec.seq().encoded;
    for (i, &rb) in r.iter().enumerate() {
        // BAM 4-bit codes: A=1 C=2 G=4 T=8, high nibble first.
        let want = match rb {
            b'A' => 1,
            b'C' => 2,
            b'G' => 4,
            b'T' => 8,
            _ => continue,
        };
        if (raw[i / 2] >> (4 * (1 - (i & 1)))) & 0xf != want {
            return false;
        }
    }
    true
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

/// Records read, filtered and BAQ'd together before they are pushed in order.
/// Large enough that most reads find same-shape partners for the SIMD lanes
/// (chr19 at the default skip: ~92% of 4-lane slots filled at 8,192).
const READ_BATCH: usize = 8192;

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

/// Filters, and the reference to run BAQ against (`None`: no BAQ for this
/// read). Returns `None` when the read is not pushed at all.
fn accept_record<'r>(
    rec: &Record,
    params: &ExtractParams,
    baq_ref: Option<&'r BaqRef>,
    whitelist: Option<&HashSet<String>>,
) -> Option<Option<&'r [u8]>> {
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
    let Some(r) = baq_ref.filter(|_| params.compute_baq && !excluded) else {
        return Some(None);
    };
    let seq: &[u8] = if baq_is_ob(flags) { &r.ob } else { &r.ot };
    Some(params.baq_mode.wants(rec, seq).then_some(seq))
}

/// Process one genomic window `[start, end)` with column-wise `bam_plp` pileup.
///
/// `baq_ref` must be the full contig starting at position 0 when `compute_baq`
/// is true (same contract as pysam `faidx_fetch_seq(..., 0, MAX_POS)`).
pub fn process_window(
    bam: &mut IndexedReader,
    fasta: &FastFaiReader,
    window: &Window,
    whitelist: Option<&HashSet<String>>,
    params: &ExtractParams,
    baq_ref: Option<&BaqRef>,
) -> Result<WindowResult> {
    let profile = params.profile;
    let t_window = Instant::now();
    let mut timing = WindowTiming::default();
    let pad = 2u64;
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
        htslib::bam_plp_set_maxcnt(plp, params.max_depth as i32);
        if params.ignore_overlaps && bam_plp_init_overlaps(plp) < 0 {
            anyhow::bail!("bam_plp_init_overlaps failed");
        }
    }

    let mut stats = CallStats::default();
    let mut intern = BarcodeIntern::default();
    let mut cells: Vec<CellMaps> = Vec::new();
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
    let mut pool: Vec<Record> = (0..READ_BATCH).map(|_| Record::new()).collect();
    let mut keep = vec![false; READ_BATCH];
    let mut jobs: Vec<BaqJob> = Vec::with_capacity(READ_BATCH);
    let mut ws = baq_simd::Workspace::new();
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
    while !eof {
        lap(profile);
        jobs.clear();
        let mut n = 0;
        while n < READ_BATCH {
            match bam.read(&mut pool[n]) {
                None => {
                    eof = true;
                    break;
                }
                Some(Err(e)) => anyhow::bail!("BAM read: {e}"),
                Some(Ok(())) => {
                    timing.records += 1;
                    let acc = accept_record(&pool[n], params, baq_ref, whitelist);
                    keep[n] = acc.is_some();
                    if let Some(Some(seq)) = acc {
                        jobs.push(BaqJob {
                            b: pool[n].inner_mut(),
                            reference: seq.as_ptr() as *const c_char,
                            ref_len: seq.len() as i64,
                        });
                    }
                    n += 1;
                }
            }
        }
        timing.read_s += lap(profile);
        timing.baq_reads += jobs.len() as u64;
        ws.run(&jobs, params.baq_width);
        timing.baq_s += lap(profile);
        for (rec, _) in pool[..n].iter().zip(&keep[..n]).filter(|(_, k)| **k) {
            timing.pushed += 1;
            if unsafe { htslib::bam_plp_push(plp, rec.inner()) } < 0 {
                anyhow::bail!("bam_plp_push failed");
            }
            drain_ready(plp, &mut consume_column)?;
        }
        plp_s += lap(profile);
    }
    if profile {
        lap(true);
    }
    if unsafe { htslib::bam_plp_push(plp, std::ptr::null()) } < 0 {
        anyhow::bail!("bam_plp_push(EOF) failed");
    }
    drain_ready(plp, &mut consume_column)?;
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

/// Full-contig BAQ references for one contig (each starts at coordinate 0):
/// `ot` for top-strand reads, `ob` for bottom-strand reads. With
/// `BaqMode::Htslib` both are the same plain sequence.
pub struct BaqRef {
    pub ot: Arc<[u8]>,
    pub ob: Arc<[u8]>,
}

impl BaqRef {
    pub(crate) fn new(seq: Vec<u8>, mode: BaqMode) -> Self {
        match mode {
            BaqMode::Htslib => {
                let s: Arc<[u8]> = seq.into();
                Self { ot: Arc::clone(&s), ob: s }
            }
            BaqMode::Taps | BaqMode::TapsFull | BaqMode::TapsPartial => {
                // FastFaiReader::fetch returns uppercase.
                let mask = |from: u8| -> Arc<[u8]> {
                    seq.iter().map(|&b| if b == from { b'N' } else { b }).collect()
                };
                Self { ot: mask(b'C'), ob: mask(b'G') }
            }
        }
    }
}

/// Full-contig BAQ references (`BaqRef`), shared by all workers. A contig is
/// fetched once, by the first window that needs it, and dropped when its last
/// window releases it.
pub struct BaqRefs<'a> {
    fasta: &'a FastFaiReader,
    slots: FxHashMap<String, BaqSlot>,
}

struct BaqSlot {
    seq: Mutex<Option<Arc<BaqRef>>>,
    remaining: AtomicUsize,
}

impl<'a> BaqRefs<'a> {
    pub fn new(fasta: &'a FastFaiReader, windows: &[Window]) -> Self {
        let mut slots: FxHashMap<String, BaqSlot> = FxHashMap::default();
        for w in windows {
            slots
                .entry(w.contig.clone())
                .or_insert_with(|| BaqSlot {
                    seq: Mutex::new(None),
                    remaining: AtomicUsize::new(0),
                })
                .remaining
                .fetch_add(1, Ordering::Relaxed);
        }
        Self { fasta, slots }
    }

    /// Sequence for `contig`. Pair every call with `release`.
    pub fn acquire(&self, contig: &str, mode: BaqMode) -> Result<Arc<BaqRef>> {
        let slot = self
            .slots
            .get(contig)
            .ok_or_else(|| anyhow::anyhow!("contig {contig} not in window plan"))?;
        let mut guard = slot.seq.lock().expect("BAQ ref mutex");
        if let Some(seq) = guard.as_ref() {
            return Ok(Arc::clone(seq));
        }
        let len = self
            .fasta
            .reference_length(contig)
            .ok_or_else(|| anyhow::anyhow!("contig {contig} missing in FASTA"))?;
        let seq = Arc::new(BaqRef::new(
            self.fasta
                .fetch(contig, 0, len)
                .with_context(|| format!("FASTA BAQ fetch {contig}"))?,
            mode,
        ));
        *guard = Some(Arc::clone(&seq));
        Ok(seq)
    }

    /// One window on `contig` is done. The last one drops the shared copy.
    pub fn release(&self, contig: &str) {
        if let Some(slot) = self.slots.get(contig) {
            if slot.remaining.fetch_sub(1, Ordering::AcqRel) == 1 {
                *slot.seq.lock().expect("BAQ ref mutex") = None;
            }
        }
    }
}

/// `process_window` with the BAQ reference taken from (and returned to) `baq_refs`.
pub fn process_window_with_baq(
    bam: &mut IndexedReader,
    fasta: &FastFaiReader,
    baq_refs: &BaqRefs,
    window: &Window,
    whitelist: Option<&HashSet<String>>,
    params: &ExtractParams,
) -> Result<WindowResult> {
    let baq_ref = if params.compute_baq {
        Some(baq_refs.acquire(&window.contig, params.baq_mode)?)
    } else {
        None
    };
    let r = process_window(bam, fasta, window, whitelist, params, baq_ref.as_deref());
    drop(baq_ref);
    baq_refs.release(&window.contig);
    r
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
    // Full-contig reference for sam_prob_realn (pysam fetches 0..MAX_POS per tid).
    let baq_refs = BaqRefs::new(&fasta, windows);
    for (i, w) in windows.iter().enumerate() {
        let r = process_window_with_baq(&mut bam, &fasta, &baq_refs, w, whitelist, params)?;
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
    let baq_refs = BaqRefs::new(&fasta, windows);
    let mut n_files = 0usize;
    for (i, w) in windows.iter().enumerate() {
        let mut p = params.clone();
        p.accumulate = true;
        let r = process_window_with_baq(&mut bam, &fasta, &baq_refs, w, whitelist, &p)?;
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
    fn taps_baq_ref_masks_convertible_base_per_strand() {
        let r = BaqRef::new(b"ACGTNCCGG".to_vec(), BaqMode::Taps);
        assert_eq!(&*r.ot, b"ANGTNNNGG");
        assert_eq!(&*r.ob, b"ACNTNCCNN");
        let h = BaqRef::new(b"ACGT".to_vec(), BaqMode::Htslib);
        assert_eq!(&*h.ot, b"ACGT");
        assert!(Arc::ptr_eq(&h.ot, &h.ob));
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
        // Masked OT reference: C -> N.
        let masked = BaqRef::new(b"AACGTTCGAA".to_vec(), BaqMode::Taps).ot;
        let perfect = record("4M", b"ATGT", 1); // ref ANGT: T vs masked C is fine
        let mismatch = record("4M", b"AAGA", 1);
        let read_n = record("4M", b"ATNT", 1); // N against reference G
        let clipped = record("1S3M", b"ATGT", 2);
        for (rec, taps, partial) in [(&perfect, false, false), (&mismatch, true, false), (&read_n, true, false), (&clipped, true, true)] {
            assert_eq!(BaqMode::Taps.wants(rec, &masked), taps, "{:?}", rec.seq().as_bytes());
            assert_eq!(BaqMode::TapsPartial.wants(rec, &masked), partial, "{:?}", rec.seq().as_bytes());
            assert!(BaqMode::TapsFull.wants(rec, &masked));
        }
        // Running off the reference end: run BAQ.
        assert!(BaqMode::Taps.wants(&record("4M", b"AAAA", 8), &masked));
    }
}
