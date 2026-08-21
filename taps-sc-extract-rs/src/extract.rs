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
//!   (flag 3). The reference must start at contig coordinate 0.
//! - `ignore_overlaps`: `bam_plp_init_overlaps` (zeros the lower-quality
//!   overlapping mate base; pysam then drops it via `min_base_quality`)
//!
//! Python's current extractor sets `compute_baq=True` but never passes
//! `fastafile=` to `AlignmentFile.pileup()`, so BAQ is a no-op there. Rust
//! actually applies BAQ when `compute_baq` is true because we always have FASTA.

use crate::accumulate::{ensure_cell, merge_window_cells, BarcodeIntern, CellMaps};
use crate::barcode::barcode_from_qname_bytes;
use crate::calling::{call_mctot, classify_strand};
use crate::context::{classify_trinucleotide, Context, TriContext};
use crate::fasta::FastFaiReader;
use crate::window::Window;
use anyhow::{Context as AnyhowContext, Result};
use rust_htslib::bam::pileup::Alignment;
use rust_htslib::bam::{IndexedReader, Read, Record};
use rust_htslib::htslib;
use std::collections::HashSet;
use std::os::raw::c_char;
use std::path::Path;

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
}

impl WindowResult {
    pub fn empty() -> Self {
        Self {
            stats: CallStats::default(),
            intern: BarcodeIntern::default(),
            cells: Vec::new(),
        }
    }

    pub fn n_cells(&self) -> usize {
        self.intern.len()
    }

    /// Merge counts and interned barcode names (not per-pos maps).
    pub fn absorb_counts(&mut self, other: &WindowResult) {
        self.stats.merge(&other.stats);
        for name in other.intern.names() {
            self.intern.intern_str(name);
        }
    }

    /// Merge counts and per-pos cell maps (barcode-string intern remap).
    pub fn absorb_maps(&mut self, other: &WindowResult) {
        self.stats.merge(&other.stats);
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
    pub decomp_threads: usize,
    /// Fill per-barcode (pos, t, c) maps. Off for stats-only (saves RAM).
    pub accumulate: bool,
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
            decomp_threads: 1,
            accumulate: false,
        }
    }
}

/// pysam / samtools default: skip unmapped, secondary, QC fail, duplicate.
const FLAG_FILTER: u16 = (htslib::BAM_FUNMAP
    | htslib::BAM_FSECONDARY
    | htslib::BAM_FQCFAIL
    | htslib::BAM_FDUP) as u16;

/// BAQ_APPLY | BAQ_EXTEND (pysam `compute_baq=True`, `redo_baq=False`).
const BAQ_FLAGS: i32 = 3;

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

fn accept_record(rec: &mut Record, params: &ExtractParams, baq_ptr: *const c_char, baq_len: i64) -> bool {
    let flags = rec.flags();
    if flags & FLAG_FILTER != 0 {
        return false;
    }
    if params.compute_baq && !baq_ptr.is_null() {
        let _ = unsafe { htslib::sam_prob_realn(rec.inner_mut(), baq_ptr, baq_len, BAQ_FLAGS) };
    }
    if rec.mapq() < params.min_mapq {
        return false;
    }
    if params.ignore_orphans
        && (flags & htslib::BAM_FPAIRED as u16) != 0
        && (flags & htslib::BAM_FPROPER_PAIR as u16) == 0
    {
        return false;
    }
    true
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
    baq_ref: Option<&[u8]>,
) -> Result<WindowResult> {
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

    let (baq_ptr, baq_len) = match baq_ref {
        Some(s) if params.compute_baq => (s.as_ptr() as *const c_char, s.len() as i64),
        _ => (std::ptr::null(), 0),
    };

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

    let mut mru_ptr: *const u8 = std::ptr::null();
    let mut mru_len: usize = 0;
    let mut mru_id: Option<u32> = None;

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

            let qname = rec.qname();
            let qptr = qname.as_ptr();
            let qlen = qname.len();
            let id_opt = if qptr == mru_ptr && qlen == mru_len {
                mru_id
            } else {
                let bc_bytes = barcode_from_qname_bytes(qname);
                let valid = if let Some(wl) = whitelist {
                    let bc = std::str::from_utf8(bc_bytes).unwrap_or("");
                    wl.contains(bc)
                } else {
                    true
                };
                let id = if valid {
                    Some(intern.intern_bytes(bc_bytes))
                } else {
                    None
                };
                mru_ptr = qptr;
                mru_len = qlen;
                mru_id = id;
                id
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

    let mut rec = Record::new();
    loop {
        match bam.read(&mut rec) {
            None => break,
            Some(Err(e)) => anyhow::bail!("BAM read: {e}"),
            Some(Ok(())) => {
                if !accept_record(&mut rec, params, baq_ptr, baq_len) {
                    continue;
                }
                if unsafe { htslib::bam_plp_push(plp, rec.inner()) } < 0 {
                    anyhow::bail!("bam_plp_push failed");
                }
                drain_ready(plp, &mut consume_column)?;
            }
        }
    }
    if unsafe { htslib::bam_plp_push(plp, std::ptr::null()) } < 0 {
        anyhow::bail!("bam_plp_push(EOF) failed");
    }
    drain_ready(plp, &mut consume_column)?;

    Ok(WindowResult {
        stats,
        intern,
        cells,
    })
}

/// Keep one full-contig sequence for `sam_prob_realn` (must start at coordinate 0).
pub fn refresh_baq_cache<'a>(
    fasta: &FastFaiReader,
    cache: &'a mut Option<(String, Vec<u8>)>,
    contig: &str,
) -> Result<Option<&'a [u8]>> {
    let need = cache
        .as_ref()
        .map(|(c, _)| c.as_str() != contig)
        .unwrap_or(true);
    if need {
        let len = fasta
            .reference_length(contig)
            .ok_or_else(|| anyhow::anyhow!("contig {contig} missing in FASTA"))?;
        let seq = fasta
            .fetch(contig, 0, len)
            .with_context(|| format!("FASTA BAQ fetch {contig}"))?;
        *cache = Some((contig.to_string(), seq));
    }
    Ok(cache.as_ref().map(|(_, s)| s.as_slice()))
}

pub fn open_bam(path: &Path, decomp_threads: usize) -> Result<IndexedReader> {
    let mut bam = IndexedReader::from_path(path)
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
    let mut total = WindowResult {
        stats: CallStats::default(),
        intern: BarcodeIntern::default(),
        cells: Vec::new(),
    };
    // Full-contig cache for sam_prob_realn (pysam fetches 0..MAX_POS per tid).
    let mut baq_cache: Option<(String, Vec<u8>)> = None;
    for (i, w) in windows.iter().enumerate() {
        let baq_ref = if params.compute_baq {
            let need_fetch = baq_cache
                .as_ref()
                .map(|(c, _)| c.as_str() != w.contig)
                .unwrap_or(true);
            if need_fetch {
                let len = fasta.reference_length(&w.contig).ok_or_else(|| {
                    anyhow::anyhow!("contig {} missing in FASTA", w.contig)
                })?;
                let seq = fasta
                    .fetch(&w.contig, 0, len)
                    .with_context(|| format!("FASTA BAQ fetch {}", w.contig))?;
                baq_cache = Some((w.contig.clone(), seq));
            }
            baq_cache.as_ref().map(|(_, s)| s.as_slice())
        } else {
            None
        };
        let r = process_window(&mut bam, &fasta, w, whitelist, params, baq_ref)?;
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
    let mut total = WindowResult {
        stats: CallStats::default(),
        intern: BarcodeIntern::default(),
        cells: Vec::new(),
    };
    let mut baq_cache: Option<(String, Vec<u8>)> = None;
    let mut n_files = 0usize;
    for (i, w) in windows.iter().enumerate() {
        let baq_ref = if params.compute_baq {
            let need_fetch = baq_cache
                .as_ref()
                .map(|(c, _)| c.as_str() != w.contig)
                .unwrap_or(true);
            if need_fetch {
                let len = fasta.reference_length(&w.contig).ok_or_else(|| {
                    anyhow::anyhow!("contig {} missing in FASTA", w.contig)
                })?;
                let seq = fasta
                    .fetch(&w.contig, 0, len)
                    .with_context(|| format!("FASTA BAQ fetch {}", w.contig))?;
                baq_cache = Some((w.contig.clone(), seq));
            }
            baq_cache.as_ref().map(|(_, s)| s.as_slice())
        } else {
            None
        };
        let mut p = params.clone();
        p.accumulate = true;
        let r = process_window(&mut bam, &fasta, w, whitelist, &p, baq_ref)?;
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
