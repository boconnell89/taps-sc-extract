//! Alignment input: BAM or CRAM with an index, built on demand.
//!
//! `prepare` runs once on the main thread before any reader is opened:
//! - plain SAM (`.sam`, `.sam.gz`) is converted to a temporary BAM (it must be
//!   coordinate-sorted; windowed extraction needs an index);
//! - a missing `.bai`/`.csi`/`.crai` is built next to the file, or in the temp
//!   directory when that location is not writable;
//! - CRAM decodes against the `-f` FASTA.
//!
//! Every reader then opens through `open_indexed`, which applies the chosen
//! index path and CRAM reference.

use anyhow::{Context, Result};
use rust_htslib::bam::{self, index, IndexedReader, Read};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::time::Instant;

#[derive(Clone, Debug, Default)]
struct OpenOptions {
    /// The prepared file these options belong to.
    path: PathBuf,
    /// Index file to use instead of htslib's default lookup next to the file.
    index: Option<PathBuf>,
    /// FASTA for CRAM decoding.
    reference: Option<PathBuf>,
}

/// Open options per prepared file (one entry per input; several with multi-BAM runs).
static OPEN: Mutex<Vec<OpenOptions>> = Mutex::new(Vec::new());

fn options(path: &Path) -> Option<OpenOptions> {
    OPEN.lock().ok()?.iter().find(|o| o.path == path).cloned()
}

/// Index file chosen by `prepare` for `path` (None: htslib's default lookup).
pub fn index_path(path: &Path) -> Option<PathBuf> {
    options(path).and_then(|o| o.index)
}

/// True for CRAM input (its index is tied to the open file, so it cannot be shared).
pub fn is_cram_path(path: &Path) -> bool {
    is_cram(path)
}

/// The alignment file to read, plus anything that must live as long as it.
pub struct PreparedInput {
    pub path: PathBuf,
    /// Seconds spent converting SAM and building an index (0 if neither).
    pub prep_s: f64,
    _temp: Option<tempfile::TempDir>,
}

fn lower_name(p: &Path) -> String {
    p.file_name().map(|n| n.to_string_lossy().to_ascii_lowercase()).unwrap_or_default()
}

fn is_cram(p: &Path) -> bool {
    lower_name(p).ends_with(".cram")
}

fn is_sam(p: &Path) -> bool {
    let n = lower_name(p);
    n.ends_with(".sam") || n.ends_with(".sam.gz")
}

/// Open an indexed BAM/CRAM with the index and reference chosen by `prepare`.
pub fn open_indexed(path: &Path) -> Result<IndexedReader> {
    let opts = options(path);
    let opts = opts.as_ref();
    let mut r = match opts.and_then(|o| o.index.as_ref()) {
        Some(idx) => IndexedReader::from_path_and_index(path, idx),
        None => IndexedReader::from_path(path),
    }
    .map_err(|e| anyhow::anyhow!("open {}: {e}", path.display()))?;
    if let Some(fa) = opts.and_then(|o| o.reference.as_ref()) {
        r.set_reference(fa)
            .map_err(|e| anyhow::anyhow!("CRAM reference {}: {e}", fa.display()))?;
    }
    Ok(r)
}

/// BAI unless a contig is too long for it (2^29 - 1 bp), then CSI.
fn index_type(path: &Path) -> Result<index::Type> {
    let r = bam::Reader::from_path(path).map_err(|e| anyhow::anyhow!("open {}: {e}", path.display()))?;
    let h = r.header();
    let longest = (0..h.target_count()).filter_map(|t| h.target_len(t)).max().unwrap_or(0);
    Ok(if longest >= (1u64 << 29) { index::Type::Csi(14) } else { index::Type::Bai })
}

/// Copy a coordinate-sorted SAM into `out` as BAM.
fn sam_to_bam(sam: &Path, out: &Path, threads: usize) -> Result<()> {
    let mut r = bam::Reader::from_path(sam).map_err(|e| anyhow::anyhow!("open SAM {}: {e}", sam.display()))?;
    let header = bam::Header::from_template(r.header());
    let mut w = bam::Writer::from_path(out, &header, bam::Format::Bam)
        .map_err(|e| anyhow::anyhow!("create {}: {e}", out.display()))?;
    if threads > 1 {
        let _ = w.set_threads(threads);
    }
    let mut rec = bam::Record::new();
    let mut last = (-1i32, -1i64);
    while let Some(res) = r.read(&mut rec) {
        res.map_err(|e| anyhow::anyhow!("read SAM {}: {e}", sam.display()))?;
        let key = (if rec.tid() < 0 { i32::MAX } else { rec.tid() }, rec.pos());
        if key < last {
            anyhow::bail!(
                "{} is not coordinate-sorted. Sort it first: samtools sort -o sorted.bam {}",
                sam.display(),
                sam.display()
            );
        }
        last = key;
        w.write(&rec).map_err(|e| anyhow::anyhow!("write {}: {e}", out.display()))?;
    }
    Ok(())
}

fn temp_dir(temp: &mut Option<tempfile::TempDir>) -> Result<PathBuf> {
    if temp.is_none() {
        *temp = Some(
            tempfile::Builder::new()
                .prefix("taps-rs-input-")
                .tempdir()
                .context("temp dir for alignment input")?,
        );
    }
    Ok(temp.as_ref().expect("set above").path().to_path_buf())
}

/// Make `path` readable through `open_indexed`; see the module docs.
pub fn prepare(path: &Path, fasta: &Path, threads: usize) -> Result<PreparedInput> {
    let t0 = Instant::now();
    let mut temp: Option<tempfile::TempDir> = None;
    let mut path = path.to_path_buf();
    if !path.exists() {
        anyhow::bail!("alignment file {} does not exist", path.display());
    }
    if is_sam(&path) {
        let dir = temp_dir(&mut temp)?;
        let out = dir.join("input.bam");
        eprintln!("Converting SAM {} to a temporary BAM in {} ...", path.display(), dir.display());
        sam_to_bam(&path, &out, threads)?;
        path = out;
    }
    let reference = is_cram(&path).then(|| fasta.to_path_buf());
    let mut opts = OpenOptions { path: path.clone(), index: None, reference };

    // Probe quietly: htslib logs an error for a missing index, which is handled below.
    let level = unsafe { rust_htslib::htslib::hts_get_log_level() };
    unsafe { rust_htslib::htslib::hts_set_log_level(rust_htslib::htslib::htsLogLevel_HTS_LOG_OFF) };
    let has_index = match &opts.reference {
        Some(fa) => IndexedReader::from_path(&path).and_then(|mut r| r.set_reference(fa)).is_ok(),
        None => IndexedReader::from_path(&path).is_ok(),
    };
    unsafe { rust_htslib::htslib::hts_set_log_level(level) };
    if !has_index {
        let csi = matches!(index_type(&path)?, index::Type::Csi(_));
        let ty = || if csi { index::Type::Csi(14) } else { index::Type::Bai };
        let ext = match (is_cram(&path), csi) {
            (true, _) => "crai",
            (false, false) => "bai",
            (false, true) => "csi",
        };
        let t = Instant::now();
        eprintln!("No index for {}; building .{ext} ...", path.display());
        let n = threads.clamp(1, 64) as u32;
        if index::build(&path, None, ty(), n).is_err() {
            let dir = temp_dir(&mut temp)?;
            let name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
            let idx = dir.join(format!("{name}.{ext}"));
            index::build(&path, Some(&idx), ty(), n).map_err(|e| {
                anyhow::anyhow!("build index for {} (beside it and in {}): {e}", path.display(), dir.display())
            })?;
            eprintln!("  (directory not writable; index written to {})", idx.display());
            opts.index = Some(idx);
        }
        eprintln!("Index built in {:.1}s", t.elapsed().as_secs_f64());
    }
    if let Ok(mut v) = OPEN.lock() {
        v.retain(|o| o.path != path);
        v.push(opts);
    }
    // Fail here, with the real error, if the file still cannot be opened.
    open_indexed(&path)?;
    Ok(PreparedInput {
        prep_s: if has_index && temp.is_none() { 0.0 } else { t0.elapsed().as_secs_f64() },
        path,
        _temp: temp,
    })
}
