//! Amethyst HDF5 writer: `/CG|CH/<barcode>/1` compound datasets and `master.h5`
//! relative ExternalLinks.
//!
//! hdf5-metno serializes every HDF5 call behind one process-wide lock, so
//! threads in one process cannot write shards in parallel. With several shards
//! the default is one `write-shard` child process per shard (each has its own
//! lock); `WriterMode::Thread` keeps the in-process pool. Chunk size is
//! `min(len, 65536)`.

use crate::shard_io::{
    list_chunk_files, read_chunk_file, shard_dir, write_chunk_to, ChunkFile,
};
use anyhow::{Context, Result};
use hdf5_metno::types::FixedAscii;
use hdf5_metno::{File, H5Type};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::path::{Path, PathBuf};

pub const AMETHYST_VERSION: &str = "amethyst2.0.0";

#[derive(H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct MethRec {
    pub chr: FixedAscii<10>,
    pub pos: i64,
    pub pct: f64,
    pub t: i64,
    pub c: i64,
}

impl MethRec {
    pub fn new(contig: &str, pos: u32, t: u32, c: u32) -> Self {
        Self::with_chr(contig_s10(contig), pos, t, c)
    }

    /// Same as `new` with the contig name already converted (once per window).
    #[inline]
    pub fn with_chr(chr: FixedAscii<10>, pos: u32, t: u32, c: u32) -> Self {
        let tot = t + c;
        let pct = if tot == 0 {
            0.0
        } else {
            100.0 * f64::from(c) / f64::from(tot)
        };
        Self {
            chr,
            pos: i64::from(pos),
            pct,
            t: i64::from(t),
            c: i64::from(c),
        }
    }
}

fn contig_s10(s: &str) -> FixedAscii<10> {
    let mut buf = [0u8; 10];
    let b = s.as_bytes();
    let n = b.len().min(10);
    buf[..n].copy_from_slice(&b[..n]);
    FixedAscii::<10>::from_ascii(&buf).unwrap_or_else(|_| {
        FixedAscii::<10>::from_ascii(b"chrN\0\0\0\0\0\0").expect("literal")
    })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum H5Compression {
    Gzip(u8),
    GzipShuffle(u8),
    Lzf,
    BloscLz4(u8, bool),
    BloscZstd(u8, bool),
    None,
}

impl H5Compression {
    pub fn parse(s: &str) -> Result<Self> {
        match s.to_ascii_lowercase().as_str() {
            "gzip" | "gzip-fast" | "gzip1" => Ok(Self::Gzip(1)),
            "gzip-shuffle" | "shuffle" => Ok(Self::GzipShuffle(1)),
            "gzip6" => Ok(Self::Gzip(6)),
            "lzf" => Ok(Self::Lzf),
            "blosc" | "blosc-lz4" => Ok(Self::BloscLz4(9, true)),
            "blosc-zstd" => Ok(Self::BloscZstd(5, true)),
            "none" => Ok(Self::None),
            other => anyhow::bail!("unsupported --compression {other} (options: gzip, gzip-shuffle, gzip6, lzf, blosc, blosc-zstd, none)"),
        }
    }
}

fn write_version(file: &File) -> Result<()> {
    let meta = file.create_group("metadata")?;
    let v = AMETHYST_VERSION.as_bytes();
    meta.new_dataset::<u8>()
        .shape(v.len())
        .create("version")?
        .write(v)
        .map_err(|e| anyhow::anyhow!("write metadata/version: {e}"))?;
    Ok(())
}

fn write_cell_dataset(
    file: &File,
    context: &str,
    barcode: &str,
    recs: &[MethRec],
    compression: H5Compression,
) -> Result<()> {
    if recs.is_empty() {
        return Ok(());
    }
    create_cell_dataset(file, context, barcode, recs.len(), compression)?
        .write(recs)
        .map_err(|e| anyhow::anyhow!("write {context}/{barcode}/1: {e}"))?;
    Ok(())
}

/// HDF5 chunk length used for every cell dataset.
const H5_CHUNK: usize = 65536;
/// Records converted and written per call (a multiple of `H5_CHUNK`, so every
/// slice starts on a chunk boundary). Bounds the `MethRec` buffer and HDF5's
/// type-conversion and filter buffers to about 12 MB, whatever the cell size.
const WRITE_SLICE: usize = 4 * H5_CHUNK;

/// Where one shard's write time went (seconds, summed over its cells).
#[derive(Clone, Debug, Default)]
pub struct ShardWriteStats {
    pub cells: usize,
    /// Group and dataset creation (HDF5 metadata).
    pub create_s: f64,
    /// Expanding 12-byte sites into 48-byte `MethRec`s.
    pub fill_s: f64,
    /// `write_slice`: filter (compression) plus file I/O, under the HDF5 lock.
    pub h5_write_s: f64,
}

/// Write one cell's sites from window chunks (`(contig, sites)` in chunk_id
/// order) in `WRITE_SLICE`-record pieces. Same dataset as `write_cell_dataset`
/// over the concatenated records.
fn write_cell_dataset_sliced(
    file: &File,
    context: &str,
    barcode: &str,
    parts: &[(&str, &[(u32, u32, u32)])],
    compression: H5Compression,
    st: &mut ShardWriteStats,
) -> Result<()> {
    let n: usize = parts.iter().map(|(_, s)| s.len()).sum();
    if n == 0 {
        return Ok(());
    }
    let t_create = std::time::Instant::now();
    let ds = create_cell_dataset(file, context, barcode, n, compression)?;
    st.create_s += t_create.elapsed().as_secs_f64();
    let t_cell = std::time::Instant::now();
    let mut write_s = 0.0f64;
    let mut buf: Vec<MethRec> = Vec::with_capacity(n.min(WRITE_SLICE));
    let mut start = 0usize;
    let mut flush = |buf: &mut Vec<MethRec>, start: &mut usize| -> Result<()> {
        let t = std::time::Instant::now();
        let end = *start + buf.len();
        ds.write_slice(buf.as_slice(), *start..end)
            .map_err(|e| anyhow::anyhow!("write {context}/{barcode}/1 [{}..{end}]: {e}", *start))?;
        *start = end;
        buf.clear();
        write_s += t.elapsed().as_secs_f64();
        Ok(())
    };
    for &(contig, sites) in parts {
        let chr = contig_s10(contig);
        for &(pos, t, c) in sites {
            buf.push(MethRec::with_chr(chr, pos, t, c));
            if buf.len() == WRITE_SLICE {
                flush(&mut buf, &mut start)?;
            }
        }
    }
    if !buf.is_empty() {
        flush(&mut buf, &mut start)?;
    }
    st.h5_write_s += write_s;
    st.fill_s += (t_cell.elapsed().as_secs_f64() - write_s).max(0.0);
    Ok(())
}

fn create_cell_dataset(
    file: &File,
    context: &str,
    barcode: &str,
    n: usize,
    compression: H5Compression,
) -> Result<hdf5_metno::Dataset> {
    let g = if file.link_exists(context) {
        file.group(context)?
    } else {
        file.create_group(context)?
    };
    let bc = g.create_group(barcode)?;
    let chunk = n.min(H5_CHUNK);
    let mut b = bc.new_dataset::<MethRec>().shape(n).chunk(chunk);
    b = match compression {
        H5Compression::Gzip(level) => b.deflate(level),
        H5Compression::GzipShuffle(level) => b.shuffle().deflate(level),
        H5Compression::Lzf => b.lzf(),
        H5Compression::BloscLz4(level, shuffle) => b.blosc(hdf5_metno::filters::Blosc::LZ4, level, shuffle),
        H5Compression::BloscZstd(level, shuffle) => b.blosc(hdf5_metno::filters::Blosc::ZStd, level, shuffle),
        H5Compression::None => b,
    };
    Ok(b.create("1")?)
}

/// Concatenate window chunks (already in chunk_id order) into Amethyst records.
pub fn records_from_chunks(chunks: &[ChunkFile]) -> FxHashMap<String, (Vec<MethRec>, Vec<MethRec>)> {
    let mut out: FxHashMap<String, (Vec<MethRec>, Vec<MethRec>)> = FxHashMap::default();
    for ch in chunks {
        for cell in &ch.cells {
            let e = out.entry(cell.barcode.clone()).or_default();
            e.0.extend(
                cell.cg
                    .iter()
                    .map(|&(pos, t, c)| MethRec::new(&ch.contig, pos, t, c)),
            );
            e.1.extend(
                cell.ch
                    .iter()
                    .map(|&(pos, t, c)| MethRec::new(&ch.contig, pos, t, c)),
            );
        }
    }
    out
}

pub fn write_shard_file(
    path: &Path,
    cells: &FxHashMap<String, (Vec<MethRec>, Vec<MethRec>)>,
    compression: H5Compression,
) -> Result<usize> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).ok();
    }
    let file = File::create(path).map_err(|e| anyhow::anyhow!("create {}: {e}", path.display()))?;
    write_version(&file)?;
    let mut bcs: Vec<&String> = cells.keys().collect();
    bcs.sort();
    for bc in bcs {
        let (cg, ch) = &cells[bc];
        write_cell_dataset(&file, "CG", bc, cg, compression)?;
        write_cell_dataset(&file, "CH", bc, ch, compression)?;
    }
    file.flush().ok();
    Ok(cells.len())
}

/// Blosc's internal compression threads per writer (`--compression-threads`;
/// 0 = Blosc's default of 1). Only the blosc codecs use it; LZF and gzip run
/// single-threaded inside HDF5. Set by `extract` before writing; passed on to
/// `write-shard` children.
pub static BLOSC_THREADS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// Set Blosc's thread count for this process (no-op for 0).
pub fn apply_blosc_threads(n: usize) {
    if n > 0 {
        hdf5_metno::filters::blosc_set_nthreads(n.min(255) as u8);
    }
}

impl H5Compression {
    pub fn is_blosc(self) -> bool {
        matches!(self, Self::BloscLz4(..) | Self::BloscZstd(..))
    }
}

/// Write one shard from window chunks, one barcode at a time.
///
/// Sites stay at 12 bytes each in the chunks; at most `WRITE_SLICE` records
/// are expanded to 48-byte `MethRec`s at once. Output matches
/// `write_shard_file(records_from_chunks(..))`: barcodes sorted, CG then CH,
/// sites in chunk_id order.
pub fn write_shard_from_chunks(
    path: &Path,
    mut chunks: Vec<ChunkFile>,
    compression: H5Compression,
) -> Result<ShardWriteStats> {
    chunks.sort_by_key(|c| c.chunk_id);
    let mut by_bc: FxHashMap<&str, Vec<(u32, u32)>> = FxHashMap::default();
    for (ci, ch) in chunks.iter().enumerate() {
        for (k, cell) in ch.cells.iter().enumerate() {
            by_bc
                .entry(cell.barcode.as_str())
                .or_default()
                .push((ci as u32, k as u32));
        }
    }
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).ok();
    }
    let file = File::create(path).map_err(|e| anyhow::anyhow!("create {}: {e}", path.display()))?;
    write_version(&file)?;
    let mut st = ShardWriteStats::default();
    let mut bcs: Vec<&str> = by_bc.keys().copied().collect();
    bcs.sort_unstable();
    for bc in &bcs {
        let refs = &by_bc[bc];
        let parts = |ctx_cg: bool| -> Vec<(&str, &[(u32, u32, u32)])> {
            refs.iter()
                .map(|&(ci, k)| {
                    let ch = &chunks[ci as usize];
                    let cell = &ch.cells[k as usize];
                    let sites = if ctx_cg { &cell.cg } else { &cell.ch };
                    (ch.contig.as_str(), sites.as_slice())
                })
                .collect()
        };
        write_cell_dataset_sliced(&file, "CG", bc, &parts(true), compression, &mut st)?;
        write_cell_dataset_sliced(&file, "CH", bc, &parts(false), compression, &mut st)?;
    }
    file.flush().ok();
    st.cells = bcs.len();
    Ok(st)
}

/// How shard files are written when there is more than one shard.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WriterMode {
    /// One `write-shard` child process per shard; each has its own HDF5 lock.
    Process,
    /// In-process rayon pool; all writers share hdf5-metno's global lock.
    Thread,
}

impl WriterMode {
    pub fn parse(s: &str, n_shards: usize) -> Result<Self> {
        match s.to_ascii_lowercase().as_str() {
            "auto" => Ok(if n_shards > 1 { Self::Process } else { Self::Thread }),
            "process" => Ok(Self::Process),
            "thread" => Ok(Self::Thread),
            o => anyhow::bail!("unknown --writer-mode {o} (use auto, process, or thread)"),
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Self::Process => "process",
            Self::Thread => "thread",
        }
    }
}

/// Where a `write-shard` child gets its chunks.
enum ShardInput {
    /// Stream-mode temp bucket directories `temp_dir/shard_XXX/` for this shard.
    Dirs(Vec<PathBuf>),
    /// Memory-mode payloads, sent on the child's stdin.
    Chunks(Vec<ChunkFile>),
}

/// What one `write-shard` child reported on its summary line.
#[derive(Clone, Debug, Default)]
pub struct WriterReport {
    pub shard: usize,
    /// Written by a `write-shard` child process (else in the extract process).
    pub child: bool,
    pub cells: usize,
    pub sites: u64,
    /// Reading the shard's chunks (temp files or stdin).
    pub read_s: f64,
    /// Writing the shard file: create + fill + h5_write.
    pub write_s: f64,
    pub create_s: f64,
    pub fill_s: f64,
    pub h5_write_s: f64,
    /// Size of the finished shard file.
    pub bytes: u64,
    pub peak_rss_mb: f64,
}

impl WriterReport {
    /// Report for a shard written in this process.
    fn local(shard: usize, path: &Path, sites: u64, read_s: f64, write_s: f64, st: &ShardWriteStats) -> Self {
        Self {
            shard,
            child: false,
            cells: st.cells,
            sites,
            read_s,
            write_s,
            create_s: st.create_s,
            fill_s: st.fill_s,
            h5_write_s: st.h5_write_s,
            bytes: std::fs::metadata(path).map(|m| m.len()).unwrap_or(0),
            peak_rss_mb: 0.0,
        }
    }

    /// The `write-shard` child's summary line (parsed back by the parent).
    pub fn child_line(&self, path: &Path) -> String {
        format!(
            "HDF5 shard {:03} cells={} sites={} read={:.3}s write={:.3}s create={:.3}s fill={:.3}s h5write={:.3}s bytes={} peak_rss={:.0}MB -> {}",
            self.shard,
            self.cells,
            self.sites,
            self.read_s,
            self.write_s,
            self.create_s,
            self.fill_s,
            self.h5_write_s,
            self.bytes,
            self.peak_rss_mb,
            path.display()
        )
    }
}

fn chunk_sites(chunks: &[ChunkFile]) -> u64 {
    chunks
        .iter()
        .flat_map(|c| c.cells.iter())
        .map(|cell| (cell.cg.len() + cell.ch.len()) as u64)
        .sum()
}

/// Read and write one shard in this process, timed.
pub fn write_shard_timed(shard: usize, path: &Path, input: impl FnOnce() -> Result<Vec<ChunkFile>>, compression: H5Compression) -> Result<WriterReport> {
    let t0 = std::time::Instant::now();
    let chunks = input()?;
    let read_s = t0.elapsed().as_secs_f64();
    let sites = chunk_sites(&chunks);
    let t1 = std::time::Instant::now();
    let st = write_shard_from_chunks(path, chunks, compression)?;
    Ok(WriterReport::local(shard, path, sites, read_s, t1.elapsed().as_secs_f64(), &st))
}

/// One line on the writers: per-shard spread (the slowest shard sets the
/// write wall time), where write time went, `master.h5` time.
pub fn writer_summary_line(reports: &[WriterReport], master_s: f64) -> Option<String> {
    if reports.is_empty() {
        return None;
    }
    let med = |mut v: Vec<f64>| -> f64 {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        v[v.len() / 2]
    };
    let slow = reports
        .iter()
        .max_by(|a, b| (a.read_s + a.write_s).partial_cmp(&(b.read_s + b.write_s)).unwrap_or(std::cmp::Ordering::Equal))?;
    let sites: Vec<f64> = reports.iter().map(|r| r.sites as f64).collect();
    let max_sites = sites.iter().cloned().fold(0.0, f64::max);
    let sum = |f: fn(&WriterReport) -> f64| reports.iter().map(f).sum::<f64>();
    let bytes: u64 = reports.iter().map(|r| r.bytes).sum();
    let busy = sum(|r| r.read_s + r.write_s);
    Some(format!(
        "Shard writers:         {} shards | sites/shard median {:.0}, max {:.0} ({:.2}x) | shard time median {:.2}s, max {:.2}s (shard {:03}, {} cells) | \
         time in read {:.0}%, create {:.0}%, fill {:.0}%, h5write {:.0}% | {:.1} MB/s per shard | master.h5 {:.2}s",
        reports.len(),
        med(sites.clone()),
        max_sites,
        if med(sites.clone()) > 0.0 { max_sites / med(sites) } else { 0.0 },
        med(reports.iter().map(|r| r.read_s + r.write_s).collect()),
        slow.read_s + slow.write_s,
        slow.shard,
        slow.cells,
        100.0 * sum(|r| r.read_s) / busy.max(1e-9),
        100.0 * sum(|r| r.create_s) / busy.max(1e-9),
        100.0 * sum(|r| r.fill_s) / busy.max(1e-9),
        100.0 * sum(|r| r.h5_write_s) / busy.max(1e-9),
        bytes as f64 / 1e6 / sum(|r| r.write_s).max(1e-9),
        master_s
    ))
}

/// Peak RSS of `n_writers` concurrent children, bounded above by the sum of the
/// `n_writers` largest child peaks.
pub fn concurrent_writer_peak_mb(reports: &[WriterReport], n_writers: usize) -> f64 {
    let mut peaks: Vec<f64> = reports.iter().filter(|r| r.child).map(|r| r.peak_rss_mb).collect();
    peaks.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    peaks.iter().take(n_writers.max(1)).sum()
}

fn parse_writer_line(shard: usize, stderr: &str) -> WriterReport {
    let mut rep = WriterReport {
        shard,
        child: true,
        ..Default::default()
    };
    let secs = |v: &str| v.strip_suffix('s').and_then(|x| x.parse().ok()).unwrap_or(0.0);
    for tok in stderr.split_whitespace() {
        if let Some((k, v)) = tok.split_once('=') {
            match k {
                "cells" => rep.cells = v.parse().unwrap_or(0),
                "sites" => rep.sites = v.parse().unwrap_or(0),
                "read" => rep.read_s = secs(v),
                "write" => rep.write_s = secs(v),
                "create" => rep.create_s = secs(v),
                "fill" => rep.fill_s = secs(v),
                "h5write" => rep.h5_write_s = secs(v),
                "bytes" => rep.bytes = v.parse().unwrap_or(0),
                "peak_rss" => rep.peak_rss_mb = v.strip_suffix("MB").and_then(|x| x.parse().ok()).unwrap_or(0.0),
                _ => {}
            }
        }
    }
    rep
}

/// Run `taps-sc-extract-rs write-shard` for one shard and wait for it.
fn run_shard_writer_process(
    shard: usize,
    out: &Path,
    compression_name: &str,
    input: ShardInput,
) -> Result<WriterReport> {
    use std::io::Write;
    use std::process::{Command, Stdio};
    let exe = std::env::current_exe().context("locate taps-sc-extract-rs for write-shard")?;
    let mut cmd = Command::new(exe);
    cmd.arg("write-shard")
        .arg("--shard")
        .arg(shard.to_string())
        .arg("--out")
        .arg(out)
        .arg("--compression")
        .arg(compression_name)
        .arg("--compression-threads")
        .arg(BLOSC_THREADS.load(std::sync::atomic::Ordering::Relaxed).to_string())
        .env("HDF5_USE_FILE_LOCKING", "FALSE")
        .stdout(Stdio::null())
        .stderr(Stdio::piped());
    let chunks = match input {
        ShardInput::Dirs(dirs) => {
            for d in dirs {
                cmd.arg("--chunk-dir").arg(d);
            }
            cmd.stdin(Stdio::null());
            None
        }
        ShardInput::Chunks(c) => {
            cmd.arg("--stdin").stdin(Stdio::piped());
            Some(c)
        }
    };
    let mut child = cmd
        .spawn()
        .with_context(|| format!("spawn write-shard for shard {shard:03}"))?;
    if let Some(chunks) = chunks {
        let stdin = child.stdin.take().expect("piped stdin");
        let mut w = std::io::BufWriter::with_capacity(1 << 20, stdin);
        let sent: Result<()> = chunks.iter().try_for_each(|c| write_chunk_to(&mut w, c));
        let flushed = w.flush().map_err(anyhow::Error::from);
        drop(w);
        drop(chunks);
        if let Err(e) = sent.and(flushed) {
            let _ = child.kill();
            let _ = child.wait();
            return Err(e.context(format!("send chunks to write-shard {shard:03}")));
        }
    }
    let output = child
        .wait_with_output()
        .with_context(|| format!("wait for write-shard {shard:03}"))?;
    let stderr = String::from_utf8_lossy(&output.stderr);
    if !crate::parallel::QUIET_PROGRESS.load(std::sync::atomic::Ordering::Relaxed) || !output.status.success() {
        eprint!("{stderr}");
    }
    if !output.status.success() {
        anyhow::bail!("write-shard {shard:03} failed: {}", output.status);
    }
    Ok(parse_writer_line(shard, &stderr))
}

/// Pool of `n_writers` threads, each driving one child process at a time.
fn write_shards_with_processes(
    dir: &Path,
    inputs: Vec<(usize, ShardInput)>,
    compression_name: &str,
    n_writers: usize,
) -> Result<Vec<WriterReport>> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_writers.max(1))
        .thread_name(|i| format!("taps-h5p-{i}"))
        .build()
        .context("HDF5 writer process pool")?;
    pool.install(|| {
        inputs
            .into_par_iter()
            .map(|(s, input)| {
                let p = dir.join(format!("shard_{s:03}.h5"));
                run_shard_writer_process(s, &p, compression_name, input)
            })
            .collect()
    })
}

pub fn write_master_h5(out_dir: &Path, n_shards: usize) -> Result<std::path::PathBuf> {
    let master_path = out_dir.join("master.h5");
    let file = File::create(&master_path)
        .map_err(|e| anyhow::anyhow!("create {}: {e}", master_path.display()))?;
    write_version(&file)?;
    let cg = file.create_group("CG")?;
    let ch = file.create_group("CH")?;
    for s in 0..n_shards.max(1) {
        let shard_fn = format!("shard_{s:03}.h5");
        let shard_path = out_dir.join(&shard_fn);
        if !shard_path.exists() {
            continue;
        }
        let shard = File::open(&shard_path)
            .map_err(|e| anyhow::anyhow!("open {}: {e}", shard_path.display()))?;
        if let Ok(g) = shard.group("CG") {
            for bc in g.member_names()? {
                cg.link_external(&shard_fn, &format!("CG/{bc}"), &bc)
                    .map_err(|e| anyhow::anyhow!("CG external link {bc}: {e}"))?;
            }
        }
        if let Ok(g) = shard.group("CH") {
            for bc in g.member_names()? {
                ch.link_external(&shard_fn, &format!("CH/{bc}"), &bc)
                    .map_err(|e| anyhow::anyhow!("CH external link {bc}: {e}"))?;
            }
        }
    }
    file.flush().ok();
    Ok(master_path)
}

/// Read every chunk file of a shard's stream-mode bucket directories, then
/// delete them.
pub fn read_shard_temp_chunks(dirs: &[PathBuf]) -> Result<Vec<ChunkFile>> {
    let mut chunks = Vec::new();
    for dir in dirs {
        if dir.is_dir() {
            for p in list_chunk_files(dir)? {
                chunks.push(read_chunk_file(&p)?);
            }
            let _ = std::fs::remove_dir_all(dir);
        }
    }
    Ok(chunks)
}

/// Hash buckets per output shard. Barcodes are hashed into
/// `n_shards * BUCKETS_PER_SHARD` buckets during extraction; `pack_buckets`
/// then deals whole buckets to shards by size, so a few large cells do not
/// pile into one shard (the slowest shard sets the write time).
pub const BUCKETS_PER_SHARD: usize = 4;

pub fn bucket_count(n_shards: usize) -> usize {
    if n_shards <= 1 {
        1
    } else {
        n_shards * BUCKETS_PER_SHARD
    }
}

/// Assign buckets to `n_shards` shards, largest first to the lightest shard
/// (LPT). Each shard's bucket list is ascending. Deterministic for equal input.
pub fn pack_buckets(weights: &[u64], n_shards: usize) -> Vec<Vec<usize>> {
    let n_shards = n_shards.max(1);
    let mut order: Vec<usize> = (0..weights.len()).collect();
    order.sort_by(|&a, &b| weights[b].cmp(&weights[a]).then(a.cmp(&b)));
    let mut load = vec![0u64; n_shards];
    let mut groups: Vec<Vec<usize>> = vec![Vec::new(); n_shards];
    for b in order {
        let s = (0..n_shards).min_by_key(|&s| (load[s], s)).expect("n_shards > 0");
        load[s] += weights[b];
        groups[s].push(b);
    }
    for g in &mut groups {
        g.sort_unstable();
    }
    groups
}

/// Sites per in-memory bucket.
pub fn bucket_sites(per_bucket: &[Vec<ChunkFile>]) -> Vec<u64> {
    per_bucket.iter().map(|v| chunk_sites(v)).collect()
}

/// Bytes of chunk files per stream-mode bucket directory.
pub fn bucket_bytes(temp_dir: &Path, n_buckets: usize) -> Vec<u64> {
    (0..n_buckets)
        .map(|b| {
            std::fs::read_dir(shard_dir(temp_dir, b))
                .map(|rd| rd.flatten().filter_map(|e| e.metadata().ok()).map(|m| m.len()).sum())
                .unwrap_or(0)
        })
        .collect()
}

/// Concatenate in-memory buckets into shards.
pub fn group_buckets(per_bucket: Vec<Vec<ChunkFile>>, groups: &[Vec<usize>]) -> Vec<Vec<ChunkFile>> {
    let mut slots: Vec<Option<Vec<ChunkFile>>> = per_bucket.into_iter().map(Some).collect();
    groups
        .iter()
        .map(|g| g.iter().flat_map(|&b| slots[b].take().unwrap_or_default()).collect())
        .collect()
}

fn assemble_one_shard_from_temp(
    dirs: Vec<PathBuf>,
    shard: usize,
    shard_h5: &Path,
    compression: H5Compression,
) -> Result<WriterReport> {
    write_shard_timed(shard, shard_h5, || read_shard_temp_chunks(&dirs), compression)
}

/// `master.h5` and the seconds it took.
fn timed_master(dir: &Path, n_shards: usize) -> Result<(PathBuf, f64)> {
    let t = std::time::Instant::now();
    let p = write_master_h5(dir, n_shards)?;
    Ok((p, t.elapsed().as_secs_f64()))
}

/// Shard files written, per-shard reports, and seconds spent on `master.h5`.
pub struct HdfOutput {
    pub path: PathBuf,
    pub reports: Vec<WriterReport>,
    pub master_s: f64,
}

fn is_h5_filename(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("h5") || ext.eq_ignore_ascii_case("hdf5"))
}

/// Where HDF5 bytes go.
///
/// One shard writes `out` when it is an `.h5` file name, or `out/output.h5`
/// when `out` is a directory. Several shards write `shard_XXX.h5` plus
/// `master.h5` into a directory. An `.h5` file name in that case becomes a
/// sibling directory named by the file stem (`sample.h5` → `sample/`), so
/// `create_dir` is not asked to turn an existing file into a directory.
#[derive(Debug, PartialEq, Eq)]
pub enum OutputTarget {
    File(PathBuf),
    ShardDir(PathBuf),
}

pub fn resolve_output_target(out: &Path, n_shards: usize) -> Result<OutputTarget> {
    let n_shards = n_shards.max(1);
    if n_shards == 1 {
        if out.is_dir() {
            return Ok(OutputTarget::File(out.join("output.h5")));
        }
        if is_h5_filename(out) {
            if let Some(parent) = out.parent() {
                if !parent.as_os_str().is_empty() {
                    std::fs::create_dir_all(parent)
                        .with_context(|| format!("mkdir {}", parent.display()))?;
                }
            }
            return Ok(OutputTarget::File(out.to_path_buf()));
        }
        std::fs::create_dir_all(out).with_context(|| format!("mkdir {}", out.display()))?;
        return Ok(OutputTarget::File(out.join("output.h5")));
    }

    let dir = if out.is_dir() {
        out.to_path_buf()
    } else if is_h5_filename(out) {
        let stem = out.file_stem().unwrap_or_else(|| out.as_os_str());
        let parent = match out.parent() {
            Some(p) if !p.as_os_str().is_empty() => p,
            _ => Path::new("."),
        };
        let dir = parent.join(stem);
        if dir.is_file() {
            anyhow::bail!(
                "output {} with {n_shards} shards would write a directory at {}, but that path is a file",
                out.display(),
                dir.display()
            );
        }
        if out != dir.as_path() {
            eprintln!(
                "Multiple shards: writing shard_XXX.h5 and master.h5 to {} (output path {} names a file)",
                dir.display(),
                out.display()
            );
        }
        dir
    } else if out.exists() && out.is_file() {
        anyhow::bail!(
            "output {} is a file; with {n_shards} shards pass a directory or an .h5 path",
            out.display()
        );
    } else {
        out.to_path_buf()
    };
    std::fs::create_dir_all(&dir).with_context(|| format!("mkdir {}", dir.display()))?;
    Ok(OutputTarget::ShardDir(dir))
}

/// Bytes in the file just written, or in `shard_*.h5` + `master.h5` of a shard directory.
pub fn written_output_bytes(h5_path: &Path) -> u64 {
    let is_master = h5_path.file_name().and_then(|n| n.to_str()) == Some("master.h5");
    if !is_master {
        return std::fs::metadata(h5_path).map(|m| m.len()).unwrap_or(0);
    }
    let Some(dir) = h5_path.parent() else {
        return std::fs::metadata(h5_path).map(|m| m.len()).unwrap_or(0);
    };
    let mut total = 0u64;
    let Ok(rd) = std::fs::read_dir(dir) else {
        return 0;
    };
    for ent in rd.flatten() {
        let name = ent.file_name();
        let name = name.to_string_lossy();
        if name == "master.h5" || (name.starts_with("shard_") && name.ends_with(".h5")) {
            total += ent.metadata().map(|m| m.len()).unwrap_or(0);
        }
    }
    total
}

/// Write `shard_XXX.h5` (+ `master.h5` if more than one shard) from stream
/// temp chunks. `groups[s]` lists the bucket directories shard `s` reads.
pub fn assemble_hdf5_from_temp(
    temp_dir: &Path,
    out: &Path,
    groups: &[Vec<usize>],
    compression: H5Compression,
    compression_name: &str,
    max_writers: usize,
    mode: WriterMode,
) -> Result<HdfOutput> {
    std::env::set_var("HDF5_USE_FILE_LOCKING", "FALSE");
    let n_shards = groups.len().max(1);
    let n_writers = n_shards.min(max_writers.max(1));
    let dirs = |s: usize| -> Vec<PathBuf> {
        groups.get(s).map(|g| g.iter().map(|&b| shard_dir(temp_dir, b)).collect()).unwrap_or_default()
    };
    match resolve_output_target(out, n_shards)? {
        OutputTarget::File(path) => {
            let rep = assemble_one_shard_from_temp(dirs(0), 0, &path, compression)?;
            Ok(HdfOutput { path, reports: vec![rep], master_s: 0.0 })
        }
        OutputTarget::ShardDir(dir) => {
            let reports = match mode {
                WriterMode::Process => {
                    let inputs = (0..n_shards)
                        .map(|s| (s, ShardInput::Dirs(dirs(s))))
                        .collect();
                    write_shards_with_processes(&dir, inputs, compression_name, n_writers)?
                }
                WriterMode::Thread => {
                    let pool = rayon::ThreadPoolBuilder::new()
                        .num_threads(n_writers)
                        .thread_name(|i| format!("taps-h5-{i}"))
                        .build()
                        .context("HDF5 writer pool")?;
                    pool.install(|| {
                        (0..n_shards)
                            .into_par_iter()
                            .map(|s| {
                                let p = dir.join(format!("shard_{s:03}.h5"));
                                let rep = assemble_one_shard_from_temp(dirs(s), s, &p, compression)?;
                                eprintln!("{}", rep.child_line(&p));
                                Ok(rep)
                            })
                            .collect::<Result<Vec<_>>>()
                    })?
                }
            };
            let (path, master_s) = timed_master(&dir, n_shards)?;
            Ok(HdfOutput { path, reports, master_s })
        }
    }
}

pub fn assemble_hdf5_from_memory(
    per_shard: Vec<Vec<ChunkFile>>,
    out: &Path,
    compression: H5Compression,
    compression_name: &str,
    max_writers: usize,
    mode: WriterMode,
) -> Result<HdfOutput> {
    std::env::set_var("HDF5_USE_FILE_LOCKING", "FALSE");
    let n_shards = per_shard.len().max(1);
    let n_writers = n_shards.min(max_writers.max(1));
    match resolve_output_target(out, n_shards)? {
        OutputTarget::File(path) => {
            let chunks = per_shard.into_iter().next().unwrap_or_default();
            let rep = write_shard_timed(0, &path, || Ok(chunks), compression)?;
            Ok(HdfOutput { path, reports: vec![rep], master_s: 0.0 })
        }
        OutputTarget::ShardDir(dir) => {
            let reports = match mode {
                WriterMode::Process => {
                    let inputs = per_shard
                        .into_iter()
                        .enumerate()
                        .map(|(s, v)| (s, ShardInput::Chunks(v)))
                        .collect();
                    write_shards_with_processes(&dir, inputs, compression_name, n_writers)?
                }
                WriterMode::Thread => {
                    let pool = rayon::ThreadPoolBuilder::new()
                        .num_threads(n_writers)
                        .thread_name(|i| format!("taps-h5-{i}"))
                        .build()
                        .context("HDF5 writer pool")?;
                    pool.install(|| {
                        per_shard
                            .into_par_iter()
                            .enumerate()
                            .map(|(s, v)| {
                                let p = dir.join(format!("shard_{s:03}.h5"));
                                let rep = write_shard_timed(s, &p, || Ok(v), compression)?;
                                eprintln!("{}", rep.child_line(&p));
                                Ok(rep)
                            })
                            .collect::<Result<Vec<_>>>()
                    })?
                }
            };
            let (path, master_s) = timed_master(&dir, n_shards)?;
            Ok(HdfOutput { path, reports, master_s })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shard_io::ChunkCell;

    #[test]
    fn pack_buckets_balances_and_covers_all() {
        assert_eq!(bucket_count(1), 1);
        assert_eq!(bucket_count(8), 8 * BUCKETS_PER_SHARD);
        let w = [90u64, 10, 50, 40, 30, 30, 20, 0];
        let g = pack_buckets(&w, 3);
        let mut all: Vec<usize> = g.iter().flatten().copied().collect();
        all.sort_unstable();
        assert_eq!(all, (0..w.len()).collect::<Vec<_>>());
        let loads: Vec<u64> = g.iter().map(|b| b.iter().map(|&i| w[i]).sum()).collect();
        assert_eq!(loads.iter().sum::<u64>(), 270);
        assert!(*loads.iter().max().unwrap() <= 90, "{loads:?}");
        assert!(g.iter().all(|b| b.windows(2).all(|p| p[0] < p[1])));
        assert_eq!(pack_buckets(&w, 3), g);
        assert_eq!(pack_buckets(&[5], 1), vec![vec![0]]);
    }

    #[test]
    fn h5_filename_with_multiple_shards_uses_stem_directory() {
        let root = std::env::temp_dir().join("taps_rs_out_resolve");
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(&root).unwrap();
        let file_out = root.join("sample.h5");
        std::fs::write(&file_out, b"not a directory").unwrap();
        match resolve_output_target(&file_out, 4).unwrap() {
            OutputTarget::ShardDir(dir) => {
                assert_eq!(dir, root.join("sample"));
                assert!(dir.is_dir());
            }
            OutputTarget::File(_) => panic!("expected a shard directory"),
        }
        assert!(file_out.is_file(), "existing .h5 file must be left in place");

        let dir_out = root.join("already_dir");
        std::fs::create_dir_all(&dir_out).unwrap();
        match resolve_output_target(&dir_out, 2).unwrap() {
            OutputTarget::ShardDir(dir) => assert_eq!(dir, dir_out),
            OutputTarget::File(_) => panic!("expected the existing directory"),
        }
        match resolve_output_target(&file_out, 1).unwrap() {
            OutputTarget::File(p) => assert_eq!(p, file_out),
            OutputTarget::ShardDir(_) => panic!("one shard writes the .h5 file"),
        }
        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn s10_and_pct() {
        let r = MethRec::new("chr19", 100, 1, 1);
        assert_eq!(r.pos, 100);
        assert!((r.pct - 50.0).abs() < 1e-9);
        assert_eq!(r.t, 1);
        assert_eq!(r.c, 1);
        let s = r.chr.as_str();
        assert!(s.starts_with("chr19"));
    }

    #[test]
    fn sliced_write_matches_whole_write() {
        let dir = std::env::temp_dir().join("taps_rs_h5_sliced");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        // One cell bigger than WRITE_SLICE, split over two windows and contigs.
        let big_a: Vec<(u32, u32, u32)> = (0..200_000u32).map(|i| (i * 3 + 1, i % 4, i % 3)).collect();
        let big_b: Vec<(u32, u32, u32)> = (0..100_001u32).map(|i| (i * 5 + 2, i % 2, i % 5)).collect();
        let chunks = vec![
            ChunkFile {
                chunk_id: 1,
                contig: "chr2".into(),
                cells: vec![ChunkCell { barcode: "big".into(), cg: big_b.clone(), ch: vec![(7, 1, 0)] }],
            },
            ChunkFile {
                chunk_id: 0,
                contig: "chr1".into(),
                cells: vec![
                    ChunkCell { barcode: "big".into(), cg: big_a.clone(), ch: vec![] },
                    ChunkCell { barcode: "small".into(), cg: vec![(5, 0, 1)], ch: vec![] },
                ],
            },
        ];
        let mut ordered = chunks.clone();
        ordered.sort_by_key(|c| c.chunk_id);
        let whole = dir.join("whole.h5");
        write_shard_file(&whole, &records_from_chunks(&ordered), H5Compression::Lzf).unwrap();
        let sliced = dir.join("sliced.h5");
        assert_eq!(write_shard_from_chunks(&sliced, chunks, H5Compression::Lzf).unwrap().cells, 2);
        let (fw, fs) = (File::open(&whole).unwrap(), File::open(&sliced).unwrap());
        for name in ["CG/big/1", "CH/big/1", "CG/small/1"] {
            let a: Vec<MethRec> = fw.dataset(name).unwrap().read_raw().unwrap();
            let b: Vec<MethRec> = fs.dataset(name).unwrap().read_raw().unwrap();
            assert_eq!(a, b, "{name}");
        }
        assert_eq!(fs.dataset("CG/big/1").unwrap().size(), 300_001);
        assert!(!fs.link_exists("CH/small"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn roundtrip_single_file() {
        let dir = std::env::temp_dir().join("taps_rs_h5_rt");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let chunks = vec![ChunkFile {
            chunk_id: 0,
            contig: "chr19".into(),
            cells: vec![ChunkCell {
                barcode: "cellA".into(),
                cg: vec![(100, 1, 1), (200, 0, 2)],
                ch: vec![(150, 3, 0)],
            }],
        }];
        let cells = records_from_chunks(&chunks);
        let p = dir.join("out.h5");
        write_shard_file(&p, &cells, H5Compression::None).unwrap();
        let p2 = dir.join("out_stream.h5");
        write_shard_from_chunks(&p2, chunks.clone(), H5Compression::None).unwrap();
        let f2 = File::open(&p2).unwrap();
        for name in ["CG/cellA/1", "CH/cellA/1"] {
            let a: Vec<MethRec> = File::open(&p).unwrap().dataset(name).unwrap().read_raw().unwrap();
            let b: Vec<MethRec> = f2.dataset(name).unwrap().read_raw().unwrap();
            assert_eq!(a, b);
        }
        let f = File::open(&p).unwrap();
        let ds = f.dataset("CG/cellA/1").unwrap();
        let recs: Vec<MethRec> = ds.read_raw().unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].pos, 100);
        assert_eq!(recs[1].pos, 200);
        assert_eq!(recs[0].t, 1);
        assert_eq!(recs[0].c, 1);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
