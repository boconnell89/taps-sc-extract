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
        let tot = t + c;
        let pct = if tot == 0 {
            0.0
        } else {
            100.0 * f64::from(c) / f64::from(tot)
        };
        Self {
            chr: contig_s10(contig),
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

/// Write one cell's sites from window chunks (`(contig, sites)` in chunk_id
/// order) in `WRITE_SLICE`-record pieces. Same dataset as `write_cell_dataset`
/// over the concatenated records.
fn write_cell_dataset_sliced(
    file: &File,
    context: &str,
    barcode: &str,
    parts: &[(&str, &[(u32, u32, u32)])],
    compression: H5Compression,
) -> Result<()> {
    let n: usize = parts.iter().map(|(_, s)| s.len()).sum();
    if n == 0 {
        return Ok(());
    }
    let ds = create_cell_dataset(file, context, barcode, n, compression)?;
    let mut buf: Vec<MethRec> = Vec::with_capacity(n.min(WRITE_SLICE));
    let mut start = 0usize;
    let flush = |buf: &mut Vec<MethRec>, start: &mut usize| -> Result<()> {
        let end = *start + buf.len();
        ds.write_slice(buf.as_slice(), *start..end)
            .map_err(|e| anyhow::anyhow!("write {context}/{barcode}/1 [{}..{end}]: {e}", *start))?;
        *start = end;
        buf.clear();
        Ok(())
    };
    for &(contig, sites) in parts {
        for &(pos, t, c) in sites {
            buf.push(MethRec::new(contig, pos, t, c));
            if buf.len() == WRITE_SLICE {
                flush(&mut buf, &mut start)?;
            }
        }
    }
    if !buf.is_empty() {
        flush(&mut buf, &mut start)?;
    }
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
) -> Result<usize> {
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
        write_cell_dataset_sliced(&file, "CG", bc, &parts(true), compression)?;
        write_cell_dataset_sliced(&file, "CH", bc, &parts(false), compression)?;
    }
    file.flush().ok();
    Ok(bcs.len())
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
    /// Stream-mode temp directory `temp_dir/shard_XXX/`.
    Dir(PathBuf),
    /// Memory-mode payloads, sent on the child's stdin.
    Chunks(Vec<ChunkFile>),
}

/// What one `write-shard` child reported on its summary line.
#[derive(Clone, Debug, Default)]
pub struct WriterReport {
    pub shard: usize,
    pub sites: u64,
    pub peak_rss_mb: f64,
}

/// Peak RSS of `n_writers` concurrent children, bounded above by the sum of the
/// `n_writers` largest child peaks.
pub fn concurrent_writer_peak_mb(reports: &[WriterReport], n_writers: usize) -> f64 {
    let mut peaks: Vec<f64> = reports.iter().map(|r| r.peak_rss_mb).collect();
    peaks.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    peaks.iter().take(n_writers.max(1)).sum()
}

fn parse_writer_line(shard: usize, stderr: &str) -> WriterReport {
    let mut rep = WriterReport {
        shard,
        ..Default::default()
    };
    for tok in stderr.split_whitespace() {
        if let Some(v) = tok.strip_prefix("sites=") {
            rep.sites = v.parse().unwrap_or(0);
        } else if let Some(v) = tok.strip_prefix("peak_rss=").and_then(|v| v.strip_suffix("MB")) {
            rep.peak_rss_mb = v.parse().unwrap_or(0.0);
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
        ShardInput::Dir(d) => {
            cmd.arg("--chunk-dir").arg(d).stdin(Stdio::null());
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

/// Read every chunk file of one stream-mode shard, then delete the directory.
pub fn read_shard_temp_chunks(dir: &Path) -> Result<Vec<ChunkFile>> {
    let mut chunks = Vec::new();
    if dir.is_dir() {
        for p in list_chunk_files(dir)? {
            chunks.push(read_chunk_file(&p)?);
        }
        let _ = std::fs::remove_dir_all(dir);
    }
    Ok(chunks)
}

fn assemble_one_shard_from_temp(
    temp_dir: &Path,
    shard: usize,
    shard_h5: &Path,
    compression: H5Compression,
) -> Result<usize> {
    let chunks = read_shard_temp_chunks(&shard_dir(temp_dir, shard))?;
    write_shard_from_chunks(shard_h5, chunks, compression)
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

/// Write `shard_XXX.h5` (+ `master.h5` if `n_shards > 1`) from stream temp chunks.
pub fn assemble_hdf5_from_temp(
    temp_dir: &Path,
    out: &Path,
    n_shards: usize,
    compression: H5Compression,
    compression_name: &str,
    max_writers: usize,
    mode: WriterMode,
) -> Result<(PathBuf, Vec<WriterReport>)> {
    std::env::set_var("HDF5_USE_FILE_LOCKING", "FALSE");
    let n_shards = n_shards.max(1);
    let n_writers = n_shards.min(max_writers.max(1));
    match resolve_output_target(out, n_shards)? {
        OutputTarget::File(path) => {
            assemble_one_shard_from_temp(temp_dir, 0, &path, compression)?;
            Ok((path, Vec::new()))
        }
        OutputTarget::ShardDir(dir) => {
            let reports = match mode {
                WriterMode::Process => {
                    let inputs = (0..n_shards)
                        .map(|s| (s, ShardInput::Dir(shard_dir(temp_dir, s))))
                        .collect();
                    write_shards_with_processes(&dir, inputs, compression_name, n_writers)?
                }
                WriterMode::Thread => {
                    let pool = rayon::ThreadPoolBuilder::new()
                        .num_threads(n_writers)
                        .thread_name(|i| format!("taps-h5-{i}"))
                        .build()
                        .context("HDF5 writer pool")?;
                    let temp_dir = temp_dir.to_path_buf();
                    pool.install(|| {
                        (0..n_shards).into_par_iter().try_for_each(|s| {
                            let p = dir.join(format!("shard_{s:03}.h5"));
                            let n = assemble_one_shard_from_temp(&temp_dir, s, &p, compression)?;
                            eprintln!("HDF5 shard {s:03} cells={n} -> {}", p.display());
                            Ok::<(), anyhow::Error>(())
                        })
                    })?;
                    Vec::new()
                }
            };
            Ok((write_master_h5(&dir, n_shards)?, reports))
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
) -> Result<(PathBuf, Vec<WriterReport>)> {
    std::env::set_var("HDF5_USE_FILE_LOCKING", "FALSE");
    let n_shards = per_shard.len().max(1);
    let n_writers = n_shards.min(max_writers.max(1));
    match resolve_output_target(out, n_shards)? {
        OutputTarget::File(path) => {
            let chunks = per_shard.into_iter().next().unwrap_or_default();
            write_shard_from_chunks(&path, chunks, compression)?;
            Ok((path, Vec::new()))
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
                        per_shard.into_par_iter().enumerate().try_for_each(|(s, v)| {
                            let p = dir.join(format!("shard_{s:03}.h5"));
                            let n = write_shard_from_chunks(&p, v, compression)?;
                            eprintln!("HDF5 shard {s:03} cells={n} -> {}", p.display());
                            Ok::<(), anyhow::Error>(())
                        })
                    })?;
                    Vec::new()
                }
            };
            Ok((write_master_h5(&dir, n_shards)?, reports))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shard_io::ChunkCell;

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
        assert_eq!(write_shard_from_chunks(&sliced, chunks, H5Compression::Lzf).unwrap(), 2);
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
