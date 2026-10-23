//! `tune`: measure this machine (and optionally one dataset) and recommend flags.
//!
//! Machine calibration needs no input data. It writes a small synthetic TAPS
//! BAM + FASTA into the temp directory and runs the real extraction and shard
//! writers on it, so every host runs the same workload and results compare
//! across machines. Probes:
//!
//! - host: CPU model, logical CPUs, physical cores, NUMA nodes, RAM
//! - CPU: single-thread cost per read with and without BAQ (split into read,
//!   BAQ, pileup and calling), and scaling over worker counts
//! - memory: copy bandwidth on one thread and on all
//! - storage, for each distinct input / temp / output filesystem: sequential
//!   write and read, random 4 KiB reads, filesystem type
//! - HDF5: shard write rate in-process and with one process per shard
//!
//! With `-b/-f` it also samples the dataset (mapped reads per contig, calls and
//! cells per read on read-weighted windows) and predicts processing time, write
//! time and memory for the full run.
//!
//! `--save-profile` stores the machine part in
//! `~/.config/taps-sc-extract/machine_profile.toml`; `extract` then takes its
//! default worker and writer counts from it when the host matches.

use crate::autotune;
use crate::extract::{ExtractParams, WindowResult};
use crate::h5_out::{assemble_hdf5_from_memory, assemble_hdf5_from_temp, written_output_bytes, H5Compression, WriterMode};
use crate::parallel::{extract_memory_parallel, new_cancel_flag, QUIET_PROGRESS};
use crate::shard_io::{ChunkCell, ChunkFile};
use crate::window::{plan_windows, Window};
use anyhow::{Context, Result};
use rust_htslib::bam::{self, record::Cigar, record::CigarString, Read as _};
use std::fmt::Write as _;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::Ordering;
use std::time::Instant;

// ---------------------------------------------------------------- host

#[derive(Clone, Debug, Default)]
pub struct HostInfo {
    pub hostname: String,
    pub cpu_model: String,
    pub nproc: usize,
    pub physical_cores: usize,
    pub numa_nodes: usize,
    pub mem_total_gb: f64,
    pub mem_available_gb: f64,
}

impl HostInfo {
    /// Identity used to decide whether a saved profile applies to this host.
    pub fn host_id(&self) -> String {
        format!("{}|{}|{}", self.hostname, self.cpu_model, self.nproc)
    }
}

pub fn host_info() -> HostInfo {
    let cpuinfo = std::fs::read_to_string("/proc/cpuinfo").unwrap_or_default();
    let mut cpu_model = String::new();
    let mut cores = std::collections::HashSet::new();
    let mut phys = String::new();
    for line in cpuinfo.lines() {
        let Some((k, v)) = line.split_once(':') else {
            continue;
        };
        let (k, v) = (k.trim(), v.trim());
        match k {
            "model name" if cpu_model.is_empty() => cpu_model = v.to_string(),
            "physical id" => phys = v.to_string(),
            "core id" => {
                cores.insert((phys.clone(), v.to_string()));
            }
            _ => {}
        }
    }
    let nproc = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    let numa_nodes = std::fs::read_dir("/sys/devices/system/node")
        .map(|rd| {
            rd.flatten()
                .filter(|e| {
                    let n = e.file_name();
                    let n = n.to_string_lossy();
                    n.starts_with("node") && n[4..].chars().all(|c| c.is_ascii_digit())
                })
                .count()
        })
        .unwrap_or(1)
        .max(1);
    let meminfo = std::fs::read_to_string("/proc/meminfo").unwrap_or_default();
    let mem_kb = |key: &str| -> f64 {
        meminfo
            .lines()
            .find(|l| l.starts_with(key))
            .and_then(|l| l.split_whitespace().nth(1))
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(0.0)
    };
    let hostname = std::fs::read_to_string("/proc/sys/kernel/hostname")
        .map(|s| s.trim().to_string())
        .unwrap_or_default();
    HostInfo {
        hostname,
        cpu_model,
        nproc,
        physical_cores: if cores.is_empty() { nproc } else { cores.len() },
        numa_nodes,
        mem_total_gb: mem_kb("MemTotal:") / 1048576.0,
        mem_available_gb: mem_kb("MemAvailable:") / 1048576.0,
    }
}

// ---------------------------------------------------------------- synthetic data

/// xorshift64*: deterministic across hosts, no dependency.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
    fn chance(&mut self, p: f64) -> bool {
        ((self.next() >> 11) as f64 / (1u64 << 53) as f64) < p
    }
}

pub struct Synthetic {
    pub bam: PathBuf,
    pub fasta: PathBuf,
    pub contigs: Vec<(String, u64)>,
    pub reads: u64,
}

const SYN_READ_LEN: usize = 150;
const SYN_INSERT: i64 = 320;

/// Write a synthetic reference and a coordinate-sorted, indexed paired-end TAPS
/// BAM: OT pairs (99/147) carry C->T and OB pairs (83/163) G->A at methylated
/// sites (CpG ~70%, CH ~1%), 0.3% sequencing errors, 3% of reads with a 2 bp
/// deletion (exercises BAQ), barcodes `CELLnnnn:` in QNAME.
pub fn make_synthetic(dir: &Path, n_contigs: usize, contig_len: u64, n_pairs: u64, n_cells: u64) -> Result<Synthetic> {
    std::fs::create_dir_all(dir)?;
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let bases = [b'A', b'C', b'G', b'T'];
    let fasta = dir.join("tune_ref.fa");
    let mut refs: Vec<(String, Vec<u8>)> = Vec::new();
    for i in 0..n_contigs {
        let seq: Vec<u8> = (0..contig_len).map(|_| bases[rng.below(4) as usize]).collect();
        refs.push((format!("tchr{}", i + 1), seq));
    }
    {
        let mut f = std::io::BufWriter::new(std::fs::File::create(&fasta)?);
        let mut fai = String::new();
        let mut offset = 0u64;
        for (name, seq) in &refs {
            let header = format!(">{name}\n");
            f.write_all(header.as_bytes())?;
            offset += header.len() as u64;
            let _ = writeln!(fai, "{name}\t{}\t{offset}\t60\t61", seq.len());
            for line in seq.chunks(60) {
                f.write_all(line)?;
                f.write_all(b"\n")?;
                offset += line.len() as u64 + 1;
            }
        }
        f.flush()?;
        std::fs::write(dir.join("tune_ref.fa.fai"), fai)?;
    }

    // Pairs: (tid, pos) sorted later.
    struct Rec {
        tid: i32,
        pos: i64,
        mpos: i64,
        flags: u16,
        qname: Vec<u8>,
        cigar: CigarString,
        seq: Vec<u8>,
        tlen: i64,
    }
    let mut recs: Vec<Rec> = Vec::with_capacity(2 * n_pairs as usize);
    let pairs_per_contig = n_pairs / n_contigs as u64;
    for (tid, (_, rseq)) in refs.iter().enumerate() {
        let span = rseq.len() as i64 - SYN_INSERT - 4;
        for k in 0..pairs_per_contig {
            let start = rng.below(span as u64) as i64;
            let ot = rng.chance(0.5);
            let cell = rng.below(n_cells);
            let qname = format!("CELL{cell:05}:{tid}:{k}").into_bytes();
            let mate_start = start + SYN_INSERT - SYN_READ_LEN as i64;
            let make = |pos: i64, rng: &mut Rng| -> (CigarString, Vec<u8>) {
                let del = rng.chance(0.03);
                let p = pos as usize;
                let ref_len = if del { SYN_READ_LEN + 2 } else { SYN_READ_LEN };
                let window = &rseq[p..p + ref_len];
                let mut seq: Vec<u8> = Vec::with_capacity(SYN_READ_LEN);
                let mut cigar = Vec::new();
                if del {
                    let cut = 40 + rng.below(70) as usize;
                    seq.extend_from_slice(&window[..cut]);
                    seq.extend_from_slice(&window[cut + 2..]);
                    cigar.push(Cigar::Match(cut as u32));
                    cigar.push(Cigar::Del(2));
                    cigar.push(Cigar::Match((SYN_READ_LEN - cut) as u32));
                } else {
                    seq.extend_from_slice(window);
                    cigar.push(Cigar::Match(SYN_READ_LEN as u32));
                }
                for i in 0..seq.len() {
                    let b = seq[i];
                    let next = seq.get(i + 1).copied().unwrap_or(b'A');
                    let prev = if i > 0 { seq[i - 1] } else { b'A' };
                    if ot && b == b'C' {
                        let p_meth = if next == b'G' { 0.7 } else { 0.01 };
                        if rng.chance(p_meth) {
                            seq[i] = b'T';
                        }
                    } else if !ot && b == b'G' {
                        let p_meth = if prev == b'C' { 0.7 } else { 0.01 };
                        if rng.chance(p_meth) {
                            seq[i] = b'A';
                        }
                    }
                    if rng.chance(0.003) {
                        seq[i] = bases[rng.below(4) as usize];
                    }
                }
                (CigarString(cigar), seq)
            };
            let (c1, s1) = make(start, &mut rng);
            let (c2, s2) = make(mate_start, &mut rng);
            // OT: R1 forward (99) + R2 reverse (147). OB: R1 reverse (83) + R2 forward (163).
            let (f_left, f_right) = if ot { (99u16, 147u16) } else { (163u16, 83u16) };
            recs.push(Rec { tid: tid as i32, pos: start, mpos: mate_start, flags: f_left, qname: qname.clone(), cigar: c1, seq: s1, tlen: SYN_INSERT });
            recs.push(Rec { tid: tid as i32, pos: mate_start, mpos: start, flags: f_right, qname, cigar: c2, seq: s2, tlen: -SYN_INSERT });
        }
    }
    recs.sort_by_key(|r| (r.tid, r.pos));

    let mut header = bam::Header::new();
    for (name, seq) in &refs {
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", name).push_tag(b"LN", seq.len());
        header.push_record(&sq);
    }
    let bam_path = dir.join("tune_reads.bam");
    {
        let mut w = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam)
            .map_err(|e| anyhow::anyhow!("create synthetic BAM: {e}"))?;
        let qual = vec![37u8; SYN_READ_LEN];
        let mut rec = bam::Record::new();
        for r in &recs {
            rec.set(&r.qname, Some(&r.cigar), &r.seq, &qual);
            rec.set_tid(r.tid);
            rec.set_pos(r.pos);
            rec.set_mtid(r.tid);
            rec.set_mpos(r.mpos);
            rec.set_flags(r.flags);
            rec.set_mapq(60);
            rec.set_insert_size(r.tlen);
            w.write(&rec).map_err(|e| anyhow::anyhow!("write synthetic BAM: {e}"))?;
        }
    }
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1)
        .map_err(|e| anyhow::anyhow!("index synthetic BAM: {e}"))?;
    Ok(Synthetic {
        bam: bam_path,
        fasta,
        contigs: refs.iter().map(|(n, s)| (n.clone(), s.len() as u64)).collect(),
        reads: recs.len() as u64,
    })
}

// ---------------------------------------------------------------- CPU

#[derive(Clone, Debug, Default)]
pub struct CpuRun {
    pub workers: usize,
    pub baq: bool,
    pub secs: f64,
    pub reads_per_s: f64,
}

/// Machine probes time the BAQ HMM on every read (`TapsFull`): its cost per
/// read is a property of the CPU. How many reads skip it depends on the data,
/// which only the dataset sample (`-b`) measures.
fn tune_params(baq: bool, profile: bool) -> ExtractParams {
    ExtractParams {
        compute_baq: baq,
        baq_mode: crate::extract::BaqMode::TapsFull,
        decomp_threads: 0,
        accumulate: true,
        profile,
        ..ExtractParams::default()
    }
}

/// Windows for `workers` threads over the synthetic genome: ~6 per worker,
/// without the 2 Mb floor that `auto_chunk_size_bp` applies to real genomes.
fn synthetic_windows(fai: &crate::fasta::FastFaiReader, contigs: &[&str], total_bp: u64, workers: usize) -> Vec<Window> {
    let chunk = (total_bp / (6 * workers.max(1) as u64)).max(100_000);
    plan_windows(fai, contigs, chunk)
}

fn run_extract(bam: &Path, fasta: &Path, windows: &[Window], params: &ExtractParams, workers: usize, shards: usize) -> Result<(WindowResult, Vec<Vec<ChunkFile>>, f64)> {
    let cancel = new_cancel_flag();
    let t = Instant::now();
    let (r, mem) = extract_memory_parallel(&[bam.to_path_buf()], fasta, windows, None, params, shards, workers, &cancel, None)?;
    Ok((r, mem, t.elapsed().as_secs_f64()))
}

/// Worker counts to try: 1, 2, 4, ... up to nproc, plus physical cores and nproc.
fn worker_ladder(h: &HostInfo) -> Vec<usize> {
    let mut v = vec![1usize];
    let mut w = 2;
    while w < h.nproc {
        v.push(w);
        w *= 2;
    }
    v.push(h.physical_cores);
    v.push(h.nproc);
    v.sort_unstable();
    v.dedup();
    v
}

// ---------------------------------------------------------------- memory

/// Copy bandwidth in GB/s (bytes read + written) with `threads` threads.
fn mem_bandwidth(threads: usize) -> f64 {
    const PER_THREAD: usize = 32 << 20;
    const REPS: usize = 8;
    let t = Instant::now();
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| {
                let src = vec![1u8; PER_THREAD];
                let mut dst = vec![0u8; PER_THREAD];
                for _ in 0..REPS {
                    dst.copy_from_slice(&src);
                    std::hint::black_box(&mut dst);
                }
            });
        }
    });
    let bytes = (2 * PER_THREAD * REPS * threads) as f64;
    bytes / t.elapsed().as_secs_f64() / 1e9
}

// ---------------------------------------------------------------- storage

#[derive(Clone, Debug, Default)]
pub struct FsProbe {
    pub role: String,
    pub dir: PathBuf,
    pub fs_type: String,
    pub seq_write_mb_s: f64,
    pub seq_read_mb_s: f64,
    pub rand_read_iops: f64,
    /// `posix_fadvise(DONTNEED)` succeeded before the read probes. The host
    /// (e.g. a WSL VHDX) may still cache underneath.
    pub cache_dropped: bool,
}

fn fs_type_of(dir: &Path) -> String {
    let canon = std::fs::canonicalize(dir).unwrap_or_else(|_| dir.to_path_buf());
    let mounts = std::fs::read_to_string("/proc/mounts").unwrap_or_default();
    let mut best = (0usize, String::from("?"));
    for line in mounts.lines() {
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 3 {
            continue;
        }
        let mp = f[1];
        if canon.starts_with(mp) && mp.len() >= best.0 {
            best = (mp.len(), format!("{} ({})", f[2], f[0]));
        }
    }
    best.1
}

fn dev_of(dir: &Path) -> u64 {
    use std::os::unix::fs::MetadataExt;
    std::fs::metadata(dir).map(|m| m.dev()).unwrap_or(0)
}

fn drop_cache(f: &std::fs::File) -> bool {
    use std::os::unix::io::AsRawFd;
    unsafe { libc::posix_fadvise(f.as_raw_fd(), 0, 0, libc::POSIX_FADV_DONTNEED) == 0 }
}

fn probe_storage(role: &str, dir: &Path, size_mb: usize) -> Result<FsProbe> {
    std::fs::create_dir_all(dir)?;
    let path = dir.join(format!(".taps_tune_io_{}", std::process::id()));
    let block = vec![0x5Au8; 8 << 20];
    let n_blocks = (size_mb / 8).max(1);
    let t = Instant::now();
    {
        let mut f = std::fs::File::create(&path).with_context(|| format!("create {}", path.display()))?;
        for _ in 0..n_blocks {
            f.write_all(&block)?;
        }
        f.sync_all()?;
    }
    let write_s = t.elapsed().as_secs_f64();
    let bytes = (n_blocks * block.len()) as f64;
    let mut f = std::fs::File::open(&path)?;
    let dropped = drop_cache(&f);
    let mut buf = vec![0u8; 8 << 20];
    let t = Instant::now();
    while f.read(&mut buf)? > 0 {}
    let read_s = t.elapsed().as_secs_f64();
    let _ = drop_cache(&f);
    let mut rng = Rng(0xD1B5_4A32_D192_ED03);
    let mut small = [0u8; 4096];
    let blocks4k = (bytes as u64) / 4096;
    let t = Instant::now();
    let mut n = 0u64;
    while t.elapsed().as_secs_f64() < 2.0 {
        f.seek(SeekFrom::Start(rng.below(blocks4k) * 4096))?;
        f.read_exact(&mut small)?;
        n += 1;
    }
    let iops = n as f64 / t.elapsed().as_secs_f64();
    drop(f);
    let _ = std::fs::remove_file(&path);
    Ok(FsProbe {
        role: role.to_string(),
        dir: dir.to_path_buf(),
        fs_type: fs_type_of(dir),
        seq_write_mb_s: bytes / write_s / 1e6,
        seq_read_mb_s: bytes / read_s / 1e6,
        rand_read_iops: iops,
        cache_dropped: dropped,
    })
}

// ---------------------------------------------------------------- HDF5

#[derive(Clone, Debug, Default)]
pub struct H5Run {
    pub mode: String,
    pub writers: usize,
    pub secs: f64,
    pub sites_per_s: f64,
    pub datasets_per_s: f64,
    /// Compressed output bytes per site.
    pub bytes_per_site: f64,
    /// Stream mode only: seconds to write the temp chunk files (during a real
    /// run this overlaps extraction).
    pub temp_write_s: f64,
}

/// Synthetic shard payloads: `n_cells` cells, sites per cell drawn around
/// `mean_sites`, split over `shards` shards by barcode.
fn synthetic_payload(n_cells: usize, mean_sites: u32, shards: usize) -> (Vec<Vec<ChunkFile>>, u64) {
    let mut rng = Rng(0xA076_1D64_78BD_642F);
    let mut per_shard: Vec<Vec<ChunkCell>> = vec![Vec::new(); shards];
    let mut total = 0u64;
    for c in 0..n_cells {
        let n = (mean_sites / 2 + rng.below(mean_sites as u64) as u32).max(1);
        let n_cg = (n / 25).max(1);
        let sites = |rng: &mut Rng, k: u32| -> Vec<(u32, u32, u32)> {
            let mut pos = 0u32;
            (0..k)
                .map(|_| {
                    pos += 1 + rng.below(400) as u32;
                    let meth = rng.chance(0.7) as u32;
                    (pos, 1 - meth, meth)
                })
                .collect()
        };
        let cg = sites(&mut rng, n_cg);
        let ch = sites(&mut rng, n);
        total += (cg.len() + ch.len()) as u64;
        per_shard[c % shards].push(ChunkCell { barcode: format!("CELL{c:05}"), cg, ch });
    }
    let payload = per_shard
        .into_iter()
        .map(|cells| vec![ChunkFile { chunk_id: 0, contig: "tchr1".into(), cells }])
        .collect();
    (payload, total)
}

/// Write a synthetic payload with `writers` writers. `stream_temp` = Some(dir):
/// stream mode (chunk files written to `dir` first, then read back by the
/// writers); None: memory mode (payload piped to the writers).
fn h5_run(
    dir: &Path,
    mode: WriterMode,
    writers: usize,
    n_cells: usize,
    mean_sites: u32,
    stream_temp: Option<&Path>,
) -> Result<H5Run> {
    let shards = writers.max(2);
    let (payload, sites) = synthetic_payload(n_cells, mean_sites, shards);
    let label = format!("{}{}", mode.label(), if stream_temp.is_some() { "-stream" } else { "" });
    let out = dir.join(format!("h5_{label}_{writers}_{n_cells}"));
    let _ = std::fs::remove_dir_all(&out);
    std::fs::create_dir_all(&out)?;
    let mut temp_write_s = 0.0;
    let t;
    let res = match stream_temp {
        Some(tmp) => {
            let t0 = Instant::now();
            let tmp = tmp.join(format!("h5_stream_{}", std::process::id()));
            crate::shard_io::prepare_shard_dirs(&tmp, shards)?;
            for (s, chunks) in payload.iter().enumerate() {
                for c in chunks {
                    let p = crate::shard_io::chunk_path(&tmp, s, c.chunk_id);
                    let mut w = std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(&p)?);
                    crate::shard_io::write_chunk_to(&mut w, c)?;
                    w.flush()?;
                    w.get_ref().sync_all()?;
                }
            }
            temp_write_s = t0.elapsed().as_secs_f64();
            drop(payload);
            let groups: Vec<Vec<usize>> = (0..shards).map(|s| vec![s]).collect();
            t = Instant::now();
            let r = assemble_hdf5_from_temp(&tmp, &out, &groups, H5Compression::Lzf, "lzf", writers, mode);
            let _ = std::fs::remove_dir_all(&tmp);
            r?
        }
        None => {
            t = Instant::now();
            assemble_hdf5_from_memory(payload, &out, H5Compression::Lzf, "lzf", writers, mode)?
        }
    };
    let secs = t.elapsed().as_secs_f64();
    let bytes = written_output_bytes(&res.path);
    let _ = std::fs::remove_dir_all(&out);
    Ok(H5Run {
        mode: label,
        writers,
        secs,
        sites_per_s: sites as f64 / secs,
        datasets_per_s: (2 * n_cells) as f64 / secs,
        bytes_per_site: bytes as f64 / sites.max(1) as f64,
        temp_write_s,
    })
}

// ---------------------------------------------------------------- profile file

pub fn profile_path() -> Option<PathBuf> {
    let base = std::env::var_os("XDG_CONFIG_HOME")
        .map(PathBuf::from)
        .or_else(|| std::env::var_os("HOME").map(|h| PathBuf::from(h).join(".config")))?;
    Some(base.join("taps-sc-extract").join("machine_profile.toml"))
}

/// Defaults `extract` takes from a saved profile for this host.
#[derive(Clone, Debug, Default)]
pub struct SavedProfile {
    pub workers: Option<usize>,
    pub writers: Option<usize>,
    pub decomp_threads: Option<usize>,
    /// Fastest BAQ kernel measured on this host ("avx512", "avx2", "off").
    pub baq_simd: Option<String>,
}

fn simd_name(w: crate::baq_simd::Width) -> &'static str {
    match w.0 {
        8 => "avx512",
        4 => "avx2",
        _ => "off",
    }
}

/// Load the saved machine profile if it was made on this host.
pub fn load_profile_for_host() -> Option<(SavedProfile, PathBuf)> {
    let path = profile_path()?;
    let text = std::fs::read_to_string(&path).ok()?;
    let mut kv = std::collections::HashMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.starts_with('#') || line.starts_with('[') {
            continue;
        }
        if let Some((k, v)) = line.split_once('=') {
            kv.insert(k.trim().to_string(), v.trim().trim_matches('"').to_string());
        }
    }
    if kv.get("host_id")? != &host_info().host_id() {
        eprintln!(
            "note: {} was made on another host; ignoring it (rerun `taps-sc-extract --test --save-profile`)",
            path.display()
        );
        return None;
    }
    let num = |k: &str| kv.get(k).and_then(|v| v.parse::<usize>().ok());
    Some((
        SavedProfile {
            workers: num("recommended_workers"),
            writers: num("recommended_writers"),
            decomp_threads: num("recommended_decomp_threads"),
            baq_simd: kv.get("recommended_baq_simd").cloned(),
        },
        path,
    ))
}

// ---------------------------------------------------------------- dataset sampling

#[derive(Clone, Debug, Default)]
pub struct DatasetSample {
    pub mapped_reads: u64,
    pub sampled_reads: u64,
    pub sample_secs: f64,
    pub workers: usize,
    pub sites_per_read: f64,
    pub calls_per_read: f64,
    pub cells: usize,
    /// Largest cell's share of the sample's sites.
    pub max_cell_share: f64,
    pub baq_share: f64,
    /// Share of sampled reads that ran the BAQ HMM under the default `--baq taps`.
    pub baq_read_share: f64,
    /// Worker-thread seconds per read with all workers busy.
    pub thread_s_per_read: f64,
}

/// Mapped reads per contig from the BAM index.
fn mapped_per_contig(bam_path: &Path, contigs: &[&str]) -> Result<Vec<(String, u64, u64)>> {
    let mut bam = crate::input::open_indexed(bam_path)
        .map_err(|e| anyhow::anyhow!("open {}: {e}", bam_path.display()))?;
    let names: Vec<String> = bam
        .header()
        .target_names()
        .into_iter()
        .map(|n| String::from_utf8_lossy(n).into_owned())
        .collect();
    let stats = bam.index_stats().map_err(|e| anyhow::anyhow!("index stats: {e}"))?;
    let mut out = Vec::new();
    for (tid, len, mapped, _) in stats {
        if tid < 0 {
            continue;
        }
        if let Some(n) = names.get(tid as usize) {
            if contigs.contains(&n.as_str()) && mapped > 0 {
                out.push((n.clone(), len, mapped));
            }
        }
    }
    Ok(out)
}

fn sample_dataset(bam_path: &Path, fasta: &Path, workers: usize, target_reads: u64) -> Result<DatasetSample> {
    if crate::chemistry::detect(bam_path, false)?.with_ys > 0 {
        anyhow::bail!("--test dataset sampling supports TAPS only for now; this BAM carries BSBolt bisulfite tags");
    }
    let primary = crate::window::bam_primary_contigs(bam_path)?;
    let contigs: Vec<&str> = primary.iter().map(|s| s.as_str()).collect();
    let per = mapped_per_contig(bam_path, &contigs)?;
    let mapped: u64 = per.iter().map(|x| x.2).sum();
    anyhow::ensure!(mapped > 0, "no mapped reads on the primary chromosomes");
    // Read-weighted sample: windows spread evenly through each contig, sized so
    // the total is about `target_reads`.
    let frac = (target_reads as f64 / mapped as f64).min(1.0);
    // At least 4 windows per worker overall, so every worker stays busy and the
    // per-read cost is measured under full contention.
    let min_windows = 4 * workers as u64;
    let reads_per_window = ((target_reads / min_windows.max(1)).max(2_000)) as f64;
    let mut windows = Vec::new();
    let mut id = 0u32;
    for (name, len, m) in &per {
        let want_bp = (*len as f64 * frac) as u64;
        let n_win = ((*m as f64 * frac) / reads_per_window).ceil().max(1.0) as u64;
        let win_bp = (want_bp / n_win).max(50_000);
        for k in 0..n_win {
            let start = (len / (n_win + 1)) * (k + 1);
            let end = (start + win_bp).min(*len);
            if end > start {
                windows.push(Window { chunk_id: id, contig: name.clone(), start, end });
                id += 1;
            }
        }
    }
    let params = ExtractParams { baq_mode: crate::extract::BaqMode::Taps, ..tune_params(true, true) };
    let (r, mem, secs) = run_extract(bam_path, fasta, &windows, &params, workers, 1)?;
    anyhow::ensure!(r.timing.records > 0, "dataset sample read no records");
    let reads = r.timing.records.max(1);
    let calls = r.stats.cg_c + r.stats.cg_t + r.stats.ch_c + r.stats.ch_t;
    let site_count: u64 = mem
        .iter()
        .flatten()
        .flat_map(|c| c.cells.iter())
        .map(|c| (c.cg.len() + c.ch.len()) as u64)
        .sum();
    let mut per_cell: std::collections::HashMap<&str, u64> = std::collections::HashMap::new();
    for c in mem.iter().flatten().flat_map(|c| c.cells.iter()) {
        *per_cell.entry(c.barcode.as_str()).or_default() += (c.cg.len() + c.ch.len()) as u64;
    }
    let max_cell = per_cell.values().copied().max().unwrap_or(0);
    let t = &r.timing;
    Ok(DatasetSample {
        max_cell_share: max_cell as f64 / site_count.max(1) as f64,
        mapped_reads: mapped,
        sampled_reads: reads,
        sample_secs: secs,
        workers,
        sites_per_read: site_count as f64 / reads as f64,
        calls_per_read: calls as f64 / reads as f64,
        cells: r.n_cells(),
        baq_share: if t.total_s > 0.0 { t.baq_s / t.total_s } else { 0.0 },
        baq_read_share: t.baq_reads as f64 / reads as f64,
        thread_s_per_read: t.total_s / reads as f64,
    })
}

// ---------------------------------------------------------------- driver

pub struct TuneArgs {
    pub temp_dir: PathBuf,
    pub out_dir: PathBuf,
    pub input_dir: PathBuf,
    pub bam: Option<PathBuf>,
    pub fasta: Option<PathBuf>,
    pub quick: bool,
    pub save_profile: bool,
    pub json: Option<PathBuf>,
}

pub fn run_tune(a: &TuneArgs) -> Result<()> {
    QUIET_PROGRESS.store(true, Ordering::Relaxed);
    let host = host_info();
    eprintln!("== taps-sc-extract tune ==");
    eprintln!(
        "Host: {}  {}  {} logical CPUs / {} physical cores, {} NUMA node(s), {:.0} GB RAM ({:.0} GB available)",
        host.hostname, host.cpu_model, host.nproc, host.physical_cores, host.numa_nodes, host.mem_total_gb, host.mem_available_gb
    );
    let build = std::env::current_exe()
        .ok()
        .and_then(|p| p.parent().map(|d| d.join("BUILD_INFO")))
        .and_then(|p| std::fs::read_to_string(p).ok())
        .and_then(|s| s.lines().find(|l| l.starts_with("target_cpu=")).map(|l| l.to_string()))
        .unwrap_or_else(|| "target_cpu=unknown (not built by install.sh)".into());
    eprintln!("Build: {build}");

    let work = a.temp_dir.join(format!("taps_tune_{}", std::process::id()));
    std::fs::create_dir_all(&work)?;
    let result = run_probes(a, &host, &work, &build);
    let _ = std::fs::remove_dir_all(&work);
    QUIET_PROGRESS.store(false, Ordering::Relaxed);
    result
}

fn run_probes(a: &TuneArgs, host: &HostInfo, work: &Path, build: &str) -> Result<()> {
    // Synthetic genome: 4 contigs, sized so the all-threads scaling run lasts a
    // few seconds (~40k BAQ reads/s per core). Worker counts below 8 run on the
    // first contig only, so the one-thread run stays short.
    let secs_at_full = if a.quick { 1.5 } else { 4.0 };
    let est_rate = 40_000.0 * host.physical_cores as f64 * 1.2;
    let pairs = ((est_rate * secs_at_full / 2.0) as u64).clamp(50_000, 4_000_000);
    let contig_len = (pairs * 2 * SYN_READ_LEN as u64 / 4 / 3).max(2_000_000);
    eprintln!("\n[1/5] Synthetic TAPS BAM ({} read pairs on 4 x {} Mb) ...", pairs, contig_len / 1_000_000);
    let syn = make_synthetic(&work.join("syn"), 4, contig_len, pairs, 500)?;
    let fai = crate::fasta::FastFaiReader::open(&syn.fasta)?;
    let names: Vec<&str> = syn.contigs.iter().map(|(n, _)| n.as_str()).collect();
    let total_bp: u64 = syn.contigs.iter().map(|c| c.1).sum();

    // --- CPU: single-thread split, then scaling.
    eprintln!("[2/5] CPU: per-read cost and worker scaling ...");
    let small_names = &names[..1];
    let small_bp = syn.contigs[0].1;
    let w1 = synthetic_windows(&fai, small_names, small_bp, 1);
    let (r_baq, _, s_baq) = run_extract(&syn.bam, &syn.fasta, &w1, &tune_params(true, true), 1, 1)?;
    let (r_nobaq, _, s_nobaq) = run_extract(&syn.bam, &syn.fasta, &w1, &tune_params(false, true), 1, 1)?;
    let reads = r_baq.timing.records.max(1) as f64;
    let us = |s: f64| s / reads * 1e6;
    let t = &r_baq.timing;
    let stage_us = [
        ("read+inflate", us(t.read_s)),
        ("BAQ", us(t.baq_s)),
        ("pileup", us(t.pileup_s)),
        ("calling", us(t.call_s)),
    ];
    let mut scaling: Vec<CpuRun> = Vec::new();
    for w in worker_ladder(host) {
        let wins = if w < 8 {
            synthetic_windows(&fai, small_names, small_bp, w)
        } else {
            synthetic_windows(&fai, &names, total_bp, w)
        };
        let (r, _, secs) = run_extract(&syn.bam, &syn.fasta, &wins, &tune_params(true, false), w, 1)?;
        let n = r.timing.records.max(1) as f64;
        scaling.push(CpuRun { workers: w, baq: true, secs, reads_per_s: n / secs });
    }
    let base = scaling[0].reads_per_s;
    // Knee: smallest worker count within 5% of the best throughput.
    let best = scaling.iter().map(|r| r.reads_per_s).fold(0.0, f64::max);
    let knee = scaling
        .iter()
        .find(|r| r.reads_per_s >= 0.95 * best)
        .map(|r| r.workers)
        .unwrap_or(host.nproc);
    let smt_gain = {
        let at = |w: usize| scaling.iter().find(|r| r.workers == w).map(|r| r.reads_per_s);
        match (at(host.physical_cores), at(host.nproc)) {
            (Some(p), Some(l)) if host.nproc > host.physical_cores => Some(l / p),
            _ => None,
        }
    };
    let inflate_share = t.read_s / t.total_s.max(1e-9);
    // BAQ off at the knee: how much of the scaled-up throughput BAQ costs.
    let knee_nobaq = {
        let wins = if knee < 8 {
            synthetic_windows(&fai, small_names, small_bp, knee)
        } else {
            synthetic_windows(&fai, &names, total_bp, knee)
        };
        let (r, _, secs) = run_extract(&syn.bam, &syn.fasta, &wins, &tune_params(false, false), knee, 1)?;
        r.timing.records.max(1) as f64 / secs
    };
    // BAQ kernels at the knee: wider vectors can lower the clock (AVX-512 on
    // some Intel parts), so time each one this CPU runs and keep the fastest.
    let mut kernels: Vec<(crate::baq_simd::Width, f64)> = Vec::new();
    {
        use crate::baq_simd::Width;
        let wins = if knee < 8 {
            synthetic_windows(&fai, small_names, small_bp, knee)
        } else {
            synthetic_windows(&fai, &names, total_bp, knee)
        };
        for w in [Width::HTSLIB, Width::AVX2, Width::AVX512] {
            if !w.supported() {
                continue;
            }
            let p = ExtractParams { baq_width: w, ..tune_params(true, true) };
            let (r, _, _) = run_extract(&syn.bam, &syn.fasta, &wins, &p, knee, 1)?;
            // BAQ thread-seconds per read: the kernel's own cost under full load.
            kernels.push((w, r.timing.records.max(1) as f64 / r.timing.baq_s.max(1e-9)));
        }
    }
    let best_kernel = kernels
        .iter()
        .copied()
        .fold(None::<(crate::baq_simd::Width, f64)>, |b, k| match b {
            Some(x) if x.1 >= k.1 => Some(x),
            _ => Some(k),
        })
        .map(|k| k.0)
        .unwrap_or(crate::baq_simd::Width::HTSLIB);
    // GPU BAQ, if a GPU is usable: synthetic processing rate at the knee with
    // the HMM on every read, GPU vs the best CPU kernel. Reported only; the
    // default stays the CPU (--baq-device gpu opts in).
    let mut gpu_probe: Option<(String, f64, f64)> = None;
    if let Ok(svc) = crate::baq_gpu::start(crate::baq_gpu::GpuBackend::Auto, "", 1024) {
        let wins = if knee < 8 {
            synthetic_windows(&fai, small_names, small_bp, knee)
        } else {
            synthetic_windows(&fai, &names, total_bp, knee)
        };
        let cpu_p = ExtractParams { baq_width: best_kernel, ..tune_params(true, false) };
        let gpu_p = ExtractParams { baq_device: crate::baq_gpu::BaqDevice::Gpu, ..cpu_p.clone() };
        let (rc, _, sc) = run_extract(&syn.bam, &syn.fasta, &wins, &cpu_p, knee, 1)?;
        let (rg, _, sg) = run_extract(&syn.bam, &syn.fasta, &wins, &gpu_p, knee, 1)?;
        gpu_probe = Some((
            format!("{} ({:?})", svc.device, svc.backend),
            rc.timing.records as f64 / sc,
            rg.timing.records as f64 / sg,
        ));
    }

    // --- Memory bandwidth.
    eprintln!("[3/5] Memory bandwidth ...");
    let bw1 = mem_bandwidth(1);
    let bwn = mem_bandwidth(host.nproc);

    // --- Storage.
    eprintln!("[4/5] Storage (input / temp / output) ...");
    let size_mb = if a.quick { 256 } else { 1024 };
    let mut fs: Vec<FsProbe> = Vec::new();
    let mut seen_dev = Vec::new();
    for (role, dir) in [("input", &a.input_dir), ("temp", &a.temp_dir), ("output", &a.out_dir)] {
        let dev = dev_of(dir);
        if let Some((_, r)) = seen_dev.iter().find(|(d, _)| *d == dev) {
            let mut p: FsProbe = fs.iter().find(|f: &&FsProbe| &f.role == r).cloned().unwrap_or_default();
            p.role = format!("{role} (same filesystem as {r})");
            p.dir = dir.to_path_buf();
            fs.push(p);
            continue;
        }
        match probe_storage(role, dir, size_mb) {
            Ok(p) => {
                seen_dev.push((dev, role.to_string()));
                fs.push(p);
            }
            Err(e) => eprintln!("  {role}: {} not probed: {e}", dir.display()),
        }
    }

    // --- HDF5 writers.
    eprintln!("[5/5] HDF5 shard writers ...");
    let (n_cells, mean_sites) = if a.quick { (1000usize, 30_000u32) } else { (2000usize, 50_000u32) };
    let h5_dir = a.out_dir.join(format!(".taps_tune_h5_{}", std::process::id()));
    let mut h5: Vec<H5Run> = Vec::new();
    // Up to every logical CPU: on a 64-thread server the old ladder stopped at
    // the 32 physical cores, and "fastest" then just meant "largest tried".
    let proc_ladder: Vec<usize> = {
        let mut v = vec![1, 4, host.physical_cores.min(16), host.physical_cores, host.nproc];
        v.sort_unstable();
        v.dedup();
        v
    };
    // Each writer gets real work (>= 5M sites, 2M with --quick), so the probe
    // measures writing rather than process start-up, within a quarter of the
    // memory budget.
    let per_writer_sites: u64 = if a.quick { 2_000_000 } else { 5_000_000 };
    let budget_sites = (0.25 * autotune::system_memory_budget_gb(None, None).0 * 1e9 / 12.0) as u64;
    let cells_for = |w: usize| -> usize {
        let want = (w as u64 * per_writer_sites).div_ceil(mean_sites as u64) as usize;
        let most = (budget_sites / mean_sites as u64).max(n_cells as u64) as usize;
        want.clamp(n_cells, most)
    };
    h5.push(h5_run(&h5_dir, WriterMode::Thread, 4, n_cells, mean_sites, None)?);
    for &w in &proc_ladder {
        h5.push(h5_run(&h5_dir, WriterMode::Process, w, cells_for(w), mean_sites, None)?);
    }
    let best_h5 = h5
        .iter()
        .filter(|r| r.mode == "process")
        .fold(None::<&H5Run>, |b, r| match b {
            Some(x) if x.sites_per_s >= r.sites_per_s * 1.05 => Some(x),
            Some(x) if r.sites_per_s < x.sites_per_s * 1.05 && r.writers > x.writers => Some(x),
            _ => Some(r),
        })
        .cloned()
        .unwrap_or_default();
    // Same payload at the best writer count through stream mode (temp files on
    // the temp filesystem), and a few very large cells (low cell count, deep
    // coverage: each cell is one writer's serial work).
    let bw = best_h5.writers.max(1);
    let h5_stream = h5_run(&h5_dir, WriterMode::Process, bw, cells_for(bw), mean_sites, Some(&a.temp_dir))?;
    let (dense_cells, dense_sites) = if a.quick { (2 * bw, 300_000u32) } else { (2 * bw, 1_000_000u32) };
    let h5_dense = h5_run(&h5_dir, WriterMode::Process, bw, dense_cells, dense_sites, None)?;
    let _ = std::fs::remove_dir_all(&h5_dir);

    // --- Recommendations.
    let rec_workers = knee;
    let rec_decomp = if inflate_share > 0.15 { 2 } else { 0 };
    // The profile's writer count caps the run's planner. When the largest count
    // tried was still the fastest, no limit was found: cap at the CPU count and
    // let the planner fit writers to memory per run.
    let top = proc_ladder.iter().copied().max().unwrap_or(1);
    let rec_writers = if best_h5.writers >= top { host.nproc.max(1) } else { best_h5.writers.max(1) };

    // --- Optional dataset sample.
    let dataset = match (&a.bam, &a.fasta) {
        (Some(b), Some(f)) => {
            eprintln!("[+] Dataset sample on {} ...", b.display());
            Some(sample_dataset(b, f, rec_workers, if a.quick { 300_000 } else { 1_500_000 })?)
        }
        _ => None,
    };

    // --- Report.
    let mut rpt = String::new();
    let _ = writeln!(rpt, "\n================ Machine profile ================");
    let _ = writeln!(rpt, "CPU, one worker, synthetic TAPS reads (us per read):");
    for (k, v) in stage_us {
        let _ = writeln!(rpt, "  {k:<13} {v:7.1}");
    }
    let _ = writeln!(rpt, "  total BAQ on  {:7.1}   BAQ off {:7.1}", us(s_baq), s_nobaq / r_nobaq.timing.records.max(1) as f64 * 1e6);
    let _ = writeln!(rpt, "Worker scaling (BAQ on):");
    for r in &scaling {
        let _ = writeln!(
            rpt,
            "  W={:<4} {:>9.0} reads/s   speedup {:5.2}x   efficiency {:5.1}%",
            r.workers,
            r.reads_per_s,
            r.reads_per_s / base,
            100.0 * r.reads_per_s / base / r.workers as f64
        );
    }
    if let Some(g) = smt_gain {
        let _ = writeln!(rpt, "  SMT: {} logical vs {} physical gives {:.2}x", host.nproc, host.physical_cores, g);
    }
    let knee_baq = scaling.iter().find(|r| r.workers == knee).map(|r| r.reads_per_s).unwrap_or(best);
    let _ = writeln!(
        rpt,
        "  W={knee} BAQ off: {knee_nobaq:.0} reads/s ({:.1}x BAQ on)",
        knee_nobaq / knee_baq.max(1.0)
    );
    if !kernels.is_empty() {
        let base_k = kernels[0].1.max(1e-9);
        let _ = writeln!(rpt, "BAQ kernels at W={knee} (BAQ reads per thread-second, HMM on every read):");
        for (w, rate) in &kernels {
            let _ = writeln!(
                rpt,
                "  {:32} {:>9.0}  {:.2}x{}",
                w.label(),
                rate,
                rate / base_k,
                if *w == best_kernel { "  <- fastest (--baq-simd auto uses it with --save-profile)" } else { "" }
            );
        }
    }
    match &gpu_probe {
        Some((dev, cpu, gpu)) => {
            let _ = writeln!(
                rpt,
                "GPU BAQ at W={knee} (HMM on every read): {dev}: {gpu:.0} reads/s vs CPU {cpu:.0} reads/s ({:+.0}%){}",
                (gpu / cpu - 1.0) * 100.0,
                if gpu > cpu { " -> --baq-device gpu is faster here" } else { " -> keep --baq-device cpu (default)" }
            );
        }
        None => {
            let _ = writeln!(rpt, "GPU BAQ: no usable GPU (backends built: {})", crate::baq_gpu::built_backends().join(", "));
        }
    }
    let _ = writeln!(rpt, "Memory copy bandwidth: {bw1:.1} GB/s (1 thread), {bwn:.1} GB/s ({} threads)", host.nproc);
    let _ = writeln!(rpt, "Storage ({} MB probe):", size_mb);
    for f in &fs {
        // A read several times faster than the synced write was served from a
        // cache the probe could not drop (e.g. the WSL host caching the VHDX).
        let cached = !f.cache_dropped || f.seq_read_mb_s > 3.0 * f.seq_write_mb_s;
        let _ = writeln!(
            rpt,
            "  {:<34} {:>7.0} MB/s write {:>7.0} MB/s read {:>8.0} IOPS 4K{}  [{}]",
            f.role,
            f.seq_write_mb_s,
            f.seq_read_mb_s,
            f.rand_read_iops,
            if cached { " (reads likely cached)" } else { "" },
            f.fs_type
        );
    }
    let _ = writeln!(
        rpt,
        "HDF5 shard writes (LZF, cells of ~{} sites; at least {:.0}M sites per process writer):",
        mean_sites,
        per_writer_sites as f64 / 1e6
    );
    for r in &h5 {
        let _ = writeln!(
            rpt,
            "  {:<7} x{:<3} {:6.2} s  {:>6.1} M sites/s  {:>7.0} datasets/s  {:>5.1} B/site",
            r.mode, r.writers, r.secs, r.sites_per_s / 1e6, r.datasets_per_s, r.bytes_per_site
        );
    }
    let _ = writeln!(
        rpt,
        "  stream x{:<3} {:6.2} s  {:>6.1} M sites/s  (+{:.2} s writing temp chunk files; memory mode {:.2} s)",
        h5_stream.writers,
        h5_stream.secs,
        h5_stream.sites_per_s / 1e6,
        h5_stream.temp_write_s,
        best_h5.secs
    );
    let _ = writeln!(
        rpt,
        "  {} large cells x ~{} sites, x{}: {:6.2} s  {:>6.1} M sites/s  ({:.1} M sites/s per writer on one cell)",
        dense_cells,
        dense_sites,
        h5_dense.writers,
        h5_dense.secs,
        h5_dense.sites_per_s / 1e6,
        h5_dense.sites_per_s / h5_dense.writers.max(1) as f64 / 1e6
    );

    let _ = writeln!(rpt, "\nLimiting factors on this machine:");
    let _ = writeln!(
        rpt,
        "  Extraction is CPU-bound (inflate {:.0}% of worker time; BAQ {:.0}%).{}",
        100.0 * inflate_share,
        100.0 * t.baq_s / t.total_s.max(1e-9),
        if inflate_share > 0.15 { " Inflate is large: try --decomp-threads 2." } else { "" }
    );
    let cpu_rate = scaling.iter().find(|r| r.workers == rec_workers).map(|r| r.reads_per_s).unwrap_or(best);
    if let Some(inp) = fs.iter().find(|f| f.role.starts_with("input")) {
        // ~70 bytes of BAM per read on disk (typical compressed TAPS paired-end).
        let need_mb_s = cpu_rate * 70.0 / 1e6;
        let _ = writeln!(
            rpt,
            "  Input read: extraction consumes ~{need_mb_s:.0} MB/s of BAM at W={rec_workers}; input storage reads {:.0} MB/s sequential, {:.0} IOPS random{}",
            inp.seq_read_mb_s,
            inp.rand_read_iops,
            if inp.seq_read_mb_s < 2.0 * need_mb_s || inp.rand_read_iops < 2000.0 {
                " -> storage can limit; copy the BAM to faster local storage"
            } else {
                " -> not limiting"
            }
        );
    }
    let _ = writeln!(
        rpt,
        "  HDF5 write: best {} process writers ({:.1} M sites/s); in-process writers share one HDF5 lock ({:.1} M sites/s)",
        rec_writers,
        best_h5.sites_per_s / 1e6,
        h5.iter().find(|r| r.mode == "thread").map(|r| r.sites_per_s / 1e6).unwrap_or(0.0)
    );
    if host.numa_nodes > 1 {
        let _ = writeln!(rpt, "  {} NUMA nodes: compare one run under `numactl --interleave=all`", host.numa_nodes);
    }

    let budget = autotune::system_memory_budget_gb(None, None).0;
    if let Some(d) = &dataset {
        let rate = d.sampled_reads as f64 / d.sample_secs;
        // Worker CPU per read (measured with all workers busy) over the worker
        // count. Wall time on the short sample is dominated by startup and is
        // not used. The 1.25 (window setup, BAQ contig loads, run tail) and the
        // 1.5 on the synthetic HDF5 rate (real cells vary in size) are
        // calibrated on one full run: 5base_ucbTn5_sp2, 53.5M reads, 7,351
        // cells, 32 threads, 2026-10-15 (41 s processing, 14 s write).
        let proc_s = d.mapped_reads as f64 * d.thread_s_per_read / d.workers as f64 * 1.25;
        let sites = d.sites_per_read * d.mapped_reads as f64;
        let write_s = sites / best_h5.sites_per_s.max(1.0) * 1.5;
        let plan = autotune::plan_memory(&autotune::MemoryInputs {
            mode: "auto",
            budget_gb: budget,
            mapped_reads: d.mapped_reads,
            site_share: 1.0,
            cell_count: d.cells,
            genome_mb: 0,
            n_shards: rec_writers.max(16),
            n_workers: rec_workers,
            cli_writers: 0,
            profile_writers: Some(rec_writers),
            in_process: false,
            streamable: true,
        });
        let mode = plan.mode.as_str();
        let _ = writeln!(rpt, "\n================ Dataset ================");
        let _ = writeln!(
            rpt,
            "Mapped reads: {}; sampled {} reads in {:.1} s at W={} ({:.0} reads/s); {} cells in sample",
            d.mapped_reads, d.sampled_reads, d.sample_secs, d.workers, rate, d.cells
        );
        let _ = writeln!(
            rpt,
            "Per read: {:.2} calls, {:.2} new sites; BAQ HMM ran on {:.0}% of reads (--baq taps), {:.0}% of worker time",
            d.calls_per_read,
            d.sites_per_read,
            100.0 * d.baq_read_share,
            100.0 * d.baq_share
        );
        let _ = writeln!(
            rpt,
            "Predicted: processing ~{proc_s:.0} s, HDF5 write ~{write_s:.0} s, ~{:.0}M sites; --memory-mode {mode} with {} writers, est. peak ~{:.0} GB (budget {budget:.0} GB)",
            sites / 1e6,
            plan.writers,
            plan.peak_gb
        );
        let _ = writeln!(
            rpt,
            "  ({:.1} worker-us per read under full load; overhead factors fitted to one run of one dataset and not yet checked on another)",
            d.thread_s_per_read * 1e6
        );
        // One cell is written by one writer, so the largest cell sets a floor
        // on write time that more shards cannot lower.
        let per_writer = h5_dense.sites_per_s / h5_dense.writers.max(1) as f64;
        let biggest = d.max_cell_share * sites;
        let _ = writeln!(
            rpt,
            "Cells: largest holds {:.1}% of sampled sites (~{:.0}M sites) -> one writer needs ~{:.1} s for it alone; {}",
            100.0 * d.max_cell_share,
            biggest / 1e6,
            biggest / per_writer.max(1.0),
            if d.cells < 2 * rec_writers {
                format!("only {} cells: use --shards <= {} (more shards cannot help)", d.cells, (d.cells / 2).max(1))
            } else {
                format!("--shards {} (writers are balanced by site count)", rec_writers.max(16))
            }
        );
        if let Some(o) = fs.iter().find(|f| f.role.starts_with("output")) {
            let out_mb_s = sites * best_h5.bytes_per_site / write_s.max(1e-9) / 1e6;
            let _ = writeln!(
                rpt,
                "Output: ~{:.1} GB of HDF5 (LZF) at ~{out_mb_s:.0} MB/s; output storage writes {:.0} MB/s sequential ({:.0}% of it){}",
                sites * best_h5.bytes_per_site / 1e9,
                o.seq_write_mb_s,
                100.0 * out_mb_s / o.seq_write_mb_s.max(1.0),
                if out_mb_s > 0.7 * o.seq_write_mb_s { " -> storage limits the write" } else { " -> CPU-bound, not storage" }
            );
        }
    }

    let _ = writeln!(
        rpt,
        "\nRecommended: -t {rec_workers} --decomp-threads {rec_decomp}; leave --shards and --max-writer-threads unset. \
         Each run fits as many writers as memory allows (up to {rec_writers}, saved with --save-profile) and splits \
         the output into more files when that lets more writers run at once. Giving --max-writer-threads makes it \
         a fixed count, used even if it does not fit."
    );
    eprint!("{rpt}");

    // --- Save.
    let mut toml = String::new();
    let _ = writeln!(toml, "# taps-sc-extract machine profile, written by `tune --save-profile`");
    let _ = writeln!(toml, "created = \"{}\"", chrono_like_now());
    let _ = writeln!(toml, "host_id = \"{}\"", host.host_id());
    let _ = writeln!(toml, "cpu_model = \"{}\"", host.cpu_model);
    let _ = writeln!(toml, "nproc = {}", host.nproc);
    let _ = writeln!(toml, "physical_cores = {}", host.physical_cores);
    let _ = writeln!(toml, "numa_nodes = {}", host.numa_nodes);
    let _ = writeln!(toml, "mem_total_gb = {:.1}", host.mem_total_gb);
    let _ = writeln!(toml, "build = \"{build}\"");
    let _ = writeln!(toml, "us_per_read_baq_1t = {:.2}", us(s_baq));
    for (k, v) in stage_us {
        let _ = writeln!(toml, "us_per_read_{} = {:.2}", k.replace('+', "_"), v);
    }
    for r in &scaling {
        let _ = writeln!(toml, "reads_per_s_w{} = {:.0}", r.workers, r.reads_per_s);
    }
    let _ = writeln!(toml, "mem_bw_1t_gb_s = {bw1:.1}");
    let _ = writeln!(toml, "mem_bw_all_gb_s = {bwn:.1}");
    for f in &fs {
        let key = f.role.split_whitespace().next().unwrap_or("fs");
        let _ = writeln!(toml, "fs_{key}_type = \"{}\"", f.fs_type);
        let _ = writeln!(toml, "fs_{key}_seq_write_mb_s = {:.0}", f.seq_write_mb_s);
        let _ = writeln!(toml, "fs_{key}_seq_read_mb_s = {:.0}", f.seq_read_mb_s);
        let _ = writeln!(toml, "fs_{key}_rand_read_iops = {:.0}", f.rand_read_iops);
    }
    for r in &h5 {
        let _ = writeln!(toml, "h5_{}_{}_sites_per_s = {:.0}", r.mode, r.writers, r.sites_per_s);
    }
    let _ = writeln!(toml, "h5_stream_{}_sites_per_s = {:.0}", h5_stream.writers, h5_stream.sites_per_s);
    let _ = writeln!(toml, "h5_stream_temp_write_s = {:.2}", h5_stream.temp_write_s);
    let _ = writeln!(toml, "h5_large_cells_{}_sites_per_s = {:.0}", h5_dense.writers, h5_dense.sites_per_s);
    let _ = writeln!(toml, "h5_bytes_per_site_lzf = {:.2}", best_h5.bytes_per_site);
    let _ = writeln!(toml, "reads_per_s_w{knee}_nobaq = {knee_nobaq:.0}");
    let _ = writeln!(toml, "recommended_workers = {rec_workers}");
    let _ = writeln!(toml, "recommended_decomp_threads = {rec_decomp}");
    let _ = writeln!(toml, "recommended_writers = {rec_writers}");
    for (w, rate) in &kernels {
        let _ = writeln!(toml, "baq_kernel_{}_reads_per_thread_s = {rate:.0}", w.0);
    }
    let _ = writeln!(toml, "recommended_baq_simd = \"{}\"", simd_name(best_kernel));
    if let Some((dev, cpu, gpu)) = &gpu_probe {
        let _ = writeln!(toml, "gpu_device = \"{}\"", dev.replace('"', "'"));
        let _ = writeln!(toml, "gpu_baq_reads_per_s = {gpu:.0}");
        let _ = writeln!(toml, "cpu_baq_reads_per_s = {cpu:.0}");
        let _ = writeln!(toml, "gpu_baq_faster = {}", gpu > cpu);
    }

    if let Some(j) = &a.json {
        std::fs::write(j, toml.clone() + "\n# report\n" + &rpt.lines().map(|l| format!("# {l}\n")).collect::<String>())
            .with_context(|| format!("write {}", j.display()))?;
        eprintln!("Wrote {}", j.display());
    }
    if a.save_profile {
        let p = profile_path().ok_or_else(|| anyhow::anyhow!("no HOME for profile"))?;
        if let Some(d) = p.parent() {
            std::fs::create_dir_all(d)?;
        }
        std::fs::write(&p, &toml).with_context(|| format!("write {}", p.display()))?;
        eprintln!("Saved machine profile: {} (extract uses it for default -t / writers on this host)", p.display());
    }
    Ok(())
}

fn chrono_like_now() -> String {
    std::process::Command::new("date")
        .arg("-Is")
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn synthetic_bam_extracts_calls() {
        let dir = std::env::temp_dir().join("taps_rs_tune_syn");
        let _ = std::fs::remove_dir_all(&dir);
        let syn = make_synthetic(&dir, 1, 200_000, 2_000, 20).unwrap();
        assert_eq!(syn.reads, 4_000);
        let mut r = bam::Reader::from_path(&syn.bam).unwrap();
        assert_eq!(r.records().count(), 4_000);
        let fai = crate::fasta::FastFaiReader::open(&syn.fasta).unwrap();
        let w = synthetic_windows(&fai, &["tchr1"], 200_000, 1);
        let (res, _, _) = run_extract(&syn.bam, &syn.fasta, &w, &tune_params(true, false), 1, 1).unwrap();
        assert!(res.stats.cg_c + res.stats.cg_t > 1_000, "CG calls");
        assert!(res.stats.cg_pct() > 50.0 && res.stats.ch_pct() < 5.0);
        assert_eq!(res.n_cells(), 20);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
