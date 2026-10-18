//! taps-sc-extract-rs: Rust core for TAPS single-cell methylation extraction.

#![allow(dead_code)]

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod accumulate;
mod autotune;
mod baq_gpu;
mod baq_simd;
mod barcode;
mod calling;
mod context;
mod extract;
mod extract_noodles;
mod fasta;
mod h5_out;
mod input;
mod parallel;
mod profile;
mod readcall;
mod shard_io;
mod tune;
mod window;

use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;
use std::time::Instant;

use crate::extract::{estimate_bam_traffic, BamTraffic, BaqMode, ContextSel, ExtractParams};

/// "taps", "htslib" or "off".
fn baq_label(p: &ExtractParams) -> &'static str {
    if p.compute_baq {
        p.baq_mode.label()
    } else {
        "off"
    }
}

/// `--baq-simd auto` takes the kernel `--test --save-profile` measured fastest
/// on this host, when it is still supported; otherwise the widest one.
fn resolve_baq_width(choice: baq_simd::BaqSimd, saved: Option<&tune::SavedProfile>) -> Result<baq_simd::Width> {
    use baq_simd::{BaqSimd, Width};
    if choice == BaqSimd::Auto {
        let from_profile = match saved.and_then(|p| p.baq_simd.as_deref()) {
            Some("avx512") => Some(Width::AVX512),
            Some("avx2") => Some(Width::AVX2),
            Some("off") => Some(Width::HTSLIB),
            _ => None,
        };
        if let Some(w) = from_profile.filter(|w| w.supported()) {
            return Ok(w);
        }
    }
    Width::resolve(choice)
}

/// Kernel label, or "-" when BAQ is off.
fn baq_kernel_label(p: &ExtractParams) -> &'static str {
    if p.compute_baq {
        p.baq_width.label()
    } else {
        "-"
    }
}
use crate::extract_noodles::extract_regions_stats_noodles;
use crate::fasta::FastFaiReader;
use crate::h5_out::{
    assemble_hdf5_from_memory, assemble_hdf5_from_temp, read_shard_temp_chunks,
    written_output_bytes, H5Compression, WriterMode,
};
use crate::shard_io::ChunkFile;
use crate::parallel::{
    extract_memory_parallel, extract_stats_parallel, extract_stream_parallel, install_ctrlc,
    new_cancel_flag,
};
use crate::window::{plan_windows, CANONICAL_CONTIGS};

#[derive(Parser)]
#[command(name = "taps-sc-extract-rs", version, about = "Rust TAPS extractor (WIP)")]
struct Cli {
    #[command(subcommand)]
    cmd: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Stats-only extraction for pileup-parity spikes (PR0/PR3).
    Stats {
        #[arg(short = 'b', long)]
        bam: PathBuf,
        #[arg(short = 'f', long)]
        fasta: PathBuf,
        /// Comma-separated contigs (default: canonical mm10).
        #[arg(short = 'c', long)]
        chroms: Option<String>,
        /// Restrict to a single 0-based window start (bp) on the first contig.
        #[arg(long)]
        start: Option<u64>,
        /// Window end (bp, exclusive). Requires --start.
        #[arg(long)]
        end: Option<u64>,
        /// Genomic window size in megabases (Rust default: 100).
        #[arg(long, default_value_t = 100)]
        chunk_size_mb: u32,
        /// Window workers (Rayon). 0 = auto (`min(nproc, 32)`).
        #[arg(short = 't', long = "workers", default_value_t = 0)]
        workers: usize,
        /// BGZF threads per worker. Omitted: more threads when fewer workers.
        #[arg(long)]
        decomp_threads: Option<usize>,
        #[arg(long, default_value_t = 20)]
        min_baseq: u8,
        #[arg(long, default_value_t = 0)]
        min_mapq: u8,
        #[arg(long, default_value_t = 250)]
        max_depth: u32,
        /// Do not clip overlapping mate bases (`bam_mplp_init_overlaps`).
        #[arg(long, default_value_t = false)]
        no_overlap_clip: bool,
        /// Keep paired reads that are not a proper pair.
        #[arg(long, default_value_t = false)]
        no_ignore_orphans: bool,
        /// Disable BAQ (`sam_prob_realn`). Default is on (Python's compute_baq
        /// flag is currently a no-op; we do not copy that).
        #[arg(long, default_value_t = false)]
        no_baq: bool,
        /// BAQ mode. The `taps*` modes mask the methylation-convertible base (C
        /// on OT reads, G on OB reads) so conversions are not penalised.
        /// `taps` skips reads that already match that reference exactly;
        /// `taps-full` runs every read; `taps-partial` only reads with an indel
        /// or clip. `htslib` uses the plain reference (under-calls methylation).
        #[arg(long, value_enum, default_value_t = BaqMode::Taps)]
        baq: BaqMode,
        /// BAQ kernel: auto (widest this CPU runs), avx512 (8 reads per
        /// vector), avx2 (4), or off (htslib, one read at a time). All give
        /// identical results (`baq-selftest`).
        #[arg(long, value_enum, default_value_t = baq_simd::BaqSimd::Auto)]
        baq_simd: baq_simd::BaqSimd,
        /// Methylation contexts to extract: both (default), CG or CH. With one
        /// context the other's pileup columns are skipped and its HDF5 group
        /// (/CG or /CH) is not written.
        #[arg(long, value_enum, default_value_t = ContextSel::Both, ignore_case = true)]
        context: ContextSel,
        /// How calls are made: `reads` (default; one read at a time, replaying
        /// htslib's pileup depth cap and mate-overlap rules) or `pileup`
        /// (htslib column pileup). Both give the same calls.
        #[arg(long, value_enum, default_value_t = extract::Caller::Reads)]
        caller: extract::Caller,
        /// Drop this many 5' sequenced bases of read 1 before calling.
        #[arg(long, default_value_t = 0)]
        trim_r1: u32,
        /// Drop this many 5' sequenced bases of read 2 before calling.
        #[arg(long, default_value_t = 10)]
        trim_r2: u32,
        #[arg(short = 'w', long)]
        whitelist: Option<PathBuf>,
        /// Complexity file (rank, barcode, total, unique reads, %); keeps barcodes
        /// whose unique-read count (column 4) is within --min/--max-complexity.
        #[arg(long)]
        complexity: Option<PathBuf>,
        /// Keep cells with at least this many unique reads (needs --complexity).
        #[arg(long)]
        min_complexity: Option<f64>,
        /// Keep cells with at most this many unique reads (needs --complexity).
        #[arg(long)]
        max_complexity: Option<f64>,
        /// Pileup engine: rust-htslib (htslib) or noodles (pure Rust CIGAR walk).
        #[arg(long, default_value = "htslib")]
        pileup: String,
        /// Write per-barcode CG/CH site TSV (barcode, ctx, pos, t, c) for parity diffs.
        #[arg(long)]
        dump_sites: Option<PathBuf>,
    },
    /// Single-thread window extract to compact shard temp files (PR3/PR4).
    Extract {
        #[arg(short = 'b', long)]
        bam: PathBuf,
        #[arg(short = 'f', long)]
        fasta: PathBuf,
        /// Output `.h5` file (`--shards 1`) or directory of shards + `master.h5`.
        #[arg(short = 'o', long)]
        out: PathBuf,
        #[arg(short = 'c', long)]
        chroms: Option<String>,
        #[arg(long)]
        start: Option<u64>,
        #[arg(long)]
        end: Option<u64>,
        /// Genomic window size in megabases. 0 = auto: about 6 windows per worker,
        /// clamped to 2-100 Mb.
        #[arg(long, default_value_t = 0)]
        chunk_size_mb: u32,
        /// Number of HDF5 shard files (0 = auto based on cell count: 1/8/16/32).
        #[arg(long, default_value_t = 0)]
        shards: usize,
        /// Window workers (Rayon). 0 = auto sized from CPU and memory budget.
        #[arg(short = 't', long = "workers", default_value_t = 0)]
        workers: usize,
        /// BGZF threads per worker. Omitted: more threads when fewer workers.
        #[arg(long)]
        decomp_threads: Option<usize>,
        /// Keep stream-mode temp chunk directory.
        #[arg(long, default_value_t = false)]
        keep_temp: bool,
        /// stream = temp chunks then HDF5; memory = keep in RAM; auto = heuristic from budget & cells.
        #[arg(long, default_value = "auto")]
        memory_mode: String,
        /// Memory budget in GB for the whole process tree (overrides --max-memory-fraction).
        #[arg(long)]
        max_memory_gb: Option<f64>,
        /// Share of MemAvailable the run may use (default 0.6).
        #[arg(long)]
        max_memory_fraction: Option<f64>,
        /// Optional expected cell count (auto-detected from whitelist if omitted).
        #[arg(long)]
        expected_cells: Option<usize>,
        /// HDF5 compression (default: blosc, LZ4 level 9: as fast as lzf on the
        /// test data with 2.7x smaller files; R needs rhdf5filters for blosc and lzf).
        #[arg(long, default_value = "blosc")]
        compression: String,
        /// Blosc internal threads per shard writer (blosc / blosc-zstd only).
        /// Default 1: more threads measured no gain (2026-10-15), because most
        /// cell datasets are a single chunk too small for Blosc to split.
        #[arg(long)]
        compression_threads: Option<usize>,
        /// Maximum concurrent shard writers (0 = auto sized from shards, cores, RAM).
        #[arg(long, default_value_t = 0)]
        max_writer_threads: usize,
        /// Shard writers: process (one child per shard, each with its own HDF5
        /// lock), thread (in-process pool sharing one lock), or auto (process
        /// when there is more than one shard).
        #[arg(long, default_value = "auto")]
        writer_mode: String,
        #[arg(long, default_value_t = 20)]
        min_baseq: u8,
        #[arg(long, default_value_t = 0)]
        min_mapq: u8,
        #[arg(long, default_value_t = 250)]
        max_depth: u32,
        #[arg(long, default_value_t = false)]
        no_overlap_clip: bool,
        #[arg(long, default_value_t = false)]
        no_ignore_orphans: bool,
        /// Disable BAQ. Default is on (do not copy Python's current no-op).
        #[arg(long, default_value_t = false)]
        no_baq: bool,
        /// BAQ mode: taps (default), taps-full, taps-partial or htslib.
        #[arg(long, value_enum, default_value_t = BaqMode::Taps)]
        baq: BaqMode,
        /// BAQ kernel: auto (widest this CPU runs), avx512 (8 reads per
        /// vector), avx2 (4), or off (htslib, one read at a time). All give
        /// identical results (`baq-selftest`).
        #[arg(long, value_enum, default_value_t = baq_simd::BaqSimd::Auto)]
        baq_simd: baq_simd::BaqSimd,
        /// Where the BAQ HMM runs: cpu (default, SIMD kernels) or gpu
        /// (identical results; see --gpu-backend). Reads too long or with a
        /// band too wide for the GPU kernel stay on the CPU.
        #[arg(long, value_enum, default_value_t = baq_gpu::BaqDevice::Cpu)]
        baq_device: baq_gpu::BaqDevice,
        /// GPU backend for --baq-device gpu: auto (CUDA, then HIP, then
        /// OpenCL), opencl (Intel Arc/Data Center, any fp64 device), cuda, hip.
        #[arg(long, value_enum, default_value_t = baq_gpu::GpuBackend::Auto)]
        gpu_backend: baq_gpu::GpuBackend,
        /// Use the first GPU whose name contains this text.
        #[arg(long, default_value = "")]
        gpu_device: String,
        /// GPU memory for HMM matrices per launch, in MB.
        #[arg(long, default_value_t = 1024)]
        gpu_mem_mb: usize,
        /// Methylation contexts to extract: both (default), CG or CH. With one
        /// context the other's pileup columns are skipped and its HDF5 group
        /// (/CG or /CH) is not written.
        #[arg(long, value_enum, default_value_t = ContextSel::Both, ignore_case = true)]
        context: ContextSel,
        /// How calls are made: `reads` (default; one read at a time, replaying
        /// htslib's pileup depth cap and mate-overlap rules) or `pileup`
        /// (htslib column pileup). Both give the same calls.
        #[arg(long, value_enum, default_value_t = extract::Caller::Reads)]
        caller: extract::Caller,
        /// Drop this many 5' sequenced bases of read 1 before calling.
        #[arg(long, default_value_t = 0)]
        trim_r1: u32,
        /// Drop this many 5' sequenced bases of read 2 before calling.
        #[arg(long, default_value_t = 10)]
        trim_r2: u32,
        #[arg(short = 'w', long)]
        whitelist: Option<PathBuf>,
        /// Complexity file (rank, barcode, total, unique reads, %); keeps barcodes
        /// whose unique-read count (column 4) is within --min/--max-complexity.
        #[arg(long)]
        complexity: Option<PathBuf>,
        /// Keep cells with at least this many unique reads (needs --complexity).
        #[arg(long)]
        min_complexity: Option<f64>,
        /// Keep cells with at most this many unique reads (needs --complexity).
        #[arg(long)]
        max_complexity: Option<f64>,
        /// Write per-phase and per-window timings (and writer reports) as JSON.
        #[arg(long)]
        profile_json: Option<PathBuf>,
    },
    /// Measure this machine (and optionally a dataset) and recommend flags.
    Tune {
        /// Where the synthetic BAM, temp files and storage probe go.
        #[arg(long, default_value = ".")]
        temp_dir: PathBuf,
        /// Output filesystem to probe (HDF5 writes land here).
        #[arg(long, default_value = ".")]
        out_dir: PathBuf,
        /// Input filesystem to probe (where BAMs live). Defaults to the BAM's directory.
        #[arg(long)]
        input_dir: Option<PathBuf>,
        /// Optional dataset to sample and predict for.
        #[arg(short = 'b', long)]
        bam: Option<PathBuf>,
        #[arg(short = 'f', long, requires = "bam")]
        fasta: Option<PathBuf>,
        /// Smaller probes (~1 minute).
        #[arg(long, default_value_t = false)]
        quick: bool,
        /// Save the machine profile for `extract` defaults on this host.
        #[arg(long, default_value_t = false)]
        save_profile: bool,
        /// Also write the profile and report to this file.
        #[arg(long)]
        report: Option<PathBuf>,
    },
    /// Check the SIMD BAQ kernels against htslib, read by read (QUAL must be
    /// byte-identical). Synthetic reads always; add -b/-f for real reads.
    BaqSelftest {
        #[arg(short = 'b', long)]
        bam: Option<PathBuf>,
        #[arg(short = 'f', long, requires = "bam")]
        fasta: Option<PathBuf>,
        /// Contig to read from the BAM.
        #[arg(short = 'c', long, default_value = "chr19")]
        contig: String,
        #[arg(long, default_value_t = 2_000_000)]
        max_reads: usize,
        #[arg(long, default_value_t = 200_000)]
        synthetic_reads: usize,
        /// Test the GPU path (HMM on the GPU) instead of the CPU kernels.
        #[arg(long, default_value_t = false)]
        gpu: bool,
        #[arg(long, value_enum, default_value_t = baq_gpu::GpuBackend::Auto)]
        gpu_backend: baq_gpu::GpuBackend,
        #[arg(long, default_value = "")]
        gpu_device: String,
    },
    /// Internal: write one shard file from window chunks (spawned by `extract`).
    #[command(hide = true)]
    WriteShard {
        #[arg(long)]
        shard: usize,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, default_value = "lzf")]
        compression: String,
        /// Blosc internal threads (0 = Blosc default).
        #[arg(long, default_value_t = 0)]
        compression_threads: usize,
        /// Stream-mode temp bucket directory with chunk_*.bin files (repeatable:
        /// one per bucket this shard writes).
        #[arg(long, conflicts_with = "stdin")]
        chunk_dir: Vec<PathBuf>,
        /// Read TAPSCK01 chunk records from stdin until EOF.
        #[arg(long, default_value_t = false)]
        stdin: bool,
    },
}

fn main() -> Result<()> {
    let cancel = new_cancel_flag();
    let _ = install_ctrlc(cancel.clone());
    let cli = Cli::parse();
    match cli.cmd {
        Commands::Tune {
            temp_dir,
            out_dir,
            input_dir,
            bam,
            fasta,
            quick,
            save_profile,
            report,
        } => {
            let input_dir = input_dir
                .or_else(|| bam.as_ref().and_then(|b| b.parent().map(|p| p.to_path_buf())))
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or_else(|| PathBuf::from("."));
            let prepared = match (&bam, &fasta) {
                (Some(b), Some(f)) => {
                    let nproc = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
                    Some(input::prepare(b, f, nproc)?)
                }
                _ => None,
            };
            let bam = prepared.as_ref().map(|p| p.path.clone()).or(bam);
            tune::run_tune(&tune::TuneArgs {
                temp_dir,
                out_dir,
                input_dir,
                bam,
                fasta,
                quick,
                save_profile,
                json: report,
            })?;
        }
        Commands::BaqSelftest { bam, fasta, contig, max_reads, synthetic_reads, gpu, gpu_backend, gpu_device } => {
            use baq_simd::Width;
            // htslib: the masked fallback (one copied span per read) must match too.
            let mut widths = vec![Width::HTSLIB, Width::AVX2, Width::AVX2_X8, Width::AVX512];
            widths.retain(|w| w.supported());
            if gpu {
                let svc = baq_gpu::start(gpu_backend, &gpu_device, 1024)?;
                eprintln!("BAQ self-test on GPU {} ({:?}); CPU kernel {} for reads the GPU does not take", svc.device, svc.backend, Width::cpu_best().label());
                widths = vec![Width::cpu_best()];
            }
            eprintln!(
                "BAQ self-test: CPU kernel {}; testing {}",
                Width::cpu_best().label(),
                if widths.is_empty() { "none (no AVX2)".to_string() } else { widths.iter().map(|w| w.label()).collect::<Vec<_>>().join(", ") }
            );
            let mut bad = 0usize;
            for w in &widths {
                let t = Instant::now();
                let r = baq_simd::selftest_synthetic(*w, synthetic_reads, 0x5eed, gpu)?;
                eprintln!("  synthetic  {:32} {:>9} reads, {:>9} via {}, {} differ ({:.1} s)", if gpu { "gpu (+cpu)" } else { w.label() }, r.reads, r.simd_reads, if gpu { "GPU" } else { "SIMD" }, r.mismatched, t.elapsed().as_secs_f64());
                bad += r.mismatched;
                if let (Some(b), Some(f)) = (&bam, &fasta) {
                    let t = Instant::now();
                    let r = baq_simd::selftest_bam(b, f, &contig, max_reads, *w, gpu)?;
                    eprintln!("  {contig:<10} {:32} {:>9} reads, {:>9} via {}, {} differ ({:.1} s)", if gpu { "gpu (+cpu)" } else { w.label() }, r.reads, r.simd_reads, if gpu { "GPU" } else { "SIMD" }, r.mismatched, t.elapsed().as_secs_f64());
                    bad += r.mismatched;
                }
            }
            anyhow::ensure!(bad == 0, "SIMD BAQ differs from htslib on {bad} read(s); use --baq-simd off");
            eprintln!("BAQ self-test passed");
        }
        Commands::WriteShard {
            shard,
            out,
            compression,
            compression_threads,
            chunk_dir,
            stdin,
        } => {
            std::env::set_var("HDF5_USE_FILE_LOCKING", "FALSE");
            let codec = H5Compression::parse(&compression)?;
            if codec.is_blosc() {
                crate::h5_out::apply_blosc_threads(compression_threads);
            }
            let input = || -> Result<Vec<shard_io::ChunkFile>> {
                match (chunk_dir.is_empty(), stdin) {
                    (false, false) => read_shard_temp_chunks(&chunk_dir),
                    (true, true) => {
                        let mut r = std::io::BufReader::with_capacity(1 << 20, std::io::stdin().lock());
                        shard_io::read_all_chunks_from(&mut r)
                    }
                    _ => anyhow::bail!("write-shard needs exactly one of --chunk-dir or --stdin"),
                }
            };
            let mut rep = crate::h5_out::write_shard_timed(shard, &out, input, codec)?;
            rep.peak_rss_mb = autotune::process_peak_rss_mb();
            eprintln!("{}", rep.child_line(&out));
        }
        Commands::Stats {
            bam,
            fasta,
            chroms,
            start,
            end,
            chunk_size_mb,
            workers,
            decomp_threads,
            min_baseq,
            min_mapq,
            max_depth,
            no_overlap_clip,
            no_ignore_orphans,
            no_baq,
            baq,
            baq_simd,
            context,
            caller,
            trim_r1,
            trim_r2,
            whitelist,
            complexity,
            min_complexity,
            max_complexity,
            pileup,
            dump_sites,
        } => {
            let fai = FastFaiReader::open(&fasta)?;
            let nproc = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
            let prepared = input::prepare(&bam, &fasta, nproc)?;
            let bam = prepared.path.clone();
            let contig_list: Vec<&str> = match &chroms {
                Some(s) => s.split(',').map(|x| x.trim()).filter(|x| !x.is_empty()).collect(),
                None => CANONICAL_CONTIGS.to_vec(),
            };
            let mut windows = plan_windows(&fai, &contig_list, u64::from(chunk_size_mb) * 1_000_000);
            if let (Some(s), Some(e)) = (start, end) {
                let contig = contig_list
                    .first()
                    .copied()
                    .ok_or_else(|| anyhow::anyhow!("no contigs"))?;
                windows = vec![window::Window {
                    chunk_id: 0,
                    contig: contig.to_string(),
                    start: s,
                    end: e,
                }];
            }
            let wl = barcode::load_cell_filter(
                whitelist.as_deref(),
                complexity.as_deref(),
                min_complexity,
                max_complexity,
            )?;
            let ignore_overlaps = !no_overlap_clip;
            let ignore_orphans = !no_ignore_orphans;
            let compute_baq = !no_baq;
            let n_workers = autotune::resolve_workers(workers);
            let decomp = decomp_threads.unwrap_or_else(|| autotune::default_decomp_threads(n_workers));
            let params = ExtractParams {
                min_baseq,
                min_mapq,
                max_depth,
                ignore_overlaps,
                ignore_orphans,
                compute_baq,
                baq_mode: baq,
                baq_width: baq_simd::Width::resolve(baq_simd)?,
                baq_device: baq_gpu::BaqDevice::Cpu,
                context,
                caller,
                decomp_threads: decomp,
                accumulate: dump_sites.is_some(),
                trim_r1,
                trim_r2,
                profile: false,
            };
            eprintln!(
                "taps-sc-extract-rs stats: pileup={pileup} workers={n_workers} decomp={decomp} baq={} overlaps={ignore_overlaps} orphans={ignore_orphans} trim_r1={trim_r1} trim_r2={trim_r2} {} window(s), contig(s)={:?}, chunk_size_mb={chunk_size_mb}",
                baq_label(&params),
                windows.len(),
                contig_list
            );
            let t0 = std::time::Instant::now();
            let r = match pileup.as_str() {
                "noodles" => extract_regions_stats_noodles(&bam, &fasta, &windows, wl.as_ref(), &params)?,
                "htslib" | "rust-htslib" => extract_stats_parallel(
                    &bam,
                    &fasta,
                    &windows,
                    wl.as_ref(),
                    &params,
                    n_workers,
                    &cancel,
                )?,
                other => anyhow::bail!("unknown --pileup {other} (use htslib or noodles)"),
            };
            let dt = t0.elapsed();
            if let Some(p) = dump_sites {
                accumulate::write_sites_tsv(&p, &r.intern, &r.cells)?;
                eprintln!("wrote site TSV {}", p.display());
            }
            println!(
                "cells={} CG={:.4}% ({}/{}) CH={:.4}% ({}/{}) CpG={:.4}% CHG={:.4}% CHH={:.4}% elapsed_s={:.2}",
                r.n_cells(),
                r.stats.cg_pct(),
                r.stats.cg_c,
                r.stats.cg_c + r.stats.cg_t,
                r.stats.ch_pct(),
                r.stats.ch_c,
                r.stats.ch_c + r.stats.ch_t,
                {
                    let t = r.stats.cpg_c + r.stats.cpg_t;
                    if t == 0 { 0.0 } else { 100.0 * r.stats.cpg_c as f64 / t as f64 }
                },
                {
                    let t = r.stats.chg_c + r.stats.chg_t;
                    if t == 0 { 0.0 } else { 100.0 * r.stats.chg_c as f64 / t as f64 }
                },
                {
                    let t = r.stats.chh_c + r.stats.chh_t;
                    if t == 0 { 0.0 } else { 100.0 * r.stats.chh_c as f64 / t as f64 }
                },
                dt.as_secs_f64()
            );
        }
        Commands::Extract {
            bam,
            fasta,
            out,
            chroms,
            start,
            end,
            chunk_size_mb,
            shards,
            workers,
            decomp_threads,
            keep_temp,
            memory_mode,
            max_memory_gb,
            max_memory_fraction,
            expected_cells,
            compression,
            compression_threads,
            max_writer_threads,
            writer_mode,
            min_baseq,
            min_mapq,
            max_depth,
            no_overlap_clip,
            no_ignore_orphans,
            no_baq,
            baq,
            baq_simd,
            baq_device,
            gpu_backend,
            gpu_device,
            gpu_mem_mb,
            context,
            caller,
            trim_r1,
            trim_r2,
            whitelist,
            complexity,
            min_complexity,
            max_complexity,
            profile_json,
        } => {
            let fai = FastFaiReader::open(&fasta)?;
            let nproc = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
            let prepared = input::prepare(&bam, &fasta, nproc)?;
            let bam = prepared.path.clone();
            let contig_list: Vec<&str> = match &chroms {
                Some(s) => s.split(',').map(|x| x.trim()).filter(|x| !x.is_empty()).collect(),
                None => CANONICAL_CONTIGS.to_vec(),
            };
            let wl = barcode::load_cell_filter(
                whitelist.as_deref(),
                complexity.as_deref(),
                min_complexity,
                max_complexity,
            )?;

            // PR7 Auto-tuning & Phase 1 Writer Auto-scaling
            let (budget_gb, budget_src) = autotune::system_memory_budget_gb(max_memory_gb, max_memory_fraction);
            let cell_source = autotune::CellCountSource::resolve(
                expected_cells,
                wl.as_ref().map(|w| w.len()),
                autotune::DEFAULT_ESTIMATED_CELLS,
            );
            let saved = tune::load_profile_for_host();
            let (n_workers, workers_src) = match (&saved, workers) {
                (Some((p, _)), 0) if p.workers.is_some() => (p.workers.unwrap(), "machine profile"),
                _ => autotune::resolve_workers_with_budget(workers, budget_gb),
            };
            let selected_bp: u64 = contig_list
                .iter()
                .filter_map(|c| fai.reference_length(c))
                .sum();
            let (chunk_bp, chunk_src) = if chunk_size_mb == 0 {
                (autotune::auto_chunk_size_bp(selected_bp, n_workers), "auto")
            } else {
                (u64::from(chunk_size_mb) * 1_000_000, "explicit")
            };
            let mut windows = plan_windows(&fai, &contig_list, chunk_bp);
            if let (Some(s), Some(e)) = (start, end) {
                let contig = contig_list
                    .first()
                    .copied()
                    .ok_or_else(|| anyhow::anyhow!("no contigs"))?;
                windows = vec![window::Window {
                    chunk_id: 0,
                    contig: contig.to_string(),
                    start: s,
                    end: e,
                }];
            }
            let decomp = decomp_threads
                .or_else(|| saved.as_ref().and_then(|(p, _)| p.decomp_threads))
                .unwrap_or_else(|| autotune::default_decomp_threads(n_workers));
            let (n_shards, shards_src) = autotune::resolve_shards(shards, cell_source.count());
            let writer_mode = WriterMode::parse(&writer_mode, n_shards)?;
            let traffic = match estimate_bam_traffic(&bam, &contig_list) {
                Ok(t) => t,
                Err(e) => {
                    eprintln!("warning: BAM read estimate unavailable: {e}");
                    let bytes = std::fs::metadata(&bam).map(|m| m.len()).unwrap_or(0);
                    BamTraffic {
                        mapped_reads: 0,
                        estimated_bytes: bytes,
                    }
                }
            };
            let mem_plan = autotune::plan_memory(&autotune::MemoryInputs {
                mode: &memory_mode,
                budget_gb,
                mapped_reads: traffic.mapped_reads,
                site_share: context.site_share(),
                cell_count: cell_source.count(),
                genome_mb: selected_bp / 1_000_000,
                n_shards,
                n_workers,
                cli_writers: max_writer_threads,
                profile_writers: if writer_mode == WriterMode::Process {
                    saved.as_ref().and_then(|(p, _)| p.writers)
                } else {
                    None
                },
                in_process: writer_mode == WriterMode::Thread,
            });
            let (writer_threads, writer_src) = (mem_plan.writers, mem_plan.writers_src);
            let (resolved_memory_mode, mode_src) = (mem_plan.mode.clone(), mem_plan.mode_src);
            let compression_name = compression.clone();
            let compression = H5Compression::parse(&compression)?;
            let blosc_threads = if compression.is_blosc() {
                compression_threads.unwrap_or(1).max(1)
            } else {
                0
            };
            crate::h5_out::BLOSC_THREADS.store(blosc_threads, std::sync::atomic::Ordering::Relaxed);
            if writer_mode == WriterMode::Thread || n_shards == 1 {
                crate::h5_out::apply_blosc_threads(blosc_threads);
            }
            let params = ExtractParams {
                min_baseq,
                min_mapq,
                max_depth,
                ignore_overlaps: !no_overlap_clip,
                ignore_orphans: !no_ignore_orphans,
                compute_baq: !no_baq,
                baq_mode: baq,
                baq_width: resolve_baq_width(baq_simd, saved.as_ref().map(|(p, _)| p))?,
                baq_device,
                context,
                caller,
                decomp_threads: decomp,
                accumulate: true,
                trim_r1,
                trim_r2,
                profile: profile_json.is_some(),
            };
            eprintln!("======================================================================");
            eprintln!("taps-sc-extract-rs v{}", env!("CARGO_PKG_VERSION"));
            eprintln!("Configuration:");
            eprintln!("  Memory budget:   {:.1} GB ({})", budget_gb, budget_src);
            if let Some((_, path)) = &saved {
                eprintln!("  Machine profile: {}", path.display());
            }
            eprintln!("  Expected cells:  {} ({})", cell_source.count(), cell_source.label());
            eprintln!("  Workers:         {} ({})", n_workers, workers_src);
            eprintln!("  Decomp threads:  {} per worker", decomp);
            eprintln!("  Shards:          {} ({})", n_shards, shards_src);
            eprintln!("  Writers:         {} {} ({})", writer_threads, writer_mode.label(), writer_src);
            eprintln!("  Memory mode:     {} ({})", resolved_memory_mode, mode_src);
            eprintln!(
                "  Memory plan:     ~{:.0}M sites; extraction ~{:.1} GB, writers {} x ~{:.1} GB; est. peak ~{:.1} GB of {:.1} GB",
                mem_plan.est_sites / 1e6,
                mem_plan.extract_gb,
                writer_threads,
                mem_plan.writer_each_gb,
                mem_plan.peak_gb,
                budget_gb
            );
            if mem_plan.peak_gb > budget_gb {
                eprintln!("  warning: estimated peak exceeds the memory budget; lower --max-writer-threads or use --memory-mode stream");
            }
            eprintln!("  Chunk size:      {:.1} Mb ({chunk_src})", chunk_bp as f64 / 1e6);
            if blosc_threads > 0 {
                eprintln!("  Compression:     {compression_name} ({blosc_threads} Blosc thread(s) per writer)");
            } else {
                eprintln!("  Compression:     {compression_name}");
            }
            let gpu_label = if params.compute_baq && params.baq_device == baq_gpu::BaqDevice::Gpu {
                let svc = baq_gpu::start(gpu_backend, &gpu_device, gpu_mem_mb)?;
                Some(format!("GPU {} ({:?})", svc.device, svc.backend))
            } else {
                None
            };
            match &gpu_label {
                Some(g) => eprintln!("  BAQ:             {} (HMM on {g}; CPU kernel {} for the rest)", baq_label(&params), baq_kernel_label(&params)),
                None => eprintln!("  BAQ:             {} (kernel: {})", baq_label(&params), baq_kernel_label(&params)),
            }
            eprintln!("  Contexts:        {}", params.context.label());
            eprintln!("  5' trim:         R1={} bp, R2={} bp", params.trim_r1, params.trim_r2);
            eprintln!("  Windows:         {} window(s) across {:?}", windows.len(), contig_list);
            eprintln!("  Output:          {}", out.display());
            eprintln!("======================================================================");
            let n_buckets = crate::h5_out::bucket_count(n_shards);
            let cpu_start = profile::cpu_times();
            let t_process = Instant::now();
            let (r, pending) = match resolved_memory_mode.as_str() {
                "memory" => {
                    let (r, mem) = extract_memory_parallel(
                        &bam,
                        &fasta,
                        &windows,
                        wl.as_ref(),
                        &params,
                        n_buckets,
                        n_workers,
                        &cancel,
                    )?;
                    (r, PendingOutput::Memory(mem))
                }
                "stream" => {
                    let tmp = tempfile::Builder::new()
                        .prefix("taps-rs-")
                        .tempdir()
                        .map_err(|e| anyhow::anyhow!("tempdir: {e}"))?;
                    let r = extract_stream_parallel(
                        &bam,
                        &fasta,
                        &windows,
                        wl.as_ref(),
                        &params,
                        tmp.path(),
                        n_buckets,
                        n_workers,
                        &cancel,
                    )?;
                    (r, PendingOutput::Stream(tmp))
                }
                other => anyhow::bail!("unknown --memory-mode {other} (use stream or memory)"),
            };
            let process_s = t_process.elapsed().as_secs_f64();
            let cpu_after_process = profile::cpu_times();
            eprintln!(
                "Extraction workers finished in {process_s:.1}s. Writing {n_shards} shard(s) with {writer_threads} {} writer(s).",
                writer_mode.label()
            );
            let t_write = Instant::now();
            let crate::h5_out::HdfOutput { path: h5_path, reports: writer_reports, master_s } = match pending {
                PendingOutput::Memory(mem) => {
                    let groups = crate::h5_out::pack_buckets(&crate::h5_out::bucket_sites(&mem), n_shards);
                    assemble_hdf5_from_memory(
                        crate::h5_out::group_buckets(mem, &groups),
                        &out,
                        compression,
                        &compression_name,
                        writer_threads,
                        writer_mode,
                    )?
                }
                PendingOutput::Stream(tmp) => {
                    let groups = crate::h5_out::pack_buckets(
                        &crate::h5_out::bucket_bytes(tmp.path(), n_buckets),
                        n_shards,
                    );
                    let res = assemble_hdf5_from_temp(
                        tmp.path(),
                        &out,
                        &groups,
                        compression,
                        &compression_name,
                        writer_threads,
                        writer_mode,
                    )?;
                    if keep_temp {
                        let kept = tmp.keep();
                        eprintln!("kept temp {}", kept.display());
                    }
                    res
                }
            };
            let write_s = t_write.elapsed().as_secs_f64();
            let cpu_end = profile::cpu_times();
            let out_bytes = written_output_bytes(&h5_path);
            let total_s = process_s + write_s;
            print_extraction_summary(
                &r,
                &traffic,
                total_s,
                process_s,
                write_s,
                out_bytes,
                n_shards,
                resolved_memory_mode == "stream",
                crate::autotune::process_peak_rss_mb(),
                crate::h5_out::concurrent_writer_peak_mb(&writer_reports, writer_threads),
                writer_reports.iter().filter(|w| w.child).count(),
                params.context,
            );
            if let Some(line) = crate::h5_out::writer_summary_line(&writer_reports, master_s) {
                eprintln!("{line}");
            }
            if wl.is_some() && r.n_cells() == 0 {
                eprintln!(
                    "WARNING: no reads matched the barcode filter (--whitelist / --complexity). \
                     Check that its barcodes use the same format as the QNAME prefix before ':'."
                );
            }
            if let Some(path) = &profile_json {
                let rp = profile::RunProfile {
                    config: vec![
                        ("version", env!("CARGO_PKG_VERSION").to_string()),
                        ("workers", n_workers.to_string()),
                        ("decomp_threads", decomp.to_string()),
                        ("chunk_bp", chunk_bp.to_string()),
                        ("windows", windows.len().to_string()),
                        ("shards", n_shards.to_string()),
                        ("writer_mode", writer_mode.label().to_string()),
                        ("writers", writer_threads.to_string()),
                        ("memory_mode", resolved_memory_mode.clone()),
                        ("compression", compression_name.clone()),
                        ("baq", baq_label(&params).to_string()),
                        ("baq_kernel", baq_kernel_label(&params).to_string()),
                        ("baq_device", gpu_label.clone().unwrap_or_else(|| "cpu".to_string())),
                        ("context", params.context.label().to_string()),
                        ("bam", bam.display().to_string()),
                        ("nproc", std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1).to_string()),
                    ],
                    process_s,
                    write_s,
                    cpu_start,
                    cpu_after_process,
                    cpu_end,
                    workers: n_workers,
                    writers: writer_threads,
                    result: &r,
                    windows: &windows,
                    writer_reports: &writer_reports,
                    master_s,
                    peak_rss_mb: crate::autotune::process_peak_rss_mb(),
                    writer_peak_mb: crate::h5_out::concurrent_writer_peak_mb(&writer_reports, writer_threads),
                };
                profile::write_profile_json(path, &rp)?;
                eprintln!("Profile JSON:          {}", path.display());
            }
            let read_mb_s = mb_per_s(traffic.estimated_bytes, process_s);
            let write_mb_s = mb_per_s(out_bytes, write_s);
            println!(
                "cells={} CG={:.4}% ({}/{}) CH={:.4}% ({}/{}) elapsed_s={:.2} process_s={:.2} write_s={:.2} read_mb_s={:.1} write_mb_s={:.1} h5={}",
                r.n_cells(),
                r.stats.cg_pct(),
                r.stats.cg_c,
                r.stats.cg_c + r.stats.cg_t,
                r.stats.ch_pct(),
                r.stats.ch_c,
                r.stats.ch_c + r.stats.ch_t,
                total_s,
                process_s,
                write_s,
                read_mb_s,
                write_mb_s,
                h5_path.display()
            );
        }
    }
    Ok(())
}

enum PendingOutput {
    Memory(Vec<Vec<ChunkFile>>),
    Stream(tempfile::TempDir),
}

fn mb_per_s(bytes: u64, seconds: f64) -> f64 {
    if seconds <= 0.0 {
        0.0
    } else {
        (bytes as f64) / 1_000_000.0 / seconds
    }
}

fn fmt_u64(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// Python `:,.Nf` grouping for the summary throughput lines.
fn fmt_grouped_f64(n: f64, decimals: usize) -> String {
    let neg = n.is_sign_negative();
    let n = if n.is_finite() { n.abs() } else { 0.0 };
    let scale = 10u64.pow(decimals as u32);
    let scaled = (n * scale as f64).round() as u64;
    let ip = if decimals == 0 { scaled } else { scaled / scale };
    let mut body = fmt_u64(ip);
    if decimals > 0 {
        let frac = scaled % scale;
        body.push('.');
        body.push_str(&format!("{frac:0width$}", width = decimals));
    }
    if neg {
        format!("-{body}")
    } else {
        body
    }
}

fn print_extraction_summary(
    r: &crate::extract::WindowResult,
    traffic: &BamTraffic,
    total_s: f64,
    process_s: f64,
    write_s: f64,
    out_bytes: u64,
    n_shards: usize,
    stream: bool,
    peak_rss_mb: f64,
    writer_peak_mb: f64,
    n_writer_procs: usize,
    context: ContextSel,
) {
    let cg_c = r.stats.cg_c;
    let cg_t = r.stats.cg_t;
    let chg_c = r.stats.chg_c;
    let chg_t = r.stats.chg_t;
    let chh_c = r.stats.chh_c;
    let chh_t = r.stats.chh_t;
    let ch_c = r.stats.ch_c;
    let ch_t = r.stats.ch_t;
    let cg_all = cg_c + cg_t;
    let chg_all = chg_c + chg_t;
    let chh_all = chh_c + chh_t;
    let ch_all = ch_c + ch_t;
    let pct = |c: u64, tot: u64| if tot == 0 { 0.0 } else { 100.0 * c as f64 / tot as f64 };
    let total_calls = cg_all + ch_all;
    let total_mod = cg_c + ch_c;
    let reads = traffic.mapped_reads;
    let reads_per_sec = if total_s > 0.0 { reads as f64 / total_s } else { 0.0 };
    let reads_per_min = reads_per_sec * 60.0;
    let cells = r.n_cells();
    let cells_per_min = if total_s > 0.0 { cells as f64 / (total_s / 60.0) } else { 0.0 };
    let read_mb = traffic.estimated_bytes as f64 / 1_000_000.0;
    let write_mb = out_bytes as f64 / 1_000_000.0;
    let mode = if stream {
        "Disk-backed temp files"
    } else {
        "In-memory IPC"
    };

    eprintln!("======================================================================");
    eprintln!("PARALLEL EXTRACTION SUMMARY & GROUND TRUTH COMPARISON");
    eprintln!("======================================================================");
    eprintln!("Total Wall Time:       {total_s:.1}s ({:.2} min)", total_s / 60.0);
    eprintln!("Processing Time:       {process_s:.1}s ({:.2} min)", process_s / 60.0);
    eprintln!("HDF5 Write Time:       {write_s:.1}s ({:.2} min)", write_s / 60.0);
    eprintln!(
        "Read Throughput:       {:.1} MB/s (est. {:.1} MB BAM)",
        mb_per_s(traffic.estimated_bytes, process_s),
        read_mb
    );
    eprintln!(
        "Write Throughput:      {:.1} MB/s ({write_mb:.1} MB HDF5)",
        mb_per_s(out_bytes, write_s)
    );
    eprintln!("Total Mapped Reads:    {}", fmt_u64(reads));
    eprintln!(
        "Throughput (Reads):    {} reads/min ({} reads/s)",
        fmt_grouped_f64(reads_per_min, 0),
        fmt_grouped_f64(reads_per_sec, 0)
    );
    eprintln!("Total Unique Barcodes: {}", fmt_u64(cells as u64));
    eprintln!(
        "Throughput (Cells):    {} cells/min",
        fmt_grouped_f64(cells_per_min, 1)
    );
    eprintln!("Peak Memory (RSS):     {peak_rss_mb:.1} MB");
    if n_writer_procs > 0 {
        eprintln!(
            "Writer Peak RSS:       {writer_peak_mb:.1} MB ({n_writer_procs} child processes; concurrent upper bound)"
        );
    }
    eprintln!("Shards Written:        {n_shards}");
    eprintln!("Streaming Mode:        {mode}");
    eprintln!("----------------------------------------------------------------------");
    let line = |label: &str, want: bool, c: u64, tot: u64| {
        if want {
            eprintln!("{label}{:.3}% ({} / {})", pct(c, tot), fmt_u64(c), fmt_u64(tot));
        } else {
            eprintln!("{label}not extracted (--context {})", if context == ContextSel::Cg { "CG" } else { "CH" });
        }
    };
    let cg_on = context.wants(crate::context::Context::Cg);
    let ch_on = context.wants(crate::context::Context::Ch);
    line("CpG Methylation:       ", cg_on, cg_c, cg_all);
    line("CHG Methylation:       ", ch_on, chg_c, chg_all);
    line("CHH Methylation:       ", ch_on, chh_c, chh_all);
    line("CH (Combined):         ", ch_on, ch_c, ch_all);
    line("Pooled (All C/G sites):", true, total_mod, total_calls);
    eprintln!("======================================================================");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn summary_numbers_match_python_grouping() {
        assert_eq!(fmt_u64(0), "0");
        assert_eq!(fmt_u64(999), "999");
        assert_eq!(fmt_u64(1000), "1,000");
        assert_eq!(fmt_u64(74_900_000), "74,900,000");
        assert_eq!(fmt_grouped_f64(15_390_411.4, 0), "15,390,411");
        assert_eq!(fmt_grouped_f64(1234.56, 1), "1,234.6");
        assert_eq!(mb_per_s(0, 0.0), 0.0);
        assert!((mb_per_s(2_000_000, 2.0) - 1.0).abs() < 1e-9);
    }
}
