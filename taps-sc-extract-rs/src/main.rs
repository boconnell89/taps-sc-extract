//! taps-sc-extract-rs: Rust core for TAPS single-cell methylation extraction.

#![allow(dead_code)]

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod accumulate;
mod autotune;
mod baq_gpu;
mod baq_simd;
mod barcode;
mod bisulfite;
mod calling;
mod chemistry;
mod context;
mod extract;
mod extract_noodles;
mod fasta;
mod h5_out;
mod input;
mod inputs;
mod mbias;
mod multibam;
mod parallel;
mod profile;
mod readcall;
mod shard_io;
mod shard_stream;
mod sitemask;
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
use crate::window::plan_windows;

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
        /// Comma-separated contigs (default: the BAM's primary chromosomes: autosomes, X, Y).
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
        /// Maximum reads per position, all cells together (htslib cap); 0 = unlimited.
        #[arg(long, default_value_t = extract::DEFAULT_MAX_DEPTH)]
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
        /// Input BAM(s): one or more paths, or @FILE naming a list (one path
        /// per line) or a sample sheet (TSV with a `bam` header column and
        /// optional complexity, min_complexity, max_complexity, whitelist).
        /// A barcode must belong to one input; reads are merged per window.
        #[arg(short = 'b', long, num_args = 1.., required = true)]
        bam: Vec<PathBuf>,
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
        /// Directory for stream-mode temp chunk files (default: $TMPDIR, else /tmp).
        #[arg(long)]
        temp_dir: Option<PathBuf>,
        /// stream = temp chunks then HDF5; memory = keep in RAM; auto = heuristic from budget & cells.
        #[arg(long, default_value = "auto")]
        memory_mode: String,
        /// Memory budget in GB for the whole process tree (overrides --max-memory-fraction).
        #[arg(long)]
        max_memory_gb: Option<f64>,
        /// Share of MemAvailable the run may use (default 0.85).
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
        /// Conversion chemistry: auto (from the BAM: BSBolt YS tags = bisulfite,
        /// paired-end without them = taps), taps (mC -> T, paired-end) or
        /// bisulfite (unmethylated C -> T, single-end BSBolt). The two use
        /// separate strand and call rules and different defaults.
        #[arg(long, value_enum, default_value_t = chemistry::ChemistryArg::Auto)]
        chemistry: chemistry::ChemistryArg,
        /// Bisulfite strand: auto (flag, checked against every read's BSBolt
        /// YS tag; stops at the first disagreement) or flag (no YS check).
        #[arg(long, value_enum, default_value_t = chemistry::StrandSource::Auto)]
        strand_source: chemistry::StrandSource,
        /// Minimum base quality (default 20).
        #[arg(long)]
        min_baseq: Option<u8>,
        /// Minimum mapping quality (default 10: bwa gives 0 to reads placed at random among equal hits).
        #[arg(long)]
        min_mapq: Option<u8>,
        /// htslib pileup cap, all cells together (default 0 = unlimited). A cap
        /// thins deep positions by read order; --site-depth-limit drops them.
        #[arg(long)]
        max_depth: Option<u32>,
        /// Most observations one cell keeps at a site: a number (0 = off) or
        /// auto (default) = --molecules-per-cell / 2 (a site is called on one
        /// strand: one molecule per haplotype), x 2 for TAPS with
        /// --no-overlap-clip. Extra observations are scaled away keeping the
        /// t:c ratio. Relax it (raise --molecules-per-cell, or 0) when real
        /// extra copies are expected: copy-number gains, segmental
        /// duplications, aneuploidy (e.g. cancer), replicating cells.
        #[arg(long)]
        cell_site_cap: Option<crate::extract::DepthArg>,
        /// Drop every call at positions covered by more kept-cell reads than
        /// this: a number (0 = off) or auto (default) = --molecules-per-cell x
        /// reads per molecule (taps 2, bisulfite 1) x cells kept by
        /// --complexity (with its cutoff) / --whitelist. Off without either.
        #[arg(long)]
        site_depth_limit: Option<crate::extract::DepthArg>,
        /// Molecules one cell can have at a position (both strands): 4 = 2
        /// strands x 2 haplotypes; 8 allows a replicated (G2) genome. Sets
        /// --cell-site-cap auto and --site-depth-limit auto.
        #[arg(long, default_value_t = crate::extract::MOLECULES_PER_CELL)]
        molecules_per_cell: u32,
        /// BED of regions whose calls are dropped for every cell (e.g. the
        /// ENCODE blacklist). The pileup is unchanged.
        #[arg(long)]
        blacklist: Option<PathBuf>,
        #[arg(long, default_value_t = false)]
        no_overlap_clip: bool,
        #[arg(long, default_value_t = false)]
        no_ignore_orphans: bool,
        /// Disable BAQ. Default: on for taps, off for bisulfite.
        #[arg(long, default_value_t = false)]
        no_baq: bool,
        /// BAQ mode: taps (default when BAQ is on), taps-full, taps-partial or
        /// htslib (taps only). Setting it turns BAQ on for bisulfite.
        #[arg(long, value_enum)]
        baq: Option<BaqMode>,
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
        /// Drop this many bases from the sequenced 5' end of read 1 (or of
        /// every single-end read; alias --trim-5p). Default 0.
        #[arg(long, alias = "trim-5p")]
        trim_r1: Option<u32>,
        /// Drop this many bases from the sequenced 5' end of read 2 (default:
        /// taps 10; no effect on single-end reads).
        #[arg(long)]
        trim_r2: Option<u32>,
        /// Drop this many bases from the sequenced 3' end of every read
        /// (default: taps 0, bisulfite 10 for the Tn5 gap fill-in).
        #[arg(long, alias = "trim-3prime")]
        trim_3p: Option<u32>,
        /// auto: measure M-bias on a sample of reads and pick the trims for
        /// ends not set explicitly (writes <bam>.mbias.tsv/.png/.pdf next to
        /// the output); manual: explicit values or chemistry defaults.
        /// Default: auto for bisulfite, manual for taps.
        #[arg(long, value_enum)]
        trim: Option<TrimMode>,
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

/// `--trim`: automatic (M-bias) or manual trims.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
enum TrimMode {
    Auto,
    Manual,
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
                        if let Ok(f) = std::env::var("TAPS_DUMP_STDIN") {
                            if shard == 0 {
                                let mut buf = Vec::new();
                                std::io::Read::read_to_end(&mut std::io::stdin().lock(), &mut buf)?;
                                std::fs::write(&f, &buf)?;
                                return shard_io::read_all_chunks_from(&mut buf.as_slice());
                            }
                        }
                        if let Ok(f) = std::env::var("TAPS_READ_FILE") {
                            let t = std::time::Instant::now();
                            let buf = std::fs::read(&f)?;
                            let t1 = t.elapsed().as_secs_f64();
                            let v = shard_io::read_all_chunks_from(&mut buf.as_slice());
                            eprintln!("file read {t1:.3}s parse {:.3}s", t.elapsed().as_secs_f64() - t1);
                            return v;
                        }
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
            if chemistry::detect(&bam, false)?.with_ys > 0 {
                anyhow::bail!("`stats` supports TAPS only; this BAM carries BSBolt bisulfite tags (use `extract`)");
            }
            // Default: the BAM's primary chromosomes (autosomes + sex chromosomes).
            let primary = match &chroms {
                Some(_) => Vec::new(),
                None => {
                    let v = crate::window::bam_primary_contigs(&bam)?;
                    eprintln!("Contigs: {} from the BAM header (--chroms to choose)", crate::window::describe_primary(&v));
                    v
                }
            };
            let contig_list: Vec<&str> = match &chroms {
                Some(s) => s.split(',').map(|x| x.trim()).filter(|x| !x.is_empty()).collect(),
                None => primary.iter().map(|s| s.as_str()).collect(),
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
                site_depth_limit: 0,
                blacklist: None,
                cell_site_cap: 0,
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
                chemistry: crate::chemistry::Chemistry::Taps,
                trim_3p: 0,
                verify_ys: false,
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
                    std::slice::from_ref(&bam),
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
            temp_dir,
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
            cell_site_cap,
            site_depth_limit,
            molecules_per_cell,
            blacklist,
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
            chemistry,
            strand_source,
            trim_r1,
            trim_r2,
            trim_3p,
            trim,
            whitelist,
            complexity,
            min_complexity,
            max_complexity,
            profile_json,
        } => {
            let fai = FastFaiReader::open(&fasta)?;
            let nproc = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
            let spec = crate::inputs::parse_bam_args(&bam)?;
            let prepared = spec
                .bams
                .iter()
                .map(|b| input::prepare(b, &fasta, nproc))
                .collect::<Result<Vec<_>>>()?;
            let bams: Vec<PathBuf> = prepared.iter().map(|p| p.path.clone()).collect();
            let bam = bams[0].clone();
            if bams.len() > 1 {
                crate::inputs::check_same_contigs(&bams)?;
                eprintln!("Inputs:          {} BAMs, merged per window", bams.len());
            }
            // Default: the BAM's primary chromosomes (autosomes + sex chromosomes).
            let primary = match &chroms {
                Some(_) => Vec::new(),
                None => {
                    let v = crate::window::bam_primary_contigs(&bam)?;
                    eprintln!("Contigs: {} from the BAM header (--chroms to choose)", crate::window::describe_primary(&v));
                    v
                }
            };
            let contig_list: Vec<&str> = match &chroms {
                Some(s) => s.split(',').map(|x| x.trim()).filter(|x| !x.is_empty()).collect(),
                None => primary.iter().map(|s| s.as_str()).collect(),
            };
            // Chemistry from the BAM, then its defaults for every unset option.
            let detected = chemistry::detect(&bam, strand_source == chemistry::StrandSource::Auto)?;
            let chem = chemistry::resolve(chemistry, &detected, strand_source)?;
            // Every input must resolve to the same chemistry: runs never mix them.
            for b in &bams[1..] {
                let d = chemistry::detect(b, strand_source == chemistry::StrandSource::Auto)?;
                let c = chemistry::resolve(chemistry, &d, strand_source)
                    .map_err(|e| anyhow::anyhow!("chemistry of {}: {e}", b.display()))?;
                if c != chem {
                    anyhow::bail!(
                        "{} is {} but {} is {}: one run takes one chemistry",
                        bam.display(),
                        chem.label(),
                        b.display(),
                        c.label()
                    );
                }
            }
            chemistry::check_reference(&fai, &bam, &contig_list)?;
            let dflt = chem.defaults();
            if chem == chemistry::Chemistry::Bisulfite && baq == Some(BaqMode::Htslib) {
                anyhow::bail!("--baq htslib treats every unmethylated C as a mismatch in bisulfite data; use a taps BAQ mode or --no-baq");
            }
            if chem == chemistry::Chemistry::Bisulfite && trim_r2.is_some() {
                eprintln!("warning: --trim-r2 has no effect on single-end reads; use --trim-5p / --trim-3p");
            }
            let min_baseq = min_baseq.unwrap_or(dflt.min_baseq);
            let min_mapq = min_mapq.unwrap_or(dflt.min_mapq);
            // Cells kept by --complexity (with its cutoff; every barcode in
            // the file without one) and/or --whitelist; with several inputs,
            // the union of each input's own filter.
            let cell_filter = crate::inputs::load_cell_filters(
                &spec,
                whitelist.as_deref(),
                complexity.as_deref(),
                min_complexity,
                max_complexity,
            )?;
            crate::inputs::check_barcodes_match(&bams, &cell_filter.listed)?;
            if bams.len() > 1 {
                for (b, (kept, listed)) in bams.iter().zip(&cell_filter.per_file) {
                    let f = |v: &Option<usize>| v.map_or_else(|| "all".to_string(), |n| n.to_string());
                    eprintln!("  {}: {} cells kept of {} listed", b.display(), f(kept), f(listed));
                }
            }
            let wl = cell_filter.keep.clone();
            let max_depth = max_depth.unwrap_or(dflt.max_depth);
            if molecules_per_cell == 0 {
                anyhow::bail!("--molecules-per-cell must be at least 1");
            }
            let (site_depth_limit, site_depth_src) = match (site_depth_limit, wl.as_ref().map(|w| w.len())) {
                (Some(crate::extract::DepthArg::Fixed(n)), _) => (n, "explicit".to_string()),
                (_, Some(n)) => (
                    crate::extract::auto_site_depth_limit(n, molecules_per_cell, chem.reads_per_molecule()),
                    format!(
                        "auto: {molecules_per_cell} molecules x {} read(s) per molecule x {n} cells",
                        chem.reads_per_molecule()
                    ),
                ),
                (Some(crate::extract::DepthArg::Auto), None) => {
                    anyhow::bail!("--site-depth-limit auto needs --complexity or --whitelist to count the cells; give a number instead")
                }
                (None, None) => (0, "auto; off without --complexity or --whitelist".to_string()),
            };
            let reads_per_call = if chem == chemistry::Chemistry::Taps && no_overlap_clip { 2 } else { 1 };
            let (cell_site_cap, cell_cap_src) = match cell_site_cap {
                Some(crate::extract::DepthArg::Fixed(n)) => (n, "explicit".to_string()),
                _ => (
                    crate::extract::auto_cell_site_cap(molecules_per_cell, reads_per_call),
                    format!("auto: {molecules_per_cell} molecules / 2 strands x {reads_per_call} call(s) per molecule"),
                ),
            };
            let blacklist_set = match blacklist.as_deref() {
                Some(p) => Some(std::sync::Arc::new(crate::sitemask::Blacklist::load(p)?)),
                None => None,
            };
            // --trim auto (M-bias) is the default for both chemistries.
            let trim_mode = trim.unwrap_or(TrimMode::Auto);
            let (set_r1, set_r2, set_3p) = (trim_r1, trim_r2, trim_3p);
            let mut trim_r1 = trim_r1.unwrap_or(dflt.trim_r1);
            let mut trim_r2 = trim_r2.unwrap_or(dflt.trim_r2);
            let mut trim_3p = trim_3p.unwrap_or(dflt.trim_3p);
            let no_baq = no_baq || (baq.is_none() && !dflt.baq);
            let baq = baq.unwrap_or(BaqMode::Taps);
            crate::h5_out::set_run_attrs(&[
                ("chemistry", chem.label().to_string()),
                ("bsbolt_version", detected.bsbolt.clone().unwrap_or_default()),
                ("trim_r1", trim_r1.to_string()),
                ("trim_r2", trim_r2.to_string()),
                ("trim_3p", trim_3p.to_string()),
                ("trim_mode", "manual".to_string()),
                ("min_mapq", min_mapq.to_string()),
                ("min_baseq", min_baseq.to_string()),
                ("max_depth", max_depth.to_string()),
                ("site_depth_limit", site_depth_limit.to_string()),
                ("cell_site_cap", cell_site_cap.to_string()),
                ("inputs", bams.iter().map(|b| b.display().to_string()).collect::<Vec<_>>().join(",")),
                ("blacklist", blacklist.as_ref().map(|p| p.display().to_string()).unwrap_or_default()),
                ("baq", if no_baq { "off".to_string() } else { format!("{baq:?}").to_lowercase() }),
                ("extractor", format!("taps-sc-extract-rs {}", env!("CARGO_PKG_VERSION"))),
            ]);

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
            // Several inputs: every worker reads every file, so per-file
            // decompression threads would multiply; inflate on the worker.
            let decomp = decomp_threads
                .or_else(|| saved.as_ref().and_then(|(p, _)| p.decomp_threads))
                .unwrap_or_else(|| if bams.len() > 1 { 0 } else { autotune::default_decomp_threads(n_workers) });
            // Summed over every input.
            let mut traffic = BamTraffic { mapped_reads: 0, estimated_bytes: 0 };
            for b in &bams {
                match estimate_bam_traffic(b, &contig_list) {
                    Ok(t) => {
                        traffic.mapped_reads += t.mapped_reads;
                        traffic.estimated_bytes += t.estimated_bytes;
                    }
                    Err(e) => {
                        eprintln!("warning: BAM read estimate unavailable for {}: {e}", b.display());
                        traffic.estimated_bytes += std::fs::metadata(b).map(|m| m.len()).unwrap_or(0);
                    }
                }
            }
            let est_sites = autotune::estimate_sites(
                traffic.mapped_reads,
                context.site_share(),
                cell_source.count(),
                selected_bp / 1_000_000,
            );
            let (n_shards, shards_src) = autotune::resolve_shards(
                shards,
                cell_source.count(),
                est_sites,
                nproc,
                crate::h5_out::is_h5_filename(&out) && !out.is_dir(),
            );
            let writer_mode = WriterMode::parse(&writer_mode, n_shards)?;
            let streams_enabled = std::env::var("TAPS_STREAM_WRITERS").map_or(true, |v| v != "0");
            if trim_mode == TrimMode::Auto {
                // --trim auto: M-bias on a fixed sample, trims for ends not set.
                let t0 = std::time::Instant::now();
                use rust_htslib::bam::Read as _;
                // Windows from every input's mapped reads together; each input
                // is counted over the same windows and the tables added.
                let mut sel: Vec<(String, u64, u64)> = Vec::new();
                for b in &bams {
                    let mut rd = crate::input::open_indexed(b)?;
                    let stats = rd.index_stats().map_err(|e| anyhow::anyhow!("index stats {}: {e}", b.display()))?;
                    let names: Vec<String> =
                        rd.header().target_names().iter().map(|n| String::from_utf8_lossy(n).into_owned()).collect();
                    for c in &contig_list {
                        let Some(tid) = names.iter().position(|n| n == c) else { continue };
                        let Some((_, len, mapped, _)) = stats.iter().find(|s| s.0 == tid as i64) else { continue };
                        match sel.iter_mut().find(|e| e.0 == *c) {
                            Some(e) => e.2 += *mapped,
                            None => sel.push((c.to_string(), *len, *mapped)),
                        }
                    }
                }
                let windows = mbias::sample_windows(&sel, start, end);
                let mut counts = mbias::count(&bam, &fai, &windows, chem, min_mapq, min_baseq)?;
                for b in &bams[1..] {
                    counts.add(&mbias::count(b, &fai, &windows, chem, min_mapq, min_baseq)?);
                }
                let mates: &[u8] = if chem == chemistry::Chemistry::Bisulfite { &[1] } else { &[1, 2] };
                let chosen = mbias::choose(&counts, chem, mates);
                if set_r1.is_none() {
                    trim_r1 = chosen.trim_r1;
                }
                if set_r2.is_none() && chem == chemistry::Chemistry::Taps {
                    trim_r2 = chosen.trim_r2;
                }
                if set_3p.is_none() {
                    trim_3p = chosen.trim_3p;
                }
                let prefix = mbias::output_prefix(&out, n_shards, spec.name_path());
                if let Some(dir) = prefix.parent() {
                    std::fs::create_dir_all(dir)?;
                }
                let tsv = PathBuf::from(format!("{}.mbias.tsv", prefix.display()));
                mbias::write_tsv(&tsv, &counts, trim_r1, trim_r2, trim_3p)?;
                eprintln!(
                    "Auto trim: 5' {trim_r1} bp, 3' {trim_3p} bp from {} sampled reads in {:.1}s (M-bias: {}.mbias.tsv/.png/.pdf)",
                    counts.reads,
                    t0.elapsed().as_secs_f64(),
                    prefix.display()
                );
                for n in &chosen.notes {
                    eprintln!("warning: auto trim: {n}");
                }
                mbias::plot(&tsv, &prefix.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default());
                crate::h5_out::set_run_attrs(&[
                    ("chemistry", chem.label().to_string()),
                    ("bsbolt_version", detected.bsbolt.clone().unwrap_or_default()),
                    ("trim_r1", trim_r1.to_string()),
                    ("trim_r2", trim_r2.to_string()),
                    ("trim_3p", trim_3p.to_string()),
                    ("trim_mode", "auto".to_string()),
                    ("min_mapq", min_mapq.to_string()),
                    ("min_baseq", min_baseq.to_string()),
                    ("max_depth", max_depth.to_string()),
                    ("site_depth_limit", site_depth_limit.to_string()),
                ("cell_site_cap", cell_site_cap.to_string()),
                ("inputs", bams.iter().map(|b| b.display().to_string()).collect::<Vec<_>>().join(",")),
                    ("blacklist", blacklist.as_ref().map(|p| p.display().to_string()).unwrap_or_default()),
                    ("baq", if no_baq { "off".to_string() } else { format!("{baq:?}").to_lowercase() }),
                    ("extractor", format!("taps-sc-extract-rs {}", env!("CARGO_PKG_VERSION"))),
                ]);
            }
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
                // Same conditions as the ShardStreams start below.
                streamable: writer_mode == WriterMode::Process && n_shards > 1 && streams_enabled,
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
                site_depth_limit,
                blacklist: blacklist_set.clone(),
                cell_site_cap,
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
                chemistry: chem,
                trim_3p,
                verify_ys: chem == chemistry::Chemistry::Bisulfite && strand_source == chemistry::StrandSource::Auto,
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
                eprintln!(
                    "  WARNING: the estimated peak exceeds the memory budget and the run may run out of memory; \
                     explicit settings are followed as given (lower --max-writer-threads, raise --shards, or use --memory-mode stream)"
                );
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
            eprintln!(
                "  Chemistry:       {} ({}{})",
                params.chemistry.label(),
                match (chemistry, params.chemistry) {
                    (chemistry::ChemistryArg::Auto, _) => "detected from the BAM",
                    _ => "explicit",
                },
                detected.bsbolt.as_ref().map(|v| format!("; BSBolt {v}")).unwrap_or_default()
            );
            eprintln!(
                "  Read filters:    MAPQ >= {}, base quality >= {}, pileup cap {}{}",
                params.min_mapq,
                params.min_baseq,
                if params.max_depth == 0 { "unlimited".to_string() } else { params.max_depth.to_string() },
                if params.chemistry == chemistry::Chemistry::Bisulfite { ", supplementary dropped" } else { "" }
            );
            if params.chemistry == chemistry::Chemistry::Bisulfite {
                eprintln!("  Trims:           5' {} bp, 3' {} bp (sequenced ends)", params.trim_r1, params.trim_3p);
            } else {
                eprintln!("  5' trim:         R1={} bp, R2={} bp; 3' {} bp", params.trim_r1, params.trim_r2, params.trim_3p);
            }
            eprintln!(
                "  Site filters:    per-cell cap {} ({cell_cap_src}); depth limit {} ({site_depth_src}); blacklist {}",
                if cell_site_cap == 0 { "off".to_string() } else { cell_site_cap.to_string() },
                if site_depth_limit == 0 { "off".to_string() } else { site_depth_limit.to_string() },
                match (&blacklist, &blacklist_set) {
                    (Some(p), Some(b)) => format!("{} ({} regions, {:.2} Mb)", p.display(), b.n_regions(), b.total_bp() as f64 / 1e6),
                    _ => "none".to_string(),
                }
            );
            eprintln!("  Windows:         {} window(s) across {:?}", windows.len(), contig_list);
            eprintln!("  Output:          {}", out.display());
            eprintln!("======================================================================");
            let n_buckets = crate::h5_out::bucket_count(n_shards);
            let cpu_start = profile::cpu_times();
            let t_process = Instant::now();
            // Memory mode with one writer process per shard: start the writers
            // now and pipe each window's payloads as it finishes, so the
            // transfer overlaps extraction (TAPS_STREAM_WRITERS=0 turns it off).
            let streams = if resolved_memory_mode == "memory"
                && writer_mode == crate::h5_out::WriterMode::Process
                && n_shards > 1
                && writer_threads >= n_shards
                && streams_enabled
            {
                match crate::h5_out::resolve_output_target(&out, n_shards)? {
                    crate::h5_out::OutputTarget::ShardDir(dir) => Some(crate::shard_stream::ShardStreams::start(
                        &dir,
                        n_shards,
                        n_buckets,
                        windows.len(),
                        &compression_name,
                    )?),
                    crate::h5_out::OutputTarget::File(_) => None,
                }
            } else {
                None
            };
            let (r, pending) = match resolved_memory_mode.as_str() {
                "memory" => {
                    let (r, mem) = extract_memory_parallel(
                        &bams,
                        &fasta,
                        &windows,
                        wl.as_ref(),
                        &params,
                        n_buckets,
                        n_workers,
                        &cancel,
                        streams.as_ref(),
                    )?;
                    match streams {
                        Some(st) => (r, PendingOutput::Streamed(st)),
                        None => (r, PendingOutput::Memory(mem)),
                    }
                }
                "stream" => {
                    let mut builder = tempfile::Builder::new();
                    builder.prefix("taps-rs-");
                    let tmp = match &temp_dir {
                        Some(d) => {
                            std::fs::create_dir_all(d).map_err(|e| anyhow::anyhow!("--temp-dir {}: {e}", d.display()))?;
                            builder.tempdir_in(d)
                        }
                        None => builder.tempdir(),
                    }
                    .map_err(|e| anyhow::anyhow!("tempdir: {e}"))?;
                    let r = extract_stream_parallel(
                        &bams,
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
            // Re-plan the write phase from the measured shard payloads and the
            // memory free now (writers fed during extraction are already running).
            let write_plan = {
                let weights = match &pending {
                    PendingOutput::Streamed(_) => None,
                    PendingOutput::Memory(mem) => {
                        Some(crate::h5_out::bucket_sites(mem).iter().map(|s| s * 12).collect::<Vec<u64>>())
                    }
                    PendingOutput::Stream(tmp) => Some(crate::h5_out::bucket_bytes(tmp.path(), n_buckets)),
                };
                weights.map(|w| {
                    let (cap, explicit) = if max_writer_threads > 0 {
                        (max_writer_threads, true)
                    } else if let Some(p) = saved.as_ref().and_then(|(p, _)| p.writers) {
                        (p, false)
                    } else {
                        (nproc, false)
                    };
                    let cap = if writer_mode == WriterMode::Thread { cap.min(autotune::MAX_AUTO_WRITERS) } else { cap };
                    let budget = autotune::write_budget_bytes(max_memory_gb, max_memory_fraction);
                    autotune::plan_write(&w, n_shards, shards_src == "explicit", cap, explicit, budget, nproc)
                })
            };
            if let Some(p) = &write_plan {
                for n in &p.notes {
                    eprintln!("note: {n}");
                }
            }
            let writer_threads = write_plan.as_ref().map_or(writer_threads, |p| p.writers);
            let n_files = write_plan.as_ref().map_or(n_shards, |p| p.groups.len());
            eprintln!(
                "Extraction workers finished in {process_s:.1}s. Writing {n_files} shard(s) with {writer_threads} {} writer(s){}.",
                writer_mode.label(),
                write_plan.as_ref().map_or(String::new(), |p| format!(" (~{:.1} GB each)", p.writer_rss / 1e9))
            );
            let t_write = Instant::now();
            let crate::h5_out::HdfOutput { path: h5_path, reports: writer_reports, master_s } = match pending {
                PendingOutput::Streamed(st) => {
                    let (dir, reports) = st.finish()?;
                    let t_master = Instant::now();
                    let path = crate::h5_out::write_master_h5(&dir, n_shards)?;
                    crate::h5_out::HdfOutput { path, reports, master_s: t_master.elapsed().as_secs_f64() }
                }
                PendingOutput::Memory(mem) => {
                    let groups = write_plan.as_ref().expect("planned for memory output").groups.clone();
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
                    let groups = write_plan.as_ref().expect("planned for stream output").groups.clone();
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
                n_files,
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
                        ("shards", n_files.to_string()),
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
    /// Memory mode with writers fed during extraction (`shard_stream`).
    Streamed(crate::shard_stream::ShardStreams),
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
    let t = &r.timing;
    if t.site_depth_bp > 0 || t.blacklist_sites > 0 || t.cap_sites > 0 {
        eprintln!(
            "Site filters:          per-cell cap trimmed {} calls at {} cell-sites; depth limit dropped {} calls at {} cell-sites ({} bp over the limit); blacklist dropped {} calls at {} cell-sites",
            t.cap_calls, t.cap_sites, t.site_depth_calls, t.site_depth_sites, t.site_depth_bp, t.blacklist_calls, t.blacklist_sites
        );
    }
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
