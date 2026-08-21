//! taps-sc-extract-rs: Rust core for TAPS single-cell methylation extraction.

#![allow(dead_code)]

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod accumulate;
mod autotune;
mod barcode;
mod calling;
mod context;
mod extract;
mod extract_noodles;
mod fasta;
mod h5_out;
mod parallel;
mod shard_io;
mod window;

use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

use crate::extract::ExtractParams;
use crate::extract_noodles::extract_regions_stats_noodles;
use crate::fasta::FastFaiReader;
use crate::h5_out::{assemble_hdf5_from_memory, assemble_hdf5_from_temp, H5Compression};
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
        #[arg(long, default_value_t = 10)]
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
        #[arg(short = 'w', long)]
        whitelist: Option<PathBuf>,
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
        #[arg(long, default_value_t = 10)]
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
        /// Optional memory budget in GB (defaults to conservative 0.6 * MemAvailable).
        #[arg(long)]
        max_memory_gb: Option<f64>,
        /// Optional expected cell count (auto-detected from whitelist if omitted).
        #[arg(long)]
        expected_cells: Option<usize>,
        #[arg(long, default_value = "gzip")]
        compression: String,
        #[arg(long, default_value_t = 6)]
        max_writer_threads: usize,
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
        #[arg(short = 'w', long)]
        whitelist: Option<PathBuf>,
    },
}

fn main() -> Result<()> {
    let cancel = new_cancel_flag();
    let _ = install_ctrlc(cancel.clone());
    let cli = Cli::parse();
    match cli.cmd {
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
            whitelist,
            pileup,
            dump_sites,
        } => {
            let fai = FastFaiReader::open(&fasta)?;
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
            let wl = match whitelist {
                Some(p) => Some(barcode::parse_annot(&p)?),
                None => None,
            };
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
                decomp_threads: decomp,
                accumulate: dump_sites.is_some(),
            };
            eprintln!(
                "taps-sc-extract-rs stats: pileup={pileup} workers={n_workers} decomp={decomp} baq={compute_baq} overlaps={ignore_overlaps} orphans={ignore_orphans} {} window(s), contig(s)={:?}, chunk_size_mb={chunk_size_mb}",
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
            expected_cells,
            compression,
            max_writer_threads,
            min_baseq,
            min_mapq,
            max_depth,
            no_overlap_clip,
            no_ignore_orphans,
            no_baq,
            whitelist,
        } => {
            let fai = FastFaiReader::open(&fasta)?;
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
            let wl = match whitelist {
                Some(p) => Some(barcode::parse_annot(&p)?),
                None => None,
            };

            // PR7 Auto-tuning
            let (budget_gb, budget_src) = autotune::system_memory_budget_gb(max_memory_gb);
            let cell_source = autotune::CellCountSource::resolve(
                expected_cells,
                wl.as_ref().map(|w| w.len()),
                autotune::DEFAULT_ESTIMATED_CELLS,
            );
            let (n_workers, workers_src) = autotune::resolve_workers_with_budget(workers, budget_gb);
            let decomp = decomp_threads.unwrap_or_else(|| autotune::default_decomp_threads(n_workers));
            let (n_shards, shards_src) = autotune::resolve_shards(shards, cell_source.count());
            let total_genome_mb = windows.len() as u64 * u64::from(chunk_size_mb);
            let (resolved_memory_mode, mode_src) = autotune::resolve_memory_mode(
                &memory_mode,
                budget_gb,
                cell_source.count(),
                total_genome_mb,
            );
            let compression = H5Compression::parse(&compression)?;
            let params = ExtractParams {
                min_baseq,
                min_mapq,
                max_depth,
                ignore_overlaps: !no_overlap_clip,
                ignore_orphans: !no_ignore_orphans,
                compute_baq: !no_baq,
                decomp_threads: decomp,
                accumulate: true,
            };
            eprintln!("======================================================================");
            eprintln!("taps-sc-extract-rs v{}", env!("CARGO_PKG_VERSION"));
            eprintln!("Configuration:");
            eprintln!("  Memory budget:   {:.1} GB ({})", budget_gb, budget_src);
            eprintln!("  Expected cells:  {} ({})", cell_source.count(), cell_source.label());
            eprintln!("  Workers:         {} ({})", n_workers, workers_src);
            eprintln!("  Decomp threads:  {} per worker", decomp);
            eprintln!("  Shards:          {} ({})", n_shards, shards_src);
            eprintln!("  Memory mode:     {} ({})", resolved_memory_mode, mode_src);
            eprintln!("  BAQ:             {}", params.compute_baq);
            eprintln!("  Windows:         {} window(s) across {:?}", windows.len(), contig_list);
            eprintln!("  Output:          {}", out.display());
            eprintln!("======================================================================");
            let t0 = std::time::Instant::now();
            let (r, h5_path) = match resolved_memory_mode.as_str() {
                "memory" => {
                    let (r, mem) = extract_memory_parallel(
                        &bam,
                        &fasta,
                        &windows,
                        wl.as_ref(),
                        &params,
                        n_shards,
                        n_workers,
                        &cancel,
                    )?;
                    eprintln!("assembling HDF5 from memory (writers≤{max_writer_threads})…");
                    let p = assemble_hdf5_from_memory(&mem, &out, compression, max_writer_threads)?;
                    (r, p)
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
                        n_shards,
                        n_workers,
                        &cancel,
                    )?;
                    eprintln!("assembling HDF5 from temp chunks (writers≤{max_writer_threads})…");
                    let p = assemble_hdf5_from_temp(
                        tmp.path(),
                        &out,
                        n_shards,
                        compression,
                        max_writer_threads,
                    )?;
                    if keep_temp {
                        let kept = tmp.keep();
                        eprintln!("kept temp {}", kept.display());
                    }
                    (r, p)
                }
                other => anyhow::bail!("unknown --memory-mode {other} (use stream or memory)"),
            };
            let dt = t0.elapsed();
            println!(
                "cells={} CG={:.4}% ({}/{}) CH={:.4}% ({}/{}) elapsed_s={:.2} h5={}",
                r.n_cells(),
                r.stats.cg_pct(),
                r.stats.cg_c,
                r.stats.cg_c + r.stats.cg_t,
                r.stats.ch_pct(),
                r.stats.ch_c,
                r.stats.ch_c + r.stats.ch_t,
                dt.as_secs_f64(),
                h5_path.display()
            );
        }
    }
    Ok(())
}
