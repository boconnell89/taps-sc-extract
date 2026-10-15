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
use std::time::Instant;

use crate::extract::{estimate_bam_traffic, BamTraffic, ExtractParams};
use crate::extract_noodles::extract_regions_stats_noodles;
use crate::fasta::FastFaiReader;
use crate::h5_out::{
    assemble_hdf5_from_memory, assemble_hdf5_from_temp, written_output_bytes, H5Compression,
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
        /// Drop this many 5' sequenced bases of read 1 before calling.
        #[arg(long, default_value_t = 0)]
        trim_r1: u32,
        /// Drop this many 5' sequenced bases of read 2 before calling.
        #[arg(long, default_value_t = 10)]
        trim_r2: u32,
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
        /// Genomic window size in megabases (Rust default: 100).
        #[arg(long, default_value_t = 100)]
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
        /// HDF5 compression (Rust default: lzf).
        #[arg(long, default_value = "lzf")]
        compression: String,
        /// Maximum parallel shard-writer threads (0 = auto sized from shards, cores, RAM).
        #[arg(long, default_value_t = 0)]
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
        /// Drop this many 5' sequenced bases of read 1 before calling.
        #[arg(long, default_value_t = 0)]
        trim_r1: u32,
        /// Drop this many 5' sequenced bases of read 2 before calling.
        #[arg(long, default_value_t = 10)]
        trim_r2: u32,
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
            trim_r1,
            trim_r2,
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
                trim_r1,
                trim_r2,
            };
            eprintln!(
                "taps-sc-extract-rs stats: pileup={pileup} workers={n_workers} decomp={decomp} baq={compute_baq} overlaps={ignore_overlaps} orphans={ignore_orphans} trim_r1={trim_r1} trim_r2={trim_r2} {} window(s), contig(s)={:?}, chunk_size_mb={chunk_size_mb}",
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
            trim_r1,
            trim_r2,
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

            // PR7 Auto-tuning & Phase 1 Writer Auto-scaling
            let (budget_gb, budget_src) = autotune::system_memory_budget_gb(max_memory_gb);
            let cell_source = autotune::CellCountSource::resolve(
                expected_cells,
                wl.as_ref().map(|w| w.len()),
                autotune::DEFAULT_ESTIMATED_CELLS,
            );
            let (n_workers, workers_src) = autotune::resolve_workers_with_budget(workers, budget_gb);
            let decomp = decomp_threads.unwrap_or_else(|| autotune::default_decomp_threads(n_workers));
            let (n_shards, shards_src) = autotune::resolve_shards(shards, cell_source.count());
            let (writer_threads, writer_src) = autotune::resolve_writer_threads(max_writer_threads, n_shards, budget_gb);
            let total_genome_mb = windows.len() as u64 * u64::from(chunk_size_mb);
            let (resolved_memory_mode, mode_src) = autotune::resolve_memory_mode(
                &memory_mode,
                budget_gb,
                cell_source.count(),
                total_genome_mb,
            );
            let compression_name = compression.clone();
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
                trim_r1,
                trim_r2,
            };
            eprintln!("======================================================================");
            eprintln!("taps-sc-extract-rs v{}", env!("CARGO_PKG_VERSION"));
            eprintln!("Configuration:");
            eprintln!("  Memory budget:   {:.1} GB ({})", budget_gb, budget_src);
            eprintln!("  Expected cells:  {} ({})", cell_source.count(), cell_source.label());
            eprintln!("  Workers:         {} ({})", n_workers, workers_src);
            eprintln!("  Decomp threads:  {} per worker", decomp);
            eprintln!("  Shards:          {} ({})", n_shards, shards_src);
            eprintln!("  Writer threads:  {} ({})", writer_threads, writer_src);
            eprintln!("  Memory mode:     {} ({})", resolved_memory_mode, mode_src);
            eprintln!("  Chunk size:      {chunk_size_mb} Mb");
            eprintln!("  Compression:     {compression_name}");
            eprintln!("  BAQ:             {}", params.compute_baq);
            eprintln!("  5' trim:         R1={} bp, R2={} bp", params.trim_r1, params.trim_r2);
            eprintln!("  Windows:         {} window(s) across {:?}", windows.len(), contig_list);
            eprintln!("  Output:          {}", out.display());
            eprintln!("======================================================================");
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
            let t_process = Instant::now();
            let (r, pending) = match resolved_memory_mode.as_str() {
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
                        n_shards,
                        n_workers,
                        &cancel,
                    )?;
                    (r, PendingOutput::Stream(tmp))
                }
                other => anyhow::bail!("unknown --memory-mode {other} (use stream or memory)"),
            };
            let process_s = t_process.elapsed().as_secs_f64();
            eprintln!(
                "Extraction workers finished in {process_s:.1}s. Writing {n_shards} shard(s) with {writer_threads} writer(s)."
            );
            let t_write = Instant::now();
            let h5_path = match pending {
                PendingOutput::Memory(mem) => {
                    assemble_hdf5_from_memory(&mem, &out, compression, writer_threads)?
                }
                PendingOutput::Stream(tmp) => {
                    let p = assemble_hdf5_from_temp(
                        tmp.path(),
                        &out,
                        n_shards,
                        compression,
                        writer_threads,
                    )?;
                    if keep_temp {
                        let kept = tmp.keep();
                        eprintln!("kept temp {}", kept.display());
                    }
                    p
                }
            };
            let write_s = t_write.elapsed().as_secs_f64();
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
            );
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
    eprintln!("Shards Written:        {n_shards}");
    eprintln!("Streaming Mode:        {mode}");
    eprintln!("----------------------------------------------------------------------");
    eprintln!(
        "CpG Methylation:       {:.3}% ({c} / {tot})",
        pct(cg_c, cg_all),
        c = fmt_u64(cg_c),
        tot = fmt_u64(cg_all)
    );
    eprintln!(
        "CHG Methylation:       {:.3}% ({c} / {tot})",
        pct(chg_c, chg_all),
        c = fmt_u64(chg_c),
        tot = fmt_u64(chg_all)
    );
    eprintln!(
        "CHH Methylation:       {:.3}% ({c} / {tot})",
        pct(chh_c, chh_all),
        c = fmt_u64(chh_c),
        tot = fmt_u64(chh_all)
    );
    eprintln!(
        "CH (Combined):         {:.3}% ({c} / {tot})",
        pct(ch_c, ch_all),
        c = fmt_u64(ch_c),
        tot = fmt_u64(ch_all)
    );
    eprintln!(
        "Pooled (All C/G sites):{:.3}% ({c} / {tot})",
        pct(total_mod, total_calls),
        c = fmt_u64(total_mod),
        tot = fmt_u64(total_calls)
    );
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
