//! Rayon window pool.
//!
//! Thread-safety:
//! - Each pool thread owns a private `IndexedReader` (thread-local). Never `Send`ed.
//! - `FastFaiReader` is `Arc`-shared; each `fetch` opens its own `File`.
//! - Interners and cell maps are per-window, then either written to a unique
//!   `shard_XXX/chunk_YYYYYY.bin` or merged by barcode string after join.
//! - Cancellation is `Arc<AtomicBool>`; workers only load it.

use crate::extract::{open_bam, process_window, ExtractParams, WindowResult};
use crate::fasta::FastFaiReader;
use crate::shard_io::{prepare_shard_dirs, write_window_shards};
use crate::window::Window;
use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rust_htslib::bam::IndexedReader;
use std::cell::RefCell;
use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

pub type CancelFlag = Arc<AtomicBool>;

/// Suppress per-window progress lines (set by `tune`, which runs many short extractions).
pub static QUIET_PROGRESS: AtomicBool = AtomicBool::new(false);

pub fn new_cancel_flag() -> CancelFlag {
    Arc::new(AtomicBool::new(false))
}

pub fn install_ctrlc(flag: CancelFlag) -> Result<()> {
    ctrlc::set_handler(move || {
        flag.store(true, Ordering::SeqCst);
        eprintln!("taps-sc-extract-rs: SIGINT, cancelling workers…");
    })
    .context("ctrlc handler")?;
    Ok(())
}

fn check_cancel(flag: &CancelFlag) -> Result<()> {
    if flag.load(Ordering::Relaxed) {
        bail!("cancelled");
    }
    Ok(())
}

struct RunClock {
    start: Instant,
    last_log: Mutex<Instant>,
    peak_rss: AtomicU64,
}

impl RunClock {
    fn new() -> Self {
        let now = Instant::now();
        Self {
            start: now,
            last_log: Mutex::new(now),
            peak_rss: AtomicU64::new(0),
        }
    }
}

fn fmt_usize(n: usize) -> String {
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

/// Python-style chunk progress: percent, chunk rate, elapsed, and ETA.
fn log_extraction_progress(
    clock: &RunClock,
    n: usize,
    k: usize,
    n_cells: usize,
    cg_c: u64,
    cg_t: u64,
    ch_c: u64,
    ch_t: u64,
) {
    if QUIET_PROGRESS.load(Ordering::Relaxed) {
        return;
    }
    let now = Instant::now();
    let mut last = clock.last_log.lock().expect("progress clock");
    let due = k == n || k % 10 == 0 || now.duration_since(*last) >= Duration::from_secs(15);
    if !due {
        return;
    }
    *last = now;
    drop(last);

    let elapsed = clock.start.elapsed().as_secs_f64();
    let pct = if n == 0 { 100.0 } else { 100.0 * k as f64 / n as f64 };
    let chunk_rate = if elapsed > 0.0 { k as f64 / elapsed * 60.0 } else { 0.0 };
    let eta_min = if k > 0 && elapsed > 0.0 {
        ((n - k) as f64 / (k as f64 / elapsed)) / 60.0
    } else {
        0.0
    };
    let cg_tot = cg_c + cg_t;
    let cg_pct = if cg_tot == 0 { 0.0 } else { 100.0 * cg_c as f64 / cg_tot as f64 };
    let ch_tot = ch_c + ch_t;
    let ch_pct = if ch_tot == 0 { 0.0 } else { 100.0 * ch_c as f64 / ch_tot as f64 };
    let rss = crate::autotune::process_rss_mb();
    let rss_bytes = (rss * 1024.0 * 1024.0) as u64;
    let prev = clock.peak_rss.fetch_max(rss_bytes, Ordering::Relaxed);
    let peak_mb = prev.max(rss_bytes) as f64 / (1024.0 * 1024.0);
    eprintln!(
        "Extraction: [{k:3}/{n} chunks | {pct:5.1}%] | Rate: {chunk_rate:4.1} chk/min | Elapsed: {elapsed:5.1}s | ETA: {eta_min:4.1}m | Cells: {} | CG: {cg_pct:.1}% | CH: {ch_pct:.2}% | RSS: {rss:.0} MB (Peak: {peak_mb:.0} MB)",
        fmt_usize(n_cells)
    );
}

thread_local! {
    static TLS_BAM: RefCell<Option<(PathBuf, usize, IndexedReader)>> = const { RefCell::new(None) };
}

fn with_thread_bam<R>(
    bam_path: &Path,
    decomp: usize,
    f: impl FnOnce(&mut IndexedReader) -> Result<R>,
) -> Result<R> {
    TLS_BAM.with(|bam_cell| {
        let mut slot = bam_cell.borrow_mut();
        let need = match slot.as_ref() {
            None => true,
            Some((p, d, _)) => p.as_path() != bam_path || *d != decomp,
        };
        if need {
            let bam = open_bam(bam_path, decomp)?;
            *slot = Some((bam_path.to_path_buf(), decomp, bam));
        }
        f(&mut slot.as_mut().unwrap().2)
    })
}

fn process_one(
    bam_path: &Path,
    fasta: &FastFaiReader,
    window: &Window,
    whitelist: Option<&HashSet<String>>,
    params: &ExtractParams,
) -> Result<WindowResult> {
    with_thread_bam(bam_path, params.decomp_threads, |bam| {
        process_window(bam, fasta, window, whitelist, params)
    })
}

fn install_pool<T>(n_workers: usize, f: impl FnOnce() -> T + Send) -> Result<T>
where
    T: Send,
{
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_workers.max(1))
        .thread_name(|i| format!("taps-win-{i}"))
        .build()
        .context("rayon thread pool")?;
    Ok(pool.install(f))
}

/// Parallel stats (optional per-pos maps if `params.accumulate`).
pub fn extract_stats_parallel(
    bam_path: &Path,
    fasta_path: &Path,
    windows: &[Window],
    whitelist: Option<&HashSet<String>>,
    params: &ExtractParams,
    n_workers: usize,
    cancel: &CancelFlag,
) -> Result<WindowResult> {
    let fasta = Arc::new(FastFaiReader::open(fasta_path)?);
    let wl = whitelist.map(|w| Arc::new(w.clone()));
    let bam_path = bam_path.to_path_buf();
    let params = params.clone();
    let n = windows.len();
    let done = AtomicUsize::new(0);
    let total = Mutex::new(WindowResult::empty());
    let clock = RunClock::new();

    install_pool(n_workers, || {
        windows.par_iter().try_for_each(|w| -> Result<()> {
            check_cancel(cancel)?;
            let r = process_one(
                &bam_path,
                fasta.as_ref(),
                w,
                wl.as_ref().map(|a| a.as_ref()),
                &params,
            )?;
            let snap = {
                let mut t = total.lock().expect("stats mutex");
                if params.accumulate {
                    t.absorb_maps(&r);
                } else {
                    t.absorb_counts(&r);
                if params.profile {
                    t.window_timings.push((w.chunk_id, r.timing.clone()));
                }
                }
                (t.n_cells(), t.stats.cg_c, t.stats.cg_t, t.stats.ch_c, t.stats.ch_t)
            };
            let k = done.fetch_add(1, Ordering::Relaxed) + 1;
            log_extraction_progress(&clock, n, k, snap.0, snap.1, snap.2, snap.3, snap.4);
            Ok(())
        })
    })??;

    check_cancel(cancel)?;
    Ok(total.into_inner().expect("stats mutex"))
}

/// Parallel stream extract: each window writes unique `chunk_YYYYYY.bin` then drops maps.
pub fn extract_stream_parallel(
    bam_path: &Path,
    fasta_path: &Path,
    windows: &[Window],
    whitelist: Option<&HashSet<String>>,
    params: &ExtractParams,
    temp_dir: &Path,
    n_shards: usize,
    n_workers: usize,
    cancel: &CancelFlag,
) -> Result<WindowResult> {
    prepare_shard_dirs(temp_dir, n_shards)?;
    let fasta = Arc::new(FastFaiReader::open(fasta_path)?);
    let wl = whitelist.map(|w| Arc::new(w.clone()));
    let bam_path = bam_path.to_path_buf();
    let temp_dir = temp_dir.to_path_buf();
    let mut params = params.clone();
    params.accumulate = true;
    let n = windows.len();
    let done = AtomicUsize::new(0);
    let n_files = AtomicUsize::new(0);
    let total = Mutex::new(WindowResult::empty());
    let clock = RunClock::new();

    install_pool(n_workers, || {
        windows.par_iter().try_for_each(|w| -> Result<()> {
            check_cancel(cancel)?;
            let r = process_one(
                &bam_path,
                fasta.as_ref(),
                w,
                wl.as_ref().map(|a| a.as_ref()),
                &params,
            )?;
            let nf = write_window_shards(&temp_dir, w, &r.intern, &r.cells, n_shards)?;
            n_files.fetch_add(nf, Ordering::Relaxed);
            let snap = {
                let mut t = total.lock().expect("stream mutex");
                t.absorb_counts(&r);
                if params.profile {
                    t.window_timings.push((w.chunk_id, r.timing.clone()));
                }
                (t.n_cells(), t.stats.cg_c, t.stats.cg_t, t.stats.ch_c, t.stats.ch_t)
            };
            let k = done.fetch_add(1, Ordering::Relaxed) + 1;
            log_extraction_progress(&clock, n, k, snap.0, snap.1, snap.2, snap.3, snap.4);
            Ok(())
        })
    })??;

    check_cancel(cancel)?;
    eprintln!(
        "Extraction workers finished. Temp chunk files: {}",
        n_files.load(Ordering::Relaxed)
    );
    Ok(total.into_inner().expect("stream mutex"))
}

/// In-memory extract: keep per-shard window payloads, no temp files.
pub fn extract_memory_parallel(
    bam_path: &Path,
    fasta_path: &Path,
    windows: &[Window],
    whitelist: Option<&HashSet<String>>,
    params: &ExtractParams,
    n_shards: usize,
    n_workers: usize,
    cancel: &CancelFlag,
) -> Result<(WindowResult, Vec<Vec<crate::shard_io::ChunkFile>>)> {
    let fasta = Arc::new(FastFaiReader::open(fasta_path)?);
    let wl = whitelist.map(|w| Arc::new(w.clone()));
    let bam_path = bam_path.to_path_buf();
    let mut params = params.clone();
    params.accumulate = true;
    let n = windows.len();
    let n_shards = n_shards.max(1);
    let done = AtomicUsize::new(0);
    let total = Mutex::new(WindowResult::empty());
    let mem: Mutex<Vec<Vec<crate::shard_io::ChunkFile>>> =
        Mutex::new((0..n_shards).map(|_| Vec::new()).collect());
    let clock = RunClock::new();

    install_pool(n_workers, || {
        windows.par_iter().try_for_each(|w| -> Result<()> {
            check_cancel(cancel)?;
            let r = process_one(
                &bam_path,
                fasta.as_ref(),
                w,
                wl.as_ref().map(|a| a.as_ref()),
                &params,
            )?;
            let payloads = crate::shard_io::window_shard_payloads(w, &r.intern, &r.cells, n_shards);
            {
                let mut m = mem.lock().expect("mem mutex");
                for (s, cf) in payloads {
                    m[s].push(cf);
                }
            }
            let snap = {
                let mut t = total.lock().expect("mem stats");
                t.absorb_counts(&r);
                if params.profile {
                    t.window_timings.push((w.chunk_id, r.timing.clone()));
                }
                (t.n_cells(), t.stats.cg_c, t.stats.cg_t, t.stats.ch_c, t.stats.ch_t)
            };
            let k = done.fetch_add(1, Ordering::Relaxed) + 1;
            log_extraction_progress(&clock, n, k, snap.0, snap.1, snap.2, snap.3, snap.4);
            Ok(())
        })
    })??;

    check_cancel(cancel)?;
    Ok((
        total.into_inner().expect("mem stats"),
        mem.into_inner().expect("mem mutex"),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_send<T: Send>() {}
    fn assert_sync<T: Sync>() {}

    #[test]
    fn send_sync_audit() {
        assert_send::<ExtractParams>();
        assert_send::<Window>();
        assert_send::<WindowResult>();
        assert_send::<FastFaiReader>();
        assert_send::<CancelFlag>();
        assert_sync::<FastFaiReader>();
        assert_sync::<ExtractParams>();
        assert_sync::<Window>();
        assert_sync::<CancelFlag>();
        // IndexedReader stays thread-local (raw htslib pointers, not shared).
    }

    #[test]
    fn cancel_flag_trips_check() {
        let f = new_cancel_flag();
        assert!(check_cancel(&f).is_ok());
        f.store(true, Ordering::SeqCst);
        assert!(check_cancel(&f).is_err());
    }
}
