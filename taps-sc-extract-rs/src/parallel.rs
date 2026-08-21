//! Rayon window pool.
//!
//! Thread-safety:
//! - Each pool thread owns a private `IndexedReader` (thread-local). Never `Send`ed.
//! - `FastFaiReader` is `Arc`-shared; each `fetch` opens its own `File`.
//! - Interners and cell maps are per-window, then either written to a unique
//!   `shard_XXX/chunk_YYYYYY.bin` or merged by barcode string after join.
//! - Cancellation is `Arc<AtomicBool>`; workers only load it.

use crate::extract::{
    open_bam, process_window, refresh_baq_cache, ExtractParams, WindowResult,
};
use crate::fasta::FastFaiReader;
use crate::shard_io::{prepare_shard_dirs, write_window_shards};
use crate::window::Window;
use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rust_htslib::bam::IndexedReader;
use std::cell::RefCell;
use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

pub type CancelFlag = Arc<AtomicBool>;

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

thread_local! {
    static TLS_BAM: RefCell<Option<(PathBuf, usize, IndexedReader)>> = const { RefCell::new(None) };
    static TLS_BAQ: RefCell<Option<(String, Vec<u8>)>> = const { RefCell::new(None) };
}

fn with_thread_bam<R>(
    bam_path: &Path,
    decomp: usize,
    f: impl FnOnce(&mut IndexedReader, &mut Option<(String, Vec<u8>)>) -> Result<R>,
) -> Result<R> {
    TLS_BAM.with(|bam_cell| {
        TLS_BAQ.with(|baq_cell| {
            {
                let mut slot = bam_cell.borrow_mut();
                let need = match slot.as_ref() {
                    None => true,
                    Some((p, d, _)) => p.as_path() != bam_path || *d != decomp,
                };
                if need {
                    let bam = open_bam(bam_path, decomp)?;
                    *slot = Some((bam_path.to_path_buf(), decomp, bam));
                    *baq_cell.borrow_mut() = None;
                }
            }
            let mut bam_slot = bam_cell.borrow_mut();
            let bam = &mut bam_slot.as_mut().unwrap().2;
            let mut baq_slot = baq_cell.borrow_mut();
            f(bam, &mut baq_slot)
        })
    })
}

fn process_one(
    bam_path: &Path,
    fasta: &FastFaiReader,
    window: &Window,
    whitelist: Option<&HashSet<String>>,
    params: &ExtractParams,
) -> Result<WindowResult> {
    with_thread_bam(bam_path, params.decomp_threads, |bam, cache| {
        let baq_ref = if params.compute_baq {
            refresh_baq_cache(fasta, cache, &window.contig)?
        } else {
            None
        };
        process_window(bam, fasta, window, whitelist, params, baq_ref)
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
            {
                let mut t = total.lock().expect("stats mutex");
                if params.accumulate {
                    t.absorb_maps(&r);
                } else {
                    t.absorb_counts(&r);
                }
            }
            let k = done.fetch_add(1, Ordering::Relaxed) + 1;
            if k % 5 == 0 || k == n {
                let t = total.lock().expect("stats mutex");
                eprintln!(
                    "Rust extract: [{k}/{n}] cells={} CG={:.2}% CH={:.3}%",
                    t.n_cells(),
                    t.stats.cg_pct(),
                    t.stats.ch_pct()
                );
            }
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
            {
                let mut t = total.lock().expect("stream mutex");
                t.absorb_counts(&r);
            }
            let k = done.fetch_add(1, Ordering::Relaxed) + 1;
            if k % 5 == 0 || k == n {
                let t = total.lock().expect("stream mutex");
                eprintln!(
                    "Rust extract: [{k}/{n}] cells={} files={} CG={:.2}% CH={:.3}%",
                    t.n_cells(),
                    n_files.load(Ordering::Relaxed),
                    t.stats.cg_pct(),
                    t.stats.ch_pct()
                );
            }
            Ok(())
        })
    })??;

    check_cancel(cancel)?;
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
            {
                let mut t = total.lock().expect("mem stats");
                t.absorb_counts(&r);
            }
            let k = done.fetch_add(1, Ordering::Relaxed) + 1;
            if k % 5 == 0 || k == n {
                let t = total.lock().expect("mem stats");
                eprintln!(
                    "Rust extract: [{k}/{n}] cells={} CG={:.2}% CH={:.3}% (memory)",
                    t.n_cells(),
                    t.stats.cg_pct(),
                    t.stats.ch_pct()
                );
            }
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
