//! Memory mode with process writers, overlapped with extraction.
//!
//! Every shard's `write-shard` child is started before extraction, and each
//! window's payloads are piped to it as soon as the window finishes, so the
//! children receive and decode their shard while the workers are still
//! running. A cell's dataset still needs every window, so the HDF5 writing
//! itself starts once extraction ends; what moves out of the write phase is
//! the transfer (parent encode, pipe, child decode).
//!
//! Shard balance: extraction hashes cells into buckets
//! (`h5_out::bucket_count`). The non-streamed path deals buckets to shards by
//! their final site counts (`pack_buckets`). Here buckets are dealt once
//! `assign_after` windows are done, by their site counts so far; payloads
//! that arrive earlier are held until then. A cell stays in one bucket, so it
//! still goes to exactly one shard, and each child sorts its chunks by
//! `chunk_id` before writing, as before.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc;
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;

use anyhow::{Context, Result};

use crate::h5_out::{chunk_pipe, finish_shard_writer, pack_buckets, spawn_shard_writer, WriterReport};
use crate::shard_io::{write_chunk_to, ChunkFile};

struct State {
    /// Bucket -> shard, once assigned.
    assign: Option<Vec<usize>>,
    /// Payloads held per bucket until assignment, and sites per bucket.
    held: Vec<Vec<ChunkFile>>,
    sites: Vec<u64>,
    windows_done: usize,
}

pub struct ShardStreams {
    state: Mutex<State>,
    senders: Vec<mpsc::Sender<ChunkFile>>,
    threads: Vec<JoinHandle<Result<WriterReport>>>,
    n_shards: usize,
    assign_after: usize,
    dir: PathBuf,
    /// Set when dropped without `finish` (extraction failed or was
    /// cancelled): the children are killed instead of writing partial shards.
    abort: Arc<AtomicBool>,
}

/// Nice level for the transfer during extraction (`TAPS_STREAM_NICE`, default
/// 10): the extraction workers get the CPU first, and senders and children
/// catch up on spare cycles (and at full speed once the workers are done).
/// `tid` 0 is the calling thread.
fn lower_priority(tid: i32) {
    let nice: i32 = std::env::var("TAPS_STREAM_NICE").ok().and_then(|v| v.parse().ok()).unwrap_or(10);
    if nice > 0 {
        let id = if tid == 0 { unsafe { libc::syscall(libc::SYS_gettid) as i32 } } else { tid };
        unsafe { libc::setpriority(libc::PRIO_PROCESS, id as libc::id_t, nice) };
    }
}

fn sites_of(c: &ChunkFile) -> u64 {
    c.cells.iter().map(|x| (x.cg.len() + x.ch.len()) as u64).sum()
}

impl ShardStreams {
    /// Start one `write-shard` child per shard, writing `dir/shard_XXX.h5`.
    pub fn start(dir: &Path, n_shards: usize, n_buckets: usize, n_windows: usize, compression_name: &str) -> Result<Self> {
        std::fs::create_dir_all(dir).with_context(|| format!("mkdir {}", dir.display()))?;
        let abort = Arc::new(AtomicBool::new(false));
        let mut senders = Vec::with_capacity(n_shards);
        let mut threads = Vec::with_capacity(n_shards);
        for s in 0..n_shards {
            let path = dir.join(format!("shard_{s:03}.h5"));
            let mut child = spawn_shard_writer(s, &path, compression_name, &[], true)?;
            lower_priority(child.id() as i32);
            let (tx, rx) = mpsc::channel::<ChunkFile>();
            senders.push(tx);
            let abort = abort.clone();
            let h = std::thread::Builder::new()
                .name(format!("taps-h5s-{s}"))
                .spawn(move || -> Result<WriterReport> {
                    lower_priority(0);
                    let mut w = chunk_pipe(&mut child);
                    let sent: Result<()> = (|| {
                        for c in rx {
                            write_chunk_to(&mut w, &c)?;
                        }
                        w.flush()?;
                        Ok(())
                    })();
                    if abort.load(Ordering::SeqCst) {
                        let _ = child.kill();
                        let _ = child.wait();
                        anyhow::bail!("write-shard {s:03} aborted");
                    }
                    drop(w);
                    if let Err(e) = sent {
                        let _ = child.kill();
                        let _ = child.wait();
                        return Err(e.context(format!("send chunks to write-shard {s:03}")));
                    }
                    finish_shard_writer(s, child)
                })
                .context("start shard sender thread")?;
            threads.push(h);
        }
        Ok(ShardStreams {
            state: Mutex::new(State {
                assign: None,
                held: (0..n_buckets).map(|_| Vec::new()).collect(),
                sites: vec![0; n_buckets],
                windows_done: 0,
            }),
            senders,
            threads,
            n_shards,
            assign_after: (n_windows / 5).max(1),
            dir: dir.to_path_buf(),
            abort,
        })
    }

    /// One finished window's payloads (bucket, chunk).
    pub fn push_window(&self, payloads: Vec<(usize, ChunkFile)>) {
        let mut st = self.state.lock().expect("shard stream state");
        match &st.assign {
            Some(a) => {
                for (b, c) in payloads {
                    // A closed channel means the sender failed; finish() reports it.
                    let _ = self.senders[a[b]].send(c);
                }
            }
            None => {
                for (b, c) in payloads {
                    st.sites[b] += sites_of(&c);
                    st.held[b].push(c);
                }
            }
        }
        st.windows_done += 1;
        if st.assign.is_none() && st.windows_done >= self.assign_after {
            self.assign(&mut st);
        }
    }

    /// Deal buckets to shards by their sites so far and release held payloads.
    fn assign(&self, st: &mut State) {
        let groups = pack_buckets(&st.sites, self.n_shards);
        let mut a = vec![0usize; st.sites.len()];
        for (s, g) in groups.iter().enumerate() {
            for &b in g {
                a[b] = s;
            }
        }
        for (b, held) in st.held.iter_mut().enumerate() {
            for c in held.drain(..) {
                let _ = self.senders[a[b]].send(c);
            }
        }
        st.assign = Some(a);
    }

    /// Close every pipe (after extraction) and wait for the writers.
    pub fn finish(mut self) -> Result<(PathBuf, Vec<WriterReport>)> {
        {
            let mut st = self.state.lock().expect("shard stream state");
            if st.assign.is_none() {
                self.assign(&mut st);
            }
        }
        self.senders.clear();
        let mut reports = Vec::with_capacity(self.threads.len());
        for h in std::mem::take(&mut self.threads) {
            reports.push(h.join().map_err(|_| anyhow::anyhow!("shard sender thread panicked"))??);
        }
        Ok((self.dir.clone(), reports))
    }
}

impl Drop for ShardStreams {
    fn drop(&mut self) {
        if self.threads.is_empty() {
            return; // finished
        }
        self.abort.store(true, Ordering::SeqCst);
        self.senders.clear();
        for h in std::mem::take(&mut self.threads) {
            let _ = h.join();
        }
    }
}
