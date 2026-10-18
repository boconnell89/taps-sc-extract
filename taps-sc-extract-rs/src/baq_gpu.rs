//! GPU BAQ (`--baq-device gpu`; the default is the CPU).
//!
//! The HMM of TAPS-aware BAQ runs on the GPU; everything around it stays on
//! the CPU in `csrc/baq_simd.c` (span, band, TAPS mask, skip rule, BAQ_EXTEND,
//! quality cap). Results are identical to the CPU path read for read: the
//! device kernel (`csrc/gpu/baq_hmm_body.h`) does htslib's double arithmetic in
//! the same order, gets HMM parameters computed on the host, and returns a
//! quality only when it cannot differ from glibc's; otherwise the host
//! recomputes it (`baq-selftest --gpu` checks this).
//!
//! Backends: OpenCL (always built; Intel Arc / Data Center GPUs, also AMD,
//! NVIDIA, PoCL), CUDA (built when `nvcc` is present) and HIP (when `hipcc` is
//! present). One service thread owns the device; workers hand it their
//! prepared reads and wait, and it pools requests from all workers into large
//! launches sorted by read length.

use crate::baq_simd::{BaqJob, Width, Workspace};
use anyhow::Result;
use std::os::raw::{c_char, c_int, c_void};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{mpsc, Mutex, OnceLock};
use std::time::{Duration, Instant};

/// `--baq-device`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum BaqDevice {
    /// SIMD kernels on the CPU (default).
    Cpu,
    /// The GPU, through `--gpu-backend`.
    Gpu,
}

/// `--gpu-backend`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum GpuBackend {
    /// CUDA, then HIP, then OpenCL: the first that finds a usable GPU.
    Auto,
    /// OpenCL (Intel Arc / Data Center GPUs; any device with fp64).
    Opencl,
    /// NVIDIA CUDA (when built with nvcc).
    Cuda,
    /// AMD HIP (when built with hipcc).
    Hip,
}

/// Per-read device input; layout matches `taps_gpu_read` (128 bytes).
#[repr(C)]
#[derive(Clone, Copy, Default)]
pub struct GpuRead {
    m: [f64; 9],
    s_m: f64,
    b_m: f64,
    b_i: f64,
    l_ref: i32,
    l_query: i32,
    bw: i32,
    i_dim: i32,
    ref_off: u32,
    base_off: u32,
    pad: [i32; 2],
}
const _: () = assert!(std::mem::size_of::<GpuRead>() == 128);

#[repr(C)]
#[derive(Default)]
struct Sizes {
    n_dev: i64,
    ref_bytes: i64,
    bases: i64,
    n_cpu: i64,
}

extern "C" {
    fn taps_baq_gpu_prepare(ws: *mut c_void, jobs: *const BaqJob, n: c_int, sz: *mut Sizes, cpu_jobs: *mut c_int);
    fn taps_baq_gpu_export(ws: *mut c_void, jobs: *const BaqJob, rd: *mut GpuRead, refb: *mut u8, qry: *mut u8, qual: *mut u8);
    fn taps_baq_gpu_finish(ws: *mut c_void, jobs: *const BaqJob, state: *const i16, q: *const u8, maxv: *const f64);
    fn taps_baq_gpu_q2p(out: *mut f32);
    fn taps_ocl_prof(out: *mut f64);
    fn taps_ocl_open(src: *const c_char, want: *const c_char, name: *mut c_char, name_len: c_int, f_budget: usize) -> *mut c_void;
    fn taps_ocl_run(
        ctx: *mut c_void, rd: *const GpuRead, n: c_int, refb: *const u8, nref: usize, qry: *const u8, qual: *const u8,
        nbase: usize, q2p: *const f32, state: *mut i16, q: *mut u8, maxv: *mut f64, nflag: *mut c_int,
    ) -> c_int;
}

#[cfg(taps_cuda)]
extern "C" {
    fn taps_cuda_open(want: *const c_char, name: *mut c_char, name_len: c_int, f_budget: usize) -> *mut c_void;
    fn taps_cuda_run(
        ctx: *mut c_void, rd: *const GpuRead, n: c_int, refb: *const u8, nref: usize, qry: *const u8, qual: *const u8,
        nbase: usize, q2p: *const f32, state: *mut i16, q: *mut u8, maxv: *mut f64, nflag: *mut c_int,
    ) -> c_int;
}

#[cfg(taps_hip)]
extern "C" {
    fn taps_hip_open(want: *const c_char, name: *mut c_char, name_len: c_int, f_budget: usize) -> *mut c_void;
    fn taps_hip_run(
        ctx: *mut c_void, rd: *const GpuRead, n: c_int, refb: *const u8, nref: usize, qry: *const u8, qual: *const u8,
        nbase: usize, q2p: *const f32, state: *mut i16, q: *mut u8, maxv: *mut f64, nflag: *mut c_int,
    ) -> c_int;
}

/// OpenCL program: prelude + shared kernel body + entry point.
const CL_SRC: &str = concat!(
    include_str!("../csrc/gpu/baq_hmm_cl_prefix.cl"),
    include_str!("../csrc/gpu/baq_hmm_body.h"),
    include_str!("../csrc/gpu/baq_hmm_cl_kernel.cl"),
    "\0"
);

/// Backends compiled into this binary.
pub fn built_backends() -> Vec<&'static str> {
    let mut v = vec!["opencl"];
    if cfg!(taps_cuda) {
        v.push("cuda");
    }
    if cfg!(taps_hip) {
        v.push("hip");
    }
    v
}

struct Engine {
    kind: GpuBackend,
    ctx: *mut c_void,
}

// The device context is used only by the service thread that opened it.
unsafe impl Send for Engine {}

impl Engine {
    fn open(kind: GpuBackend, want: &str, f_budget: usize) -> Result<(Engine, String), String> {
        let want_c = std::ffi::CString::new(want).unwrap_or_default();
        let mut name = vec![0 as c_char; 512];
        let ctx = unsafe {
            match kind {
                GpuBackend::Opencl => taps_ocl_open(
                    CL_SRC.as_ptr() as *const c_char,
                    want_c.as_ptr(),
                    name.as_mut_ptr(),
                    name.len() as c_int,
                    f_budget,
                ),
                #[cfg(taps_cuda)]
                GpuBackend::Cuda => taps_cuda_open(want_c.as_ptr(), name.as_mut_ptr(), name.len() as c_int, f_budget),
                #[cfg(taps_hip)]
                GpuBackend::Hip => taps_hip_open(want_c.as_ptr(), name.as_mut_ptr(), name.len() as c_int, f_budget),
                _ => return Err(format!("{kind:?} backend not built into this binary (built: {})", built_backends().join(", "))),
            }
        };
        let msg = unsafe { std::ffi::CStr::from_ptr(name.as_ptr()) }.to_string_lossy().into_owned();
        if ctx.is_null() {
            Err(msg)
        } else {
            Ok((Engine { kind, ctx }, msg))
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn run(&mut self, rd: &[GpuRead], refb: &[u8], qry: &[u8], qual: &[u8], q2p: &[f32; 256], state: &mut [i16], q: &mut [u8], maxv: &mut [f64]) -> Result<i32, i32> {
        let mut nflag: c_int = 0;
        let args = (
            self.ctx,
            rd.as_ptr(),
            rd.len() as c_int,
            refb.as_ptr(),
            refb.len(),
            qry.as_ptr(),
            qual.as_ptr(),
            qry.len(),
            q2p.as_ptr(),
            state.as_mut_ptr(),
            q.as_mut_ptr(),
            maxv.as_mut_ptr(),
            &mut nflag as *mut c_int,
        );
        let rc = unsafe {
            match self.kind {
                GpuBackend::Opencl => taps_ocl_run(args.0, args.1, args.2, args.3, args.4, args.5, args.6, args.7, args.8, args.9, args.10, args.11, args.12),
                #[cfg(taps_cuda)]
                GpuBackend::Cuda => taps_cuda_run(args.0, args.1, args.2, args.3, args.4, args.5, args.6, args.7, args.8, args.9, args.10, args.11, args.12),
                #[cfg(taps_hip)]
                GpuBackend::Hip => taps_hip_run(args.0, args.1, args.2, args.3, args.4, args.5, args.6, args.7, args.8, args.9, args.10, args.11, args.12),
                _ => -1,
            }
        };
        if rc == 0 {
            Ok(nflag)
        } else {
            Err(rc)
        }
    }
}

/// One worker's prepared reads (offsets from 0) and where to send the result.
struct Request {
    rd: Vec<GpuRead>,
    refb: Vec<u8>,
    qry: Vec<u8>,
    qual: Vec<u8>,
    reply: mpsc::SyncSender<Result<Reply, String>>,
}

struct Reply {
    state: Vec<i16>,
    q: Vec<u8>,
    /// Posterior maxima, filled only when some quality was flagged.
    maxv: Vec<f64>,
}

/// The running GPU service.
pub struct Service {
    tx: Mutex<mpsc::Sender<Request>>,
    pub device: String,
    pub backend: GpuBackend,
}

static SERVICE: OnceLock<Service> = OnceLock::new();
static WARNED: AtomicBool = AtomicBool::new(false);

/// Reads pooled per launch before the service stops waiting for more.
const TARGET_READS: usize = 65_536;
/// How long the service waits for more requests after the first one.
const POOL_WAIT: Duration = Duration::from_micros(500);

/// Start the GPU service (once per process). `want`: device name substring.
pub fn start(backend: GpuBackend, want: &str, mem_mb: usize) -> Result<&'static Service> {
    if let Some(s) = SERVICE.get() {
        return Ok(s);
    }
    let order: Vec<GpuBackend> = match backend {
        GpuBackend::Auto => vec![GpuBackend::Cuda, GpuBackend::Hip, GpuBackend::Opencl],
        b => vec![b],
    };
    let (tx, rx) = mpsc::channel::<Request>();
    let (ready_tx, ready_rx) = mpsc::channel::<Result<(GpuBackend, String), String>>();
    let want = want.to_string();
    let f_budget = mem_mb.max(64) * 1024 * 1024;
    std::thread::Builder::new()
        .name("taps-gpu-baq".into())
        .spawn(move || {
            let mut errors = Vec::new();
            let mut engine = None;
            for kind in order {
                match Engine::open(kind, &want, f_budget) {
                    Ok((e, name)) => {
                        engine = Some((e, name));
                        break;
                    }
                    Err(msg) => errors.push(format!("{kind:?}: {msg}")),
                }
            }
            let Some((mut engine, name)) = engine else {
                let _ = ready_tx.send(Err(errors.join("; ")));
                return;
            };
            let _ = ready_tx.send(Ok((engine.kind, name)));
            serve(&mut engine, rx);
        })?;
    let (kind, device) = ready_rx
        .recv()
        .map_err(|_| anyhow::anyhow!("GPU service thread exited"))?
        .map_err(|e| anyhow::anyhow!("no usable GPU for BAQ ({e})"))?;
    let _ = SERVICE.set(Service { tx: Mutex::new(tx), device, backend: kind });
    Ok(SERVICE.get().expect("set above"))
}

/// A pooled launch: requests concatenated (each keeps a contiguous base
/// range), reads sorted for the kernel.
struct Packed {
    rd: Vec<GpuRead>,
    refb: Vec<u8>,
    qry: Vec<u8>,
    qual: Vec<u8>,
    ranges: Vec<(usize, usize)>,
    replies: Vec<mpsc::SyncSender<Result<Reply, String>>>,
    pack_s: f64,
}

/// Packer thread: pools requests from all workers and concatenates them,
/// while the launcher thread runs the previous batch on the device.
fn pack_loop(rx: mpsc::Receiver<Request>, tx: mpsc::SyncSender<Packed>) {
    while let Ok(first) = rx.recv() {
        let mut batch = vec![first];
        let mut n = batch[0].rd.len();
        let deadline = Instant::now() + POOL_WAIT;
        while n < TARGET_READS {
            match rx.recv_timeout(deadline.saturating_duration_since(Instant::now())) {
                Ok(r) => {
                    n += r.rd.len();
                    batch.push(r);
                }
                Err(_) => break,
            }
        }
        let t0 = Instant::now();
        let (mut rd, mut refb, mut qry, mut qual) = (Vec::with_capacity(n), Vec::new(), Vec::new(), Vec::new());
        let mut ranges = Vec::with_capacity(batch.len());
        let mut replies = Vec::with_capacity(batch.len());
        for r in batch {
            let (ro, bo) = (refb.len() as u32, qry.len() as u32);
            rd.extend(r.rd.iter().map(|g| GpuRead { ref_off: g.ref_off + ro, base_off: g.base_off + bo, ..*g }));
            ranges.push((bo as usize, r.qry.len()));
            refb.extend_from_slice(&r.refb);
            qry.extend_from_slice(&r.qry);
            qual.extend_from_slice(&r.qual);
            replies.push(r.reply);
        }
        // Band-7 reads (the kernel's fast path) first, then by length: a SIMD
        // group that mixes fast and generic reads runs both paths.
        rd.sort_by_key(|g| (!(g.bw == 7 && g.i_dim == 51), g.l_query, g.i_dim));
        let pack_s = t0.elapsed().as_secs_f64();
        if tx.send(Packed { rd, refb, qry, qual, ranges, replies, pack_s }).is_err() {
            break;
        }
    }
}

/// Launcher (the thread that opened the device): runs packed batches and
/// returns each request's slice of the results.
fn serve(engine: &mut Engine, rx: mpsc::Receiver<Request>) {
    let mut q2p = [0f32; 256];
    unsafe { taps_baq_gpu_q2p(q2p.as_mut_ptr()) };
    let (ptx, prx) = mpsc::sync_channel::<Packed>(1);
    if std::thread::Builder::new().name("taps-gpu-pack".into()).spawn(move || pack_loop(rx, ptx)).is_err() {
        return;
    }
    // TAPS_GPU_PROFILE=N: cumulative service timings every N launches.
    let profile_every: u64 = std::env::var("TAPS_GPU_PROFILE").ok().and_then(|v| v.parse().ok()).unwrap_or(0);
    let (mut launches, mut reads, mut fast_reads) = (0u64, 0u64, 0u64);
    let (mut t_wait, mut t_pack, mut t_run, mut t_split) = (0f64, 0f64, 0f64, 0f64);
    let mut idle = Instant::now();
    // Result buffers reused across launches (no per-launch allocation).
    let (mut state, mut q, mut maxv): (Vec<i16>, Vec<u8>, Vec<f64>) = (Vec::new(), Vec::new(), Vec::new());
    for p in prx {
        let t_start = Instant::now();
        t_wait += t_start.duration_since(idle).as_secs_f64();
        t_pack += p.pack_s;
        let nb = p.qry.len();
        if state.len() < nb {
            state.resize(nb, 0);
            q.resize(nb, 0);
            maxv.resize(nb, 0.);
        }
        let res = engine.run(&p.rd, &p.refb, &p.qry, &p.qual, &q2p, &mut state[..nb], &mut q[..nb], &mut maxv[..nb]);
        let t_ran = Instant::now();
        launches += 1;
        reads += p.rd.len() as u64;
        if profile_every > 0 {
            fast_reads += p.rd.iter().filter(|g| g.bw == 7 && g.i_dim == 51).count() as u64;
        }
        for (reply, (b0, len)) in p.replies.into_iter().zip(p.ranges) {
            let out = match res {
                Ok(nflag) => Ok(Reply {
                    state: state[b0..b0 + len].to_vec(),
                    q: q[b0..b0 + len].to_vec(),
                    maxv: if nflag > 0 { maxv[b0..b0 + len].to_vec() } else { Vec::new() },
                }),
                Err(rc) => Err(format!("GPU BAQ launch failed (code {rc})")),
            };
            let _ = reply.send(out);
        }
        idle = Instant::now();
        t_run += t_ran.duration_since(t_start).as_secs_f64();
        t_split += idle.duration_since(t_ran).as_secs_f64();
        if profile_every > 0 && launches % profile_every == 0 {
            let mut o = [0f64; 3];
            if engine.kind == GpuBackend::Opencl {
                unsafe { taps_ocl_prof(o.as_mut_ptr()) };
            }
            eprintln!(
                "GPU BAQ profile: {launches} launches, {reads} reads ({:.0}/launch, {:.1}% band-7 fast path); launcher s: waiting {t_wait:.1} run {t_run:.2} split {t_split:.2}; packer {t_pack:.1} s (overlapped); run = upload {:.2} kernel {:.2} download {:.2}; {:.0} reads/s while running",
                reads as f64 / launches as f64,
                100.0 * fast_reads as f64 / reads.max(1) as f64,
                o[0],
                o[1],
                o[2],
                reads as f64 / t_run.max(1e-9)
            );
        }
    }
}

/// A batch handed to the GPU service, to be finished by `complete`.
pub struct Pending {
    reply: Option<mpsc::Receiver<Result<Reply, String>>>,
    n_dev: usize,
    cpu: Vec<c_int>,
    n_cpu: usize,
}

/// Prepare `jobs` (CPU part of BAQ before the HMM) and send the device's
/// share to the GPU service; returns at once. Nothing in the reads changes
/// until `complete`, which must get the same `ws` and `jobs`.
pub fn submit(ws: &mut Workspace, jobs: &[BaqJob]) -> Pending {
    let mut p = Pending { reply: None, n_dev: 0, cpu: Vec::new(), n_cpu: 0 };
    if jobs.is_empty() {
        return p;
    }
    let svc = SERVICE.get().expect("GPU BAQ service not started");
    let mut sz = Sizes::default();
    p.cpu = vec![0 as c_int; jobs.len()];
    unsafe { taps_baq_gpu_prepare(ws.raw(), jobs.as_ptr(), jobs.len() as c_int, &mut sz, p.cpu.as_mut_ptr()) };
    p.n_cpu = sz.n_cpu as usize;
    p.n_dev = sz.n_dev as usize;
    if p.n_dev > 0 {
        let mut rd = vec![GpuRead::default(); p.n_dev];
        let mut refb = vec![0u8; sz.ref_bytes as usize];
        let mut qry = vec![0u8; sz.bases as usize];
        let mut qual = vec![0u8; sz.bases as usize];
        unsafe {
            taps_baq_gpu_export(ws.raw(), jobs.as_ptr(), rd.as_mut_ptr(), refb.as_mut_ptr(), qry.as_mut_ptr(), qual.as_mut_ptr())
        };
        let (reply_tx, reply_rx) = mpsc::sync_channel(1);
        let sent = svc
            .tx
            .lock()
            .map_err(|_| ())
            .and_then(|tx| tx.send(Request { rd, refb, qry, qual, reply: reply_tx }).map_err(|_| ()));
        if sent.is_ok() {
            p.reply = Some(reply_rx);
        }
    }
    p
}

/// Wait for a submitted batch and finish it: the GPU's states and qualities,
/// then BAQ_EXTEND and the quality cap on the CPU; reads the device did not
/// take run on the CPU kernel at `width`. If the device failed, the whole
/// batch runs on the CPU. Returns the number of reads whose HMM ran on the GPU.
pub fn complete(p: Pending, ws: &mut Workspace, jobs: &[BaqJob], width: Width) -> usize {
    if jobs.is_empty() {
        return 0;
    }
    let mut on_gpu = 0;
    if p.n_dev > 0 {
        let reply = p.reply.as_ref().and_then(|rx| rx.recv().ok());
        match reply {
            Some(Ok(r)) => {
                let maxv = if r.maxv.is_empty() { std::ptr::null() } else { r.maxv.as_ptr() };
                unsafe { taps_baq_gpu_finish(ws.raw(), jobs.as_ptr(), r.state.as_ptr(), r.q.as_ptr(), maxv) };
                on_gpu = p.n_dev;
            }
            other => {
                if !WARNED.swap(true, Ordering::Relaxed) {
                    let why = match other {
                        Some(Err(e)) => e,
                        _ => "GPU service unavailable".to_string(),
                    };
                    eprintln!("warning: {why}; BAQ continues on the CPU");
                }
                // Nothing was modified yet: run the whole batch on the CPU.
                return ws.run(jobs, width);
            }
        }
    }
    if p.n_cpu > 0 {
        let sub: Vec<BaqJob> = p.cpu[..p.n_cpu].iter().map(|&i| jobs[i as usize]).collect();
        ws.run(&sub, width);
    }
    on_gpu
}

/// `submit` then `complete` (no overlap).
pub fn run_batch(ws: &mut Workspace, jobs: &[BaqJob], width: Width) -> usize {
    let p = submit(ws, jobs);
    complete(p, ws, jobs, width)
}

/// Batches each worker keeps in flight with `--baq-device gpu`, so it calls
/// one batch while the GPU runs BAQ on the next (TAPS_GPU_DEPTH, default 4).
pub fn depth() -> usize {
    std::env::var("TAPS_GPU_DEPTH").ok().and_then(|v| v.parse().ok()).filter(|&d| d >= 1).unwrap_or(4)
}
