# GPU acceleration plan (v3, 2026-10-17)

Developed by two rounds of adversarial review: Claude drafts, Grok 4.6 (xhigh; measured on this machine's GPU)
and Gemini 3.8 Flash (high; source + literature review), each round cross-examining the other. Working files,
prototype code and raw measurements: `scratch/gpu-plan/` (`plan-v1.md`, `plan-v2.md`, `grok/review-r1.md`,
`grok/review-r2.md`, `gemini/review-r1.md`, `gemini/review-r2.md`; Grok's prototypes and `results*.txt` are in the
review worktree, copied to `scratch/gpu-plan/grok/`).

User decisions: server GPU unknown; output exact by default, opt-in fast mode allowed if per-cell differences are
measured; platforms NVIDIA, AMD and Intel.

## 1. Verdict

- **An exact GPU BAQ is possible.** A port of htslib `probaln_glocal` (one read per GPU thread, reads interleaved
  in memory, posterior pass fused into the backward pass, float d/e, no FMA contraction, IEEE division, sequential
  band sums) matched the CPU on 1,048,576 synthetic reads (104M bases): 0 differences in alignment state, posterior
  maximum, or quality. The same kernel in OpenCL C was also exact.
- **On this workstation it loses, and should stay off.** RTX 3060 (PCIe 3.0 x4), best exact path including
  transfers: 0.91M alignment-model reads/s. The CPU's 32 AVX2 threads do ~2.7M/s. FP32 (inexact) is only 2.2x the
  FP64 kernel (2.38M/s), still below the CPU, and changed the alignment state of 12.7% of reads (state feeds the
  quality: `finish()` zeroes BAQ where the state disagrees). Overlap cannot hide a co-processor slower than the
  work it absorbs.
- **Only strong-FP64 GPUs can win in exact mode**: NVIDIA A100/H100/H200, AMD Instinct MI210/MI250X/MI300X, Intel
  Data Center GPU Max. Scaled from the 3060 measurement these would reach roughly 20-60M HMM reads/s (estimates,
  unverified), against ~2.7M on this CPU and an estimated ~9M on a 96-core AVX-512 server. Inference GPUs (T4, L4,
  A10) and consumer cards are estimated at or below a large CPU; fast mode might rescue some (L4/L40S), unmeasured.
- **Drop all non-BAQ GPU work.** HDF5 compression on GPU cannot emit the Blosc container (filter 32001) Amethyst
  readers need, and writing is already faster than processing; BGZF inflate costs more PCIe time than CPU inflate;
  GPU pileup/calling would be a second pileup engine (overlap clipping, depth cap). Both reviewers agree.

So the plan is a portable, exact, measured-before-enabled BAQ co-processor for servers with the right GPU, and
nothing on this workstation by default.

## 2. Gates

| Gate | Question | How | Stops the project if |
|---|---|---|---|
| G0 | What GPU and PCIe link do the servers have? | `nvidia-smi` / `rocm-smi` / `xpu-smi`; run the standalone kernel benchmark (Phase 1) there | no GPU, or an inference/consumer GPU and fast mode is not wanted |
| G1 | Is the production kernel exact on real data? | `baq-selftest --gpu`: QUAL byte-identical to htslib on all chr19 reads, 200k synthetic, then the whole test BAM | any difference |
| G2 | Is the device faster than the host CPU's BAQ? | standalone benchmark on the target: transfers included, realistic read-length mix | device < 2x the host's CPU BAQ rate |
| G3 | Does the whole run get faster? | integrated run on the target, reads/s end to end, identical output per cell | no measured end-to-end win |

## 3. Phase 1 — standalone exact kernel + benchmark (cheap; enables G0 to G2 anywhere)

- **One kernel body, many backends.** The HMM lives in one plain-C header (like `csrc/baq_kernel.h`) compiled into
  each backend, so exactness rules are written once: CUDA (NVIDIA), HIP from the same source (AMD; `hipcc`, kernel
  unchanged), and OpenCL C (Intel compute GPUs via the NEO runtime; also proven exact on NVIDIA through PoCL).
  Required flags: CUDA `-fmad=false --prec-div=true`; HIP `-ffp-contract=off`; OpenCL `#pragma OPENCL FP_CONTRACT
  OFF` with fp64. Rejected: SYCL/AdaptiveCpp and Kokkos (heavy toolchains at install); OpenCL as the only backend
  (no OpenCL on NVIDIA's WSL driver; PoCL-CUDA needs CUDA's libdevice and has gaps).
- **Layout**: one read per work-item, per-read arrays interleaved (coalesced); posterior fused into the backward
  pass (two rows of b; ~41 KB of f per read).
- **Variable-shape launches (new, untested).** The CPU groups reads by (length, reference span, band) because SIMD
  lanes run in lockstep. On a GPU this fragments real batches into tiny groups: a realistic 80-150 bp mix gave a
  largest group of 418 reads and ~8k reads/s. GPU threads can follow their own loop bounds, so the kernel should
  take per-read shapes in one launch (sorted by length so warps diverge little, f padded per launch). This is the
  main open technical question; if it fails, only the dominant read length (70% at 99 bp here) goes to the GPU and
  every other shape stays on AVX2. A tiny GPU launch is never better than AVX2.
- **Quality on the device, exactly**: compute q with the device log; flag any base whose value is within 1e-9 of
  an integer (q's rounding point) and recompute it on the host with glibc, which the CPU path uses. Measured:
  device log vs glibc gave 0 q differences on 20M inputs; the flag caught every near-boundary input and never
  fired on 104M real-shape HMM bases. (glibc itself differs in the last bit between FMA and non-FMA CPUs on
  1,509 of 20M inputs, with 0 q differences.) fdlibm is not a glibc substitute (9% of results differ in bits).
- **finish() on the device for all CIGARs** (soft clips, insertions, deletions; only all-match was tested), so the
  device returns final qualities only (99 bytes per read instead of 1,188). A read with a flagged base is
  re-finished on the host. Prepare (CIGAR span, band, masked translation) stays on the CPU (~300 B/read upload).
- **Benchmark command** (`taps-sc-extract-rs gpu-bench` or a script): exactness self-test plus HMM reads/s with
  transfers, for the dominant shape and a realistic mix, warm and cold clocks, PCIe gen/width. This is what runs on
  the server for G0/G2, before any integration work.

## 4. Phase 2 — integration (only after G2 passes on a real target)

- One device service thread per GPU; workers submit prepared reads into shared queues; launches when enough reads
  are pooled across workers (one worker's 8,192-record batch yields only ~2,300 HMM reads, too few).
- Workers must not stall: a worker's batch goes to the device only when the device's expected completion beats
  running it locally on AVX2; on a slow device this means never.
- Overlapping calling with device BAQ needs reads to be BAQ'd ahead of the pileup. Two designs to prototype and
  compare: a decoupled BAQ-ahead stage (readers prepare and queue reads across windows; pileup workers consume
  finished batches in file order), or two windows in flight per worker. Records must still enter the pileup in
  file order with final qualities.
- `--baq-device auto|gpu|cpu`: auto enables the GPU only if `--test` measured an end-to-end win on this host
  (saved in the machine profile). On this workstation auto = cpu.
- Build: `install.sh` builds a backend only when its compiler is present (`nvcc`, `hipcc`, OpenCL headers + ICD),
  via the existing `cc` build script; otherwise today's CPU build. The binary falls back to the CPU at run time.

## 5. Phase 3 — opt-in fast mode

`--baq-precision fast`: FP32 HMM on the device (and possibly an 8-lane AVX2 float kernel on the CPU). Ships only
with a per-cell CG/CH difference table on chr19 and the full test BAM, and a stated budget (to be agreed). Never
default. It does not rescue this RTX 3060 (FP32 2.38M/s < CPU 2.7M/s).

## 6. Not GPU, but found during the review

- Fusing the posterior pass into the backward pass in the CPU AVX2 kernel is exact and halves its working set
  (331 KB to 175 KB per 4-read group) but measured only +3% single-threaded and within noise at 32 threads. Low
  priority.
- After BAQ, the largest CPU cost is the calling callback (36% of worker time). That is the next CPU target, and
  not a GPU problem.

## 7. Open questions

1. Server GPU model(s) and PCIe link (G0). Decides everything.
2. Fast-mode budget: how many differing CG cell-sites per million are acceptable?
3. AMD and Intel targets: which exact models, so HIP and OpenCL/NEO can be tested on real hardware (neither was
   available here)?
4. Is a persistent GPU context acceptable during `extract` (keeps GPU clocks up between launches)?
