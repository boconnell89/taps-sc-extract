# GPU acceleration plan (v4, 2026-10-17)

Developed by two rounds of adversarial review: Claude drafts, Grok 4.6 (xhigh; measured on this machine's GPU)
and Gemini 3.8 Flash (high; source + literature review), each round cross-examining the other. Working files,
prototype code and raw measurements: `scratch/gpu-plan/` (`plan-v1.md`, `plan-v2.md`, `grok/review-r1.md`,
`grok/review-r2.md`, `gemini/review-r1.md`, `gemini/review-r2.md`; Grok's prototypes and `results*.txt` are in the
review worktree, copied to `scratch/gpu-plan/grok/`).

User decisions: server GPU **Intel Arc Pro B70** (Xe2 "Battlemage", 32 Xe-cores, 32 GB, 22.9 TFLOPS FP32);
CUDA and AMD backends must also exist for future machines; output exact by default, opt-in fast mode allowed if
per-cell differences are measured; **BAQ runs on the CPU by default, `--baq-device gpu` opts in**.

## 0. Intel Arc Pro B70 (the target)
- Intel publishes no FP64 rate for the B70. Xe2 discrete GPUs do have hardware FP64 (Alchemist did not): the Arc
  B580 measured as a strong double-precision device in clpeak, between AMD RX 7900 GRE and XT and ahead of the
  RTX 4070 Ti Super (Phoronix). If Xe2's FP64 rate is 1/16 of FP32 (unconfirmed), the B70 has ~1.4 TFLOPS FP64,
  ~16x this RTX 3060. Scaled from the 3060 measurement that is on the order of 10-18M exact HMM reads/s (estimate,
  unverified), against an estimated ~5M on a 64-thread AVX2 server. **First task on the server: `clinfo` (fp64
  support) and the standalone kernel benchmark.**
- Intel's compute runtime (NEO) exposes OpenCL 3.0 and Level Zero on Linux. The OpenCL C kernel (proven exact here
  through PoCL-CUDA) is therefore the Intel backend; SYCL is not needed.

## 1. Verdict

- **An exact GPU BAQ is possible.** A port of htslib `probaln_glocal` (one read per GPU thread, reads interleaved
  in memory, posterior pass fused into the backward pass, float d/e, no FMA contraction, IEEE division, sequential
  band sums) matched the CPU on 1,048,576 synthetic reads (104M bases): 0 differences in alignment state, posterior
  maximum, or quality. The same kernel in OpenCL C was also exact.
- **TAPS-aware BAQ is what was tested.** The GPU kernel computes the same HMM the CPU's TAPS mode does: Grok's
  synthetic reads were generated against a TAPS-masked reference (C -> N for top-strand reads) with SNPs, indels
  and extra N, and the CPU wrapper keeps doing what is TAPS-specific: the mask (applied while translating the
  reference), the strand choice, and the skip of reads that already match. Not yet tested: real reads (chr19,
  gate G1) and the `taps-partial` / `htslib` modes (same kernel, different inputs).
- **On this workstation it loses, and should stay off.** RTX 3060 (PCIe 3.0 x4), best exact path including
  transfers: 0.91M alignment-model reads/s. The CPU's 32 AVX2 threads do ~2.7M/s. FP32 (inexact) is only 2.2x the
  FP64 kernel (2.38M/s), still below the CPU, and changed the alignment state of 12.7% of reads (state feeds the
  quality: `finish()` zeroes BAQ where the state disagrees). Overlap cannot hide a co-processor slower than the
  work it absorbs.
- **Only strong-FP64 GPUs can win in exact mode**: NVIDIA A100/H100/H200, AMD Instinct MI210/MI250X/MI300X, Intel
  Data Center GPU Max. Scaled from the 3060 measurement these would reach roughly 20-60M HMM reads/s (estimates,
  unverified), against ~2.7M on this CPU and an estimated ~9M on a 96-core AVX-512 server. Inference GPUs (T4, L4,
  A10) and consumer cards are estimated at or below a large CPU; fast mode might rescue some (L4/L40S), unmeasured.
- **Drop GPU HDF5 compression and BGZF inflate.** HDF5 compression on GPU cannot emit the Blosc container (filter
  32001) Amethyst readers need, and writing is already faster than processing; BGZF inflate costs more PCIe time
  than CPU inflate. Both reviewers agree.
- **Pileup/calling: not as htslib's column pileup, but yes as read-centric calling (Phase 4).** See section 6.

So the plan is a portable, exact BAQ co-processor (OpenCL for Intel, CUDA for NVIDIA, HIP for AMD) that is off
unless `--baq-device gpu` is given, followed by read-centric calling that can run BAQ and calling together on the
GPU.

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
- `--baq-device cpu|gpu` with **cpu the default** (user decision). `gpu` uses the first fp64-capable device found
  (OpenCL, CUDA or HIP, whichever backends were built); `--test` reports whether it would help on this host.
  A later `auto` may enable it from the saved machine profile.
- Build: `install.sh` builds a backend only when its compiler is present (`nvcc`, `hipcc`, OpenCL headers + ICD),
  via the existing `cc` build script; otherwise today's CPU build. The binary falls back to the CPU at run time.

## 5. Phase 3 — opt-in fast mode

`--baq-precision fast`: FP32 HMM on the device (and possibly an 8-lane AVX2 float kernel on the CPU). Ships only
with a per-cell CG/CH difference table on chr19 and the full test BAM, and a stated budget (to be agreed). Never
default. It does not rescue this RTX 3060 (FP32 2.38M/s < CPU 2.7M/s).

## 6. Phase 4 — read-centric calling (CPU first, then GPU)

Why not htslib's pileup on the GPU: `bam_plp` walks the genome column by column, keeps a linked list of the reads
covering each column, and decides per column; that is pointer chasing with a data-dependent stream, the opposite
of what GPUs do well. But the result it produces can be computed per read:
- Each read's calls are independent of other reads except for (a) mate-overlap clipping, which only involves the
  two mates of a pair (htslib 1.21 `tweak_overlap_quality`, patched in here), and (b) the `--max-depth` cap, which
  depends on how many reads start at the same position in file order. Both can be resolved before calling: pair
  mates on the CPU and give each pair to one GPU thread; compute the depth-cap drop list from read starts.
- Then one thread per read (or pair) walks its CIGAR, applies trims, base quality, strand, context and barcode,
  and emits (cell, context, position, t/c) records; a GPU radix sort + reduce merges them into per-cell sorted
  site vectors, which is exactly what the HDF5 writers consume.
- Fused with GPU BAQ, a read goes to the device once and only the calls come back: no quality round trip.
- Costs: PCIe volume (~250 B/read up; ~16 B per cell-site down, ~13 GB for CG+CH on the test BAM: ~0.3 s on PCIe
  5.0 x16, ~4 s on this x4 link); CG-only output is 30x smaller. And exact parity with htslib's pileup (overlap
  clipping, depth cap, deletion/skip handling, Rust-equals-Python per cell) must be re-proven.

Order: build read-centric calling **on the CPU first**, proven equal per cell to the current pileup on chr19 and the
full test BAM. Calling plus pileup is 48% of worker time today (more with `--no-baq`), and per-column overhead is
much of it, so this may be a large CPU win by itself. Then port the same algorithm to the GPU behind
`--call-device gpu`, fused with GPU BAQ.

## 7. Not GPU, but found during the review

- Fusing the posterior pass into the backward pass in the CPU AVX2 kernel is exact and halves its working set
  (331 KB to 175 KB per 4-read group) but measured only +3% single-threaded and within noise at 32 threads. Low
  priority.
- After BAQ, the largest CPU cost is the calling callback (36% of worker time). That is the next CPU target, and
  not a GPU problem.

## 8. Open questions

1. B70: fp64 support and rate (`clinfo`), PCIe link, Linux driver/compute-runtime version on the server (G0).
2. Fast-mode budget: how many differing CG cell-sites per million are acceptable?
3. Which AMD GPU to test the HIP build on (none available here; it would ship compiled but untested).
4. Is a persistent GPU context acceptable during `extract` (keeps GPU clocks up between launches)?
