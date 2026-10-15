# Workstation profiling snapshot — 2026-10-15

As of **2:00am on 10-15-26**. The Rust engine is the uncommitted working tree on top of `026c52e`, built natively with `./install.sh`. It includes the Rust upgrade steps 1, 2, 4a and 5, the Phase 5 change to sorted per-cell vectors, and the barcode-cache fix. Same dataset and machine as [`nvme-profiling-2026-10-14.md`](nvme-profiling-2026-10-14.md):
- **Dataset:** `5base_ucbTn5_sp2.srt.bam`, copied to WSL ext4: 53,461,257 mapped records, 7,351 cells, mm10 canonical contigs.
- **Defaults:** `--trim-r1 0 --trim-r2 10 --min-baseq 20`, LZF, BAQ on.
- **Machine:** AMD Ryzen 9 5900XT, 16 cores / 32 threads, 1 NUMA node, ~60 GB RAM, WSL2. Inputs, temp and output are all on `/dev/sdd` ext4 (VHDX on NVMe). The BAM is page-cached.

`--no-baq` is a test setting; its calls differ from BAQ on. Process = extraction; write = HDF5 shards + `master.h5`. "Writer peak" is the summed peak RSS of the concurrent `write-shard` child processes (an upper bound).

## Before and after (same BAM, `-t 32 --shards 16`)

| | 2026-10-14 (`026c52e`) | 2026-10-15 |
|---|---:|---:|
| Process | 63.5 s | **40.2–42.1 s** |
| HDF5 write | 31.5 s | **11.4–16.9 s** |
| Wall | 95.0 s | **52.7–57.7 s** |
| Reads/s (wall) | 562,565 | ~930,000–1,010,000 |
| Peak RSS | 29.7 GB (extract process) | 11.3–11.5 GB extract + ~10.2 GB writers (memory mode); 4.9 GB + ~10.4 GB (stream mode) |
| `--no-baq`, process | 13.8 s (10 Mb windows) | **7.3 s** |

Cells, CG and CH totals are identical across every BAQ-on row: 21,499,669 / 30,341,601 and 9,625,744 / 814,768,999. Per-cell, the barcode fix moved 0.01% of calls to their correct cell (see `rust-implementation-log.md`).

## Final matrix

| Run | Process | Write | Wall | Extract peak | Writer peak |
|---|---:|---:|---:|---:|---:|
| T=32 memory (a) | 40.2 s | 16.9 s | 57.1 s | 11.5 GB | 10.1 GB |
| T=32 memory (b) | 42.1 s | 15.6 s | 57.7 s | 11.3 GB | 10.2 GB |
| T=32 stream | 41.3 s | 11.4 s | 52.7 s | 4.9 GB | 10.4 GB |
| T=16 memory | 56.7 s | 15.2 s | 71.9 s | 10.5 GB | 10.0 GB |
| T=16 stream | 57.0 s | 10.4 s | 67.4 s | 4.7 GB | 10.2 GB |
| T=32 memory `--no-baq` | 7.3 s | 14.1 s | 21.4 s | 10.7 GB | 10.1 GB |
| T=32, **1 shard** (auto for <10k cells) | 41.0 s | **61.9 s** | 102.9 s | 9.1 GB | — |

Write time ranged 10.4–16.9 s across identical 16-shard runs. Three back-to-back repeats with nothing else running gave 13.3, 14.1 and 14.7 s. Never time while compiling or running other jobs.

## What each step changed (T=32, memory, 16 shards; process seconds unless noted)

| Step | Change | Effect |
|---|---|---|
| 0 | 10 Mb windows instead of 100 Mb (37 windows) | 63.5 → 54.1 s; calls unchanged |
| 1+2 | libdeflate, MAPQ/orphan before BAQ, shared BAQ contigs, decomp 0, auto windows (14.2 Mb, 203 windows) | → 50–55 s |
| 4a | One `write-shard` process per shard; cell-at-a-time sliced writes | Write 31.5–41 s → 11–17 s; extract peak 30 → 11.5 GB; writer peak 28 → 10 GB |
| 5 | Native build (`-march=native`, also htslib's BAQ) | ~7% at T=32, ~2% at T=16 |
| fix | Barcode cache keyed on `bam1_t.id` | → 45.0 s (fewer barcode parses) |
| P5 | Sorted per-cell `Vec` instead of `FxHashMap` per call | → 40.4 s; `--no-baq` 16.6 → 10.0 s (7.3 s final) |
| rejected | BAQ only on callable reads | No measurable gain (99% of reads are callable) and 2 CH calls changed |
| skipped | Shared decompression pool (Step 3) | Inflate is 1–6% of worker time |

## Where the time goes (`--profile-json`, T=32, before the P5 change)

Worker-thread seconds per stage, summed over workers:

| | read+inflate | BAQ | pileup | calling |
|---|---:|---:|---:|---:|
| BAQ on | 39 (2.4%) | 1,081 (67%) | 152 (9.5%) | 327 (20%) |
| `--no-baq` | 28 (6%) | 2 | 131 (27%) | 310 (65%) |
| `--no-baq`, after P5 | 27 | 2 | 133 | **126** |

- Extraction CPU utilisation was 0.98 with BAQ on and 1.00 with `--no-baq`.
- The window tail is small: with BAQ on, ideal 50.2 s against 54.7 s actual.
- The slowest window was chr9:0–14 Mb at 16.9 s.

## Limiting factors now

- **BAQ on (production):** CPU-bound in htslib's `sam_prob_realn` (about 70% of worker time). Two levers remain: fewer HMM calls, which is not bit-exact, or a faster HMM. The pileup and the calling callback are now about equal and together cost under 30%.
- **HDF5 write:** CPU-bound in the 16 writer processes (utilisation ~0.8); the disk is not the cap. One shard is 4–5× slower, so use `--shards` ≥ the number of writers you want.
- **Thread scaling on this machine:** T=16 → T=32 gives 1.4× processing (SMT). `taps-sc-extract --test` measures the same curve on synthetic data: W=16 → W=32 gives 1.40×.
- **Memory:** stream mode has the smallest tree (about 5 GB + 10 GB of writers) and is not slower than memory mode on local NVMe.
- **BAM read:** not limiting from local NVMe (about 90 MB/s of BAM consumed). The page-cached BAM does not test cold reads. `--test` reports when the read probe was served from cache.

## `--test` against the real run

`taps-sc-extract --test-quick -b … -f …` on this machine (about 40 s):

| | `--test` | Real run |
|---|---|---|
| Sites | 709M | 701.9M |
| Memory, memory mode | 23 GB | 21.7 GB (11.5 + 10.2) |
| Processing | 33 s before calibration | 41 s |
| Write | 9 s before calibration | 14 s |

The processing and write gaps were then fitted as factors of 1.25 and 1.5 to this same run (noted in `tune.rs`). The predictions match this dataset by construction and have not yet been checked against a run on another dataset.

Recommended flags matched the best rows: `-t 32 --decomp-threads 0 --writer-mode process --max-writer-threads 16`.

Logs: `/home/oconnelb/taps_io_test/logs/fin_*.log`, `step0_*.log`, `prof_*.json`, `wr_rep*.log`.
