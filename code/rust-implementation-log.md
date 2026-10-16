# Rust implementation log

Source of truth for the design: [`rust-acceleration-plan.md`](rust-acceleration-plan.md).

## 2026-08-20 — Plan saved; crate skeleton started

- Copied the approved acceleration plan into the repo.
- Plan updates from review:
  - `--decomp-threads` is user-settable at run time. If omitted, use more BGZF threads when fewer workers (`4` if ≤4 workers, `2` if 5–7, `1` if ≥8).
  - `--chunk-size-mb` default stays **10**. The user may override. No automatic shrink at high worker counts.
- Scaffolded `taps-sc-extract-rs/` with calling, context, barcode, FASTA `.fai` reader, window planner, and autotune helpers (PR1–PR2 scope). Python `taps_sc_extract/` is unchanged.
- PR0 (rust-htslib vs Python chr19 pileup) not run yet: needs a rust-htslib/htslib build on this host.

## 2026-08-20 — rust-htslib builds; PR0 region spike

- `rust-htslib` 0.49 compiles (bundled htslib + bindgen). Needs conda `clang`/`libclang` (`LIBCLANG_PATH`).
- `taps-sc-extract-rs stats` does column-wise window extract (MAPQ, baseq, OT/OB flags, approximate mate-overlap).
- BAQ (`compute_baq`) not applied yet.

### Spike results (PR0) — `taps_ucbTn5_sp2.srt.bam` chr19

| Region | Engine | Cells | CG c/tot | CH c/tot |
|---|---|---|---|---|
| 5.0–5.1 Mb | Python (BAQ on) | 192 | 1120/2856 | 120/41340 |
| 5.0–5.1 Mb | Rust | 192 | 1120/2858 | 120/41351 |
| 5–6 Mb | Python | 545 | 12523/23842 | 1132/376631 |
| 5–6 Mb | Rust | 545 | 12523/23840 | 1132/376649 |

- **Methylated counts (`c`) and cell counts match exactly.**
- Unmethylated (`t`) differs by <0.1% (BAQ + overlap-clip approximation).
- **Go:** continue with rust-htslib pileup; add `sam_prob_realn` BAQ next to close the `t` gap.
- noodles: not needed; rust-htslib is the pileup engine.

## Spike results (PR0)

- rust-htslib vs Python call counts: **GO** (exact `c` and cells; tiny `t` delta without BAQ)
- noodles vs Python: **same as rust-htslib** (exact `c` and cells; same tiny `t` delta). Noodles CIGAR walk matches rust-htslib bit-for-bit on the tested windows; debug build ~2.5× slower. Incorporated as `--pileup noodles`. Default production engine remains rust-htslib (faster; BAQ via `sam_prob_realn` is available there).

## 2026-08-20 — rust-htslib `bam_mplp` + BAQ (`sam_prob_realn`)

Replaced `IndexedReader::pileup()` (`bam_plp`, no overlaps, no BAQ) with a samtools-stepper `bam_mplp`:

- flag filter UNMAP|SECONDARY|QCFAIL|DUP
- `ignore_orphans` (paired but not proper)
- `bam_mplp_init_overlaps` (not the earlier per-QNAME higher-qual clip)
- optional `sam_prob_realn(b, full_contig_from_0, len, 3)` (`BAQ_APPLY|BAQ_EXTEND`)
- CLI: `--no-baq`, `--no-overlap-clip`, `--no-ignore-orphans`
- Python spike script: `--fasta-baq` passes `pysam.FastaFile` so BAQ actually runs

### Important: Python `compute_baq=True` is currently a no-op

`taps_sc_extract` sets `compute_baq=True` but never passes `fastafile=` to `bam.pileup()`. pysam only calls `sam_prob_realn` when a FASTA is attached. Confirmed: Python with `--fasta-baq` changes counts; `--fasta-baq --no-baq` matches production.

For HDF5 A/B vs **today’s Python**, run Rust `--no-baq`. Rust default BAQ-on matches pysam *with* FASTA, not current production. Do **not** silently change Python to pass `fastafile` without an explicit decision (it would shift chr19 CG ~64.37% → ~63.75%).

### Parity vs Python — `taps_ucbTn5_sp2.srt.bam` chr19

| Region | Mode | Engine | Cells | CG c/tot | CH c/tot |
|---|---|---|---|---|---|
| 5.0–5.1 Mb | no BAQ | Python production | 192 | 1120/2856 | 120/41340 |
| 5.0–5.1 Mb | no BAQ | Rust `--no-baq` | 192 | **1120/2856** | 120/41343 |
| 5.0–5.1 Mb | BAQ | Python `--fasta-baq` | 192 | 1079/2815 | 114/41181 |
| 5.0–5.1 Mb | BAQ | Rust default | 192 | **1079/2815** | 114/41183 |
| 5–6 Mb | no BAQ | Python | 545 | 12523/23842 | 1132/376631 |
| 5–6 Mb | no BAQ | Rust `--no-baq` | 545 | **12523/23842** | 1132/376636 |
| 5–6 Mb | BAQ | Python `--fasta-baq` | 545 | 12135/23428 | 1043/375208 |
| 5–6 Mb | BAQ | Rust default | 545 | **12135/23428** | 1043/375211 |
| full chr19, 7×10 Mb | no BAQ | Python (sum of windows) | — | 579151/899696 | 63693/21378196 |
| full chr19, 7×10 Mb | no BAQ | Rust `--no-baq` | **5781** | 579179/899731 | 63705/21378519 |

- Cells match on every window.
- CG is exact on the 100 kb and 1 Mb spikes; full chr19 CG tot **+35 / 899696 (+0.004%)**.
- Residual CH tot **+323 / 21.4M (+0.0015%)** on full chr19 (~2–5 extra `t` per 1 Mb). Same direction with BAQ on. Not leftmost-read ownership (column-wise windows).
- Debug wall time, full chr19: Python ~120 s vs Rust `--no-baq` ~31 s vs Rust BAQ-on ~45 s (single thread, unoptimized).

noodles on 5.0–5.1 Mb `--no-baq` is still the old CIGAR-walk counts (1120/2858, 120/41351): no `bam_mplp_init_overlaps`, no BAQ.

### Next

- Residual CH `t` (3 sites on 100 kb) still open; overlap_push not pairing some mates. Not a Python bug — real htslib clip we want. BAQ-on is the intended default (Python’s `compute_baq=True` without `fastafile` is a no-op; we do **not** copy that).
- PR6: HDF5 from compact temp chunks.

## 2026-08-20 — PR4: Rayon window pool

- `parallel.rs`: Rayon pool, **thread-local `IndexedReader`** (never sent across threads), `Arc<FastFaiReader>` (each fetch opens its own `File`).
- Interners/maps stay per-window; shard files unique by `(shard, chunk_id)`.
- `Arc<AtomicBool>` cancel; SIGINT via `ctrlc` sets it; workers load-only.
- CLI: `-t/--workers` (0 = `min(nproc, 32)`). `--decomp-threads` omitted → 4/2/1 by worker count. BAQ **on** by default.
- Send/Sync audit test: params, Window, FastFaiReader, CancelFlag. `IndexedReader` is thread-local.

### Measured (debug, BAQ on, `taps_ucbTn5_sp2.srt.bam` chr19)

| Run | Workers | Cells | CG c/tot | CH c/tot | elapsed |
|---|---|---|---|---|---|
| 5–6 Mb | 1 | 545 | 12135/23428 | 1043/375211 | 1.56 s |
| 5–6 Mb | 4 (1 window) | 545 | same | same | 1.52 s |
| full chr19, 7×10 Mb | 8 | 5781 | 561520/880793 | 59075/21283380 | **10.0 s** |
| full chr19 serial (prior) | 1 | 5781 | same | same | 45.4 s |

8-way Rayon matches serial BAQ-on counts exactly (~4.5× wall-clock on 7 windows). Extract `--workers 4 --shards 8` writes 8 chunk files.

## 2026-08-20 — PR5/PR6: memory mode + Amethyst HDF5

- `hdf5-metno` 0.14, gzip deflate level 1 or none. `HDF5_DIR` for conda libhdf5.
- `--memory-mode stream` (temp chunks, then assemble) or `memory` (keep window maps).
- Writer pool cap 6. Single-pass `create` with `chunk = min(len, 65536)`.
- `n_shards==1`: one `.h5` file. `n_shards>1`: `shard_XXX.h5` + `master.h5` relative ExternalLinks (`CG/<bc>` → `shard_XXX.h5:/CG/<bc>`).
- TempDir for stream chunks; `--keep-temp` retains it. Dropped after assemble otherwise.
- Dtype via h5py: `('chr','S10'),('pos','<i8'),('pct','<f8'),('t','<i8'),('c','<i8')`, `t` before `c`, 1-based pos, `metadata/version = amethyst2.0.0`.

### Smoke (BAQ on, chr19 5.0–5.1 Mb)

Stream `--shards 2`: 192 cells, CG 1079/2815, CH 114/41183; `master.h5` + two shards; ExternalLinks resolve in h5py. Memory `--shards 1` same counts, single file.

## 2026-08-20 — PR3: intern, site maps, compact temps, residual isolated

- `accumulate.rs`: barcode intern (no per-base `String`), per-cell CG/CH `pos→(t,c)` (1-based).
- `stats --dump-sites TSV` and Python `--dump-sites` for parity diffs.
- `extract -o dir --shards N --no-baq` writes `shard_XXX/chunk_YYYYYY.bin` (`TAPSCK01`). 100 kb chr19 → 8/8 shards, 192 cells.
- Hot loop: intern by QNAME bytes; maps only when `accumulate` (stats-only skips them).

### Residual CH `t` on chr19 5.0–5.1 Mb (`--no-baq`)

Site TSV: **41247/41247 sites identical keys**; **3 value mismatches**, all extra Rust `t` (+3 total):

| barcode | ctx | pos (1-based) | Python (t,c) | Rust (t,c) |
|---|---|---|---|---|
| ATTGGCTCAATGATCCTGGCCTCGGTCAAT | CH | 5063986 | (3,0) | (4,0) |
| ATTGGCTCAATGATCCTGGCCTCGGTCAAT | CH | 5063987 | (3,0) | (4,0) |
| CTGGCTTAGTTTCCATTCTTTGGCCGCAAT | CH | 5073239 | (1,0) | (2,0) |

The third site is a proper pair (flags 99/147, same QNAME). Python `bam_mplp_init_overlaps` sums quals (24+40=64) and zeros the mate. Rust still sees both bases at 24 and 40. `bam_plp_init_overlaps` returns success and the overlap hash is non-NULL; `overlap_push` is still not pairing these mates (same start pos 5073215, isize ±78). Not leftmost-read ownership. Next: make `overlap_push` actually tweak this pair (or match pysam’s mpileup path bit-for-bit).

## 2026-08-20 — Full chr19 extract timings + Amethyst check

BAM `260731/taps_ucbTn5_sp2.srt.bam`, FASTA mm10, `-c chr19`, 7×10 Mb windows, BAQ on, `--workers 24 --shards 8 --compression gzip --max-writer-threads 6`. Release binary.

| Mode | elapsed_s (tool) | `/usr/bin/time` wall | user/sys | max RSS | output |
|---|---|---|---|---|---|
| `stream` | **13.21 s** | 13.84 s | 31.77 / 6.35 | 1.51 GiB | `/tmp/taps_chr19_stream/` 122 MB |
| `memory` | **13.24 s** | 13.27 s | 31.67 / 5.92 | 1.82 GiB | `/tmp/taps_chr19_memory/` 124 MB |

Counts identical in both modes: **cells=5781 CG=63.7516% (561520/880793) CH=0.2776% (59075/21283380)**.

HDF5 (h5py): dtype `chr S10 / pos i8 / pct f8 / t i8 / c i8` (`t` before `c`), gzip, 1-based pos, `master.h5` relative `ExternalLink`s to `shard_XXX.h5`. 5025 CG barcodes, 5781 CH barcodes (756 CH-only). Stream vs memory: same barcode sets, same CG/CH row counts (832201 / 20062257), per-cell arrays equal.

Amethyst in conda env **`amethyst_r`** only (`/home/oconnelb/miniforge3/envs/amethyst_r`, amethyst 1.0.5, rhdf5 2.54.1, R 4.5.3). Other envs have no Amethyst. On stream `master.h5`:

- `createObject` instant (5025 CG barcodes)
- `indexChr(obj, "CG", chrList="chr19")` **189.57 s** — chr19, 5025 cells, 832201 sites
- `indexChr(obj, "CH", chrList="chr19")` **194.57 s** — chr19, 5025 cells, 20031676 sites
- Printed `AMETHYST_CHR19_OK`

CH `indexChr` used the CG barcode list (Amethyst `h5paths` from `/CG`); CH-only barcodes are present in the shards but not in that object.

## 2026-08-20 — PR7 & PR8: Auto-tuning, memory flags, and Python CLI delegation

- **PR7 Auto-Tuning**:
  - Implemented `/proc/meminfo` parsing (`get_available_memory_gb()`) with conservative 0.6× multiplier (`system_memory_budget_gb`).
  - Implemented cell count auto-detection: explicit `--expected-cells` > whitelist cardinality > 10,000 default estimate.
  - Implemented `--memory-mode auto` heuristic: dynamically selects `memory` mode if `expected_cells * genome_sites * 40B * 2.5 < 0.4 * budget_gb`, otherwise selects `stream` mode.
  - Implemented automatic shard count selection (1 / 8 / 16 / 32) when `--shards 0` or omitted.
  - Fixed `open_bam` decompression thread setting (`decomp_threads >= 1` enables background BGZF threads via `bam.set_threads`).
- **PR8 Python CLI Integration & A/B Harness**:
  - Added `--engine auto|rust|python` in Python `taps_sc_extract.cli`. Default is `auto` (auto-detects `taps-sc-extract-rs` binary and delegates transparently).
  - Added automated A/B benchmark script `scripts/benchmark_ab.py`.
  - Tagged `python-baseline-20260820` on Git main branch.

### A/B Benchmark Results (`taps_ucbTn5_sp2.srt.bam`, mm10, 24 workers, 16 shards)

| Dataset / Contigs | Python Baseline | Rust Core (v0.1.0 Initial) | Rust Core (Optimized) | Total Speedup | Peak RAM (Rust) | Parity / Calls |
|---|---|---|---|---|---|---|
| `chr19` (7 windows, 2.1M reads) | 39.04 s | 13.21 s | **9.89 s** | **3.95× faster** | 1.52 GiB | 5,781 cells, 63.75% mCG, 0.278% mCH |
| `chr1–3` (56 windows, 16.6M reads) | 176.3 s | 74.02 s | **67.25 s** | **2.62× faster** | 7.93 GiB | 7,268 cells, 65.86% mCG, 0.275% mCH |
| **Whole Genome mm10** (286 windows, 74.9M reads) | **12.70 min** (762 s) | 5.68 min (340.88 s) | **4.86 min (292.12 s)** | **2.61× faster** | 31.5 GiB | **7,355 cells, 38.7M CpG calls (65.36% mCG), 1.075B CH calls (0.278% mCH)** |

### Performance Optimization Summary
1. **Thread Balancing**: `--decomp-threads 0` (synchronous BAM reader per worker) achieves highest throughput with $\ge 8$ workers because 24–32 CPU cores are already saturated without context-switch latency.
2. **`rustc_hash::FxHashMap`**: Replacing SipHash with fast non-cryptographic word hashing saves ~70 seconds of CPU core time over 1 billion pileup calls.
3. **Phase 1 Shard Writer Scaling**: Auto-scaling shard writer threads up to 16 compresses all 16 shards concurrently, cutting whole-genome tail assembly time significantly.
4. **Phase 2 Hot-Loop MRU Barcode Caching**: Avoids redundant ASCII string scanning and intern hash queries for consecutive pileup columns on the same read (>90% cache hit rate).
5. **Phase 3 Compression Filters**: Supports `gzip` (level 1 deflate), `gzip-shuffle` (byte shuffling + level 1 deflate), `gzip6`, and `none`.

## 2026-10-14 9:47pm — Workstation remeasure (`026c52e`)

Full numbers: [`nvme-profiling-2026-10-14.md`](nvme-profiling-2026-10-14.md). Dataset `5base_ucbTn5_sp2.srt.bam`, 53.5M mapped, 7,351 cells, `-t 32 --decomp-threads 0 --engine rust`, LZF.

- E: 9p HDD, BAQ on, 100 Mb, memory, 16 shards: **171.2 s** wall (80.0 process / 91.3 write), 312,209 reads/s.
- WSL NVMe ext4, same settings: **95.0 s** (63.5 / 31.5), 562,565 reads/s. `/dev/sdd` ~1% util; BAM was page-cached.
- WSL NVMe, `--no-baq --chunk-size-mb 10` (test): **45.3 s** (13.8 / 31.6), 1,178,919 reads/s. HDF5 write did not move. `--no-baq` changes CG/CH rates (more bases pass Q20); not a production default.
- Process limiter with BAQ on: `sam_prob_realn` plus 37 coarse windows. After `--no-baq` and 10 Mb windows, HDF5 (~14.7k LZF datasets) is ~70% of wall.




## 2026-10-15 — Fixed: Rust put a small share of calls under the wrong barcode

**Affects every Rust-engine HDF5 written before this fix, including `026c52e` and all `260731/` Rust extracts. Python-engine output is not affected.**

- **Cause.** The callback cached the read's barcode keyed on the QNAME pointer and length. htslib recycles pileup nodes LIFO, and `bam_copy1` reuses the data buffer. A read pushed after the cached read was freed could land at the same address, with a QNAME of the same length, and take the old read's barcode.
- **Fix.** Key the cache on `bam1_t.id`. htslib numbers every read it stores in the iterator 0, 1, 2, … and never reuses a number within a window.
- **Size** (chr19, BAQ on, against `test_chr19_rust_026c52e`):
  - 79 of 837,906 CG calls (0.009%) and 2,914 of 19,956,092 CH calls (0.015%) moved to their correct barcode.
  - Every per-position total (t and c summed over cells) is unchanged, so global and per-region methylation percentages are unaffected. Per-cell analyses saw the small misattribution.
- **Against Python** (`--no-baq`, chr19):
  - Before the fix, 86 CG and 2,756 CH positions had calls under a different cell with the same position total.
  - After it, zero.
  - A separate gap remains: 31 CG and 309 CH positions where Rust has more calls than Python. Its cause is still open (see the upgrade plan, Parity).
- **Tests.**
  - New frozen reference: `test_chr19_rust_cellfix_2026-10-15/`.
  - `test_chr19_rust_python_same_cell_per_call` now fails if any call changes barcode between engines. The frozen-reference test alone could not catch this class of bug.

## 2026-10-15 — Rust upgrade steps 0–5 and Phase 5 (working tree on `026c52e`)

Full numbers: [`nvme-profiling-2026-10-15.md`](nvme-profiling-2026-10-15.md). Work order: `taps-sc-extract-rust-upgrade-plan.md`, as revised in the approved plan.

- **Headline, T=32, 16 shards, NVMe:** wall 95.0 s → 52.7–57.7 s. Processing 63.5 → 40–42 s; HDF5 write 31.5 → 11–17 s; extraction-process peak 29.7 → 4.9 GB (stream) or 11.5 GB (memory). All calls are identical except the barcode fix above.
- **Step 0:** 10 Mb windows took 63.5 s → 54.1 s with BAQ on and 22.0 → 13.8 s with `--no-baq`. Totals were unchanged. BAQ is about 75% of processing.
- **Step 1:**
  - libdeflate on.
  - MAPQ/orphan tests run before `sam_prob_realn`.
  - BAQ contigs shared (`BaqRefs`, released after a contig's last window).
  - Rust decomp auto is 0 at 8 or more workers; the Python CLI no longer forces 1.
  - `default_workers` no longer capped at 64.
  - Lifting the 6-writer clamp made in-process writing worse (16 writers: 86 s, 47 GB), so auto in-process writers stay capped at 6 (`MAX_AUTO_WRITERS`).
- **Step 2:** `--chunk-size-mb 0` (default) = about 6 windows per worker, 2–100 Mb. Per-cell output is identical at 100, 10 and 2 Mb (chr19).
- **Step 3 (shared BGZF pool):** not built. Inflate is 1–6% of worker time, and rust-htslib's `ThreadPool` is `Arc<RefCell>` (not `Sync`).
- **Step 4a:**
  - hdf5-metno's process-wide lock serializes all in-process HDF5 calls, so each shard now gets its own `write-shard` child process.
  - Memory mode pipes `TAPSCK01` chunks to the child's stdin; stream mode passes the temp directory.
  - Cells are written one at a time in 262,144-record slices.
  - Each child reports its peak RSS, and the summary prints the writer total.
  - Step 4b (direct chunk write) is not needed.
- **Step 5:**
  - `install.sh` writes `.cargo/config.toml` from `config.toml.in`, holding the conda prefix, `target-cpu`, and a forced `CFLAGS="-O3 -march=…"`. The conda env's own `CFLAGS` (`-march=nocona -O2`) otherwise won over native.
  - About 7% gain at T=32.
- **Phase 5:** per-cell sites became sorted `Vec<(pos, t, c)>` filled in column order (no hashing). Processing 45.0 → 40.4 s; `--no-baq` 16.6 → 10.0 s. Callback cost fell 310 → 126 worker-seconds (no-BAQ).
- **BAQ-gating experiment** (`sam_prob_realn` only on 99/147/83/163 reads): no gain, because 99% of records are callable, and 2 CH calls changed. Removed.
- **`--profile-json`:** per-phase CPU and utilisation, per-stage worker time (read, BAQ, pileup, calling), per-window timings, writer memory. Adds about 5–10% to processing.
- **`--test` (`tune`):**
  - Synthetic TAPS BAM → per-read stage cost, worker scaling, memory bandwidth, storage per filesystem, HDF5 writers in-process vs per process.
  - Optional dataset sampling with time and memory prediction.
  - `--save-profile` makes the recommendations the defaults on this host.
- **Open:**
  - Auto shards is 1 below 10k cells, and one shard writes 4–5× slower (61.9 s against 11–16 s).
  - (Resolved later the same day: the parity gap was htslib 1.19.1 vs 1.21 overlap clipping. See below.)
  - The 96-core server and the NUMA check have not been run.
  - (Fixed later the same day: memory mode and writer count now account for writer-child memory. See below.)
  - `--compression-threads` is accepted by the Python CLI but ignored by Rust.
  - The `--test` prediction factors (1.25 processing, 1.5 write) were fitted to this one run and are not yet checked on another dataset.

## 2026-10-15 — Memory plan covers the writer processes

- `autotune::plan_memory` picks memory vs stream mode and the writer count together, so the extraction process plus the concurrent `write-shard` children fit the budget.
- Inputs: mapped reads from the BAM index (about 14 sites per read), memory-mode extraction at 18 B/site, stream extraction at 1 GB + 0.15 GB per worker, and each writer at 20 B/site of its shard plus 0.1 GB. All are measured on 2026-10-15, with margin added.
- New `--max-memory-fraction` (default 0.6 of MemAvailable); `--max-memory-gb` still overrides it. The plan prints in the config block.
- Full mm10 test BAM:

  | Budget | Chosen | Estimated peak | Measured |
  |---|---|---|---|
  | 34.2 GB | memory, 16 writers | 31.0 GB | 11.3 + 10.4 = 21.6 GB |
  | 0.3 → 16.8 GB | stream, 15 writers | 16.5 GB | 4.5 + 10.1 = 14.6 GB |

  Run time was unchanged in both.
- Amethyst 1.0.5 (`amethyst_r` env): `createObject` and `indexChr` on a 4-shard chr19 output from process writers ran without per-barcode errors. CG cells and sites match h5py exactly (4,977 / 792,338).

## 2026-10-15 — Parity with the Python engine: htslib 1.21 overlap clipping

- **Symptom.** With `--no-baq` on chr19, Rust had extra calls at 31 CG and 309 CH positions (343 cell-positions, nearly all +1 `t`). At every one, the read's mate overlapped the position but had a **deletion** there.
- **Cause.** pysam 0.23.3 bundles htslib 1.21, while rust-htslib 0.49 / hts-sys 2.2.1 bundle htslib 1.19.1 (still the newest on crates.io).
  - htslib 1.21's `tweak_overlap_quality` also down-weights a mate's bases where the other mate has a deletion inside the overlap: ×0.8 on one read, 0 on the other, chosen by QNAME hash.
  - 1.19.1 skipped those positions ("don't know what to do with indels"), so Rust kept a Q≥20 base that Python drops.
  - `overlap_push` is identical in both versions. The earlier equal-start and supplementary-read hypotheses were wrong.
- **Fix.**
  - `taps-sc-extract-rs/patches/hts-sys-2.2.1-overlap-deletions.patch` replaces that one function with htslib 1.21's.
  - `install.sh` downloads `hts-sys-2.2.1.crate`, checks its crates.io sha256, extracts it into `vendor/` (git-ignored) and applies the patch.
  - `Cargo.toml` `[patch.crates-io]` points hts-sys at the vendored copy. A plain `cargo build` now needs `install.sh` to have run once.
- **Result.**
  - chr19 `--no-baq`: Rust equals Python at every cell and position (856,147 CG, 20,043,609 CH calls). `test_chr19_rust_matches_python_per_cell` now requires exact equality.
  - BAQ on, chr19: −18 CG / −167 CH calls. New frozen reference `test_chr19_rust_overlapfix_2026-10-15/`.
  - BAQ on, full mm10 test BAM: CG 30,341,601 → 30,341,182 (−419), CH 814,768,999 → 814,759,288 (−9,711). Run time is unchanged.

## 2026-10-15 — `--compression-threads` in the Rust engine

- New flag, forwarded by the Python CLI. It sets Blosc's thread count in each `write-shard` child, or in-process for thread mode or one shard. It applies only to `blosc` and `blosc-zstd`; LZF and gzip run single-threaded inside HDF5.
- Measured at T=32, stream mode, full mm10 test BAM (HDF5 write seconds):

  | Codec | Shards | 1 thread | 2 threads (CPUs/writers) | 4 threads | 32 threads |
  |---|---|---|---|---|---|
  | Blosc LZ4 L9 | 16 | 9.9 | 10.4 | 10.4 | — |
  | Blosc LZ4 L9 | 1 | 54.0 | — | — | 53.2 |
  | Blosc-Zstd | 16 | 24.2 | 23.2 | — | — |
  | Blosc-Zstd | 1 | — | — | — | 70.9 |

  Output sizes: LZ4 1.4 GB, Zstd 976 MB (LZF 4.2 GB).
- **No gain**, so the default stays 1. Most cell datasets are a single small chunk, and write time is per-dataset HDF5 work plus shard parallelism, not compression throughput inside a chunk.

## 2026-10-15 — `install.sh` re-patches hts-sys when the patch changes
- `vendor/hts-sys-2.2.1/.taps_patched` now holds the sha256 of the patch it was built with. A rerun re-extracts and re-patches a clean copy whenever `patches/hts-sys-2.2.1-overlap-deletions.patch` differs, so an upgrade never builds a stale patch. Checkouts with the earlier empty marker re-patch once.
- The verified `hts-sys-2.2.1.crate` is kept in `vendor/` and reused, so re-patching needs no network. A corrupt cached crate is downloaded again and checked against the crates.io sha256.
- Checked: old marker → re-vendor, build and chr19 smoke run; unchanged patch → skip; changed patch with network blocked → re-patch from cache; corrupt cache offline → clean failure; corrupt cache online → re-download.
- README has a new "Upgrading an existing install" section.

## 2026-10-15 — Complexity cutoffs (`--complexity`, `--min-complexity`, `--max-complexity`)
- Both engines read a complexity file (`rank barcode total_reads unique_reads percent`, whitespace-separated, optional header, `.gz` accepted) and keep barcodes whose unique-read count (column 4) is within the inclusive bounds. With `--whitelist` too, a cell must pass both. A warning is printed when no read matches the filter (for example a barcode format that differs from the QNAME prefix).
- The filter acts at call time, the same place as the whitelist, so kept cells' calls equal an unfiltered run. Rust also skips BAQ for reads of filtered cells. They are still pushed into the pileup, so depth caps and overlap pairing are unchanged.
- Tests: Rust and Python unit tests (bounds, header, short rows, `.gz`, intersection); chr19 BAQ-on filtered run equals the frozen reference on kept cells; chr19 `--no-baq` filtered Rust equals Python on kept cells.

## 2026-10-15 — Write path for few, deep cells; input handling; `--test` additions

Prompted by a 64-core run (454 cells, 56.5M reads, 422M sites): 2.9 s HDF5 write with LZF, about 9 cores busy while writing, and 64 writers no faster than 32. Plan reviewed by Grok 4.6 and Gemini 3.8 Flash; claims below were checked against the source and by measurement.

**Test data.** `bigcell_450.bam`: the workstation BAM with its 7,351 barcodes merged round-robin into 450 cells, about 1.9M sites per cell on average. The largest merged cell has 37.8M sites. Mates share a QNAME, so pairs stay together. Calls differ from the original by 240 CG positions because overlap clipping picks the mate by QNAME hash. Rebuild it with:
  ```bash
  samtools view -h -@4 5base_ucbTn5_sp2.srt.bam | mawk 'BEGIN{OFS="\t"} /^@/{print;next} {i=index($1,":"); bc=substr($1,1,i-1); if(!(bc in m)){m[bc]=sprintf("CELL%03d", n%450); n++} $1=m[bc] substr($1,i); print}' | samtools view -b -@8 -o bigcell_450.bam - && samtools index bigcell_450.bam
  ```

Balanced and hash-sharded outputs are identical per cell (701.7M sites over 900 datasets), and Amethyst loads the balanced `master.h5`.

**Writer instrumentation.** Each shard reports cells, sites, read / create / fill / HDF5-write seconds, bytes and peak RSS. The run prints a `Shard writers:` line (sites per shard median and max, slowest shard, time split, `master.h5` time), and `--profile-json` has the per-shard fields and `master_h5_s`.

**Measured on the 450-cell BAM** (T=32, memory mode, LZF, 16-core workstation):

| | 32 shards | 16 shards |
|---|---:|---:|
| Before: write s, max/median sites per shard | 11.9 s, 2.60× | 13.9 s, 2.14× |
| Writer time: read / fill / HDF5 write | 37% / 11% / 51% | 30% / 15% / 55% |
| After balancing + bulk chunk I/O: write s, max/median | 12.5 s, 2.27× (one 37.8M-site cell) | 14.1 s, 1.08× |

- **Hash skew was real.** With hundreds of cells, `md5 % N` left the busiest shard at 2–2.6× the median, and the slowest shard sets the write time. Shards are now balanced: barcodes hash into 4 buckets per shard during extraction, and buckets are packed into shards by site count (memory mode) or temp-file bytes (stream mode). A cell is never split, so the largest cell is a floor; at 32 shards here that floor (one 37.8M-site cell) is what remains.
- **On this 16-core box the write is CPU-bound in total.** With every writer busy, wall time did not move. On the 64-core server, where writers were idle most of the write, balancing should shorten the write; that needs a rerun there.
- **Chunk-site encoding** now reads and writes each cell's sites in one block instead of three 4-byte calls per site. Stream temp files are read whole and written with 1 MB buffers. Writer read share fell from 37% to 22–25%.
- **Contig name** is converted once per window instead of once per site.
- **Larger HDF5 chunks: rejected.** 262,144-record chunks (same slice memory) wrote in 14.7 s LZF / 17.3 s Blosc against 13.4 / 15.4 s at 65,536, with the same sizes. Amethyst's `getGeneM` reads gene-sized slices (`h5read(start=, count=)`), and every slice decompresses whole chunks, so bigger chunks would also slow Amethyst. `H5_CHUNK` stays 65,536.
- **`--min-mapq 0` skip: not done.** The check is one `u8` compare per read.
- **Integer `pct`: rejected.** The Amethyst dtype has `pct` as `<f8`.

**`--profile-json` overhead.** It read the clock around every record stage and every pileup column; the 64-core run took 12 s with it against 7.2 s without. Now 1 in 16 records and columns are timed and scaled up: 41.5 s → 42.2 s on the workstation BAM, and stage sums still add up to window totals (1,250 of 1,251 worker-seconds). Check `/sys/devices/system/clocksource/clocksource0/current_clocksource` on servers; anything other than `tsc` makes each clock read a system call.

**Inputs (both engines).** A missing `.bai` (`.csi` if a contig is ≥ 2^29 bp) or `.crai` is built next to the file, or in a temp directory when that directory is read-only. CRAM decodes against `-f`. Plain SAM is copied to a temporary BAM and must be coordinate-sorted. Rust: `src/input.rs`, through which every reader opens. Python: `taps_sc_extract/inputs.py`. On a chr19 slice, unindexed BAM, read-only BAM, SAM and CRAM gave identical calls to the indexed BAM in both engines; an unsorted SAM fails with a `samtools sort` hint.

**`--test` additions.**
- BAQ on vs off at the recommended worker count (workstation: 6.3× more reads/s with BAQ off at W=32).
- Stream-mode write through the temp filesystem vs memory mode.
- A few very large cells at the best writer count (per-writer sites/s on one cell).
- Output bytes per site, and with `-b`: the largest cell's share of sites and its write-time floor, the shard advice that follows from the cell count, and the predicted HDF5 MB/s as a share of the output disk's sequential write speed. On `bigcell_450.bam` it predicted the largest cell at ~40M sites (actual 37.8M).

**BAQ.** Python never runs BAQ; Rust runs htslib's by default, like asTair. On chr19 htslib BAQ lowers methylated CpG calls by 3.0% and unmethylated by 0.4% (CpG 67.07% → 66.48%; full genome 71.54% → 70.86%). In TAPS a methylated C reads as T, a mismatch BAQ treats as possible misalignment. Two TAPS-aware variants were tried on chr19 at the same cost as htslib BAQ: masking the convertible base to N (C on OT reads, G on OB reads) gave 67.24%; a 3-letter version (ref and read C→T on OT, G→A on OB) gave 67.24%. Both remove the methylation-specific loss; t calls still drop ~0.8%. Not shipped: choosing the default is a user decision. The experiment is kept as `scratch/baq-taps-aware-experiment.diff` (not in git). BAQ output is the same with and without FMA contraction (`-ffp-contract=off`, chr19 identical, 16 vs 0 FMA instructions in `probaln.o`).
- Faster bit-exact BAQ: no faster drop-in exists (newest htslib still uses scalar double `probaln_glocal`). The realistic options are reusing the HMM's per-read buffers (it zeroes ~120 KB per read; a few percent to ~15% of BAQ) or running several reads in SIMD lanes with the operation order unchanged (larger gain, much more work). Neither is started.
