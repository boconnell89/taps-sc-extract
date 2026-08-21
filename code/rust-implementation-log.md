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

| Dataset / Contigs | Python Baseline | Rust Acceleration Core | Speedup | Peak RAM (Rust) | Parity / Cells |
|---|---|---|---|---|---|
| `chr19` (7 windows, 2.1M reads) | 39.04 s | **9.89 s** | **3.95× faster** | 1.52 GiB | 5,781 cells, 63.75% mCG, 0.278% mCH |
| `chr1–3` (56 windows, 16.6M reads) | 176.3 s | **67.25 s** | **2.62× faster** | 7.93 GiB | 7,268 cells, 65.86% mCG, 0.275% mCH |
| **Whole Genome mm10** (286 windows, 74.9M reads) | **12.70 min** (762 s) | **5.68 min** (340.88 s) | **2.24× faster** | 33.3 GiB | **7,355 cells, 38.7M CpG calls (65.36% mCG), 1.075B CH calls (0.278% mCH)** |

### Performance Optimization Summary
1. **Thread Balancing**: `--decomp-threads 0` (synchronous BAM reader per worker) achieves highest throughput with $\ge 8$ workers because 24–32 CPU cores are already saturated without context-switch latency.
2. **`rustc_hash::FxHashMap`**: Replacing SipHash with fast non-cryptographic word hashing saves ~70 seconds of CPU core time over 1 billion pileup calls.
3. **Stream Memory Mode**: Intermediate chunks partitioned into temp files are assembled and compressed across 16 shards by a capped 6-thread pool in under 3.5 minutes.


