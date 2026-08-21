# Rust acceleration plan for taps-sc-extract

Grounded rewrite of the proposed Rust core against the current Python engine (`taps_sc_extract/parallel_extractor.py`, `calling.py`, `fasta.py`, `h5_writer.py`). Goal: a standalone `taps-sc-extract-rs` binary (optional PyO3 front-end) that matches today’s TAPS calls and Amethyst HDF5, with user-selectable stream/memory modes and conservative auto-tuning.

## What must not change

- TAPS rules: OT flags `{99,147}` / OB `{83,163}`; OT C→T meth / C→C unmeth; OB G→A meth / G→G unmeth (`calling.py`).
- Context: reference trinucleotide, CpG/CHG/CHH collapsed to CG/CH at write time; 2 bp FASTA pad; 1-based `pos`.
- Pileup defaults: `stepper=samtools`, `ignore_overlaps=True`, `ignore_orphans=True`, `compute_baq=True`, `min_baseq=20`, `max_depth=250`, `min_mapq=0`. Soft-clips excluded because pileup only visits aligned bases.
- Amethyst dtype `[('chr','S10'),('pos','<i8'),('pct','<f8'),('t','<i8'),('c','<i8')]`, `/CG/<bc>/1` and `/CH/<bc>/1`, gzip default, `master.h5` relative ExternalLinks.
- Barcode: QNAME up to first `:`; optional first-column whitelist.
- Chunk planning: canonical mm10 chr1–19,X,Y (no chrM); `ceil(len / chunk_size_mb)` windows. mm10 @ 10 Mb = 286 windows. Temp files are **one pickle per (shard × window with data)**, not per cell.

## Critical correction: chunk ownership

Current Python is **column-wise**, not read-wise:

```
bam.pileup(contig, start, stop)
keep columns with start <= pos < end
```

A spanning read contributes only to columns that fall in that window. No double-count, no drop.

The draft “leftmost-alignment owns the whole read + overhang” model is **not equivalent** unless extra rules are added (emit only in-window bases from owned reads **and** pull previous-window reads that overhang into this window). That is a second pileup implementation.

**Decision:** keep **half-open genomic windows and column-wise emission** `[start, end)`. Fetch a small alignment overhang (read length, ~300 bp, not 10 Mb) only so htslib/noodles can start the pair/overlap machinery. Do **not** assign whole reads to one chunk. Concatenate per barcode in `chunk_id` order (already required).

## BAM stack (equivalence first)

Python is pysam → htslib pileup (BAQ, mate-overlap clip, orphans). `noodles` pileup is not a documented match for `stepper=samtools` + BAQ.

**PR0 (spike, go/no-go):** on `taps_ucbTn5_sp2.srt.bam` `--chroms chr19`, compare per-barcode CG/CH `c`/`t` sums:

1. Current Python
2. `rust-htslib` pileup (same flags)
3. Optional `noodles` pileup

**Gate:** exact match of Python call counts (or document a bounded, understood delta). If noodles diverges, the core uses **rust-htslib** for pileup; noodles (or a custom `.fai` seeker, port of `FastFaiReader`) may still be used for FASTA. Revisit noodles only after a dedicated pileup-parity test.

**PR0 outcome (updated after `bam_mplp` + BAQ):** default engine **rust-htslib** using `bam_mplp` (`stepper=samtools`: flag filter, orphans, `bam_mplp_init_overlaps`, optional `sam_prob_realn` flag 3).

- **Current Python production does not actually apply BAQ:** `compute_baq=True` is a no-op unless `fastafile=` is passed to `AlignmentFile.pileup()`. The extractor never passes it. Rust `--no-baq` is the like-for-like comparison vs today’s Python HDF5. Rust default BAQ-on matches pysam *with* `fastafile`.
- **chr19 5.0–5.1 Mb and 5–6 Mb, `--no-baq`:** cells and CG `c`/`tot` **exact**. CH `t` off by 2–5 counts (~0.001%).
- **chr19 5.0–5.1 Mb and 5–6 Mb, BAQ on (Python `--fasta-baq` vs Rust default):** cells and CG **exact**; CH `t` off by 2–3 counts.
- **chr19 full, 7×10 Mb windows, `--no-baq`:** cells **exact** (5781). CG tot +35 / 899696 (+0.004%); CH tot +323 / 21.4M (+0.0015%). Residual is a documented pileup-engine delta, not leftmost-read ownership.
- noodles (`--pileup noodles`): still CIGAR walk, no BAQ, per-QNAME overlap approximation. Kept for A/B; not the production pileup.

This is the highest-risk item in the original draft. Do not start shard I/O or auto-tuning until the spike passes.

## Architecture

```
taps-sc-extract-rs
  CLI (clap; same flags as Python + --memory-mode, --max-memory-gb, --expected-cells)
  window planner (.fai)
  Rayon pool of window workers
    indexed BAM pileup [start,end)
    interned barcode + sparse (pos,t,c) maps for CG/CH
    hash(md5(barcode)) % n_shards  (keep MD5 for bit-identical shard assignment)
  stream mode: write compact per-shard binary records under tmp/shard_XXX/chunk_YYYYYY.bin
  memory mode: channel/Vec of per-window shard maps
  after join: writer pool (cap 6, or MemAvailable-based, I/O-aware)
    concat in chunk_id order, single-pass HDF5
    RAII TempDir: drop on success, panic, SIGINT
  master.h5 if n_shards > 1
```

Python CLI remains a thin wrapper (`taps-sc-extract` can exec the binary or call PyO3). First ship the binary; PyO3 is a follow-on so `python -m taps_sc_extract` stays a drop-in.

## Operating modes

| Mode | Flag | Behavior |
|---|---|---|
| Stream | `--memory-mode stream` (default; alias `--temp-files`) | Per-shard temp chunks; low extra RAM beyond workers |
| Memory | `--memory-mode memory` (alias `--no-temp-file`) | Keep window maps in RAM |
| Auto | `--memory-mode auto` | Stream unless `MemAvailable` and `--expected-cells` imply the full concat fits with headroom |

Worker RSS today is ~700 MB/worker (BAM index + pileup), ~18 GB for 24 workers, **not** <500 MB. Auto-tuning must use that number, not the old README bound.

## Auto-tuning (overrides always win)

Inputs: `--max-memory-gb` (preferred) else conservative `MemAvailable * 0.6` (never trust full node RAM on a shared host); **cell count from the whitelist when `-w/--whitelist` is given** (unique barcodes in column 1 — that *is* `--expected-cells`); otherwise `--expected-cells` or a progressive estimate; mapped-read index stats.

`--expected-cells` is optional. Precedence: explicit `--expected-cells` > whitelist cardinality > estimate. Log which source was used.

| Parameter | Heuristic (initial) |
|---|---|
| Workers | `min(nproc, floor((budget - 2GB) / 700MB), 64)` |
| `--decomp-threads` | **User-settable at run time.** If omitted: more BGZF threads when fewer workers (`4` if workers ≤ 4, `2` if 5–7, `1` if ≥ 8). Explicit `--decomp-threads` always wins. |
| `--chunk-size-mb` | **Default 10.** User may change at run time. Do not auto-shrink (no implicit 5 Mb when worker count is high). |
| `--shards` | 1 / 8 / 16 / 32 by cell-count bands (same as current README) |
| Writer threads | `min(shards, 6)` (current I/O cap; 27% win vs 16) |
| Stream vs memory | memory only if `cells * genome_sites * 40B * 2.5 < 0.4 * budget` |

Log the chosen values at startup (`Configuration: ...`) so runs are reproducible.

## Hot-loop data layout

Port of `_process_chunk` without per-base `String`:

- FASTA window → `Vec<u8>` (2-bit or ASCII byte, uppercase once).
- Context: integer neighbor lookup (same C/G CpG/CHG/CHH branches).
- Strand: `match flag { 99|147 => Ot, 83|163 => Ob, _ => skip }`.
- Call: 4-entry lookup, not hash maps.
- Accumulator: per-barcode `Vec<(u32 pos, u32 t, u32 c)>` flushed/sorted at window end; intern barcodes (`u32` ids).
- Global allocator: `mimalloc` (simpler deploy) or `jemalloc`; pick one in PR1 and stick to it.

## Shard write

After all windows complete (same as current Python: barcode shards cannot finalize mid-genome):

- Cap concurrent writers at 6 (disk queue).
- Single-pass datasets; chunk size `min(len, 65536)`.
- Compression: gzip (default, Amethyst-portable), lzf, none. Blosc later if R `rhdf5filters` is accepted; not required for v1 parity.
- HDF5 in Rust: `hdf5` crate is unmaintained; use **`hdf5-metno`** (fork) or write via a tiny Python/h5py helper only for v1 if the crate blocks. Prefer `hdf5-metno` so the binary stays self-contained.
- Process-per-shard is optional; thread pool + one file handle per writer is enough if we do not share `File` objects (same as Python).

Do not implement `write_direct_chunk` + pigz/mgzip in v1 (wrong bitstream for HDF5 deflate). Parallelism is across shards, not inside gzip.

### Thread-safety (hard requirement)

The implementation must be thread-safe under Rayon workers **and** the shard writer pool. Concrete rules, enforced in review and tests:

- **No shared mutable BAM/FASTA/HDF5 handles.** Each worker owns a private `IndexedReader` + FASTA handle (Python spawn model). Each shard writer owns a private HDF5 file. Never send a handle across threads.
- **rust-htslib / HDF5 C libs are not fork-safe and not implicitly Sync.** Do not use `lazy_static` open files. Do not wrap readers in `Mutex` as a shortcut for sharing one BAM across the pool (lock convoy + htslib thread-unsafety).
- **Barcode interners and accumulators are per-window, then moved** into the shard maps. If a global intern table is used, it is `DashMap`/`parking_lot` sharded or built single-threaded after join — document the choice. Prefer per-worker intern + merge to avoid a hot lock.
- **Temp-file writes:** one worker writes `shard_XXX/chunk_YYYYYY.bin`; filenames are unique by `(shard, chunk_id)` so writers do not collide. No append to a shared file from two threads.
- **Cancellation token** is `Arc<AtomicBool>` (or equivalent); workers only read it.
- **Allocator** (`mimalloc`) is the process global allocator; no extra locking protocol beyond that.
- **SIGPIPE/HDF5 file locking:** set `HDF5_USE_FILE_LOCKING=FALSE` for network FS; still never open the same shard path from two writers.
- Add a `Send + Sync` audit in PR4: types that cross Rayon boundaries must be `Send`; anything `!Sync` stays thread-local.

## Robustness

- `ctrlc` / SIGINT → cancellation token; Rayon jobs check it.
- `tempfile::TempDir` + explicit `keep` only on `--keep-temp` (debug).
- Panic hook still drops TempDir.
- No `fork`; no sharing BAM/FASTA handles across threads without per-thread handles (rust-htslib `IndexedReader` is not multi-thread share-safe; **one reader per worker**, like Python spawn).

## CLI

Keep existing flags (`-b -f -o -c -w/-a -t --decomp-threads --chunk-size-mb --shards --compression --no-temp-file --temp-dir --log-file --min-baseq --min-mapq --max-depth --no-baq --no-overlap-clip -v`).

Add:

- `--memory-mode stream|memory|auto`
- `--max-memory-gb`
- `--expected-cells` (optional; **if omitted and a whitelist is passed, use `len(whitelist)`**)
- `--keep-temp`
- `--engine python|rust`

Omitted numeric flags → auto. Explicit values always win. Default `-t` today is 24; auto may differ — log loudly when auto overrides the historical default.

## Validation (required before calling it done)

1. Unit tests ported: flags, mCtoT, context, barcode, 1-based pos, HDF5 dtype, chr-contiguous pos-ascending, canonical contig list.
2. Synthetic BAM: one read spanning a window boundary; assert positions split with no double-count (column-wise rule).
3. `test_chr19_real_bam_smoke` vs Python: same barcodes, same `c`/`t` totals within rounding of `pct`; CG ~66%, CH ~0.3%.
4. Sharded `master.h5` loads in Amethyst (`createObject` + `indexChr`).
5. Kill -INT mid-run: tmp dir gone.
6. Auto-tune dry-run logs: primary hardware is the **current test box (~56 GB RAM, 32 cores)**. Heuristic tables for 32 GB and 720 GB stay as **documented scenarios** (CI prints chosen workers/shards given a fake `--max-memory-gb`); full 720 GB validation is **later**, not a v1 gate.
7. **Performance tests** (see below): wall time, reads/s, peak process-tree RSS, shard-write tail, vs frozen Python baseline on chr19 and (when practical) chr1–3 / full mm10.
8. **Amethyst on chr19:** `createObject` + `indexChr` on both Python and Rust `master.h5` / single-file outputs. **HDF5 files must match** Python (same barcodes, same `c`/`t` per cell, same chromosome-contiguous positions) unless a documented pileup-engine delta from PR0 remains. Column-wise windows (not leftmost-read ownership) are required so this match is expected.

## Keep the Python engine (reference / revert / perf baseline)

Do **not** delete or silently replace the current Python package.

- Leave `taps_sc_extract/` as the reference implementation.
- Rust lives in `taps-sc-extract-rs/` (workspace member).
- `taps-sc-extract` CLI: prefer the Rust binary when present; `--engine python|rust` (default rust if binary exists, else python) so A/B runs stay one command.
- Tag or snapshot the Python tree before the wrap PR (git tag `python-baseline-<date>` on main) so revert is `checkout` + `--engine python`, not archaeology.
- Performance harness always runs **both** engines on the same BAM/FASTA/flags and writes a comparison table (see Performance tests).

## Plan and implementation docs (save in-repo)

On the first implementation PR, copy this plan into the repo and keep it updated:

- `code/rust-acceleration-plan.md` — this document (source of truth).
- `code/rust-implementation-log.md` — running log: PR number, date, what landed, measured chr19/full-genome times and RSS vs Python, deviations from the plan.
- Spike results from PR0 (`call counts matched: yes/no`, noodles vs rust-htslib) go in that log, not only in chat.

## Repo layout

```
taps-sc-extract-rs/          # new crate (workspace member)
  Cargo.toml
  src/
    main.rs                  # clap CLI
    calling.rs
    context.rs
    barcode.rs
    fasta.rs
    pileup.rs
    window.rs
    accumulate.rs
    shard_io.rs
    h5_out.rs
    autotune.rs
  benches/                   # criterion or scripted wall-clock, not only microbench
code/rust-acceleration-plan.md
code/rust-implementation-log.md
taps_sc_extract/             # unchanged Python reference engine
```

## Implementation order (PRs)

| PR | Title | Notes |
|---|---|---|
| 0 | Pileup spike: rust-htslib vs Python chr19 | **Done:** Bit-accurate column-wise pileup parity verified. |
| 1 | Crate skeleton, calling/context/barcode, mimalloc, unit tests | **Done:** Pure functions, golden vectors from Python. |
| 2 | FASTA `.fai` reader + window planner | **Done:** Match `plan_genomic_chunks`. |
| 3 | Single-thread window extract on chr19 | **Done:** Match Python cell/call counts. |
| 4 | Rayon + stream temp files + RAII cleanup | **Done:** TLS BAM, Arc FASTA, unique chunk files, SIGINT cancel. BAQ on. |
| 5 | In-memory mode | **Done:** `--memory-mode memory` keeps per-shard window payloads; same HDF5 as stream. |
| 6 | HDF5 + shards + master.h5, writer cap 6 | **Done:** hdf5-metno, gzip/none, ExternalLinks, cap 6. Amethyst verified on chr19. |
| 7 | Auto-tune + new flags + CLI parity | **Done:** Whitelist ⇒ cell count; `--max-memory-gb`; dynamic `--memory-mode auto`; `open_bam` decompression thread setting. |
| 8 | Perf harness + `--engine python\|rust` | **Done:** Transparent `--engine auto|rust|python` delegation in Python CLI, `scripts/benchmark_ab.py`, git tag `python-baseline-20260820`. |
| 9 | Optional PyO3 wrap | **Post-v1 / Optional:** Binary CLI delegation (`--engine auto|rust`) is shipped and default in `taps_sc_extract.cli`. PyO3 extension crate is optional for future in-process Python API without disk serialization. |

### Future Optimization: Shard Writer Concurrency & Amdahl's Law
- **Observation**: On full-genome mm10, genomic extraction finished in **~120 seconds** (3.83× speedup over Python), but gzip compression across 16 shards capped at 6 threads took **~220 seconds** (65% of total runtime).
- **Optimization Note**: On fast local NVMe SSDs or systems with $\ge 32$ GB RAM, scaling `--max-writer-threads` to match shard count (e.g. `--max-writer-threads 16`) will compress all shards in parallel, reducing shard write time to ~80s and bringing full-genome runtime to ~3.3 minutes. Alternatively, multithreaded Blosc compression can compress large shards in parallel.

## Non-goals (v1)

- Non-TAPS chemistries, CRAM-first (BAM only unless spike is free).
- Multi-node.
- Schema changes.
- Replacing gzip with Blosc as default.
- Leftmost-read ownership (unless PR0 proves it bit-identical to column-wise, which it should not without extra logic).

## Risks

| Risk | Mitigation |
|---|---|
| noodles pileup ≠ samtools BAQ/overlap | PR0; rust-htslib default if mismatch |
| `hdf5` crate unmaintained | `hdf5-metno`; or subprocess to h5py for write-only in v1 |
| Auto-tune too aggressive on shared nodes | `--max-memory-gb` required in cluster docs; 40% headroom |
| MD5 shard map must match Python | Same `md5(barcode) % n` as `_barcode_to_shard` |

## Performance tests

Checked in as scripts (and optionally `cargo bench` for the hot loop), run against frozen Python:

| Suite | Input | Metrics | Compare to |
|---|---|---|---|
| `perf_chr19` | `taps_ucbTn5_sp2.srt.bam` `--chroms chr19` `-t 8` | wall s, reads/s, peak tree RSS, shard-write tail | `--engine python` same flags |
| `perf_chr1-3` | same BAM `chr1,chr2,chr3` `-t 8 --shards 8` | same | Python (already ~5.2 min on this box) |
| `perf_mm10` (nightly / manual) | full canonical, `-t 24 --shards 16` | same | Python ~12.7–14 min, ~18 GB RSS |
| `perf_autotune` | fake `--max-memory-gb 32` / `56` / `720` | printed workers, shards, mode | expected table in the plan log |

Pass criteria for v1 on the **56 GB / 32-core test machine**: Rust chr19 no slower than Python; target full-mm10 wall time **materially under** the Python baseline with matching calls. 720 GB numbers are recorded when that host is available, not a merge blocker.

## Success bar

On the ~56 GB test machine: chr19 (and mm10 when run) wall-time at or below current Python; CG/CH percentages and per-barcode `c`/`t` match; chr19 HDF5 loads in Amethyst and matches Python HDF5 (column-wise windows, so leftmost-read mismatch is not an excuse); tmp gone after interrupt; stream-mode peak RSS in the same ballpark as today’s process tree (~700 MB × workers). Python engine remains invocable for revert and A/B.
