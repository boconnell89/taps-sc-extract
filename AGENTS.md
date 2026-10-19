# AGENTS.md

> **Guide for AI Coding Assistants (Claude, Gemini, Grok, ChatGPT, etc.) working on `taps-sc-extract`.**

---

## 1. Project Overview

`taps-sc-extract` is a high-performance Python package and CLI tool for extracting single-cell DNA methylation from TAPS (mC-to-T chemistry) coordinate-sorted BAM files. It writes base-resolution methylation calls into **Amethyst-compatible HDF5 files** (`/CG/<barcode>/1` and `/CH/<barcode>/1`).

### Target Scale
- **Workstations**: 16 cores / 32 threads, ~60 GB RAM on this host. Latest snapshot: [`code/nvme-profiling-2026-10-18.md`](code/nvme-profiling-2026-10-18.md) (README [Performance](README.md#performance) has the summary and history); earlier: [`2026-10-17`](code/nvme-profiling-2026-10-17.md), [`2026-10-16`](code/nvme-profiling-2026-10-16.md), [`2026-10-15`](code/nvme-profiling-2026-10-15.md), [`2026-10-14`](code/nvme-profiling-2026-10-14.md). README performance figures are rates (mapped reads/s or reads/min) and percent changes, never clock times; raw times live in the snapshots. On a new machine run `taps-sc-extract --test` before quoting or tuning anything.
- **Production Servers**: 96 cores, 720 GB RAM (scales to >200k cells and >1B reads).

---

## 2. Codebase Architecture & File Map

```
taps-sc-extract/
├── pyproject.toml              # Build config, dependencies, CLI entrypoint (taps-sc-extract)
├── README.md                   # User documentation and performance guide
├── AGENTS.md                   # This AI developer reference
├── install.sh                  # Host-optimised build (native CPU flags), writes .cargo/config.toml
├── code/nvme-profiling-*.md    # Dated workstation timing snapshots
├── code/rust-implementation-log.md  # Dated record of Rust engine changes and measurements
├── scripts/
│   ├── benchmark_ab.py         # Automated A/B performance benchmark (Python vs Rust)
│   └── plot_performance_history.py  # README history plots (docs/img/perf_history_*.svg); update its table with the README's
├── taps_sc_extract/            # Reference Python extraction engine
│   ├── __init__.py             # Version definition
│   ├── __main__.py             # Module execution entrypoint (python -m taps_sc_extract)
│   ├── cli.py                  # CLI argument parser & engine dispatcher (--engine auto|rust|python)
│   ├── parallel_extractor.py   # High-throughput multiprocessing extraction engine
│   ├── extractor.py            # Single-process reference extractor
│   ├── h5_writer.py            # Amethyst HDF5 writer & sharded directory manager
│   ├── fasta.py                # FastFaiReader: process-safe .fai indexed FASTA reader
│   ├── calling.py              # TAPS mC-to-T lookup tables & SAM flag classification
│   ├── context.py              # Sequence context classifier (CpG, CHG, CHH, CG, CH)
│   ├── barcode.py              # Barcode extraction from QNAME/tags, whitelist & complexity-file parsing
│   └── inputs.py               # BAM/CRAM/SAM input: builds a missing index, converts SAM, CRAM reference
├── taps-sc-extract-rs/         # Rust acceleration core
│   ├── Cargo.toml              # Rust crate manifest (mimalloc, rust-htslib+libdeflate, hdf5-metno, rayon, libc)
│   ├── .cargo/config.toml.in   # Template install.sh fills (conda prefix, target CPU, forced CFLAGS incl. -ffp-contract=off)
│   ├── build.rs                # Compiles csrc/ against the vendored htslib headers (cc)
│   ├── csrc/baq_simd.c         # Batched BAQ: realn.c wrapper + grouping by HMM shape; AVX2/AVX-512 dispatch
│   ├── csrc/baq_kernel.h       # Lane-parallel probaln_glocal (included once per vector width)
│   ├── csrc/olap.c             # htslib 1.21 tweak_overlap_quality + cigar helpers, verbatim (build.rs checks against vendored sam.c)
│   ├── csrc/gpu/baq_hmm_body.h # GPU BAQ HMM, one read per work-item; shared by OpenCL, CUDA and HIP
│   ├── csrc/gpu/baq_opencl.c   # OpenCL backend (libOpenCL loaded at run time; Intel Arc/Data Center GPUs)
│   ├── csrc/gpu/baq_cuda.cu    # CUDA backend; the same file is the HIP backend with -DTAPS_HIP
│   ├── csrc/gpu/*.cl, taps_gpu.h # OpenCL prelude/entry point; host-side struct shared with the kernel
│   └── src/
│       ├── main.rs             # Rust CLI: extract, stats, tune, baq-selftest, write-shard (hidden, spawned by extract)
│       ├── baq_simd.rs         # FFI to csrc/baq_simd.c, --baq-simd widths, per-read self-tests
│       ├── baq_gpu.rs          # --baq-device gpu: backend selection, GPU service thread (pools reads across workers)
│       ├── extract.rs          # Batched read/BAQ/push loop, bam_plp column callback (--caller pileup), BAQ modes/masks/reference slices, stage timing
│       ├── readcall.rs         # Read-centric calling (--caller reads, default): replays bam_plp's depth cap and mate overlap, calls each read when htslib would free it
│       ├── extract_noodles.rs  # Slower pure-Rust reference pileup (noodles)
│       ├── parallel.rs         # Rayon window pool, thread-local readers, cancellation
│       ├── window.rs           # Genomic window planner
│       ├── accumulate.rs       # Barcode intern & per-cell sorted (pos, t, c) site vectors
│       ├── h5_out.rs           # Amethyst HDF5 writer: sliced cell datasets, process-per-shard writers, shard balancing, writer reports, master.h5
│       ├── input.rs            # BAM/CRAM/SAM input: builds a missing index, converts SAM, CRAM reference; all readers open through it
│       ├── shard_io.rs         # TAPSCK01 chunk format (temp files and writer stdin)
│       ├── shard_stream.rs     # Memory mode: writer children fed during extraction (overlapped transfer)
│       ├── autotune.rs         # Worker / window / writer / memory-mode defaults
│       ├── profile.rs          # --profile-json output
│       └── tune.rs             # --test: machine probes on a synthetic BAM, dataset sampling, saved profile
└── tests/
    └── test_taps_sc_extract.py # Pytest unit & integration test suite
```

---

## 3. Critical Invariants & Rules

When modifying or extending this codebase, you **MUST** adhere to the following rules:

### A. Dual Engine Architecture & Parity
- **Preserve Python Reference Engine**: Never remove the Python reference implementation in `taps_sc_extract/`. It serves as the baseline ground truth and fallback.
- **Unified Entrypoint**: `taps-sc-extract` CLI defaults to `--engine auto`, executing the compiled Rust binary when available and falling back to Python.

### B. Multiprocessing & Multithreading Safety
- **Python Workers**: Always use `mp.get_context("spawn")`. Never use `fork` with `pysam` or `h5py`.
- **Rust Rayon Workers**: Thread-local `IndexedReader` instances (never send BAM handles across threads); atomic `CancelFlag` for clean SIGINT cancellation.

### C. FASTA Reading Across Workers
- **Never use `pysam.FastaFile` in parallel Python workers**: `htslib`'s underlying `bgzf_read_block` is not safe when multiple processes read the same indexed FASTA.
- **Always use `FastFaiReader`**: Pure Python / Rust binary byte `seek()` and `read()` using `.fai` offsets is 100% process-safe.

### D. Memory Bounding & Shard Writer Concurrency
- **Process-tree accounting**: The summary reports the extraction process's peak RSS and the writer children's combined peak (each `write-shard` child reports its own peak). Count both when sizing memory.
- **Disk-Streaming vs. Memory Mode**: In `stream` mode, intermediate chunks are partitioned into `temp_dir/shard_XXX/chunk_YYYYYY.bin` and purged on write. In `memory` mode, per-window payloads stay in RAM and are piped to the writers.
- **HDF5 global lock**: hdf5-metno wraps every HDF5 call in one process-wide `ReentrantMutex` (`hdf5-metno-sys` `LOCK`), compression included. Threads in one process never write HDF5 in parallel. Parallel writing means one `write-shard` child process per shard (`--writer-mode process`, the default with >1 shard). In-process (`thread`) writers are auto-capped at 6.
- **Shard balancing**: Rust extraction hashes barcodes into `bucket_count(n_shards)` = 4 × shards buckets (`h5_out::BUCKETS_PER_SHARD`); after extraction `pack_buckets` deals whole buckets to shards by sites (memory) or temp bytes (stream). A cell stays in one bucket, so its records keep chunk_id order. Shard placement therefore differs from Python's `md5 % N`; compare outputs per cell, never per shard file.
- **Streamed writers (memory mode, process writers, at least one writer per shard; `TAPS_STREAM_WRITERS=0` disables)**: `shard_stream::ShardStreams` starts every `write-shard` child before extraction and pipes each window's payloads as the window finishes. Payload site vectors are moved, not cloned. Buckets are dealt to shards (`pack_buckets`) after 20% of windows, by sites so far; earlier payloads are held until then. A cell's dataset still needs every window, so HDF5 writing starts only after extraction. Senders and children run at nice 10 (`TAPS_STREAM_NICE`) so the extraction workers get the CPU first. Dropped without `finish` (error or cancel), the children are killed, never left to write partial shards. Shard placement (and so per-shard files) differs run to run; compare per cell.
- **Bounded writer memory**: `write_shard_from_chunks` writes one cell at a time in `WRITE_SLICE` (4 × 65,536-record) slices. Never build a whole shard's or a whole cell's `MethRec` vector.

### E. Per-Read State in the Pileup Callback
- **Key per-read caches on `bam1_t.id`**, the number htslib gives each read it stores in an iterator (0, 1, 2, … never reused within a window). **Never key on the record or QNAME address**: htslib recycles pileup nodes LIFO and `bam_copy1` reuses buffers, so a new read can sit at the previous read's address. A QNAME-pointer cache did this until 2026-10-15 and put ~0.01% of calls under the wrong barcode (see `code/rust-implementation-log.md`).
- **Barcode filters (`--whitelist`, `--complexity` + `--min/--max-complexity`) act at call time**, in both engines: a filtered cell's reads stay in the pileup, so depth caps and mate overlaps for kept cells are unchanged. Rust skips only their BAQ (`accept_record`), which touches no other read. Never drop filtered reads before `bam_plp_push`; that frees `--max-depth` slots and changes kept cells' calls. Complexity files: barcode in column 2, unique reads in column 4 (`barcode::parse_complexity`, `barcode.py::parse_complexity`).
- **`--context CG|CH` skips at the reference position** (`ContextSel::wants` in `readcall::call_read` for the default caller, `consume_column` for `--caller pileup`; Python `if ctx not in contexts`), after the context is classified and before any per-read work, in both engines. Never filter reads for it: the pileup, BAQ, depth caps and overlaps must stay identical, so the kept context's calls equal a `both` run's. Groups (`/CG`, `/CH`) are created only when a cell has data, in shards, single files and `master.h5`, so a skipped context has no group.
- **Rust default compression is `blosc`** (LZ4 level 9; `RUST_DEFAULT_COMPRESSION` in `cli.py`, `--compression` default in `main.rs`). Reading it from Python needs `import hdf5plugin` before `h5py` (tests do this in `_cells_of`); R needs `rhdf5filters`, as for `lzf`.
- **Read-centric calling (`--caller reads`, the default) must equal htslib's column pileup per cell** (`--caller pileup` keeps the old path for A/B). `readcall::ReadPileup` replays `bam_plp_push` + drain for each kept read in file order, with fresh state per window:
  - **Depth cap**: drop a read when it starts at the start of the last accepted read and live reads + 1 (htslib's sentinel node) > `--max-depth`. A dropped read removes its QNAME from the overlap table and does not move the position. The cap fires at the default 250 on the test BAM (805,351 reads dropped).
  - **Retiring**: a read is freed, and only then called, once a later accepted read starts past its end (`end < pos`), or at window end. Reads with no reference span are never linked.
  - **Mate overlap**: `overlap_push` rules verbatim (proper pair, mate mapped, same contig, the `isize`/`mpos >= end` exclusion, store only while the mate is ahead), then `taps_tweak_overlap_quality` (`csrc/olap.c`, htslib's code copied byte for byte; never rewrite it). Freeing any read removes its QNAME's entry. Supplementary reads are not filtered, so they take part.
  - **Batches**: reads are copied (`bam_copy1`) at push, so batch slots can be reused after their push as before.
  - **Accumulation**: per read, calls are collected in position order and merged into the cell with `CellMaps::add_run`; a mate that overlaps the previous read is merged into the tail.
  - **Checks after any change**: `test_read_caller_matches_pileup` (includes `--max-depth 1`/`3`), plus per-site dumps (`stats --dump-sites`) of both callers on chr19.
- **Per-cell sites arrive in position order** (columns ascend within a window), so `CellMaps::add` appends or adds to the last entry. Keep that invariant; `normalize()` repairs any out-of-order input.
- **Parity checks must be per cell**, not totals: `test_chr19_rust_matches_python_per_cell` requires Rust `--no-baq` to equal Python exactly on chr19.
- **BAQ is TAPS-aware by default** (`extract::BaqMode`): the HMM sees C as N for OT reads and G as N for OB reads (`baq_mask`, strand from `baq_is_ob`), so methylation conversions are not treated as mismatches. No masked reference is stored: each window fetches one unmasked slice padded by `BAQ_PAD` (16 kb), every `BaqJob` carries the slice bounds and its mask byte, and `csrc/baq_simd.c` applies the mask while translating the reference (the htslib fallback masks a per-read copy of `taps_baq_need`'s span). `baq_need` (Rust) and `taps_baq_need` (C) must agree and must cover every base `sam_prob_realn` can read; a read past the padded slice gets its own fetch. The default `--baq taps` skips the HMM for reads that are all-M and match the masked reference exactly (`BaqMode::wants`; ~72% of reads, 0.0007% of CG calls differ from `--baq taps-full`). `--baq taps-partial` runs it only for reads with an indel or clip. Any change to the skip rule needs a new frozen chr19 reference. `--baq htslib` is the plain reference (biased against methylated calls). Calling always uses the unmasked reference. The Python engine never applies BAQ, so parity with Python is only defined for `--no-baq`.
- **SIMD BAQ must equal htslib read by read.** `csrc/baq_kernel.h` is htslib 1.19.1 `probaln_glocal` with `double` replaced by a vector of 4 or 8 reads; every lane does the same IEEE operations in the same order. Keep `-ffp-contract=off` (no FMA) on both it and htslib (`.cargo/config.toml.in`, `build.rs`), never use `-ffast-math` or a vector `log`, and zero `f`/`b` per group (the forward pass reads zeroed cells outside the band). HMM parameters stay `float` (`probaln_par_t { float d, e; }`): htslib evaluates `(1 - c->d) / l_ref` etc. in float, and using doubles shifted ~1 read in a million by one quality point. Reads the kernel does not cover (BQ/ZQ tag, >1000 bp, empty span) go to `sam_prob_realn`, as does everything on non-x86 hosts (`TAPS_BAQ_NO_SIMD` builds that path on x86 for testing). After any change run `taps-sc-extract-rs baq-selftest -b BAM -f FA` (QUAL memcmp per read; `install.sh` runs it and fails the install on a mismatch). Reads are BAQ'd per 8,192-record batch, then pushed and drained one at a time in file order; never push a batch before draining, or `--max-depth` drops different reads.
- **GPU BAQ must equal the CPU read by read** (`--baq-device gpu`, default `cpu`). One kernel body (`csrc/gpu/baq_hmm_body.h`) serves every backend; keep the arithmetic identical to `csrc/baq_kernel.h` (no contraction: CUDA `-fmad=false --prec-div=true`, HIP `-ffp-contract=off`, OpenCL `FP_CONTRACT OFF` and no fast-math build options). HMM parameters (`m[]`, `sM`, `bM`, `bI`) and the quality table are computed on the host (`taps_baq_gpu_export`, `taps_baq_gpu_q2p`) with probaln.c's float expressions; the device does only double + - * / and log. Bases whose `-4.343*log(1-max)+.499` is within 1e-9 of an integer (or not finite) come back as 0xFF and are recomputed on the host with glibc. Rows of f and b are zeroed to `i_dim + 3` cells because htslib reads up to 3 cells past a row (zeros in the next row). Reads the device does not take (htslib-only, l_query > 192, i_dim > 72) and every read when the device fails go to `taps_baq_batch` after `taps_baq_gpu_finish`. Check with `baq-selftest --gpu -b BAM -f FA`; `install.sh` runs it when a GPU is usable. Batches are shared with the CPU: `baq_gpu::take_batch(slack)` sends a batch to the device only if its queue (backlog / measured rate) clears within the worker's slack (its push-and-call time per batch, `CALL_BATCH_S`, × batches already in flight); otherwise `ws.run` does it on the worker. Never route by a fixed queue cap: on a CPU-bound host (B70) a 25 ms cap sent 27% of batches to the CPU and lost 4%. With the GPU, `process_window` keeps `baq_gpu::depth()` batch slots (each its own record pool, jobs, reference extras and C workspace): `submit` starts BAQ on a new batch, `complete` finishes the oldest, which is then pushed and drained; batches are pushed strictly in file order, and a slot's records are never reused before its push.
- **htslib version matters for calls.** Rust builds hts-sys 2.2.1 (htslib 1.19.1) patched with htslib 1.21's overlap clipping (`taps-sc-extract-rs/patches/`, applied by `install.sh` into `vendor/`, re-applied whenever the patch file's checksum changes) to match pysam 0.23 (htslib 1.21). Re-check parity whenever pysam or rust-htslib changes.

### F. TAPS Chemistry & Amethyst Schema
- **TAPS Calling Logic**:
  - Original Top strand (`OT`, `+`): Reference `C` $\to$ Read `T` is **Methylated** (`c += 1`), Read `C` is **Unmethylated** (`t += 1`).
  - Original Bottom strand (`OB`, `-`): Reference `G` $\to$ Read `A` is **Methylated** (`c += 1`), Read `G` is **Unmethylated** (`t += 1`).
  - **5′ sequenced-base trim**: drop `--trim-r2` (default 10) bases of read 2 and `--trim-r1` (default 0) of read 1 before calling. Reverse-strand BAM SEQ is reverse-complemented, so the original 5′ end is at the high end of SEQ. Implemented by `keep_after_5prime_trim()` in Python and Rust.
- **Amethyst Structured Array Dtype**:
  `[('chr', 'S10'), ('pos', '<i8'), ('pct', '<f8'), ('t', '<i8'), ('c', '<i8')]`
  - `t` (unmethylated) comes **before** `c` (methylated).
  - Coordinates are **1-based** integers.
  - Chromosome names are **ASCII bytes** (`b'chr1'`).

---

## 4. How to Run the Code

### CLI Execution
```bash
# 1. Standard run (auto-selects Rust core if compiled, else Python)
taps-sc-extract -b sample.srt.bam -f ref.fa -o output.h5 -t 24

# 2. High-throughput sharded directory output with Rust backend
taps-sc-extract -b sample.srt.bam -f ref.fa -o sharded_dir/ -t 24 --shards 16 --engine rust

# 3. Explicit Python reference engine run
taps-sc-extract -b sample.srt.bam -f ref.fa -o output.h5 -t 24 --engine python

# 4. Automated A/B performance comparison
python3 scripts/benchmark_ab.py -b sample.srt.bam -f ref.fa -c chr1,chr2,chr3 -t 24 --shards 16

# 5. Tune a new machine (synthetic workload; add -b/-f/-o to sample a dataset)
taps-sc-extract --test --save-profile

# 6. Where does the time go? (per-stage and per-window timings)
taps-sc-extract -b sample.srt.bam -f ref.fa -o out/ --shards 16 --profile-json run_profile.json

# 7. Build the Rust core for this host
./install.sh --no-pip
```

### Python Programmatic API
```python
from taps_sc_extract.parallel_extractor import extract_methylation_parallel

summary = extract_methylation_parallel(
    bam_path="/path/to/sample.srt.bam",
    fasta_path="/path/to/ref.fa",
    out_h5_path="/path/to/output_dir",
    n_workers=24,
    n_shards=16,
    chunk_size_mb=10,
    use_temp_files=True,
    trim_r1=0,
    trim_r2=10,
)
print("Extraction Summary:", summary)
```

---

## 5. R & Amethyst Integration Verification

To verify that generated HDF5 files or `master.h5` files are compatible with Amethyst:

```R
library(amethyst)
library(rhdf5)

# For single file: "output.h5" | For sharded directory: "sharded_dir/master.h5"
h5_path <- "sharded_dir/master.h5"

ls_df <- h5ls(h5_path)
barcodes <- ls_df$name[ls_df$group == "/CG"]

h5_paths <- data.frame(path = h5_path, barcode = barcodes)
obj <- createObject(h5paths = h5_paths)

# Index CG and CH contexts
obj_cg <- indexChr(obj, "CG")
obj_ch <- indexChr(obj, "CH")
```

---

## 6. Running Tests

```bash
# Run all unit tests
pytest tests/ -k "not test_chr19_real_bam_smoke" -v

# Run full test suite including chr19 BAM tests (Python smoke, Rust frozen reference, Rust-vs-Python per cell)
pytest tests/ -v

# Rust unit tests
(cd taps-sc-extract-rs && cargo test --release)
```

---

## 7. Troubleshooting & Common Pitfalls

| Issue | Root Cause | Solution |
| :--- | :--- | :--- |
| `HDF5 file locking error` | Network filesystem or multi-process file open | Set `os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"` (already set in `h5_writer.py`). |
| Worker freeze / hang | `pysam` used with `fork` | Always use `mp.get_context("spawn")`. |
| Memory spike / OOM | Loading all chunks into RAM | Ensure `use_temp_files=True` (default) and use `--shards 16`. |
| `Chunk size must be <= data shape` | `chunks=(65536,)` on a cell with 10 calls | Use `chunk_size = min(len(records), 65536)`. |
| Missing `master.h5` link target | Relative paths misaligned | Use relative file names `shard_XXX.h5` inside `master.h5` so the directory is portable. |
| More writer threads, slower HDF5 write and more RAM | hdf5-metno's global lock serializes in-process writers | Use `--writer-mode process` (default with >1 shard). |
| `failed to load source for dependency hts-sys` | `vendor/hts-sys-2.2.1` missing | Run `./install.sh` once (downloads and patches hts-sys). |
| Rust calls drift from Python after a pull | `patches/` changed but the build used the old vendored copy | Rerun `./install.sh`; it re-patches `vendor/` when the patch checksum differs from `vendor/hts-sys-*/.taps_patched`. |
| SIGILL from `taps-sc-extract-rs` | Native build copied to a different CPU | Rerun `./install.sh` on that host, or build with `--portable`. |
| Timings vary 2× between identical runs | Concurrent builds or other jobs on the box | Never time while compiling; repeat runs and report the spread. |
