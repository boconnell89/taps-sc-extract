# AGENTS.md

> **Guide for AI Coding Assistants (Claude, Gemini, Grok, ChatGPT, etc.) working on `taps-sc-extract`.**

---

## 1. Project Overview

`taps-sc-extract` is a high-performance Python package and CLI tool for extracting single-cell DNA methylation from TAPS (mC-to-T chemistry) coordinate-sorted BAM files. It writes base-resolution methylation calls into **Amethyst-compatible HDF5 files** (`/CG/<barcode>/1` and `/CH/<barcode>/1`).

### Target Scale
- **Workstations**: 32 cores, 32 GB RAM (processes ~75M mapped reads across mm10 in **~12.7 minutes** with **~18 GB peak total process-tree RAM** across 24 workers).
- **Production Servers**: 96 cores, 720 GB RAM (scales to >200k cells and >1B reads).

---

## 2. Codebase Architecture & File Map

```
taps-sc-extract/
├── pyproject.toml              # Build config, dependencies, CLI entrypoint (taps-sc-extract)
├── README.md                   # User documentation and performance guide
├── AGENTS.md                   # This AI developer reference
├── scripts/
│   └── benchmark_ab.py         # Automated A/B performance benchmark (Python vs Rust)
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
│   └── barcode.py              # Barcode extraction from QNAME/tags & whitelist parsing
├── taps-sc-extract-rs/         # High-Performance Rust Acceleration Core (2.6x faster, 56% less RAM)
│   ├── Cargo.toml              # Rust crate manifest (mimalloc, rust-htslib, hdf5-metno, rayon, rustc-hash)
│   └── src/
│       ├── main.rs             # Rust CLI entrypoint (extract, stats)
│       ├── extract.rs          # Column-wise bam_mplp samtools-style pileup & BAQ
│       ├── parallel.rs         # Rayon thread pool & cancellation tokens
│       ├── accumulate.rs       # Barcode intern & FxHashMap cell position maps
│       ├── h5_out.rs           # Amethyst structured array HDF5 generator & master.h5
│       ├── shard_io.rs         # Compact binary intermediate chunk partition files
│       └── autotune.rs         # Hardware & memory budget auto-tuning heuristics
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
- **Accurate Process-Tree Accounting**: Memory is monitored across the entire process tree (`/proc/<pid>/statm` or `sys_info`).
- **Disk-Streaming vs. Memory Mode**: In `stream` mode (default), intermediate chunks are partitioned into `temp_dir/shard_XXX/chunk_YYYYYY.bin` and purged on write. In `memory` mode, maps are kept in RAM.
- **Capped Shard Writer Pool**: Sized from memory budget and capped at **6 parallel threads** by default (`--max-writer-threads 6`) to prevent disk I/O queue depth bottlenecks while maximizing parallel compression.

### E. TAPS Chemistry & Amethyst Schema
- **TAPS Calling Logic**:
  - Original Top strand (`OT`, `+`): Reference `C` $\to$ Read `T` is **Methylated** (`c += 1`), Read `C` is **Unmethylated** (`t += 1`).
  - Original Bottom strand (`OB`, `-`): Reference `G` $\to$ Read `A` is **Methylated** (`c += 1`), Read `G` is **Unmethylated** (`t += 1`).
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

# Run full test suite including BAM integration smoke test
pytest tests/ -v
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
