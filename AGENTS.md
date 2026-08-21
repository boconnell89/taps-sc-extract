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
├── taps_sc_extract/
│   ├── __init__.py             # Version definition
│   ├── __main__.py             # Module execution entrypoint (python -m taps_sc_extract)
│   ├── cli.py                  # CLI argument parser and logging setup
│   ├── parallel_extractor.py   # High-throughput multiprocessing extraction engine
│   ├── extractor.py            # Single-process reference extractor
│   ├── h5_writer.py            # Amethyst HDF5 writer & sharded directory manager
│   ├── fasta.py                # FastFaiReader: process-safe .fai indexed FASTA reader
│   ├── calling.py              # TAPS mC-to-T lookup tables & SAM flag classification
│   ├── context.py              # Sequence context classifier (CpG, CHG, CHH, CG, CH)
│   └── barcode.py              # Barcode extraction from QNAME/tags & whitelist parsing
└── tests/
    └── test_taps_sc_extract.py # Pytest unit & integration test suite
```

---

## 3. Critical Invariants & Rules

When modifying or extending this codebase, you **MUST** adhere to the following rules:

### A. Multiprocessing & File Descriptors
- **Always use `mp.get_context("spawn")`**: Never use `fork` with `pysam` or `h5py`. Forked C file descriptors cause deadlocks and corrupted reads.
- **Never pass open file handles to worker processes**: Workers must initialize their own private `pysam.AlignmentFile` and `FastFaiReader` handles in `_init_worker`.

### B. FASTA Reading Across Workers
- **Never use `pysam.FastaFile` in parallel workers**: `htslib`'s underlying `bgzf_read_block` is not safe when multiple processes read the same indexed FASTA, causing C-level memory corruption.
- **Always use [`taps_sc_extract.fasta.FastFaiReader`](file:///mnt/e/sciMET_TAPS/taps_sc_extract/fasta.py)**: It performs pure Python binary `seek()` and `read()` using `.fai` byte offsets, which is 100% process-safe.

### C. HDF5 Single-Pass Creation vs. Resizing
- **Never call `dataset.resize()` in a tight loop**: For 10,000+ cells, incremental resizing causes massive B-tree metadata churn (>20 minutes vs. <10 seconds).
- **Always use single-pass dataset creation (`create_cell_dataset`)**: Pre-aggregate calls per cell/shard and call `create_dataset('1', data=records, ...)` once with the final exact shape.
- **Chunk size constraint**: When creating a chunked dataset in `h5py`, chunk dimension cannot exceed data length: `chunk_size = min(len(records), 65536)`.

### D. Memory Bounding & Shard Writer Concurrency
- **Accurate Process-Tree Accounting**: Memory is monitored across the entire process tree (`/proc/<pid>/statm`) including all worker children. Each worker reaches a steady-state buffer (~700 MB) and plateaus.
- **Never load all genomic chunk arrays into memory at once**: In disk-streaming mode (`use_temp_files=True`), workers write chunk outputs partitioned into `temp_dir/shard_XXX/chunk_YYYYYY.bin`.
- **Capped Shard Writer Pool**: `_shard_writer_concurrency()` sizes writer threads from `MemAvailable` and caps at **6 parallel threads** to prevent disk I/O queue depth congestion while maximizing parallel CPU compression. Each shard purges its temporary chunk files immediately upon writing.

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
# 1. Standard parallel extraction
taps-sc-extract -b sample.srt.bam -f ref.fa -o output.h5 -t 24

# 2. Sharded directory output (for large datasets)
taps-sc-extract -b sample.srt.bam -f ref.fa -o sharded_dir/ -t 24 --shards 16 --log-file extraction.log

# 3. High-memory in-memory mode
taps-sc-extract -b sample.srt.bam -f ref.fa -o output.h5 -t 64 --no-temp-file
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
