# taps-sc-extract

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC_BY--NC_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

**`taps-sc-extract`** is a high-throughput, memory-bounded, parallel single-cell DNA methylation extraction tool designed for TAPS (TET-Assisted Pyridine Borane Sequencing / mC-to-T chemistry) coordinate-sorted BAM files.

It outputs base-resolution methylation datasets directly into **Amethyst** (and Facet) compatible HDF5 formats (`/CG/<barcode>/1` and `/CH/<barcode>/1`), supporting both single monolithic HDF5 files and multi-file sharded directories with portable `master.h5` index files.

---

## Key Features

- **Blistering Performance**: Rust acceleration core processes **>15.4 million mapped reads/minute (256,400+ reads/s)** (genome-wide mouse mm10 with 75M reads and 1.11B base calls extracted in **~4.86 minutes** across 24 workers, **2.61× faster** than the multiprocessing Python reference engine).
- **Predictable, Bounded Memory**: Process-tree memory tracks accurately via `/proc` and plateaus at steady state (~700 MB / worker, ~18 GB total across 24 workers) regardless of genome length.
- **Auto-Scaled Shard Writer Pool**: Dynamically sizes writer concurrency (`--max-writer-threads 0` auto) to compress all output shards concurrently without bottlenecking disk queue depth.
- **High-Speed HDF5 Compression Options**: Supports `gzip` (level 1 deflate, 100% portable), `gzip-shuffle` (byte shuffling + level 1 deflate, 25.5× compression ratio), `lzf` (near raw write speed with 87% size reduction), `blosc` (multithreaded LZ4), `blosc-zstd`, and uncompressed `none`.
- **Process-Safe Indexed FASTA Reader**: Custom byte-seeking `.fai` reader eliminates BGZF memory collisions across dozens of multiprocessing workers.
- **Single-Pass HDF5 Writer**: Bypasses HDF5 B-tree resizing overhead, writing tens of thousands of cell datasets in **seconds** rather than hours.
- **Multi-File Sharded Directories**: Partitions cell barcodes across $N$ parallel shard files (`shard_000.h5`..`shard_NNN.h5`) and creates a portable `master.h5` with relative `ExternalLink` references.
- **Native Amethyst Compatibility**: Output files can be loaded directly into R with `amethyst::createObject()` and `amethyst::indexChr()`.

---

## Installation

### Using pip
```bash
git clone https://github.com/boconnell89/taps-sc-extract.git
cd taps-sc-extract
pip install -e .
```

### Using Conda / Mamba
```bash
mamba create -n taps_extract python=3.11 pysam h5py numpy pytest -y
conda activate taps_extract
pip install -e .
```

Verify the installation:
```bash
taps-sc-extract --help
```

---

## Quickstart

### 1. Standard Parallel Run (Auto Engine: Uses Compiled Rust Core)
```bash
taps-sc-extract \
  -b /path/to/aligned_taps.srt.bam \
  -f /path/to/reference.fa \
  -o all_cells.h5 \
  -t 24 \
  --chunk-size-mb 10 \
  --log-file extraction.log
```

### 2. High-Throughput Sharded Output (Recommended for >10k Cells)
```bash
taps-sc-extract \
  -b /path/to/aligned_taps.srt.bam \
  -f /path/to/reference.fa \
  -o /path/to/sharded_output/ \
  -t 24 \
  --shards 16 \
  --log-file extraction.log
```

### 3. Maximum Speed Run with LZF or Raw Writing
```bash
taps-sc-extract \
  -b /path/to/aligned_taps.srt.bam \
  -f /path/to/reference.fa \
  -o /path/to/sharded_output/ \
  -t 24 \
  --shards 16 \
  --compression lzf \
  --log-file extraction.log
```

---

## CLI Reference

```
usage: taps-sc-extract [-h] -b BAM -f FASTA -o OUT [-c CHROMS] [-w WHITELIST]
                       [--engine {auto,rust,python}] [--memory-mode {auto,stream,memory}]
                       [--max-memory-gb MAX_MEMORY_GB] [--expected-cells EXPECTED_CELLS]
                       [-t WORKERS] [--decomp-threads DECOMP_THREADS]
                       [--chunk-size-mb CHUNK_SIZE_MB] [--shards SHARDS]
                       [--compression {gzip,gzip-shuffle,gzip6,lzf,blosc,blosc-zstd,none}]
                       [--max-writer-threads MAX_WRITER_THREADS]
                       [--no-temp-file] [--temp-dir TEMP_DIR]
                       [--log-file LOG_FILE] [--min-baseq MIN_BASEQ]
                       [--min-mapq MIN_MAPQ] [--max-depth MAX_DEPTH]
                       [--no-baq] [--no-overlap-clip] [-v] [--version]
```

### Argument Details

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `-b, --bam` | `str` | *Required* | Path to coordinate-sorted, indexed BAM (`.bam` + `.bai`). |
| `-f, --fasta` | `str` | *Required* | Path to reference FASTA (`.fa` + `.fai`). |
| `-o, --out` | `str` | *Required* | Output `.h5` file path or output directory (when `--shards > 1`). |
| `-c, --chroms` | `str` | `None` | Comma-separated contigs (e.g. `chr1,chr2` or `chr19`). Default: canonical autosomes + `chrX`, `chrY`. |
| `-w, -a, --whitelist` | `str` | `None` | Path to optional barcode whitelist or cell annotation file. |
| `--engine` | `str` | `auto` | Backend engine: `auto` (prefers Rust, fallback Python), `rust` (fastest), `python` (reference). |
| `--memory-mode` | `str` | `auto` | Memory mode: `auto` (heuristic based on budget/cells), `stream` (disk temp), `memory` (RAM). |
| `--max-memory-gb` | `float` | `None` | Optional RAM budget in GB (default: 0.6 × available memory). |
| `--expected-cells` | `int` | `None` | Optional cell count (auto-detected from `-w/--whitelist` cardinality if omitted). |
| `-t, --threads, --workers`| `int` | `24` | Number of parallel chunk worker processes (or Rayon threads in Rust). `0` = auto. |
| `--decomp-threads` | `int` | `0` | BAM BGZF decompression threads per worker process (`0` for synchronous reading). |
| `--chunk-size-mb` | `int` | `10` | Genomic window chunk size in megabases. |
| `--shards` | `int` | `1` | Number of output HDF5 shard files to write in parallel (`0` = auto: 1/8/16/32). |
| `--compression` | `str` | `gzip` | Compression: `gzip` (level 1 deflate), `gzip-shuffle`, `gzip6`, `lzf`, `blosc`, `blosc-zstd`, `none`. |
| `--max-writer-threads` | `int` | `0` | Maximum parallel shard-writer threads after extraction (`0` = auto-scales up to shard count). |
| `--no-temp-file` | `flag` | `False` | In-memory mode: keeps chunk results in RAM instead of disk (equivalent to `--memory-mode memory`). |
| `--temp-dir` | `str` | `None` | Custom directory for temporary chunk streaming (default: `/tmp`). |
| `--log-file` | `str` | `None` | Path to output log file for timestamps and performance tracking. |
| `--min-baseq` | `int` | `20` | Minimum base quality for pileup. |
| `--min-mapq` | `int` | `0` | Minimum mapping quality for alignment filtering. |
| `--max-depth` | `int` | `250` | Maximum pileup depth. |
| `--no-baq` | `flag` | `False` | Disable Base Alignment Quality (BAQ) computation. |
| `--no-overlap-clip` | `flag` | `False` | Do not ignore overlapping mate read bases. |
| `-v, --verbose` | `flag` | `False` | Enable verbose debug logging. |

---

## Performance Benchmark & Flags Guide

### Whole-Genome mm10 Performance Comparison (74.9M Reads, 1.11 Billion Calls)

| Engine / Mode | Wall Time | Elapsed Job Time | Extraction Throughput | Peak RAM | Call Parity |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Python Reference Engine** | 12.70 min (762 s) | 12.70 min | 98.3k reads/s | 18.1 GB | 7,355 cells |
| **Rust Acceleration Core** | **5.00 min (300 s)** | **4.86 min (292 s)** | **256.4k reads/s** | **31.5 GB** | **100% Match (7,355 cells)** |
| **Speedup** | **2.54× faster** | **2.61× faster** | **2.61× faster** | Bounded | Exact Bitwise Match |

---

### Compression Algorithm Trade-Offs (`--compression`)

| Algorithm | Extraction + Write Time (`chr1–3`) | Shard Output Size | Compression Ratio | R / Bioconductor Requirements |
| :--- | :---: | :---: | :---: | :--- |
| **`none`** | **46.47 s** | 11.0 GB | 1.0× | Native `rhdf5` (No plugins needed) |
| **`lzf`** | **47.88 s** *(Fastest Compressed)* | 1.4 GB | **7.8×** | `library(rhdf5filters)` |
| **`blosc`** (LZ4) | **62.82 s** | 532 MB | **20.7×** | `library(rhdf5filters)` |
| **`gzip`** *(Default)* | **66.53 s** | 728 MB | **15.1×** | **Native `rhdf5` (Zero extra packages)** |
| **`gzip-shuffle`** | **76.79 s** | 431 MB | **25.5×** | Native `rhdf5` |
| **`blosc-zstd`** | **74.71 s** | **385 MB** *(Smallest)* | **28.6×** | `library(rhdf5filters)` |

### 2. `-t, --workers` (Parallel Worker Threads/Processes)
- **Mechanism**: Splits the genome into $N$ non-overlapping genomic windows and processes them concurrently.
- **Recommendation**: Set to $0.75 \times$ to $1.0 \times$ available physical CPU cores (e.g., `24` on a 32-core workstation; `64` on a 96-core server). `0` automatically sizes workers from hardware CPU count and memory budget.
- **Scaling**: Near-linear scaling up to 64 workers.

### 3. `--decomp-threads` (BAM BGZF Decompression Threads)
- **Mechanism**: Allocates background `htslib` threads per worker for BGZF block decompression.
- **Tuning by Worker Count**:
  - **High Worker Count ($\ge 8$ workers, e.g. `-t 24` or `-t 32`)**: Use `--decomp-threads 0` or `1`. Because all CPU cores are already fully saturated processing genomic windows, spawning 24–48 additional threads causes thread oversubscription, context switching, and cache eviction.
  - **Low Worker Count ($\le 4$ workers, e.g. single-chromosome extractions)**: Use `--decomp-threads 2` to `4`. Spare CPU cores are utilized to decompress BAM blocks ahead of pileup.

### 3. `--chunk-size-mb` (Genomic Granularity & Load Balancing)
- **Mechanism**: Defines the window size of each genomic chunk (default: `10` Mb $\to$ 286 chunks for mm10).
- **Trade-offs**:
  - **Smaller chunks (e.g., `5` Mb)**: Better load balancing across many workers; prevents workers from idling at the end of chromosomes.
  - **Larger chunks (e.g., `20` Mb)**: Reduces chunk file count and IPC overhead; slightly higher memory usage per chunk.
- **Recommendation**: `10` Mb is optimal for mammalian genomes (mm10 / hg38).

### 4. `--shards` (Parallel HDF5 Writing & Scalability)
- **Mechanism**: Partitions cell barcodes deterministically across $N$ shard files (`shard_000.h5`..`shard_NNN.h5`) using `hash(bc) % N` and creates `master.h5` with relative `ExternalLink` references.
- **Benefits**:
  - Eliminates single-file HDF5 write lock contention.
  - Assembles and compresses all $N$ shards **concurrently** in parallel.
  - Reduces individual file sizes for easier downstream parallelization in R/Amethyst.
- **Recommendation**:
  - $\le 10,000$ cells: `--shards 1` (single file) or `--shards 4`.
  - $10,000$–$50,000$ cells: `--shards 8` or `--shards 16`.
  - $> 100,000$ cells: `--shards 16` or `--shards 32`.

### 5. `--max-writer-threads` (Shard Writer Concurrency & Storage Optimization)
- **Mechanism**: Controls the maximum number of parallel worker threads that assemble, compress, and write shard HDF5 files after chunk extraction finishes.
- **Why it matters for large genomes**: On whole-genome mammalian datasets (>1 billion calls across 16 shards), gzip compression in the writer phase can become the dominant bottleneck. Increasing `--max-writer-threads` from 6 to 16 on high-core NVMe systems allows all shards to compress concurrently, cutting shard assembly time from ~3.5 minutes down to ~1 minute.
- **Recommended Defaults by Storage Type**:
  - **Local Fast NVMe / PCIe SSD with $\ge 32$ GB RAM**: `--max-writer-threads 8` to `16` &mdash; Compresses all shards simultaneously for maximum throughput.
  - **Standard Workstations (Default)**: `--max-writer-threads 6` &mdash; Safe balance between parallel CPU compression and disk I/O queue depth.
  - **Network File Systems (NFS / Lustre / GPFS / SMB)**: `--max-writer-threads 2` to `4` &mdash; Minimizes concurrent metadata and lock contention across shared storage nodes.
  - **Spinning Disk (HDD) or Memory-Constrained Systems**: `--max-writer-threads 1` &mdash; Strict single-thread sequential writing with zero seek contention and minimal RAM.

### 6. `--no-temp-file` vs. Default Streaming Mode
- **Default Mode (Disk-Streaming)**:
  - Workers write compact binary chunk files to fast disk (`/tmp`) as each genomic window finishes.
  - Each shard reads only its own chunk files in coordinate order and purges them immediately upon writing.
  - **RAM Usage**: Process-tree RAM scales as $\approx 700\text{ MB / worker}$ (fixed heap buffers + BAM indices) and plateaus at steady state (~18 GB total across 24 workers), never growing linearly with genome size.
  - **Writer Pool**: Automatically sizes writer concurrency from `MemAvailable` and respects `--max-writer-threads`.
  - **Recommendation**: Always use for workstations, desktops, and standard production servers.
- **`--no-temp-file` (In-Memory Mode)**:
  - Workers pass structured array dictionaries directly over IPC without writing to disk.
  - **RAM Usage**: Accumulates all chunk arrays directly in memory across the entire run (~25–35 GB for whole-genome mammalian datasets).
  - **Recommendation**: Use on high-memory servers ($\ge 128$ GB RAM, 64–96 cores) with shared RAM to maximize I/O throughput.

### 7. `--temp-dir` (Scratch Space Location)
- **Mechanism**: Specifies the filesystem location where temporary chunk files are stored during extraction.
- **Recommendation**: Point to a fast local NVMe SSD or memory-backed filesystem (e.g. `/tmp` or `/dev/shm`). Avoid slow network filesystems (NFS/Lustre) for temporary files.

---

## TAPS Chemistry & Calling Logic

TAPS uses TET oxidation followed by pyridine borane reduction to convert methylated cytosines ($5\text{mC}$ and $5\text{hmC}$) into dihydrouracil ($\text{DHU}$), which is read by DNA polymerase as **Thymine ($\text{T}$)**. Unmethylated cytosines remain **Cytosine ($\text{C}$)**.

$$\text{Methylated } 5\text{mC} / 5\text{hmC} \xrightarrow{\text{TAPS}} \text{T}$$
$$\text{Unmethylated } \text{C} \xrightarrow{\text{TAPS}} \text{C}$$

### Strand & Base Interpretation

| Alignment Strand | Reference Base | Read Base | Interpretation | Methylation State |
| :--- | :---: | :---: | :---: | :---: |
| **OT** (Original Top, Strand `+`) | `C` | `T` | Modified ($5\text{mC}/5\text{hmC}$) | **Methylated** (`c += 1`) |
| **OT** (Original Top, Strand `+`) | `C` | `C` | Unmodified ($\text{C}$) | **Unmethylated** (`t += 1`) |
| **OB** (Original Bottom, Strand `-`) | `G` | `A` | Modified ($5\text{mC}/5\text{hmC}$) | **Methylated** (`c += 1`) |
| **OB** (Original Bottom, Strand `-`) | `G` | `G` | Unmodified ($\text{C}$) | **Unmethylated** (`t += 1`) |

*Note: In the Amethyst HDF5 output schema, `t` stores unmethylated counts and `c` stores methylated counts.*

---

## HDF5 Output Structure & Schema

The output HDF5 matches the Amethyst 1.0+ / Facet specification:

```
all_cells.h5 (or master.h5)
├── metadata/
│   └── version = "amethyst2.0.0"
├── CG/
│   ├── <barcode_1>/
│   │   └── 1  (Dataset: structured array)
│   └── <barcode_2>/
│       └── 1  (Dataset: structured array)
└── CH/
    ├── <barcode_1>/
    │   └── 1  (Dataset: structured array)
    └── <barcode_2>/
        └── 1  (Dataset: structured array)
```

### Dataset Dtype Specification
Each `/1` dataset is a coordinate-sorted NumPy structured array with dtype:
```python
dtype = [
    ('chr', 'S10'),   # Chromosome name (ASCII bytes, e.g. b'chr1')
    ('pos', '<i8'),   # 1-based genomic coordinate
    ('pct', '<f8'),   # Methylation percentage: 100.0 * c / (c + t)
    ('t',   '<i8'),   # Unmethylated read count
    ('c',   '<i8'),   # Methylated read count
]
```

---

## Amethyst R Integration Guide

### Loading a Single HDF5 File
```R
library(amethyst)
library(rhdf5)

h5_path <- "all_cells.h5"

# List barcodes under /CG
ls_df <- h5ls(h5_path)
barcodes <- ls_df$name[ls_df$group == "/CG"]

# Create Amethyst object
h5_paths <- data.frame(path = h5_path, barcode = barcodes)
obj <- createObject(h5paths = h5_paths)

# Index chromosomes
obj_cg <- indexChr(obj, "CG")
obj_ch <- indexChr(obj, "CH")
```

### Loading a Sharded Directory (`master.h5`)
```R
library(amethyst)
library(rhdf5)

master_path <- "sharded_output/master.h5"

# List barcodes from master index
ls_df <- h5ls(master_path)
barcodes <- ls_df$name[ls_df$group == "/CG"]

# Create Amethyst object pointing to master.h5
h5_paths <- data.frame(path = master_path, barcode = barcodes)
obj <- createObject(h5paths = h5_paths)

# Index chromosomes (resolves all external links transparently)
obj_cg <- indexChr(obj, "CG")
obj_ch <- indexChr(obj, "CH")
```

---

## Running Tests

Run the test suite with `pytest`:

```bash
pytest tests/ -v
```

---

## License
 
This project is licensed under the **Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)** License. See the [LICENSE](LICENSE) file for details.
 
Free for non-commercial academic, scientific, and educational research use. Commercial use requires explicit permission.
