# taps-sc-extract

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC_BY--NC_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

**`taps-sc-extract`** is a high-throughput, memory-bounded, parallel single-cell DNA methylation extraction tool designed for TAPS (TET-Assisted Pyridine Borane Sequencing / mC-to-T chemistry) coordinate-sorted BAM files.

It outputs base-resolution methylation datasets directly into **Amethyst** (and Facet) compatible HDF5 formats (`/CG/<barcode>/1` and `/CH/<barcode>/1`), supporting both single monolithic HDF5 files and multi-file sharded directories with portable `master.h5` index files.

---

## Key Features

- **Rust acceleration core**: Window-parallel extraction with a Python reference engine behind the same CLI. Current workstation timings: [`code/nvme-profiling-2026-10-15.md`](code/nvme-profiling-2026-10-15.md).
- **Bounded memory**: The run summary reports the extraction process's peak RSS and, separately, the shard-writer processes' combined peak.
- **One writer process per shard**: hdf5-metno serializes every HDF5 call in a process behind one lock. With more than one shard, the Rust engine writes each shard in its own child process (`--writer-mode process`, default), so shards compress and write in parallel.
- **Machine tuning (`--test`)**: Measures this machine on a built-in synthetic workload, optionally samples your BAM, and recommends `-t`, decompression and writer settings. `--save-profile` makes them the defaults on this host.
- **High-Speed HDF5 Compression Options**: Supports `gzip` (level 1 deflate, 100% portable), `gzip-shuffle` (byte shuffling + level 1 deflate, 25.5× compression ratio), `lzf` (near raw write speed with 87% size reduction), `blosc` (multithreaded LZ4), `blosc-zstd`, and uncompressed `none`.
- **Process-Safe Indexed FASTA Reader**: Custom byte-seeking `.fai` reader eliminates BGZF memory collisions across dozens of multiprocessing workers.
- **Single-Pass HDF5 Writer**: Bypasses HDF5 B-tree resizing overhead, writing tens of thousands of cell datasets in **seconds** rather than hours.
- **Multi-File Sharded Directories**: Partitions cell barcodes across $N$ parallel shard files (`shard_000.h5`..`shard_NNN.h5`) and creates a portable `master.h5` with relative `ExternalLink` references.
- **Native Amethyst Compatibility**: Output files can be loaded directly into R with `amethyst::createObject()` and `amethyst::indexChr()`.

---

## Installation

### Rust core, built for this machine (recommended)
```bash
git clone https://github.com/boconnell89/taps-sc-extract.git
cd taps-sc-extract
conda install -c conda-forge hdf5 clang libclang pkg-config   # build dependencies
./install.sh            # native build for this CPU, then pip install -e .
```
`install.sh` downloads hts-sys 2.2.1 and patches in htslib 1.21's read-pair overlap clipping, so Rust calls match the Python engine. It writes `taps-sc-extract-rs/.cargo/config.toml` for your conda prefix and builds with `-C target-cpu=native` (Rust) and `-O3 -march=native` (the vendored htslib, including BAQ). A native binary can fail with SIGILL on a different CPU, so **rerun `install.sh` on each host**, or use `./install.sh --portable` for an x86-64-v3 build. `./install.sh --smoke-bam BAM --smoke-fasta FA` also extracts chr19:5.0–5.1 Mb as a check.

### Python reference engine only (pip)
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

### 0. Tune this machine (once per host)
```bash
taps-sc-extract --test --save-profile                     # machine only, a few minutes
taps-sc-extract --test -b sample.srt.bam -f ref.fa -o out/   # also predict time and memory for this BAM
```
The report names the limiting factor for each phase (CPU, BAM read, HDF5 write, storage) and recommends flags. `--test-quick` takes about 40 seconds.

### 1. Standard Parallel Run (Auto Engine: Uses Compiled Rust Core)
```bash
taps-sc-extract \
  -b /path/to/aligned_taps.srt.bam \
  -f /path/to/reference.fa \
  -o all_cells.h5 \
  -t 24 \
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
                       [--no-baq] [--no-overlap-clip]
                       [--trim-r2 TRIM_R2] [--trim-r1 TRIM_R1]
                       [--writer-mode {auto,process,thread}] [--profile-json PROFILE_JSON]
                       [--test] [--test-quick] [--save-profile] [--test-report TEST_REPORT]
                       [-v] [--version]
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
| `--max-memory-gb` | `float` | `None` | Memory budget in GB for the whole run, extraction plus writer processes. Overrides `--max-memory-fraction`. |
| `--max-memory-fraction` | `float` | `0.6` | Rust: share of available memory the run may use. |
| `--expected-cells` | `int` | `None` | Optional cell count (auto-detected from `-w/--whitelist` cardinality if omitted). |
| `-t, --threads, --workers`| `int` | `24` | Number of parallel chunk worker processes (or Rayon threads in Rust). `0` = auto. |
| `--decomp-threads` | `int` | Rust: 4 / 2 / 0 by worker count; Python: 1 | BAM BGZF decompression threads per worker (`0` = inflate on the worker thread). |
| `--chunk-size-mb` | `int` | `10` Python / auto Rust | Genomic window size in megabases. Rust auto: about 6 windows per worker, 2–100 Mb. |
| `--shards` | `int` | `1` | Number of output HDF5 shard files to write in parallel (`0` = auto: 1/8/16/32). An `.h5` file name with multiple shards is written to a directory named after the file stem. |
| `--compression` | `str` | `gzip` Python / `lzf` Rust | Compression: `gzip` (level 1 deflate), `gzip-shuffle`, `gzip6`, `lzf`, `blosc`, `blosc-zstd`, `none`. |
| `--max-writer-threads` | `int` | `0` | Maximum concurrent shard writers after extraction (`0` = auto: up to shard count, cores and memory; in-process writers capped at 6). |
| `--writer-mode` | `str` | `auto` | Rust shard writers: `process` (one child per shard, default when shards > 1) or `thread` (in-process, one shared HDF5 lock). |
| `--profile-json` | `str` | `None` | Rust: write per-phase CPU use, per-stage and per-window timings and writer memory to JSON. Adds ~5–10% to processing time. |
| `--test` | `flag` | `False` | Run the tuning probes instead of extracting (see Quickstart 0). `-b/-f` optional. |
| `--test-quick` / `--save-profile` / `--test-report` | | | Shorter probes; save the machine profile for later runs on this host; also write the report to a file. |
| `--no-temp-file` | `flag` | `False` | In-memory mode: keeps chunk results in RAM instead of disk (equivalent to `--memory-mode memory`). |
| `--temp-dir` | `str` | `None` | Custom directory for temporary chunk streaming (default: `/tmp`). |
| `--log-file` | `str` | `None` | Path to output log file for timestamps and performance tracking. |
| `--min-baseq` | `int` | `20` | Minimum base quality for pileup. |
| `--min-mapq` | `int` | `0` | Minimum mapping quality for alignment filtering. |
| `--max-depth` | `int` | `250` | Maximum pileup depth. |
| `--no-baq` | `flag` | `False` | Disable Base Alignment Quality (BAQ) computation. |
| `--no-overlap-clip` | `flag` | `False` | Do not ignore overlapping mate read bases. |
| `--trim-r2` | `int` | `10` | Drop this many 5′ sequenced bases of read 2 before calling (Tn5/TAPS M-bias). |
| `--trim-r1` | `int` | `0` | Drop this many 5′ sequenced bases of read 1 before calling. |
| `-v, --verbose` | `flag` | `False` | Enable verbose debug logging. |

---

## Performance Flags Guide

Current workstation measurements (Rust engine, `5base_ucbTn5_sp2.srt.bam`, 53.5M reads, 7,351 cells, mm10, 16-core / 32-thread Ryzen 9 5900XT, WSL2 ext4 on NVMe): [`code/nvme-profiling-2026-10-15.md`](code/nvme-profiling-2026-10-15.md). The 2026-10-14 starting point is in [`code/nvme-profiling-2026-10-14.md`](code/nvme-profiling-2026-10-14.md). The run summary prints processing time, HDF5 write time, estimated read and write MB/s, and peak memory of the extraction process and the writer processes. **Throughput (Reads)** is mapped records over wall time, including HDF5 write. For a new machine, start with `taps-sc-extract --test`.

### 1. Compression (`--compression`)

Measured on chr1–3 before `026c52e` with the older in-process writers. Ratios and plugin requirements still apply; the times do not.

| Algorithm | Extraction + Write Time (`chr1–3`) | Shard Output Size | Compression Ratio | R / Bioconductor Requirements |
| :--- | :---: | :---: | :---: | :--- |
| **`none`** | 46.47 s | 11.0 GB | 1.0× | Native `rhdf5` (No plugins needed) |
| **`lzf`** *(Rust default)* | 47.88 s | 1.4 GB | **7.8×** | `library(rhdf5filters)` |
| **`blosc`** (LZ4, level 9: high compression) | 62.82 s | 532 MB | **20.7×** | `library(rhdf5filters)` |
| **`gzip`** *(Python default)* | 66.53 s | 728 MB | **15.1×** | **Native `rhdf5` (Zero extra packages)** |
| **`gzip-shuffle`** | 76.79 s | 431 MB | **25.5×** | Native `rhdf5` |
| **`blosc-zstd`** | 74.71 s | **385 MB** *(Smallest)* | **28.6×** | `library(rhdf5filters)` |

### 2. `-t, --workers`
- **Mechanism**: Rayon threads, each processing genomic windows. Extraction is CPU-bound; with BAQ on, BAQ (`sam_prob_realn`) is about 70% of worker time.
- **Recommendation**: Use all logical CPUs. On the 16-core / 32-thread test machine, 32 workers ran about 1.4× faster than 16 (SMT). `0` = auto (CPUs, bounded by the memory budget), or the saved `--test` profile.

### 3. `--decomp-threads` (BAM BGZF decompression)
- **Mechanism**: htslib threads per worker that inflate BGZF blocks ahead of the reader.
- **Recommendation**: Leave it unset. Rust picks 0 at 8 or more workers (inflate runs on the worker thread), 2 at 5–7 and 4 at 1–4. With libdeflate, inflate measured 1–6% of worker time, so extra threads only add contention.

### 4. `--chunk-size-mb` (window size)
- **Rust default: auto.** About 6 windows per worker, clamped to 2–100 Mb (mm10: ~28 Mb at 16 workers, ~14 Mb at 32, ~4.7 Mb at 96). Enough windows keep every worker busy to the end; the old fixed 100 Mb gave 37 windows, so 32 workers had a long tail and 96 workers would sit mostly idle.
- Window size does not change calls (verified per cell at 100, 10 and 2 Mb).
- The Python engine defaults to `10` Mb.

### 5. `--shards`
- **Mechanism**: Partitions barcodes across $N$ shard files by `hash(bc) % N` and writes `master.h5` with relative `ExternalLink`s.
- **Recommendation**: Depends on the filesystem. On local NVMe, 16 shards wrote 4.2 GB in 11–17 s against 62 s for one shard. On a server filesystem with fewer cells, more shards were much slower, so measure there (`--test` probes the output filesystem). The default stays 1 shard below 10,000 cells.

### 6. `--max-writer-threads` and `--writer-mode`
- **`--writer-mode process`** (default with more than one shard): one `taps-sc-extract-rs write-shard` child per shard, up to `--max-writer-threads` at a time. Each child streams its shard one cell at a time in bounded slices, so a writer needs roughly 15–60 bytes per site of its shard.
- **`--writer-mode thread`**: in-process writers. hdf5-metno takes one process-wide lock for every HDF5 call, so they do not write in parallel; auto caps them at 6.
- **Auto count**: up to min(shards, CPUs), and only as many as fit the memory budget at about 20 bytes per site of a shard each, or the saved `--test` profile (still memory-bounded). Lower it on network filesystems (NFS/Lustre: 2–4) or spinning disks (1).

### 7. `--memory-mode` / `--no-temp-file`
- **`stream`**: window results go to compact temp files under `--temp-dir`; each writer reads its own shard's files. Lowest memory (about 5 GB peak for the extraction process on the test BAM).
- **`memory`**: window results stay in RAM and are piped to the writer processes (about 12 GB on the test BAM). Not faster than `stream` on local NVMe; use it when temp storage is slow.
- **`auto`**: the Rust engine estimates sites from the BAM index's mapped-read count (about 14 per read). It picks memory mode only if the extraction process plus the writer processes fit the budget without fewer writers than stream mode would allow; otherwise stream. The writer count is bounded by the same budget. The config block prints the plan (`Memory plan: ... est. peak ~X GB of Y GB`). Set the budget with `--max-memory-fraction` (default 0.6 of available memory) or `--max-memory-gb`.

### 8. `--temp-dir`
- Point it at fast local storage. Keep inputs, temp and output off slow network or 9p mounts: on the test machine the same job took 171 s from a 9p-mounted HDD against 95 s from local NVMe (before these upgrades).

### 9. Profiling a run (`--profile-json`)
Writes per-phase wall and CPU time, utilisation, a per-stage split of worker time (read/inflate, BAQ, pileup, calling), per-window timings and writer memory. Use it to see what limits a run on a new system.

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

Run the test suite with `pytest` and `cargo`:

```bash
pytest tests/ -v                        # includes chr19 BAM tests when the data is present
(cd taps-sc-extract-rs && cargo test --release)
```

The chr19 tests check the Rust engine bit for bit against a frozen reference (`test_chr19_rust_overlapfix_2026-10-15/`), and require Rust `--no-baq` to equal the Python engine at every cell and position.

---

## License
 
This project is licensed under the **Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)** License. See the [LICENSE](LICENSE) file for details.
 
Free for non-commercial academic, scientific, and educational research use. Commercial use requires explicit permission.
