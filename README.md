# taps-sc-extract

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC_BY--NC_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

**`taps-sc-extract`** is a high-throughput, memory-bounded, parallel single-cell DNA methylation extraction tool designed for TAPS (TET-Assisted Pyridine Borane Sequencing / mC-to-T chemistry) coordinate-sorted BAM files.

It outputs base-resolution methylation datasets directly into **Amethyst** (and Facet) compatible HDF5 formats (`/CG/<barcode>/1` and `/CH/<barcode>/1`), supporting both single monolithic HDF5 files and multi-file sharded directories with portable `master.h5` index files.

---

## Key Features

- **Rust acceleration core**: Window-parallel extraction with a Python reference engine behind the same CLI. Measured speed: [Performance](#performance).
- **Bounded memory**: The run summary reports the extraction process's peak RSS and, separately, the shard-writer processes' combined peak.
- **One writer process per shard**: hdf5-metno serializes every HDF5 call in a process behind one lock. With more than one shard, the Rust engine writes each shard in its own child process (`--writer-mode process`, default), so shards compress and write in parallel.
- **Machine tuning (`--test`)**: Measures this machine on a built-in synthetic workload, optionally samples your BAM, and recommends `-t`, decompression and writer settings. `--save-profile` makes them the defaults on this host.
- **HDF5 compression options**: `blosc` (Rust default; LZ4 level 9), `lzf`, `gzip` (Python default; readable by plain `rhdf5`), `gzip-shuffle`, `gzip6`, `blosc-zstd` (smallest) and uncompressed `none`. See [Compression](#1-compression---compression).
- **Context selection (`--context CG|CH`)**: extract only CpG (or only CH). CG-only skips the CH columns, writes no `/CH` group, and runs at 6.0M reads/s end to end on the test data, against 2.45M reads/s for both contexts.
- **Process-Safe Indexed FASTA Reader**: Custom byte-seeking `.fai` reader eliminates BGZF memory collisions across dozens of multiprocessing workers.
- **Single-Pass HDF5 Writer**: Writes each cell's dataset once, in bounded slices, without HDF5 dataset resizing.
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

BAQ runs on SIMD kernels that process 8 reads per vector (AVX-512F) or 4 (AVX2). Both are always compiled, including in `--portable` builds, and the binary picks the widest the CPU runs. After building, `install.sh` runs `taps-sc-extract-rs baq-selftest`, which runs every supported kernel and htslib's scalar BAQ on the same reads and requires identical base qualities for every read; the install stops on any difference. Override with `--baq-simd avx512|avx2|off`. `taps-sc-extract --test --save-profile` times each kernel under full load and saves the fastest, since some Intel CPUs lower their clock for AVX-512.

BAQ can also run on a GPU (`--baq-device gpu`, off by default). The OpenCL backend is always built and needs, at run time, a GPU with double precision and its vendor's OpenCL runtime (Intel Arc / Data Center GPUs: Intel's compute runtime, e.g. the `intel-opencl-icd` package; check with `clinfo | grep -i fp64`). The CUDA backend is built when `nvcc` is on the PATH and the AMD HIP backend when `hipcc` is. When a usable GPU is present, `install.sh` also runs `baq-selftest --gpu` and stops on any difference. See [GPU BAQ](#11-gpu-baq---baq-device-gpu).

### Upgrading an existing install
```bash
cd taps-sc-extract
git checkout main
git pull
./install.sh            # same options as the first install (--prefix, --portable, --no-pip)
```
Rerunning `install.sh` is safe, and is required after every pull: it regenerates `.cargo/config.toml`, rebuilds the Rust core from clean, and reinstalls the Python CLI. It also checks the hts-sys patch. The patched copy in `taps-sc-extract-rs/vendor/` records the checksum of the patch it was built with, so a pull that changes `taps-sc-extract-rs/patches/` makes `install.sh` re-extract and re-patch it. The downloaded `hts-sys-2.2.1.crate` is kept in `vendor/`, so re-patching works without network. The first install on a machine needs `curl`, `patch`, `sha256sum` and access to crates.io.

Upgrading from a checkout older than 2026-10-15 (`026c52e` or earlier): `taps-sc-extract-rs/.cargo/config.toml` is no longer tracked, and `install.sh` now generates it. If you edited that file locally, `git pull` will refuse. Discard the edit with `git checkout -- taps-sc-extract-rs/.cargo/config.toml`, pull, and pass any non-default conda location as `./install.sh --prefix /path/to/env`. Outputs written by Rust builds before 2026-10-15 put about 0.01% of calls under the wrong barcode and differ slightly from the Python engine at read-pair overlaps, so re-extract them (see `code/rust-implementation-log.md`).

After installing or upgrading on a new machine, run `taps-sc-extract --test --save-profile` to measure it and save tuned defaults.

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
The report names the limiting factor for each phase (CPU, BAM read, HDF5 write, storage) and recommends flags. It probes the input, temp and output filesystems separately, measures BAQ on and off at the recommended worker count, compares memory-mode and stream-mode writing, and times a few very large cells. With `-b`, it also reports the largest cell's share of sites (one writer writes a whole cell, so that cell sets a floor on write time) and the predicted HDF5 MB/s as a share of the output disk's sequential write speed. `--test-quick` takes about 40 seconds.

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

### 3. CpG Only
```bash
taps-sc-extract \
  -b /path/to/aligned_taps.srt.bam \
  -f /path/to/reference.fa \
  -o /path/to/sharded_output/ \
  -t 24 \
  --shards 16 \
  --context CG \
  --log-file extraction.log
```
The output has `/CG` only. In Amethyst, index only that context: `indexChr(obj, "CG")` (indexing `"CH"` on it prints a missing-object error for every cell).

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
                       [--no-baq] [--baq {taps,taps-full,taps-partial,htslib}] [--baq-simd {auto,avx512,avx2,off}]
                       [--baq-device {cpu,gpu}] [--gpu-backend {auto,opencl,cuda,hip}]
                       [--gpu-device NAME] [--gpu-mem-mb MB]
                       [--context {both,CG,CH}] [--no-overlap-clip]
                       [--trim-r2 TRIM_R2] [--trim-r1 TRIM_R1]
                       [--writer-mode {auto,process,thread}] [--profile-json PROFILE_JSON]
                       [--test] [--test-quick] [--save-profile] [--test-report TEST_REPORT]
                       [--complexity COMPLEXITY] [--min-complexity N] [--max-complexity N]
                       [-v] [--version]
```

### Argument Details

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `-b, --bam` | `str` | *Required* | Coordinate-sorted alignments: BAM, CRAM (decoded against `-f`) or SAM (copied to a temporary BAM first). A missing index (`.bai`, `.csi` for contigs over 512 Mb, `.crai`) is built next to the file, or in the temp directory if that directory is read-only; the time is logged separately. CRAM decoding costs extra CPU per read. |
| `-f, --fasta` | `str` | *Required* | Path to reference FASTA (`.fa` + `.fai`). |
| `-o, --out` | `str` | *Required* | Output `.h5` file path or output directory (when `--shards > 1`). |
| `-c, --chroms` | `str` | `None` | Comma-separated contigs (e.g. `chr1,chr2` or `chr19`). Default: canonical autosomes + `chrX`, `chrY`. |
| `-w, -a, --whitelist` | `str` | `None` | Path to optional barcode whitelist or cell annotation file (first column = barcode; `.gz` accepted). |
| `--complexity` | `str` | `None` | Complexity file: `rank barcode total_reads unique_reads percent` per line (whitespace-separated, `.gz` accepted). Keeps barcodes whose unique-read count (column 4) is within `--min-complexity` / `--max-complexity`. With `-w`, a cell must pass both. |
| `--min-complexity` / `--max-complexity` | `float` | `None` | Inclusive unique-read bounds for `--complexity` (either may be omitted). Filtered cells are left out of the HDF5; kept cells' calls are identical to an unfiltered run. Rust also skips BAQ for filtered cells' reads. |
| `--engine` | `str` | `auto` | Backend engine: `auto` (prefers Rust, fallback Python), `rust` (fastest), `python` (reference). |
| `--memory-mode` | `str` | `auto` | Memory mode: `auto` (heuristic based on budget/cells), `stream` (disk temp), `memory` (RAM). |
| `--max-memory-gb` | `float` | `None` | Memory budget in GB for the whole run, extraction plus writer processes. Overrides `--max-memory-fraction`. |
| `--max-memory-fraction` | `float` | `0.6` | Rust: share of available memory the run may use. |
| `--expected-cells` | `int` | `None` | Optional cell count (auto-detected from the `-w` / `--complexity` filter's barcode count if omitted). |
| `-t, --threads, --workers`| `int` | `24` | Number of parallel chunk worker processes (or Rayon threads in Rust). `0` = auto. |
| `--decomp-threads` | `int` | Rust: 4 / 2 / 0 by worker count; Python: 1 | BAM BGZF decompression threads per worker (`0` = inflate on the worker thread). |
| `--chunk-size-mb` | `int` | `10` Python / auto Rust | Genomic window size in megabases. Rust auto: about 6 windows per worker, 2–100 Mb. |
| `--shards` | `int` | `1` | Number of output HDF5 shard files to write in parallel (`0` = auto: 1/8/16/32). An `.h5` file name with multiple shards is written to a directory named after the file stem. Rust balances shards by site count (barcodes are hashed into 4 buckets per shard, and buckets are packed by size), so which `shard_XXX.h5` holds a cell differs from the Python engine; `master.h5` links every cell either way. A single cell is never split, so the largest cell sets a floor on write time. |
| `--compression` | `str` | `gzip` Python / `blosc` Rust | Compression: `gzip` (level 1 deflate), `gzip-shuffle`, `gzip6`, `lzf`, `blosc` (LZ4 level 9), `blosc-zstd`, `none`. Reading `blosc`/`lzf` in R needs `rhdf5filters`; in Python, `import hdf5plugin` before `h5py` for `blosc`. |
| `--max-writer-threads` | `int` | `0` | Maximum concurrent shard writers after extraction (`0` = auto: up to shard count, cores and memory; in-process writers capped at 6). |
| `--writer-mode` | `str` | `auto` | Rust shard writers: `process` (one child per shard, default when shards > 1) or `thread` (in-process, one shared HDF5 lock). |
| `--profile-json` | `str` | `None` | Rust: write per-phase CPU use, per-stage and per-window timings (1 in 16 reads and pileup columns are timed and scaled up) and per-shard writer reports (read, create, fill, HDF5 write seconds, bytes, peak memory) to JSON. Adds ~2% to processing time. |
| `--test` | `flag` | `False` | Run the tuning probes instead of extracting (see Quickstart 0). `-b/-f` optional. |
| `--test-quick` / `--save-profile` / `--test-report` | | | Shorter probes; save the machine profile for later runs on this host; also write the report to a file. |
| `--no-temp-file` | `flag` | `False` | In-memory mode: keeps chunk results in RAM instead of disk (equivalent to `--memory-mode memory`). |
| `--temp-dir` | `str` | `None` | Custom directory for temporary chunk streaming (default: `/tmp`). |
| `--log-file` | `str` | `None` | Path to output log file for timestamps and performance tracking. |
| `--min-baseq` | `int` | `20` | Minimum base quality for pileup. |
| `--min-mapq` | `int` | `0` | Minimum mapping quality for alignment filtering. |
| `--max-depth` | `int` | `250` | Maximum pileup depth. |
| `--no-baq` | `flag` | `False` | Disable Base Alignment Quality (BAQ) computation. |
| `--baq` | `taps`\|`taps-full`\|`taps-partial`\|`htslib` | `taps` | BAQ mode (Rust engine). The `taps` modes realign OT reads against the reference with every C masked to N and OB reads with every G masked, so methylation conversions (all contexts) are not penalised. `taps` skips reads that already match that reference exactly (processing 5.1M vs 2.2M reads/s for `taps-full`; 0.0007% of CG calls differ); `taps-full` runs every read; `taps-partial` only reads with an indel or clip (7.8M reads/s; 0.06% of CG calls differ). `htslib` is samtools/asTair BAQ on the plain reference and under-calls methylation. The Python engine never applies BAQ. See [BAQ for TAPS](#baq-for-taps). |
| `--caller` | `reads`\|`pileup` | `reads` | How calls are made (Rust engine). `reads` calls one read at a time, replaying htslib's pileup depth cap and mate-overlap rules; `pileup` uses htslib's column pileup. Identical calls; `reads` processes 36% faster. |
| `--baq-simd` | `auto`\|`avx512`\|`avx2`\|`off` | `auto` | BAQ kernel: 8 reads per vector (AVX-512), 4 (AVX2), or htslib one read at a time. All give identical results. |
| `--baq-device` | `cpu`\|`gpu` | `cpu` | Where the BAQ HMM runs. Identical results. See [GPU BAQ](#11-gpu-baq---baq-device-gpu). |
| `--gpu-backend` | `auto`\|`opencl`\|`cuda`\|`hip` | `auto` | GPU backend: auto tries CUDA, then HIP, then OpenCL. Intel GPUs use `opencl`. |
| `--gpu-device` | `str` | `None` | Use the first GPU whose name contains this text. |
| `--gpu-mem-mb` | `int` | `1024` | GPU memory for HMM matrices per launch. |
| `--context` | `both`\|`CG`\|`CH` | `both` | Contexts to extract. `CG` or `CH` skips the other context's pileup columns and leaves its HDF5 group out. See [`--context`](#10---context-cg-or-ch-only). |
| `--no-overlap-clip` | `flag` | `False` | Do not ignore overlapping mate read bases. |
| `--trim-r2` | `int` | `10` | Drop this many 5′ sequenced bases of read 2 before calling (Tn5/TAPS M-bias). |
| `--trim-r1` | `int` | `0` | Drop this many 5′ sequenced bases of read 1 before calling. |
| `-v, --verbose` | `flag` | `False` | Enable verbose debug logging. |

---

## Performance Flags Guide

Rates in this guide are mapped reads per second (M = million) on the reference workstation and test dataset described under [Performance](#performance): Rust engine, `-t 32 --shards 16 --memory-mode memory`, defaults otherwise, measured 2026-10-17. "Processing" is extraction (BAM read, BAQ, pileup, calling), "write" is the HDF5 write phase, "end to end" is the whole run. Single runs vary by about ±4%. They will shift on other machines; for a new machine, start with `taps-sc-extract --test`. The run summary prints processing time, HDF5 write time, estimated read and write MB/s, and peak memory of the extraction process and the writer processes. **Throughput (Reads)** is mapped records over wall time, including HDF5 write.

### 1. Compression (`--compression`)

Compression runs inside the writer processes, so it changes only the HDF5 write phase. One run per codec (`blosc`: median of five), measured 2026-10-18.

| Codec | End to end | HDF5 write | End to end vs `blosc` | Output size vs `none` | R / Bioconductor requirement |
| :--- | ---: | ---: | :---: | :---: | :--- |
| **`blosc`** (LZ4, level 9) *(Rust default)* | 2.45M reads/s | 4.73M reads/s | — | 96% smaller | `library(rhdf5filters)` |
| **`lzf`** | 2.14M reads/s | 4.20M reads/s | 13% slower | 88% smaller | `library(rhdf5filters)` |
| **`gzip`** (deflate level 1; Python default) | 2.13M reads/s | 3.61M reads/s | 13% slower | 94% smaller | **Plain `rhdf5`, no extra packages** |
| **`gzip-shuffle`** (byte shuffle + deflate 1) | 1.91M reads/s | 3.11M reads/s | 22% slower | 96% smaller | Plain `rhdf5` |
| **`gzip6`** (deflate level 6) | 1.64M reads/s | 2.45M reads/s | 33% slower | 95% smaller | Plain `rhdf5` |
| **`blosc-zstd`** | 1.54M reads/s | 2.31M reads/s | 37% slower | **97% smaller** (smallest) | `library(rhdf5filters)` |
| **`none`** | 1.01M reads/s | 1.29M reads/s | 59% slower | 0% (34 GB here) | Plain `rhdf5` |

- **`blosc` is the Rust default since 2026-10-16**: it writes faster than `lzf` with files 63% smaller (1.5 GB vs 4.0 GB here), and R needs the same `rhdf5filters` package for both. Python readers need `import hdf5plugin` before opening the file with `h5py` (`lzf` is built into h5py; Blosc is not).
- **`gzip-shuffle`** is the smallest output that plain `rhdf5` reads without plugins.
- **`blosc-zstd`** gives the smallest files at about half the write rate.
- **`none`** is slowest here because 34 GB has to reach the disk; the writers are disk-bound, not compression-bound.

### 2. `-t, --workers`
- **Mechanism**: Rayon threads, each processing genomic windows. Extraction is CPU-bound. With the default BAQ, worker time splits about 55% BAQ, 23% calling, 13% BAM read and inflate, 7% pileup bookkeeping (read-centric calling).
- **Recommendation**: Use all logical CPUs. On the 16-core / 32-thread reference machine, 32 workers process 5.08M reads/s against 2.89M with 16 (SMT); end to end 2.45M vs 1.73M reads/s. `0` = auto (CPUs, bounded by the memory budget), or the saved `--test` profile.

### 3. `--decomp-threads` (BAM BGZF decompression)
- **Mechanism**: htslib threads per worker that inflate BGZF blocks ahead of the reader.
- **Recommendation**: Leave it unset. Rust picks 0 at 8 or more workers (inflate runs on the worker thread), 2 at 5–7 and 4 at 1–4. With libdeflate, inflate measured 1–6% of worker time, so extra threads only add contention.

### 4. `--chunk-size-mb` (window size)
- **Rust default: auto.** About 6 windows per worker, clamped to 2–100 Mb (mm10: ~28 Mb at 16 workers, ~14 Mb at 32, ~4.7 Mb at 96). Enough windows keep every worker busy to the end; the old fixed 100 Mb gave 37 windows, so 32 workers had a long tail and 96 workers would sit mostly idle.
- Window size does not change calls (verified per cell at 100, 10 and 2 Mb).
- The Python engine defaults to `10` Mb.

### 5. `--shards`
- **Mechanism**: Writes $N$ shard files and `master.h5` with relative `ExternalLink`s. Rust hashes barcodes into 4N buckets and packs buckets into shards by site count, so a few deep cells do not land in one shard (the slowest shard sets the write time). Python assigns `md5(bc) % N`.
- **Recommendation**: Depends on the filesystem and the cell count. On local NVMe with `blosc`, 16 shards write 4.73M reads/s against 3.61M for one shard (end to end 2.45M vs 2.09M reads/s); with the larger `lzf` output one shard was much slower (2026-10-16: 1.8M vs 4.1M reads/s write). On a server filesystem with fewer cells, more shards were much slower, so measure there (`--test` probes the output filesystem). The default stays 1 shard below 10,000 cells.
- **Few, deep cells** (hundreds of cells, ~1M sites each): one cell is one writer's serial work, so beyond about cells/2 shards, more shards add nothing. The run summary's `Shard writers:` line shows the sites per shard (median, max), the slowest shard, and where writer time went (read, create, fill, HDF5 write).

### 6. `--max-writer-threads` and `--writer-mode`
- **`--writer-mode process`** (default with more than one shard): one `taps-sc-extract-rs write-shard` child per shard, up to `--max-writer-threads` at a time. Each child streams its shard one cell at a time in bounded slices, so a writer needs roughly 15–60 bytes per site of its shard.
- **`--compression-threads`** (Blosc codecs only, default 1): Blosc threads per writer. On the test BAM, 1, 2, 4 and 32 threads wrote equally fast (within 5%) at 16 shards and at 1 shard. Most cell datasets are a single chunk too small for Blosc to split, so parallelism comes from shards.
- **`--writer-mode thread`**: in-process writers. hdf5-metno takes one process-wide lock for every HDF5 call, so they do not write in parallel; auto caps them at 6.
- **Auto count**: up to min(shards, CPUs), and only as many as fit the memory budget at about 20 bytes per site of a shard each, or the saved `--test` profile (still memory-bounded). Lower it on network filesystems (NFS/Lustre: 2–4) or spinning disks (1).

### 7. `--memory-mode` / `--no-temp-file`
- **`stream`**: window results go to compact temp files under `--temp-dir`; each writer reads its own shard's files. The extraction process peaks at 3.7 GB against 11.0 GB in memory mode on the test BAM, at about the same rate (2.38M vs 2.45M reads/s end to end).
- **`memory`**: window results stay in RAM and are piped to the writer processes. No faster than `stream` on local NVMe (within run-to-run noise); use it when temp storage is slow.
- **`auto`**: the Rust engine estimates sites from the BAM index's mapped-read count (about 14 per read). It picks memory mode only if the extraction process plus the writer processes fit the budget without fewer writers than stream mode would allow; otherwise stream. The writer count is bounded by the same budget. The config block prints the plan (`Memory plan: ... est. peak ~X GB of Y GB`). Set the budget with `--max-memory-fraction` (default 0.6 of available memory) or `--max-memory-gb`.

### 8. `--temp-dir`
- Point it at fast local storage. Keep inputs, temp and output off slow network or 9p mounts: on the reference machine the same job ran at 0.31M reads/s from a 9p-mounted HDD against 0.56M reads/s from local NVMe (measured before the 2026-10-15 upgrades).

### 9. Profiling a run (`--profile-json`)
Writes per-phase wall and CPU time, utilisation, a per-stage split of worker time (read/inflate, BAQ, pileup, calling), per-window timings and writer memory. Use it to see what limits a run on a new system.

### 10. `--context` (CG or CH only)
- **Mechanism**: each reference position's context is decided before any per-read work; with `--context CG` CH positions are skipped there, so their calls are never made, stored or written, and no `/CH` group is created (in shards, single files or `master.h5`). The kept context's calls are identical to a `both` run, cell for cell. BAQ and the pileup itself still run on every read, since any read can cover a CpG. Both engines support it.
- **Measured** (full test dataset; CH is 96% of sites):

| Run | Processing | HDF5 write | End to end | Extraction memory | Writer memory |
|---|---:|---:|---:|---:|---:|
| `both` (default) | 5.08M reads/s | 4.73M reads/s | 2.45M reads/s | 11.0 GB | 10.7 GB |
| `--context CG` | 6.40M reads/s | 94M reads/s | **5.99M reads/s** | 2.1 GB | 1.4 GB |
| `--context CH` | 4.75M reads/s | 4.64M reads/s | 2.35M reads/s | 10.6 GB | 10.3 GB |
| `--no-baq` | 8.67M reads/s | 4.72M reads/s | 3.06M reads/s | 10.9 GB | 10.7 GB |
| `--context CG --no-baq` | 16.7M reads/s | 92M reads/s | **14.1M reads/s** | 1.9 GB | 1.4 GB |

  CG-only is 144% faster end to end than both contexts (362% with `--no-baq`), with 81% less extraction memory and 87% less writer memory. CH-only saves nothing beyond noise. (Measured 2026-10-18.)
- **Amethyst**: index only the extracted context (`indexChr(obj, "CG")`); list barcodes from that group (`ls_df$group == "/CG"`, or `"/CH"` for CH-only output).

### 11. GPU BAQ (`--baq-device gpu`)
- **What runs where**: the HMM at the core of BAQ (the dense double-precision part) runs on the GPU; everything TAPS-specific stays on the CPU: the C/G mask, the strand choice, skipping reads that already match, BAQ_EXTEND and the quality cap. Reads longer than 192 bp or with a very wide band (large indels) also stay on the CPU. One service thread owns the GPU and pools reads from all workers into large launches. Each worker keeps several batches in flight (4 by default, `TAPS_GPU_DEPTH`), so it runs pileup and calling on one batch while the GPU computes BAQ for the next; batches still enter the pileup in file order.
- **Exact**: the GPU kernel is htslib's arithmetic operation for operation (no fused multiply-add, IEEE division, sums in order), with the HMM parameters computed on the CPU. A base whose quality sits within 1e-9 of a rounding point is recomputed on the CPU with glibc. `taps-sc-extract-rs baq-selftest --gpu -b BAM -f FA` checks every read: all 1,431,323 chr19 reads identical through CUDA and through OpenCL (99.98% of them on the GPU); chr19 extraction reproduces the frozen references in all BAQ modes.
- **When it helps**: on GPUs with enough double precision relative to the host CPU. On an Intel Arc Pro B70 (Battlemage host, Ryzen 7 1700X, 16 threads) it is faster than the CPU; on the reference workstation's RTX 3060 (weak FP64, PCIe 3.0 x4) it is slower, so the default stays `cpu`. Full test dataset:

| Host | BAQ | Device | Processing | End to end |
|---|---|---|---:|---:|
| Battlemage | `taps` (default) | CPU (Ryzen 7 1700X) | 1.34M reads/s (2026-10-18: 1.40M) | 0.92M reads/s (0.95M) |
| Battlemage | `taps` | **Arc Pro B70, OpenCL** (2026-10-18: overlapped batches, faster CPU-side preparation) | **2.21M reads/s (+58% vs 1.40M)** | **1.27M reads/s (+34% vs 0.95M)** |
| Battlemage | `taps-full` | CPU | 0.65M reads/s | 0.54M reads/s |
| Battlemage | `taps-full` | **Arc Pro B70, OpenCL** | **1.23M reads/s (+89%)** | **0.87M reads/s (+61%)** |

On the reference workstation (32 threads, 2026-10-18):

| BAQ | Device | Processing | End to end |
|---|---|---:|---:|
| `taps` (default) | CPU, AVX2 | 5.08M reads/s | 2.45M reads/s |
| `taps` | RTX 3060, CUDA | 2.65M reads/s | 1.74M reads/s |
| `taps-full` | CPU, AVX2 | 2.19M reads/s | 1.48M reads/s |
| `taps-full` | RTX 3060, CUDA | 0.76M reads/s | 0.66M reads/s |

The 3060 is 37% faster than its first GPU build (1.94M reads/s processing) but stays well behind this CPU. OpenCL through PoCL on the 3060 was exact earlier (1.86M reads/s on the first build) and was not re-measured.

- **Check a host**: `taps-sc-extract --test` times BAQ on the GPU against the CPU and says which is faster (reported only; `--baq-device` still defaults to `cpu`).
- **Intel Arc Pro B70**: works through Intel's compute runtime (`--gpu-backend opencl` or `auto`); `clinfo` lists double precision (`cl_khr_fp64`). Measured double-precision rate 0.70 TFLOP/s (8x the RTX 3060). The kernel runs at about 2.1M reads/s on real reads (the `--test` probe reports it), using Intel's 256-register mode automatically. With the GPU, BAQ stops being the limit on that host: `--no-baq` processes 2.71M reads/s. Workers overlap their GPU batches with calling, and the GPU is idle part of the time. What remains is CPU work around the HMM: preparing reads, BAQ_EXTEND and the quality cap, all of which stay on the CPU because they are TAPS-specific and cheap per read.

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

### BAQ for TAPS

BAQ (htslib `sam_prob_realn`) realigns each read to the reference with a profile HMM and lowers base qualities where the alignment is uncertain; bases that fall below `--min-baseq` are not called. On the plain reference, a methylated C read as T (A on OB) is a mismatch, so BAQ lowers quality at and around methylated sites and drops methylated calls preferentially.

The `taps` modes give the HMM a reference with the convertible base masked to `N`: every `C` for OT reads and every `G` for OB reads, in all contexts (CpG and CH). `N` scores as neither match nor mismatch in the HMM, so a converted or unconverted base costs nothing, while indels and other mismatches still lower quality. Only BAQ sees the masked reference; calling uses the real one.

chr19, change in calls against `--no-baq`:

| BAQ | CG methylated | CG unmethylated | CG % | CH methylated | CH unmethylated | CH % |
|---|---:|---:|---:|---:|---:|---:|
| off (`--no-baq`) | 574,228 | 281,919 | 67.07 | 59,592 | 19,984,017 | 0.297 |
| `--baq htslib` | −3.00% | −0.37% | 66.48 | −4.93% | −0.42% | 0.284 |
| `--baq taps-full` | −0.07% | −0.84% | 67.24 | −0.60% | −0.38% | 0.297 |

BAQ is the largest single cost in processing. About 72% of reads have an all-M CIGAR and match the masked reference at every A/C/G/T base; for those the HMM mostly lowers end bases that the R2 trim or `--min-baseq` drop anyway. The default `--baq taps` skips them. Full test dataset, AVX2 kernel:

| BAQ | Reads through the HMM | Processing | End to end | CpG % | CG cell-sites differing from `taps-full` |
|---|---:|---:|---:|---:|---:|
| `taps-full` (HMM on every read) | 99% | 2.19M reads/s | 1.48M reads/s | 71.615 | — |
| `taps` (default) | 28% | 5.08M reads/s | 2.45M reads/s | 71.615 | 172 of 25.9M (0.0007%) |
| `taps-partial` (indel/clip reads only) | 8.5% | 7.80M reads/s | 2.90M reads/s | 71.585 | 15,671 (0.06%) |
| `--no-baq` | 0% | 8.67M reads/s | 3.06M reads/s | 71.543 | 64,056 (0.25%) |
| `htslib` (plain reference) | 99% | 2.32M reads/s | 1.54M reads/s | 70.859 | (biased; see above) |

Measured 2026-10-18 (`taps-full`: mean of two runs, 2.02–2.35M reads/s processing).

The HMM runs on SIMD kernels (`--baq-simd`, default `auto`): 8 reads per vector with AVX-512, 4 with AVX2. Results are identical to htslib's scalar code for every read (whole test dataset, 53.7M reads, checked with `baq-selftest`). The TAPS mask is applied inside the kernel while it reads the reference, so no masked copy of any chromosome is stored. AVX2 against htslib (`--baq-simd off`) at the default BAQ: processing goes from 3.15M to 5.08M reads/s (61% faster), end to end from 1.86M to 2.45M reads/s (32% faster).

---

## Performance

### Reference workstation (2026-10-18)

| | |
|---|---|
| CPU | AMD Ryzen 9 5900XT: 16 cores / 32 threads, 3.3 GHz base, up to 4.8 GHz boost (AMD spec), AVX2 (no AVX-512) |
| Memory | ~58 GB visible to WSL2 |
| Storage | WSL2 ext4 VHDX on NVMe; input, temp and output all on it |
| Dataset | `5base_ucbTn5_sp2.srt.bam`: 53.5M mapped reads, 7,351 cells, mm10, 3.7 GB BAM |
| Build | `./install.sh` (native); BAQ kernel AVX2 |
| Settings | `-t 32 --shards 16 --memory-mode memory`, defaults otherwise: `--compression blosc`, `--baq taps`, `--caller reads`, `--context both`, auto window size (14.2 Mb, 203 windows), 16 writer processes. Measured 2026-10-18. |

Best case on this machine with the defaults, median of five runs (range 2.39–2.56M reads/s end to end):

| Phase | Rate |
|---|---:|
| Processing (BAM read, BAQ, calling) | 5.08M reads/s |
| HDF5 write (1.5 GB, 16 shards) | 4.73M reads/s |
| **End to end** | **2.45M reads/s (147M reads/min)** |
| Peak memory | 11.0 GB extraction + 10.7 GB writers (stream mode: 3.7 GB + 10.7 GB, 2.38M reads/s) |

Processing and writing now run at about the same rate, so end to end is bounded by both: neither alone can double it.

Faster options on the same machine:

| Options | End to end | Output |
|---|---:|---|
| `--baq taps-partial` | 2.90M reads/s | same contexts; 0.06% of CG calls differ |
| `--no-baq` | 3.06M reads/s | no BAQ |
| `--context CG` | 5.99M reads/s | CpG only, 0.15 GB |
| `--context CG --no-baq` | 14.1M reads/s | CpG only, no BAQ |

How speed changed across releases: [Performance History](#performance-history).

### Servers

*To be filled in: 64-thread server with an Intel Arc Pro B70 (`install.sh`, `clinfo`, `--test` including the GPU BAQ probe, `baq-selftest --gpu`, T=32/64, `--baq-device cpu` vs `gpu`, shard counts on its filesystem).*

---

## HDF5 Output Structure & Schema

The output HDF5 matches the Amethyst 1.0+ / Facet specification. With `--context CG` or `--context CH`, only that group is present.

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

The chr19 tests check the Rust engine bit for bit against a frozen reference (`test_chr19_rust_tapsskip_2026-10-16/`, default BAQ; `--baq taps-full` against `test_chr19_rust_tapsbaq_2026-10-16/`, `--baq htslib` against `test_chr19_rust_overlapfix_2026-10-15/`), and require Rust `--no-baq` to equal the Python engine at every cell and position.

---

## Performance History

End-to-end and per-phase rates on the reference workstation and test dataset under [Performance](#performance) (53.5M mapped reads), `-t 32 --shards 16`, each release's default options. "% faster" is the rate increase over the first measured build (2026-10-14, `026c52e`).

| Release | End to end | Processing | HDF5 write | Peak memory |
|---|---:|---:|---:|---:|
| First measured build (2026-10-14): htslib BAQ, 100 Mb windows, in-process writers, LZF | 0.56M reads/s (34M/min) | 0.84M reads/s | 1.70M reads/s | 29.7 GB |
| 2026-10-15: balanced windows, libdeflate, shared BAQ references, one writer process per shard, native build, per-cell sorted vectors | 0.97M reads/s (72% faster) | 1.30M (54% faster) | 3.78M (123% faster) | 21.7 GB (27% less) |
| 2026-10-16: TAPS-aware BAQ (correctness fix; same cost) | 1.00M reads/s (77% faster) | 1.31M (56% faster) | 4.18M (146% faster) | 23.2 GB (22% less) |
| 2026-10-16: skip BAQ for reads that already match the masked reference | 1.64M reads/s (192% faster) | 2.60M (208% faster) | 4.49M (165% faster) | 23.2 GB (22% less) |
| 2026-10-16: SIMD BAQ (AVX2) | 1.76M reads/s (212% faster) | 3.24M (285% faster) | 4.14M (144% faster) | 24.0 GB (19% less) |
| 2026-10-16: `blosc` default compression (output 63% smaller) | 1.95M reads/s (246% faster) | 3.41M (304% faster) | 4.55M (168% faster) | 23.9 GB (20% less) |
| 2026-10-17: TAPS mask inside the BAQ kernel | 1.85M reads/s (229% faster) | 3.29M (291% faster) | 4.29M (153% faster) | 21.2 GB (29% less); stream mode 13.5 GB (55% less) |
| 2026-10-18: branch-free BAQ skip check, table-driven BAQ preparation | 2.07M reads/s (270% faster) | 3.84M (357% faster) | 4.49M (164% faster) | 21.7 GB (27% less) |
| 2026-10-18: read-centric calling (`--caller reads`; same calls as the column pileup) (current) | **2.45M reads/s (338% faster)** | **5.08M (505% faster)** | **4.73M (178% faster)** | **21.7 GB (27% less)** |

- The two 2026-10-18 rows are medians of three and five runs (end to end 2.03–2.18M and 2.39–2.56M reads/s). Rows above the 2026-10-17 one other than the end-to-end column are unchanged.
- The 2026-10-16 `blosc` and 2026-10-17 rows run at the same rate within run-to-run spread (today's five runs: 1.78–1.93M reads/s; the 2026-10-16 runs: 1.94–1.96M). The 2026-10-17 change removed memory, not time.
- Peak memory counts the extraction process plus the writer processes.
- The first build used htslib BAQ (CpG 70.86%); later builds use TAPS-aware BAQ (CpG 71.62%). Running BAQ on every read (`--baq taps-full`), the current build processes 2.19M reads/s, 161% faster than the first build.
- Not in the table: `--context CG`, 5.99M reads/s end to end (144% faster than the current default).
- Detailed snapshots: [`code/nvme-profiling-2026-10-14.md`](code/nvme-profiling-2026-10-14.md), [`2026-10-15`](code/nvme-profiling-2026-10-15.md), [`2026-10-16`](code/nvme-profiling-2026-10-16.md), [`2026-10-17`](code/nvme-profiling-2026-10-17.md), [`2026-10-18`](code/nvme-profiling-2026-10-18.md); change log: [`code/rust-implementation-log.md`](code/rust-implementation-log.md).

---

## License
 
This project is licensed under the **Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)** License. See the [LICENSE](LICENSE) file for details.
 
Free for non-commercial academic, scientific, and educational research use. Commercial use requires explicit permission.
