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
- **HDF5 compression options**: `lzf` (default), `blosc` (LZ4 level 9: as fast as `lzf`, 2.7× smaller files), `gzip` (readable by plain `rhdf5`), `gzip-shuffle`, `gzip6`, `blosc-zstd` (smallest) and uncompressed `none`. See [Compression](#1-compression---compression).
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

BAQ runs on SIMD kernels that process 8 reads per vector (AVX-512F) or 4 (AVX2). Both are always compiled, including in `--portable` builds, and the binary picks the widest the CPU runs. After building, `install.sh` runs `taps-sc-extract-rs baq-selftest`, which runs every supported kernel and htslib's scalar BAQ on the same reads and requires identical base qualities for every read; the install stops on any difference. Override with `--baq-simd avx512|avx2|off`. `taps-sc-extract --test --save-profile` times each kernel under full load and saves the fastest, since some Intel CPUs lower their clock for AVX-512.

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
                       [--no-baq] [--baq {taps,taps-full,taps-partial,htslib}] [--no-overlap-clip]
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
| `--compression` | `str` | `gzip` Python / `lzf` Rust | Compression: `gzip` (level 1 deflate), `gzip-shuffle`, `gzip6`, `lzf`, `blosc`, `blosc-zstd`, `none`. |
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
| `--baq` | `taps`\|`taps-full`\|`taps-partial`\|`htslib` | `taps` | BAQ mode (Rust engine). The `taps` modes realign OT reads against the reference with every C masked to N and OB reads with every G masked, so methylation conversions (all contexts) are not penalised. `taps` skips reads that already match that reference exactly (~2× faster than `taps-full`, 0.0007% of CG calls differ); `taps-full` runs every read; `taps-partial` only reads with an indel or clip (~4×, 0.06%). `htslib` is samtools/asTair BAQ on the plain reference and under-calls methylation. The Python engine never applies BAQ. See [BAQ for TAPS](#baq-for-taps). |
| `--no-overlap-clip` | `flag` | `False` | Do not ignore overlapping mate read bases. |
| `--trim-r2` | `int` | `10` | Drop this many 5′ sequenced bases of read 2 before calling (Tn5/TAPS M-bias). |
| `--trim-r1` | `int` | `0` | Drop this many 5′ sequenced bases of read 1 before calling. |
| `-v, --verbose` | `flag` | `False` | Enable verbose debug logging. |

---

## Performance Flags Guide

Comparisons in this guide are relative (fold change or percent) and come from the reference workstation and test dataset described under [Performance](#performance), with the Rust engine at `-t 32 --shards 16`, LZF, default BAQ, unless noted. They will shift on other machines; for a new machine, start with `taps-sc-extract --test`. The run summary prints processing time, HDF5 write time, estimated read and write MB/s, and peak memory of the extraction process and the writer processes. **Throughput (Reads)** is mapped records over wall time, including HDF5 write.

### 1. Compression (`--compression`)

Rust engine, full test dataset, 16 shards, one run per codec (2026-10-16). Compression runs inside the writer processes, so it changes only the HDF5 write phase; processing is the same for every codec. Write-time differences under about 10% are within run-to-run noise.

| Codec | HDF5 write time vs `lzf` | Total wall vs `lzf` | Output size vs `none` | R / Bioconductor requirement |
| :--- | :---: | :---: | :---: | :--- |
| **`lzf`** *(Rust default)* | 1.00× | 1.00× | 8.6× smaller | `library(rhdf5filters)` |
| **`blosc`** (LZ4, level 9) | 0.93× | 0.97× | 23× smaller | `library(rhdf5filters)` |
| **`gzip`** (deflate level 1; Python default) | 1.13× | 1.06× | 17× smaller | **Plain `rhdf5`, no extra packages** |
| **`gzip-shuffle`** (byte shuffle + deflate 1) | 1.40× | 1.16× | 28× smaller | Plain `rhdf5` |
| **`gzip6`** (deflate level 6) | 1.73× | 1.28× | 18× smaller | Plain `rhdf5` |
| **`blosc-zstd`** | 1.92× | 1.41× | **32× smaller** (smallest) | `library(rhdf5filters)` |
| **`none`** | 3.6× | 2.2× | 1× (34 GB here) | Plain `rhdf5` |

- **`blosc` writes as fast as `lzf` and its files are 2.7× smaller**, with the same R requirement. It is the better choice when you have `rhdf5filters`.
- **`gzip-shuffle`** is the smallest output that plain `rhdf5` reads without plugins.
- **`none`** is slowest here because 34 GB has to reach the disk; the writers are disk-bound, not compression-bound.

### 2. `-t, --workers`
- **Mechanism**: Rayon threads, each processing genomic windows. Extraction is CPU-bound. With the default BAQ, worker time splits about 38% BAQ, 36% calling, 14% BAM read and inflate, 13% pileup.
- **Recommendation**: Use all logical CPUs. On the 16-core / 32-thread reference machine, 32 workers process 1.4× faster than 16 (SMT). `0` = auto (CPUs, bounded by the memory budget), or the saved `--test` profile.

### 3. `--decomp-threads` (BAM BGZF decompression)
- **Mechanism**: htslib threads per worker that inflate BGZF blocks ahead of the reader.
- **Recommendation**: Leave it unset. Rust picks 0 at 8 or more workers (inflate runs on the worker thread), 2 at 5–7 and 4 at 1–4. With libdeflate, inflate measured 1–6% of worker time, so extra threads only add contention.

### 4. `--chunk-size-mb` (window size)
- **Rust default: auto.** About 6 windows per worker, clamped to 2–100 Mb (mm10: ~28 Mb at 16 workers, ~14 Mb at 32, ~4.7 Mb at 96). Enough windows keep every worker busy to the end; the old fixed 100 Mb gave 37 windows, so 32 workers had a long tail and 96 workers would sit mostly idle.
- Window size does not change calls (verified per cell at 100, 10 and 2 Mb).
- The Python engine defaults to `10` Mb.

### 5. `--shards`
- **Mechanism**: Writes $N$ shard files and `master.h5` with relative `ExternalLink`s. Rust hashes barcodes into 4N buckets and packs buckets into shards by site count, so a few deep cells do not land in one shard (the slowest shard sets the write time). Python assigns `md5(bc) % N`.
- **Recommendation**: Depends on the filesystem and the cell count. On local NVMe, 16 shards write 2.2× faster than one shard (1.5× faster end to end). On a server filesystem with fewer cells, more shards were much slower, so measure there (`--test` probes the output filesystem). The default stays 1 shard below 10,000 cells.
- **Few, deep cells** (hundreds of cells, ~1M sites each): one cell is one writer's serial work, so beyond about cells/2 shards, more shards add nothing. The run summary's `Shard writers:` line shows the sites per shard (median, max), the slowest shard, and where writer time went (read, create, fill, HDF5 write).

### 6. `--max-writer-threads` and `--writer-mode`
- **`--writer-mode process`** (default with more than one shard): one `taps-sc-extract-rs write-shard` child per shard, up to `--max-writer-threads` at a time. Each child streams its shard one cell at a time in bounded slices, so a writer needs roughly 15–60 bytes per site of its shard.
- **`--compression-threads`** (Blosc codecs only, default 1): Blosc threads per writer. On the test BAM, 1, 2, 4 and 32 threads wrote equally fast (within 5%) at 16 shards and at 1 shard. Most cell datasets are a single chunk too small for Blosc to split, so parallelism comes from shards.
- **`--writer-mode thread`**: in-process writers. hdf5-metno takes one process-wide lock for every HDF5 call, so they do not write in parallel; auto caps them at 6.
- **Auto count**: up to min(shards, CPUs), and only as many as fit the memory budget at about 20 bytes per site of a shard each, or the saved `--test` profile (still memory-bounded). Lower it on network filesystems (NFS/Lustre: 2–4) or spinning disks (1).

### 7. `--memory-mode` / `--no-temp-file`
- **`stream`**: window results go to compact temp files under `--temp-dir`; each writer reads its own shard's files. The extraction process peaks at about 55% of memory mode's (7.3 GB vs 13.3 GB on the test BAM).
- **`memory`**: window results stay in RAM and are piped to the writer processes. No faster than `stream` on local NVMe (within run-to-run noise); use it when temp storage is slow.
- **`auto`**: the Rust engine estimates sites from the BAM index's mapped-read count (about 14 per read). It picks memory mode only if the extraction process plus the writer processes fit the budget without fewer writers than stream mode would allow; otherwise stream. The writer count is bounded by the same budget. The config block prints the plan (`Memory plan: ... est. peak ~X GB of Y GB`). Set the budget with `--max-memory-fraction` (default 0.6 of available memory) or `--max-memory-gb`.

### 8. `--temp-dir`
- Point it at fast local storage. Keep inputs, temp and output off slow network or 9p mounts: on the reference machine the same job ran 1.8× slower from a 9p-mounted HDD than from local NVMe (measured before the 2026-10-15 upgrades).

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

### BAQ for TAPS

BAQ (htslib `sam_prob_realn`) realigns each read to the reference with a profile HMM and lowers base qualities where the alignment is uncertain; bases that fall below `--min-baseq` are not called. On the plain reference, a methylated C read as T (A on OB) is a mismatch, so BAQ lowers quality at and around methylated sites and drops methylated calls preferentially.

The `taps` modes give the HMM a reference with the convertible base masked to `N`: every `C` for OT reads and every `G` for OB reads, in all contexts (CpG and CH). `N` scores as neither match nor mismatch in the HMM, so a converted or unconverted base costs nothing, while indels and other mismatches still lower quality. Only BAQ sees the masked reference; calling uses the real one.

chr19, change in calls against `--no-baq`:

| BAQ | CG methylated | CG unmethylated | CG % | CH methylated | CH unmethylated | CH % |
|---|---:|---:|---:|---:|---:|---:|
| off (`--no-baq`) | 574,228 | 281,919 | 67.07 | 59,592 | 19,984,017 | 0.297 |
| `--baq htslib` | −3.00% | −0.37% | 66.48 | −4.93% | −0.42% | 0.284 |
| `--baq taps-full` | −0.07% | −0.84% | 67.24 | −0.60% | −0.38% | 0.297 |

BAQ is the largest single cost in processing. About 72% of reads have an all-M CIGAR and match the masked reference at every A/C/G/T base; for those the HMM mostly lowers end bases that the R2 trim or `--min-baseq` drop anyway. The default `--baq taps` skips them. Full test dataset, AVX2 kernel, speed relative to the default:

| BAQ | Reads through the HMM | Processing speed | Total wall speed | CpG % | CG cell-sites differing from `taps-full` |
|---|---:|---:|---:|---:|---:|
| `taps-full` (HMM on every read) | 99% | 0.59× | 0.71× | 71.615 | — |
| `taps` (default) | 28% | 1.00× | 1.00× | 71.615 | 172 of 25.9M (0.0007%) |
| `taps-partial` (indel/clip reads only) | 8.5% | 1.5× | 1.26× | 71.585 | 15,671 (0.06%) |
| `--no-baq` | 0% | 2.0× | 1.36× | 71.543 | 64,056 (0.25%) |

`htslib` costs the same as `taps-full`.

The HMM runs on SIMD kernels (`--baq-simd`, default `auto`): 8 reads per vector with AVX-512, 4 with AVX2. Results are identical to htslib's scalar code for every read (whole test dataset, 53.7M reads, checked with `baq-selftest`). AVX2 against htslib (`--baq-simd off`) on the reference machine:

| BAQ | BAQ CPU time | Processing speed |
|---|---:|---:|
| `taps-full` | 1.9× less | 1.5–1.6× faster |
| `taps` (default) | 1.9× less | 1.2× faster |

---

## Performance

### Reference workstation (2026-10-16)

| | |
|---|---|
| CPU | AMD Ryzen 9 5900XT: 16 cores / 32 threads, 3.3 GHz base, up to 4.8 GHz boost (AMD spec), AVX2 (no AVX-512) |
| Memory | ~58 GB visible to WSL2 |
| Storage | WSL2 ext4 VHDX on NVMe; input, temp and output all on it |
| Dataset | `5base_ucbTn5_sp2.srt.bam`: 53.5M mapped reads, 7,351 cells, mm10, 3.7 GB BAM |
| Build | `./install.sh` (native), commit `a733a66`; BAQ kernel AVX2 |
| Settings | `-t 32 --shards 16 --memory-mode memory --compression lzf`, default BAQ (`--baq taps`), auto window size (14.2 Mb, 203 windows), 16 writer processes |

Best case on this machine, three runs:

| Phase | Time |
|---|---:|
| Processing (BAM read, BAQ, pileup, calling) | 16.1–17.5 s |
| HDF5 write (4.0 GB, 16 shards) | 12.2–14.0 s |
| **Total wall** | **28.3–30.5 s** (~1.8M reads/s) |
| Peak memory | 13.3 GB extraction + 10.7 GB writers (stream mode: 7.3 GB + 10.7 GB, same speed) |

Fastest options on the same machine: `--baq taps-partial` 23.6 s, `--no-baq` 21.8 s, `--compression blosc` 28.8 s with 2.7× smaller output.

### Improvement history (same machine and dataset)

| Change | Processing | HDF5 write | Total wall | Peak memory |
|---|---:|---:|---:|---:|
| Baseline, 2026-10-14 (`026c52e`; htslib BAQ, 100 Mb windows, in-process writers) | 1.0× | 1.0× | 1.0× | 29.7 GB |
| Balanced windows, libdeflate, shared BAQ references, one writer process per shard, native build, per-cell sorted vectors (2026-10-15) | 1.5× | 2.2× | 1.7× | 11.5 + 10.2 GB |
| TAPS-aware BAQ (correctness; same cost) | 1.6× | 2.5× | 1.8× | 12.5 + 10.7 GB |
| Skip BAQ for reads that already match (new default) | 3.0× | 2.4× | 2.8× | 12.5 + 10.7 GB |
| SIMD BAQ, AVX2 (current) | **3.8×** | **2.4×** | **3.2×** | 13.3 + 10.7 GB |

Baseline output used htslib BAQ (CpG 70.86%); current output uses TAPS-aware BAQ (CpG 71.62%). Running the HMM on every read (`--baq taps-full`) is 2.2× faster than the baseline's processing. Detailed snapshots: [`code/nvme-profiling-2026-10-14.md`](code/nvme-profiling-2026-10-14.md), [`code/nvme-profiling-2026-10-15.md`](code/nvme-profiling-2026-10-15.md), [`code/rust-implementation-log.md`](code/rust-implementation-log.md).

### Servers

*To be filled in: 64-thread server (`install.sh`, `--test`, T=32/64, shard counts on its filesystem, AVX-512 if available).*

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

The chr19 tests check the Rust engine bit for bit against a frozen reference (`test_chr19_rust_tapsskip_2026-10-16/`, default BAQ; `--baq taps-full` against `test_chr19_rust_tapsbaq_2026-10-16/`, `--baq htslib` against `test_chr19_rust_overlapfix_2026-10-15/`), and require Rust `--no-baq` to equal the Python engine at every cell and position.

---

## License
 
This project is licensed under the **Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)** License. See the [LICENSE](LICENSE) file for details.
 
Free for non-commercial academic, scientific, and educational research use. Commercial use requires explicit permission.
