# Troubleshooting and Error Guide

Solutions to common errors, validation failures, and runtime bottlenecks. [Back to the documentation index](index.md). Options are described in [usage.md](usage.md).

---

## Chemistry Detection and Strand Inconsistencies

### Detection Rejections

#### `cannot tell the chemistry from the BAM (... single-end reads sampled, no BSBolt YS tags); pass --chemistry taps or --chemistry bisulfite`
- **Cause**: The BAM contains single-end reads without BSBolt `YS` strand tags in the first 10,000 mapped reads, preventing automatic classification.
- **Fix**: Specify the chemistry explicitly using `--chemistry taps` (if single-end TAPS) or `--chemistry bisulfite` (if bisulfite without BSBolt tags, typically paired with `--strand-source flag`).

#### `--chemistry taps, but the BAM carries BSBolt bisulfite tags (YS); use --chemistry bisulfite`
- **Cause**: `--chemistry taps` was specified on a BAM containing BSBolt bisulfite tags.
- **Fix**: Switch to `--chemistry bisulfite` or remove the explicit flag to let `--chemistry auto` select bisulfite.

#### `--chemistry taps needs paired-end reads, but all ... sampled reads are single-end`
- **Cause**: Directional TAPS calling relies on paired-end orientation flags (99/147/83/163).
- **Fix**: Ensure the input alignment contains paired-end reads.

#### `--chemistry bisulfite supports single-end reads only, but the BAM has paired reads`
- **Cause**: Bisulfite extraction currently supports single-end directional protocols (e.g., standard BSBolt pipelines).
- **Fix**: Paired-end bisulfite is not supported yet. Align the reads single-end, or extract them another way for now.

#### `--chemistry bisulfite, but no read carries a BSBolt YS strand tag; pass --strand-source flag to take the strand from the alignment flag`
- **Cause**: Bisulfite was requested, but none of the first 10,000 reads has BSBolt's `YS` strand tag (another aligner, or tags stripped).
- **Fix**: If the library is directional and the alignment flag gives the strand (flag 0 = original top, 16 = original bottom), pass `--strand-source flag`.

#### `<bam1> is <chem1> but <bam2> is <chem2>: one run takes one chemistry`
- **Cause**: In a multi-BAM run, different input files resolved to conflicting chemistries (e.g., one TAPS and one bisulfite).
- **Fix**: Group inputs by chemistry and process each library type in separate runs.

### Strand Guard and Tag Failures

#### `read <qname> has YS:Z:<ys>: non-directional or PBAT bisulfite libraries are not supported`
- **Cause**: A read carries a BSBolt strand tag containing `G2A` (`W_G2A` or `C_G2A`), indicating a non-directional or PBAT bisulfite library.
- **Fix**: Only directional bisulfite libraries (where unmethylated C converts to T on OT/OB strands) are supported. Re-align with directional settings if appropriate.

#### `read <qname> has YS:Z:<ys> but flag <flag> (<strand>): the strand tag and the alignment disagree (use --strand-source flag to trust the flag)`
- **Cause**: The alignment reverse-strand flag (0x10) contradicts the BSBolt `YS` tag (`W_C2T` must be forward; `C_C2T` must be reverse).
- **Fix**: If an upstream processing step modified alignment flags or if your aligner writes non-standard tags, bypass tag verification by adding `--strand-source flag`.

#### `read <qname> has no BSBolt YS strand tag (use --strand-source flag to take the strand from the flag)`
- **Cause**: A bisulfite run with `--strand-source auto` encountered a read lacking a `YS` tag.
- **Fix**: Pass `--strand-source flag` to derive strand orientation directly from SAM bitwise flags.

#### `--baq htslib treats every unmethylated C as a mismatch in bisulfite data; use a taps BAQ mode or --no-baq`
- **Cause**: `--baq htslib` was requested on bisulfite data. Standard htslib BAQ evaluates against the unmodified reference, which severely penalizes unmethylated C→T conversions.
- **Fix**: Use `--no-baq` (default for bisulfite) or a TAPS BAQ mode (`--baq taps`).

---

## Reference Genome and Alignment Validation

#### `the FASTA looks bisulfite-converted on <first> (C ...%, G ...% of ... bases): pass the unconverted genome the aligner index was built from`
- **Cause**: More than 10,000 bases sampled from the middle of the first contig showed <5% C or <5% G bases, indicating a bisulfite-converted 3-letter FASTA index.
- **Fix**: Pass the native, unconverted reference genome FASTA used to construct the aligner index. Calling logic requires native genomic reference sequence to identify cytosine contexts.

#### `contig <c> is in the BAM but not in the FASTA`
- **Cause**: The BAM `@SQ` header references a chromosome or contig not present in the FASTA index (`.fai`).
- **Fix**: Ensure the reference FASTA matches the exact reference used during alignment. If processing a subset of chromosomes, specify them via `-c/--chroms` (e.g., `-c chr1,chr2`).

#### `contig <c>: FASTA length <l> differs from the BAM header (<bam_len>)`
- **Cause**: Chromosome names match, but the sequence length in `ref.fa.fai` differs from the BAM `@SQ LN` attribute.
- **Fix**: Alignments and reference FASTA must share identical builds (e.g., do not mix UCSC `chr1` with Ensembl `1`, or hg19 with GRCh38).

#### `<bam> is not coordinate-sorted. Sort it first: samtools sort -o sorted.bam <bam>`
- **Cause**: A SAM input (`.sam` or `.sam.gz`) is not coordinate-sorted. SAM is converted to a temporary BAM, which needs sorted input. A BAM is not re-checked here; htslib refuses to index an unsorted BAM.
- **Fix**: Sort the alignment file by coordinate before running:
  ```bash
  samtools sort -@ 8 -o sample.srt.bam sample.bam
  ```

#### Missing Index File (`.bai`, `.csi`, `.crai`)
- **Behavior**: If missing, an index is automatically generated alongside the alignment file.
- **Read-only directory**: the index is built in a temporary directory instead, and the log says `(directory not writable; index written to ...)`. If that fails too, the run stops with `build index for ... (beside it and in ...)`; index the file yourself (`samtools index sample.bam`).

---

## Multi-BAM and Sample Sheet Rejections

#### `-b @FILE takes one list or sheet; do not mix it with other -b paths`
- **Cause**: An `@file` argument was passed alongside other `-b` paths.
- **Fix**: Pass either `@list.txt` / `@sheet.tsv` by itself, or list multiple space-separated BAM paths directly.

#### `<file>: Excel sheets are read by the taps-sc-extract wrapper; save as TSV for taps-sc-extract-rs`
- **Cause**: A `.xlsx` file was passed directly to the Rust binary `taps-sc-extract-rs`.
- **Fix**: Run via the Python CLI wrapper `taps-sc-extract` (which converts `.xlsx` to TSV automatically) or export the sheet to TSV manually.

#### `<file>: unknown column <col> (columns: bam, complexity, min_complexity, max_complexity, whitelist)`
- **Cause**: A sample sheet header contained an unsupported column name.
- **Fix**: Ensure sample sheet column headers match allowed fields exactly: `bam`, `complexity`, `min_complexity`, `max_complexity`, `whitelist`.

#### `barcode <bc> is in the complexity files of both <bam1> and <bam2>; a barcode must belong to one input`
- **Cause**: A single-cell barcode appears in the complexity files of two different BAM files.
- **Fix**: Barcodes must be unique to a single input library. Demultiplex or rename barcodes across libraries (e.g., prefixing library IDs).

#### `some inputs have a cell filter (complexity/whitelist) and some do not; give every input one, or none`
- **Cause**: In a sample sheet, some rows specified complexity or whitelist files while others left them blank without a global fallback.
- **Fix**: Supply filter files for all rows, specify global `--complexity` / `-w` flags, or omit filters for all inputs.

#### `<bam1> and <bam2> have different @SQ contigs or lengths; inputs must share one reference`
- **Cause**: Input BAM files were aligned against different reference genomes or headers.
- **Fix**: Re-align all libraries against the same reference genome version and contig ordering.

#### `only <hit> of <n> sampled reads in <bam> have a barcode in its complexity file; is the sheet row pointing at the right file?`
- **Cause**: Fewer than 50% of the first 5,000 mapped reads in an input BAM carry barcodes found in its designated complexity file.
- **Fix**: Verify that each row in the sample sheet correctly aligns each BAM file with its corresponding complexity file.

#### `a finite --max-depth with several inputs is not supported by the Python engine (it runs one pileup per file); use the Rust engine or --max-depth 0`
- **Cause**: The Python reference engine was executed with multiple BAM inputs and a finite `--max-depth`.
- **Fix**: Use the default unlimited depth (`--max-depth 0`) or run with `--engine rust`.

---

## Site Filtering and Depth Limit Errors

#### `--site-depth-limit auto needs --complexity or --whitelist to count the cells; give a number instead`
- **Cause**: `--site-depth-limit auto` was requested, but neither `--complexity` nor `--whitelist` was supplied to establish the expected cell count.
- **Fix**: Supply a cell list via `-w` or `--complexity`, or set an explicit numeric depth limit (e.g., `--site-depth-limit 1000`), or disable it with `--site-depth-limit 0`.

#### `--molecules-per-cell must be at least 1`
- **Cause**: `--molecules-per-cell 0` was passed.
- **Fix**: Pass a value of at least 1 (default: 4 for diploid cells; 8 for replicated G2 cells or doublets).

---

## Output Reading and HDF5 Compression Plugins

### Python: `KeyError: "Unknown compression filter '32001'"`
- **Cause**: Attempting to read `blosc` (filter ID 32001) or `blosc-zstd` datasets with `h5py` without the Blosc HDF5 plugin registered.
- **Fix**: Import `hdf5plugin` (installed with taps-sc-extract) before `h5py` touches the file:
  ```python
  import hdf5plugin  # registers the Blosc filters; import first
  import h5py

  with h5py.File("output_shards/master.h5", "r") as f:
      data = f["/CG/cell_barcode/1"][:]
  ```

### R / Amethyst: Filter Errors Reading Datasets
- **Cause**: The `rhdf5` package requires `rhdf5filters` to decompress `blosc` or `lzf` datasets.
- **Fix**: Install and load `rhdf5filters` in your R session:
  ```R
  BiocManager::install("rhdf5filters")
  library(rhdf5filters)
  library(amethyst)
  library(rhdf5)
  ```
  *(Note: Standard `gzip` compression requires no additional plugins in R).*

---

## Binary Architecture, Native Builds, and Illegal Instructions

### `SIGILL` / Illegal Instruction Fault
- **Cause**: A native build compiled on one host (`./install.sh`, which uses `-C target-cpu=native` and `-march=native`) was transferred to another machine with an older CPU lacking supported vector extensions (e.g., AVX2 or AVX-512).
- **Fix**: Rebuild natively on the target host:
  ```bash
  ./install.sh
  ```
  Or create a portable x86-64-v3 binary:
  ```bash
  ./install.sh --portable
  ```

### `failed to load source for dependency hts-sys`
- **Cause**: The vendored `hts-sys` source archive is missing from `taps-sc-extract-rs/vendor/`.
- **Fix**: Run `./install.sh`, which downloads and unpacks `hts-sys-2.2.1` and applies overlap clipping patches.

### Parity Drift Between Rust and Python After Pulling Updates
- **Cause**: Updates were pulled from git that modified `taps-sc-extract-rs/patches/`, but the vendored source was not re-patched.
- **Fix**: Run `./install.sh`. The install script detects patch checksum discrepancies against `vendor/hts-sys-*/.taps_patched` and re-applies patches automatically.

### `SIMD BAQ differs from htslib on <bad> read(s); use --baq-simd off`
- **Cause**: Hardware floating-point contractions or vector compiler optimizations caused SIMD BAQ base qualities to diverge from scalar htslib glocal realignment.
- **Fix**: Pass `--baq-simd off` to force scalar htslib BAQ execution, and file an issue report with CPU model details.

---

## Memory Limits, Open Files, and Concurrency

### Out-of-Memory (OOM) Process Termination
- **Symptoms**: Process abruptly terminated by the Linux OOM killer (`Killed`).
- **Root Cause**: Running `--memory-mode memory` on large datasets with excessive concurrent shard writers.
- **Solutions**:
  1. Switch to disk-streaming mode:
     ```bash
     TMPDIR=/scratch/tmp taps-sc-extract ... --memory-mode stream
     ```
  2. Set an explicit memory ceiling:
     ```bash
     taps-sc-extract ... --max-memory-gb 32
     ```
  Stream mode writes its temporary files under `$TMPDIR` (default `/tmp`) with the Rust engine; `--temp-dir` applies only to the Python engine.
  3. Restrict concurrent writer processes:
     ```bash
     taps-sc-extract ... --max-writer-threads 4
     ```

### `Too many open files` / `EMFILE` Error
- **Cause**: In multi-BAM runs, each worker maintains open file descriptors for every input BAM (`workers × inputs` descriptors). If the shell descriptor limit is low (e.g., 1024), descriptor exhaustion occurs.
- **Fix**: Raise the shell file limit before executing:
  ```bash
  ulimit -n 65535
  ```
  The Rust engine raises its own soft limit up to the hard limit, so this error means the hard limit is too low: raise it (`ulimit -Hn`, or `/etc/security/limits.conf`), or lower `-t`.

### HDF5 Concurrency Degradation with Many Writer Threads
- **Symptoms**: Increasing writer threads slows down HDF5 generation while memory usage climbs.
- **Root Cause**: In-process writer threads (`--writer-mode thread`) serialize behind HDF5's global library mutex.
- **Fix**: Ensure `--writer-mode process` is enabled (the default whenever `--shards > 1`). Each shard is written by a dedicated child process with an independent HDF5 instance.

---

## `--test` needs the Rust engine

`taps-sc-extract --test` runs the Rust binary's tuner. Build it with `./install.sh`; the pip-only (Python engine) install cannot tune.

## Common Pitfalls Summary

| Issue | Root Cause | Solution |
| :--- | :--- | :--- |
| **`HDF5 file locking error`** | Network filesystem (NFS/Lustre) or concurrent multi-process access | Set `export HDF5_USE_FILE_LOCKING=FALSE` in your environment (set automatically in writer subprocesses). |
| **Worker Freeze / Deadlock** | Python engine using `fork` with `pysam` or `h5py` | The Python engine always uses `spawn`; only custom scripts that fork after opening BAM/HDF5 files hit this. |
| **Memory Spike / OOM** | Holding all chunk results in RAM | Use `--memory-mode stream`, set `--max-memory-gb`, or lower `--max-writer-threads`. |
| **Broken `master.h5` Links** | Shard files moved or separated from `master.h5` | `master.h5` uses relative `ExternalLink` targets (`shard_XXX.h5`). Move the entire shard directory together without renaming individual files. |
| **Slow Multi-Threaded Writing** | In-process threads contending on HDF5 mutex | Use `--writer-mode process` (default when `--shards > 1`). In-process threads are capped at 6. |
| **Build Dependency Error** | Missing vendored `hts-sys` crate | Run `./install.sh` to download and patch dependencies. |
| **`SIGILL` Crash** | Host-native CPU build moved to an older CPU | Rebuild on the target host using `./install.sh`, or build with `./install.sh --portable`. |
| **Execution Timing Variance** | Concurrent compilation or background I/O | Do not run builds during benchmarks; repeat measurements and report spreads. |
