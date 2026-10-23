# Extraction Modes and Use Cases

Worked examples and copy-pasteable command templates for single-cell DNA methylation extraction workflows. [Back to the documentation index](index.md). Every option is in [usage.md](usage.md); what the settings do to the calls is in [settings.md](settings.md).

---

## Overview of Common Invariants

Across all extraction modes:
- **Automatic Engine Dispatch**: `taps-sc-extract` defaults to `--engine auto`, using the compiled Rust acceleration binary (`taps-sc-extract-rs`) when available and falling back to the multiprocessing Python engine.
- **Automatic Resource Tuning**: Window workers (`-t`), BGZF decompression threads (`--decomp-threads`), genomic window sizes (`--chunk-size-mb`), output shards (`--shards`), and concurrent shard writers (`--max-writer-threads`) auto-tune based on CPU cores, available memory, and dataset scale unless explicitly specified.
- **Safety Defaults**: The pileup is unlimited (`--max-depth 0`). The per-cell cap (`--cell-site-cap auto`, `--molecules-per-cell` / 2 = 2 observations per cell per site) is always on. The position drop (`--site-depth-limit auto`) is on only when a cell list is given (`--complexity` or `--whitelist`); otherwise it is off. See [settings.md](settings.md#depth-and-site-filters).
- **Default Compression**: The Rust engine writes Amethyst-compatible HDF5 datasets using `blosc` (LZ4 level 9); the Python engine defaults to `gzip` (level 1 deflate).

---

## TAPS Paired-End (Default Protocol)

Use this for standard single-cell TAPS (TET-Assisted Pyridine Borane Sequencing) libraries where 5mC/5hmC is chemically converted to dihydrouracil and read as Thymine ($5\text{mC} \to \text{T}$), while unmodified Cytosines remain Cytosine ($C \to C$).

```bash
taps-sc-extract \
  -b sample.bam \
  -f ref.fa \
  -o output_shards/
```

### What Happens
- **Chemistry Detection**: Resolves automatically to `taps` from paired-end alignments lacking bisulfite strand tags.
- **TAPS-Aware BAQ**: Enables `--baq taps` (TAPS-aware profile HMM). OT reads realign against a C-masked reference and OB reads against a G-masked reference so that true methylation calls are not penalized as mismatches. Reads with an all-M alignment that match the masked reference exactly skip the HMM (~72% of reads); 0.0007% of CG calls differ from running BAQ on every read. See [settings.md](settings.md#baq).
- **Trims & Filters**: `--trim auto` measures M-bias per mate and end and trims what it finds (typically read 2's first 9 bp, the Tn5 gap); `--trim manual` uses the fixed read 2 5′ trim of 10 bp. The M-bias table and plots are written next to the output. Minimum base quality is 20 (`--min-baseq 20`) and minimum MAPQ is 10 (`--min-mapq 10`).
- **Overlaps**: Automatically clips mate overlaps to avoid double-counting read pairs.

---

## Bisulfite / EM-seq Single-End (BSBolt)

Use this for directional, single-end bisulfite or EM-seq libraries aligned with BSBolt, where unmethylated Cytosines convert to Thymine ($C \to T$) while methylated Cytosines remain protected ($5\text{mC} \to C$).

```bash
taps-sc-extract \
  -b sample_bsbolt.bam \
  -f ref.fa \
  -o output_shards/ \
  --chemistry bisulfite
```

### What Happens
- **Calling Inversion**: Uses the bisulfite call rules ([chemistry.md](chemistry.md#bisulfite-bsbolt-single-end)). Read C on OT (flag 0) and Read G on OB (flag 16) represent methylated calls.
- **BSBolt Verification**: Inspects read `YS` strand tags (`W_C2T` for OT, `C_C2T` for OB). The run halts immediately if tags disagree with alignment flags or if non-directional/PBAT tags (`*_G2A`) appear. Use `--strand-source flag` only if using an aligner that does not write `YS` tags.
- **Unconverted Genome Check**: Verifies that `ref.fa` is the native reference, rejecting C→T or G→A converted index FASTA files.
- **Automatic M-bias Trimming (`--trim auto`)**: Measures methylation by distance from sequenced ends across ~300,000 sampled reads, identifies where CpG levels stabilize, and trims ends accordingly. TSV, PNG, and PDF diagnostic plots are written to `<output_dir>/<bam_prefix>.mbias.*`. Manual fallback defaults to 5′ 0 bp and 3′ 10 bp (Tn5 end-repair fill-in).
- **Filters**: Enforces `--min-mapq 10` (BSBolt marks non-unique alignments as MAPQ 0), drops supplementary reads (flag `0x800`), and disables BAQ (`--no-baq`) by default.

---

## Multiple Input BAMs (Paths, Lists, and Sample Sheets)

Use this when an experiment spans multiple libraries or index batches (e.g., separate Tn5 index files) that should be merged into a single Amethyst output. Every cell barcode must belong to exactly one input.

### Option A: Direct Command-Line Paths
```bash
taps-sc-extract \
  -b lib-A.bam lib-B.bam lib-C.bam \
  -f ref.fa \
  -o output_shards/
```

### Option B: Path List File (`@list.txt`)
```bash
taps-sc-extract \
  -b @inputs.txt \
  -f ref.fa \
  -o output_shards/
```
`inputs.txt` format (paths relative to `inputs.txt` or absolute; `#` comments and blank lines ignored):
```text
# Batch 1 libraries
lib-A.bam
lib-B.bam
lib-C.bam
```

### Option C: Sample Sheet with Per-File Cutoffs (`@sheet.tsv` or `@sheet.xlsx`)
```bash
taps-sc-extract \
  -b @samples.tsv \
  -f ref.fa \
  -o output_shards/
```
`samples.tsv` format (header must contain `bam`; optional columns: `complexity`, `min_complexity`, `max_complexity`, `whitelist`):
| bam | complexity | min_complexity | max_complexity |
| :--- | :--- | :--- | :--- |
| `lib-A.bam` | `lib-A.complexity.txt` | 5000 | 50000 |
| `lib-B.bam` | `lib-B.complexity.txt` | 7500 | 75000 |

> [!NOTE]
> `.xlsx` workbooks are supported when using the `taps-sc-extract` wrapper CLI (requires `openpyxl`). The wrapper converts the sheet to a temporary TSV before passing it to the engine.

### Operational Behavior
- Reads from all files are coordinate-merged per genomic window via a k-way merge.
- Contigs, lengths, and chemistry must be identical across all inputs.
- Decompression threads per worker default to 0 (`--decomp-threads 0`) with multiple inputs to prevent CPU thread contention.

---

## Filtering Cells (Whitelist and Complexity Files)

Use this to extract only validated single cells, discarding barcodes associated with debris, low-depth droplets, or doublets.

### Whitelist Only
```bash
taps-sc-extract \
  -b sample.bam \
  -f ref.fa \
  -o output_shards/ \
  -w whitelist.txt
```
`whitelist.txt`: one barcode per line (first whitespace-delimited column; `.gz` supported).

### Complexity File with Unique Read Thresholds
```bash
taps-sc-extract \
  -b sample.bam \
  -f ref.fa \
  -o output_shards/ \
  --complexity cells.complexity.txt \
  --min-complexity 5000 \
  --max-complexity 50000
```
`cells.complexity.txt`: five whitespace-separated columns (`rank`, `barcode`, `total_reads`, `unique_reads`, `percent`). Cutoffs filter against column 4 (`unique_reads`).

### Combining Whitelist and Complexity
```bash
taps-sc-extract \
  -b sample.bam \
  -f ref.fa \
  -o output_shards/ \
  -w whitelist.txt \
  --complexity cells.complexity.txt \
  --min-complexity 5000
```
Both filters are applied (intersection).

### Pileup Invariant
Cell filtering occurs strictly at call time:
- Reads from non-whitelisted barcodes remain in the pileup traversal so that depth caps and mate-overlap tracking for valid cells are unaffected.
- The Rust engine skips BAQ computation for non-whitelisted reads to maximize throughput.
- Kept cell counts feed directly into `--site-depth-limit auto`.

---

## Context Selection (CG-Only or CH-Only)

Use this when investigating only CpG methylation or non-CpG (CHG/CHH) methylation.

### Extracting Only CpG (`--context CG`)
```bash
taps-sc-extract \
  -b sample.bam \
  -f ref.fa \
  -o output_shards/ \
  --context CG
```

### Extracting Only Non-CpG (`--context CH`)
```bash
taps-sc-extract \
  -b sample.bam \
  -f ref.fa \
  -o output_shards/ \
  --context CH
```

### Operational Differences
- Positions matching the excluded context are dropped immediately after reference trinucleotide classification, before per-read calling.
- Groups for excluded contexts (`/CH` or `/CG`) are not written to the HDF5 output.
- `--context CG` is roughly 100% faster end-to-end and reduces output file size and memory footprint by ~75–85% because CH sites constitute ~96% of genomic cytosines in mammalian genomes.
- Downstream Amethyst analysis must index only the extracted context: `indexChr(obj, "CG")`.

---

## Single File vs. Sharded Directory Output

Choose between a single monolithic `.h5` file and a multi-file sharded directory based on cell count and filesystem architecture.

### Single file
Use when downstream scripts need one file path:
```bash
taps-sc-extract \
  -b sample.bam \
  -f ref.fa \
  -o all_cells.h5
```
When `-o` specifies a file ending in `.h5` and `--shards` is omitted, `--shards 1` is selected automatically.

### Sharded directory
Give a directory and leave `--shards` unset. The Rust engine picks the shard count from the cell count and the estimated output size, capped at half the CPUs and the cell count. For example, the 7,351-cell TAPS test set gets 16 shards. Shards are written by parallel processes.
```bash
taps-sc-extract \
  -b sample.bam \
  -f ref.fa \
  -o output_shards/
```
The output directory contains `shard_000.h5` … `shard_NNN.h5` and a `master.h5` with relative links to every cell. Pass `--shards N` only to pin the count; with `-o experiment.h5 --shards 16`, the output goes to a directory named `experiment/`. Details: [output.md](output.md#single-file-or-sharded-directory).

---

## Memory Mode vs. Disk-Streaming Mode

Control how intermediate window extractions are transferred to the shard writers.

### In-Memory Streaming (`--memory-mode memory`)
Use on workstations and servers with sufficient RAM (≥32–64 GB):
```bash
taps-sc-extract \
  -b sample.bam \
  -f ref.fa \
  -o output_shards/ \
  --memory-mode memory
```
Alternatively, pass `--no-temp-file`. When running multi-process shard writers (`--writer-mode process`), the Rust engine launches writer children before extraction and streams window results over pipes as extraction progresses, reducing post-extraction write time.

### Disk-Streaming Mode (`--memory-mode stream`)
Use in memory-constrained environments or when extracting massive cell populations (>100,000 cells):
```bash
taps-sc-extract \
  -b sample.bam \
  -f ref.fa \
  -o output_shards/ \
  --memory-mode stream
```
Chunk payloads spill to temporary files and are removed as shards complete. Add `--temp-dir /scratch/tmp` to put them on fast local storage rather than a network mount (default `$TMPDIR`, else `/tmp`).

---

## CRAM and SAM Alignment Inputs

`taps-sc-extract` handles CRAM and SAM alignments natively alongside BAM.

### CRAM Input
```bash
taps-sc-extract \
  -b sample.cram \
  -f ref.fa \
  -o output_shards/
```
Decodes against `ref.fa`. If a `.crai` index is missing, it is built automatically. CRAM decoding requires additional CPU cycles during extraction.

### Coordinate-Sorted SAM Input
```bash
taps-sc-extract \
  -b sample.sam \
  -f ref.fa \
  -o output_shards/
```
SAM files must be coordinate-sorted. The input handler converts the SAM to a temporary BAM in the temporary directory, indexes it, runs extraction, and purges the temporary BAM upon completion.

---

## Python Reference Engine

Use the reference multiprocessing Python engine for validation, environments where Rust binaries cannot be compiled, or programmatic integration.

### CLI Execution
```bash
taps-sc-extract \
  -b sample.bam \
  -f ref.fa \
  -o output.h5 \
  --engine python
```
The Python engine does not auto-size workers (it uses 24 unless `-t` is given) and writes 1 shard unless `--shards` is given.

### Python Programmatic API
```python
from taps_sc_extract.parallel_extractor import extract_methylation_parallel

summary = extract_methylation_parallel(
    bam_path="sample.bam",
    fasta_path="ref.fa",
    out_h5_path="output_shards/",
    n_workers=24,
    n_shards=16,
    compression="gzip",
    trim_r1=0,
    trim_r2=10,
    min_base_quality=20,
    min_mapq=0,
    max_depth=0,
    chemistry="taps",          # the API does not detect the chemistry
    cell_site_cap=2,           # the API leaves the per-cell cap off unless set
)
print("Extracted cells:", summary["total_unique_barcodes"])
```
The API takes options as given: it does not detect the chemistry or fill per-chemistry defaults, and the per-cell cap and position drop are off unless set. The CLI does all of these.

> [!NOTE]
> The Python engine does not calculate BAQ (equivalent to Rust `--no-baq`) and defaults to `gzip` compression. Parity between Rust and Python is verified with `--no-baq`.

---

## GPU-Accelerated BAQ

Use when running on machines equipped with GPUs offering strong double-precision (FP64) compute (e.g., Intel Arc Pro B70 / Data Center GPU, NVIDIA GPUs via CUDA, or AMD GPUs via HIP).

```bash
taps-sc-extract \
  -b sample.bam \
  -f ref.fa \
  -o output_shards/ \
  --baq-device gpu \
  --gpu-backend auto
```

### Options
- `--gpu-backend {auto,opencl,cuda,hip}`: Backend driver selection (`auto` tries CUDA, then HIP, then OpenCL).
- `--gpu-device <NAME>`: Substring match to select a specific GPU in multi-GPU systems.
- `--gpu-mem-mb <MB>`: Dedicated device buffer allocation (default: 1024 MB).

### Hybrid CPU-GPU Work Sharing
The GPU driver pools reads from all worker threads. If the GPU request queue exceeds the worker's processing slack time, the worker processes the batch on the CPU using AVX SIMD kernels instead of idling. Output base qualities are bit-for-bit identical to CPU execution.

---

## Tuning a new machine

Run once on a new host to benchmark CPU SIMD, memory bandwidth, storage systems, and HDF5 writer throughput using a built-in synthetic dataset.

### Machine-Only Calibration
```bash
taps-sc-extract --test --save-profile
```
Generates diagnostic recommendations and saves optimal worker and writer thread counts to `~/.config/taps-sc-extract/machine_profile.toml`. Subsequent Rust runs load these defaults automatically.

### Calibration with Dataset Sampling
```bash
taps-sc-extract --test \
  -b sample.bam \
  -f ref.fa \
  -o output_shards/ \
  --test-report host_profile.json
```
Samples the target BAM and reference to predict total runtime, peak resident memory (RSS), and I/O saturation for the complete run.

---

## Run Profiling (`--profile-json`)

Use to audit performance, profile phase execution times, and diagnose hardware bottlenecks.

```bash
taps-sc-extract \
  -b sample.bam \
  -f ref.fa \
  -o output_shards/ \
  --profile-json run_profile.json
```

### Output Data
`run_profile.json` captures:
- Wall-clock and CPU seconds for extraction and HDF5 write phases
- Process-level CPU utilization percentages
- Stage timing distributions across BAM reading/decompression, BAQ realignment, pileup management, and calling
- Per-window execution timings
- Per-shard writer statistics (read, dataset creation, data population, HDF5 flush, and peak RSS)

---

## Direct Rust CLI Execution (`taps-sc-extract-rs`)

The compiled Rust binary can be executed directly without Python wrapper overhead.

### Direct Extraction
```bash
taps-sc-extract-rs extract \
  -b sample.bam \
  -f ref.fa \
  -o output_shards/
```
It takes the same extraction options as the wrapper ([usage.md](usage.md#rust-binary-taps-sc-extract-rs)).

### Statistics and Parity Dumps (`stats`)
Calculate methylation rates and dump site-level TSV tables without writing HDF5 datasets:
```bash
taps-sc-extract-rs stats \
  -b sample.bam \
  -f ref.fa \
  -c chr19 \
  --dump-sites chr19_sites.tsv
```

### SIMD and GPU Kernel Self-Test (`baq-selftest`)
Validate SIMD and GPU BAQ implementations against scalar htslib glocal realignment:
```bash
taps-sc-extract-rs baq-selftest \
  -b sample.bam \
  -f ref.fa \
  -c chr19
```

---

## Command Option Quick Reference

| Mode / Feature | Key Flags | Default Value |
| :--- | :--- | :--- |
| **Standard TAPS** | `--chemistry taps` | `auto` (detected) |
| **Bisulfite / EM-seq** | `--chemistry bisulfite --trim auto` | `auto` (detected) |
| **Multiple Inputs** | `-b <BAM1> <BAM2> ...` or `-b @<FILE>` | None (required) |
| **Cell Filtering** | `-w <FILE> --complexity <FILE> --min-complexity <N>` | None (all cells) |
| **Context Selection** | `--context {both,CG,CH}` | `both` |
| **Sharded Directory** | `--shards <N>` | Auto (`1` for an `.h5` output; otherwise by cells and estimated output size) |
| **Memory Mode** | `--memory-mode {auto,stream,memory}` | `auto` |
| **GPU Acceleration** | `--baq-device gpu --gpu-backend auto` | `cpu` |
| **Machine Tuning** | `--test --save-profile` | Extraction run |
| **Performance Profile** | `--profile-json <FILE>` | Disabled |
