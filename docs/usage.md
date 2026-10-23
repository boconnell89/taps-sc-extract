# CLI reference

Every option of `taps-sc-extract`, its type and its real default, plus the Rust binary's subcommands. What the settings do to the calls is in [settings.md](settings.md). Copy-paste examples are in [modes.md](modes.md). [Back to the documentation index](index.md).

`taps-sc-extract` is the usual entrypoint. It defaults to `--engine auto`: the compiled Rust core when `taps-sc-extract-rs` is on `PATH` or in `taps-sc-extract-rs/target/release/`, otherwise the Python reference engine.

Defaults in the tables below are the **effective** values for `taps-sc-extract` with the Rust engine, after chemistry detection. The wrapper leaves many flags unset on purpose so the Rust binary (and the detected chemistry) can fill them. Python-engine defaults that differ are in [Python engine](#python-engine). Flags that exist only on `taps-sc-extract-rs` are in [Rust binary](#rust-binary-taps-sc-extract-rs).

```bash
taps-sc-extract -b sample.bam -f ref.fa -o out/
taps-sc-extract --help
taps-sc-extract --version
```

`-b`, `-f` and `-o` are required except under `--test`. `--test` with a dataset needs both `-b` and `-f`.

---

## Inputs and outputs

| Option | Type | Default | Description |
|---|---|---|---|
| `-b`, `--bam` | path(s) | required | Coordinate-sorted alignments: BAM, CRAM (decoded against `-f`) or SAM / `.sam.gz` (must be coordinate-sorted; copied to a temporary BAM first). A missing index (`.bai`, `.csi` for a contig ≥ 512 Mb, `.crai` for CRAM) is built next to the file, or in a temp directory if that location is not writable. Several inputs: list paths, or one `-b @FILE` (see below). |
| `-f`, `--fasta` | path | required | Reference FASTA with a `.fai` index. Must be the **unconverted** genome the aligner used (a C→T / G→A converted FASTA is refused). |
| `-o`, `--out` | path | required | Output `.h5` file, or a directory of shards plus `master.h5`. Naming an `.h5` file keeps the automatic shard count at 1. With `--shards N` (N > 1) and an `.h5` name, the output is a directory named after the file stem (`out.h5` → `out/`). Rust with a directory and one shard writes `out/output.h5`; the Python engine writes one file at the path given, so give it an `.h5` name. Layout: [output.md](output.md). |
| `-c`, `--chroms` | string | primary chromosomes in the BAM header | Comma-separated contigs, e.g. `chr19` or `chr1,chr2`. Default: names that are a number 1–999 or X, Y, Z, W, with or without a `chr` / `Chr` / `CHR` prefix, sorted numerically then X, Y, Z, W. chrM/MT, unplaced scaffolds, alts, decoys and spike-ins are left out. The log prints what was found, e.g. `22 autosomes + chrX, chrY`. |

**Several inputs.** One run can read several coordinate-sorted, indexed BAMs and write one output. Every barcode must belong to exactly one input.

```bash
taps-sc-extract -b lib-A.bam lib-B.bam -f ref.fa -o out/
taps-sc-extract -b @inputs.txt -f ref.fa -o out/          # one path per line; # comments
taps-sc-extract -b @samples.tsv -f ref.fa -o out/         # sample sheet
taps-sc-extract -b @samples.xlsx -f ref.fa -o out/        # Excel (wrapper converts it; needs openpyxl: pip install .[xlsx])
```

- `@FILE` cannot be mixed with other `-b` paths.
- Relative paths in a list or sheet are resolved against that file's directory.
- A **list** is one path per line (`#` comments and blank lines skipped).
- A **sample sheet** is a TSV (or `.xlsx`) whose header has a `bam` column. Optional per-file columns: `complexity`, `min_complexity`, `max_complexity`, `whitelist`. A blank cell uses the matching global flag. Unknown columns are an error.
- The cells kept are the union of every file's own filter. Either every input has a cell filter, or none do.
- Each window merges every input's reads by position. With `--max-depth 0` (the default), calls equal `samtools merge` followed by a single-BAM run. Each input's index is loaded once and shared; each worker opens its own handle, so the run keeps workers × inputs files open.
- With several inputs, per-file decompression threads default to 0.
- The run stops if the inputs' `@SQ` contigs or lengths differ, their chemistries differ, a barcode is in two inputs' complexity files, or fewer than half of an input's first 5,000 mapped primary reads have a barcode in its own complexity file.
- Excel sheets are converted by the Python wrapper. `taps-sc-extract-rs extract -b @sheet.xlsx` is refused; save as TSV.
- `/metadata` records the input paths in the `inputs` attribute.

Worked examples: [modes.md](modes.md).

---

## Cells

Barcodes are the QNAME prefix before the first `:`, or before `/` if there is no colon.

| Option | Type | Default | Description |
|---|---|---|---|
| `-w`, `--whitelist` | path | none | Barcode whitelist or annotation file. First whitespace-separated column is the barcode. `#` comments skipped. `.gz` accepted. Aliases: `-a`, `--annot` (wrapper only). |
| `--complexity` | path | none | Complexity file: `rank barcode total_reads unique_reads percent` (whitespace-separated). Column 2 is the barcode, column 4 the unique-read count. A header line whose 4th column is not a number is skipped. `.gz` accepted. |
| `--min-complexity` | float | none | Inclusive lower bound on unique reads. Needs `--complexity`. |
| `--max-complexity` | float | none | Inclusive upper bound on unique reads. Needs `--complexity`. |
| `--expected-cells` | int | from the cell filter, else 10,000 | Cell count for auto shards and the memory plan. Precedence: this flag, then the number of barcodes kept by `--whitelist` / `--complexity`, then 10,000. Does not change which cells are extracted. |

With both a whitelist and a complexity file, a cell must pass both. Filtered cells are left out of the HDF5; kept cells' calls are identical to an unfiltered run (the filtered reads stay in the pileup). The Rust engine also skips BAQ for filtered cells' reads.

`--min-complexity` / `--max-complexity` without `--complexity` is an error.

---

## Read and base filters

| Option | Type | Default | Description |
|---|---|---|---|
| `--min-mapq` | int | 10 | Minimum mapping quality. `bwa mem` gives 0 to reads placed at random among equally good hits, and BSBolt to reads that are not bisulfite-unique. `--min-mapq 0` keeps them. |
| `--min-baseq` | int | 20 | Minimum base quality. Bases below this are not called. BAQ can lower qualities first (Rust engine). |
| `--no-overlap-clip` | flag | off (clipping on) | Keep overlapping mate bases. Default: htslib overlap clipping zeros the lower-quality overlapping base, which then fails `--min-baseq`. TAPS only; single-end bisulfite has no mates. With this flag, `--cell-site-cap auto` doubles (two mates can both call). |

Always on:

- Unmapped, secondary, QC-fail and duplicate reads are dropped (SAM flag bits `0x704`).
- Paired reads that are not a proper pair are dropped (`ignore_orphans`; the wrapper has no switch; the Rust binary has `--no-ignore-orphans`).
- TAPS calls only exact flags 99 / 147 (OT) and 83 / 163 (OB). Other reads can sit in the pileup.
- Bisulfite drops supplementary alignments before the pileup (`0x800`). TAPS leaves them in the pileup; they are never called.

Effects: [settings.md](settings.md#mapq-and-base-quality).

---

## Trims

Trims are in BAM SEQ coordinates (soft clips included). Reverse-strand SEQ is reverse-complemented, so the sequenced 5′ end is the high end of SEQ. One window `[lo, hi)` is applied per read: 5′ from R1/R2 (`BAM_FREAD2`) and a shared 3′.

| Option | Type | Default | Description |
|---|---|---|---|
| `--trim` | `auto` \| `manual` | `auto` (both chemistries) | `auto`: measure M-bias on a sample of reads and pick trims for ends you did not set; writes `<bam>.mbias.tsv` / `.png` / `.pdf` next to the output. `manual`: chemistry defaults or the explicit values. An explicit `--trim-5p` / `--trim-r2` / `--trim-3p` always wins for that end. |
| `--trim-r1`, `--trim-5p` | int | 0 | Drop this many sequenced 5′ bases of read 1, or of every single-end read. |
| `--trim-r2` | int | TAPS 10, bisulfite 0 | Drop this many sequenced 5′ bases of read 2. No effect on single-end reads (a warning is printed). |
| `--trim-3p`, `--trim-3prime` | int | TAPS 0, bisulfite 10 (when `--trim manual`) | Drop this many sequenced 3′ bases of every read. Bisulfite: the Tn5 9 bp gap is filled with unmethylated C at end repair. |

`--trim auto` samples about 300,000 mapped reads (every read when the selection is that small; otherwise 64 windows). It counts CpG and CH methylation by distance from each sequenced end on reads at least 40 bp, and trims each end up to where its CpG methylation returns to that end's own level at 15–19 bp (within 2 percentage points or 3 standard errors, for 3 positions in a row; cap 15 bp). Each mate and end is judged on its own, because TAPS read 2 sits about 10 points below read 1 mid-read. With too few calls an end falls back to 0 bp (TAPS read 2's 5′ end: 10 bp; bisulfite 3′: 10 bp). CH (1 point) only warns, unless CpG has too few calls. Both engines choose the same trims from the same table. Plots are drawn by `python -m taps_sc_extract.mbias`.

Effects: [settings.md](settings.md#trims-and-m-bias).

---

## BAQ

Rust engine only. The Python engine leaves BAQ off: pysam's pileup is not given a FASTA, so the HMM is skipped.

| Option | Type | Default | Description |
|---|---|---|---|
| `--baq` | `taps` \| `taps-full` \| `taps-partial` \| `htslib` | TAPS: `taps` (on). Bisulfite: off until this flag is set | How BAQ realigns. The `taps*` modes mask the convertible base (C on OT, G on OB) so methylation conversions are not treated as mismatches. `taps` skips reads that already match that reference exactly (~72% of reads). `taps-full` runs every read. `taps-partial` runs only reads with an indel or clip. `htslib` uses the plain reference (biased against methylated TAPS calls; refused for bisulfite). Setting `--baq` turns BAQ on for bisulfite. |
| `--no-baq` | flag | off | Turn BAQ off, including when `--baq` is also given. |
| `--baq-simd` | `auto` \| `avx512` \| `avx2` \| `off` | `auto` | CPU kernel: 8 reads per vector (AVX-512), 4 (AVX2), or htslib one read at a time. All give identical qualities. `auto` is the widest the CPU runs, or the fastest measured by `--test --save-profile`. |
| `--baq-device` | `cpu` \| `gpu` | `cpu` | Where the BAQ HMM runs. Results are identical. |
| `--gpu-backend` | `auto` \| `opencl` \| `cuda` \| `hip` | `auto` | For `--baq-device gpu`: `auto` tries CUDA, then HIP, then OpenCL. Intel Arc / Data Center GPUs use `opencl`. |
| `--gpu-device` | string | none | Use the first GPU whose name contains this text. |
| `--gpu-mem-mb` | int | 1024 | GPU memory for HMM matrices per launch, in MB. |

Check kernels with `taps-sc-extract-rs baq-selftest`. Effects: [settings.md](settings.md#baq). GPU rates: [performance.md](performance.md).

---

## Depth and site filters

These act after a window is called. They never change the pileup, so calls at other positions stay the same.

| Option | Type | Default | Description |
|---|---|---|---|
| `--max-depth` | int | 0 (unlimited) | htslib pileup cap, all cells together. A finite cap thins a deep position by read order (htslib drops reads that start where the last kept read started). `0` is remapped internally to unlimited (pysam's `0` would leave htslib's 8000). Deep positions are meant to be handled by `--site-depth-limit`. |
| `--site-depth-limit` | int or `auto` | `auto` | Drop every call at a position covered by more kept-cell reads than this, for every cell. Depth counts reads of kept cells that pass the read filters (flags, MAPQ, orphans), deletions included. `auto` = `--molecules-per-cell` × reads per molecule (TAPS 2, bisulfite 1) × cells kept by `--complexity` (with its cutoff; every barcode in the file without one) or `--whitelist`. Off (`0`) without a cell list. Passing `auto` without `--complexity` or `--whitelist` is an error. |
| `--cell-site-cap` | int or `auto` | `auto` | Most observations one cell keeps at an output site. `auto` = `--molecules-per-cell` / 2 (a site is called on one strand: one molecule per haplotype), × 2 for TAPS with `--no-overlap-clip`. A cell over the cap keeps its t:c ratio, scaled to the cap (methylated rounded half up). `0` = off. |
| `--molecules-per-cell` | int | 4 | Molecules one cell can have at a position, both strands: 2 strands × 2 haplotypes. `8` allows a replicated (G2) genome or a doublet. Must be ≥ 1. Sets `--cell-site-cap auto` and `--site-depth-limit auto`. |
| `--blacklist` | BED path | none | Drop every call inside these regions, for every cell. BED columns 1–3 (0-based, half-open); `track` / `browser` / `#` lines skipped; `.gz` read too. |

When to raise or turn off the caps: [settings.md](settings.md#depth-and-site-filters).

The Python API (`extract_methylation_parallel`) defaults `--cell-site-cap` to off (`0`). The CLI uses `auto`.

---

## Contexts and callers

| Option | Type | Default | Description |
|---|---|---|---|
| `--context` | `both` \| `CG` \| `CH` | `both` | Contexts to extract. `CG` or `CH` skips the other context at the reference position, after classification, and leaves its HDF5 group (`/CG` or `/CH`) out. The kept context's calls equal a `both` run, cell for cell. BAQ and the pileup still run on every read. Case-insensitive (`cg` / `ch`). |
| `--caller` | `reads` \| `pileup` | `reads` | Rust engine. `reads` calls one read at a time, replaying htslib's pileup depth cap and mate-overlap rules. `pileup` uses htslib's column pileup. Same calls; `reads` is faster. The Python engine is always a column pileup. |

Effects: [settings.md](settings.md#context) and [settings.md](settings.md#caller).

---

## Chemistry

| Option | Type | Default | Description |
|---|---|---|---|
| `--chemistry` | `auto` \| `taps` \| `bisulfite` | `auto` | Conversion chemistry. `auto` samples 10,000 mapped primary reads: any BSBolt `YS` tag → bisulfite; paired reads without `YS` → TAPS; otherwise an error. An explicit choice that contradicts the BAM is an error. One chemistry per run; several inputs must all resolve to the same one. |
| `--strand-source` | `auto` \| `flag` | `auto` | Bisulfite. `auto`: strand from the flag, and every read's `YS` must be directional (`W_C2T` / `C_C2T`) and agree with the flag (the run stops at the first disagreement). `flag`: skip the `YS` check (aligners without `YS`). Non-directional `*_G2A` / PBAT libraries are not supported. |

Per-chemistry defaults fill only options you did not set:

| Option | TAPS | Bisulfite |
|---|---|---|
| `--min-mapq` | 10 | 10 |
| `--min-baseq` | 20 | 20 |
| `--max-depth` | 0 | 0 |
| `--trim` | `auto` | `auto` |
| `--trim-r1` | 0 | 0 |
| `--trim-r2` | 10 | 0 (no read 2) |
| `--trim-3p` | 0 | 10 when `--trim manual` |
| BAQ | on (`taps`) | off |
| Supplementary | kept in the pileup | dropped before the pileup |

Calling tables: [chemistry.md](chemistry.md). Effects: [settings.md](settings.md#chemistry-and-defaults).

---

## Performance and parallelism

Tuning advice and measured rates: [performance.md](performance.md).

| Option | Type | Default | Description |
|---|---|---|---|
| `-t`, `--threads`, `--workers` | int | auto | Window workers (Rayon threads in Rust). Unset or `0`: min(logical CPUs, workers that fit the stream-mode footprint in the memory budget), or the saved `--test` profile. The wrapper alias `--threads` is not a Rust-binary flag; call `taps-sc-extract-rs` with `-t` / `--workers`. |
| `--decomp-threads` | int | 4 / 2 / 0 by worker count | BAM BGZF decompression threads per worker. `0` = inflate on the worker thread. Auto: 4 at 1–4 workers, 2 at 5–7, 0 at 8 or more (or the saved profile). Several inputs: 0. |
| `--chunk-size-mb` | int | auto (0) | Genomic window size in megabases. Auto: about 6 windows per worker, clamped to 2–100 Mb. Window size does not change calls. |
| `--shards` | int | auto (0) | Number of output HDF5 shard files. Auto: the larger of a cell-count rule (1 / 8 / 16 / 32 at ≤10k / 50k / 100k / more cells) and estimated sites / 48 million, capped at max(cell rule, nproc/2) and the cell count; after extraction the writer planner may raise it (see `--max-writer-threads`). An `.h5` `-o` keeps 1. An explicit count is kept as given. |
| `--max-writer-threads`, `--writer-threads` | int | 0 (auto) | Maximum concurrent shard writers (`--writer-threads` is a wrapper-only alias). Auto: decided after extraction from the measured shard sizes and the memory free then, as many as fit (a writer holds its whole shard, ~1.15 × its payload), up to the `--test` profile's count or the CPU count. With automatic `--shards` the output may be split into more files (whole buckets, up to 4 per requested shard) when that lets more writers run at once. One writer per shard when writers are fed during extraction. In-process (`thread`) writers are capped at 6. **An explicit value is used as given**, with a WARNING in the log if it is not expected to fit. |
| `--writer-mode` | `auto` \| `process` \| `thread` | `auto` | `process`: one `taps-sc-extract-rs write-shard` child per shard (default when shards > 1). `thread`: in-process pool sharing one HDF5 lock. `auto` = `process` when shards > 1, else `thread`. |
| `--compression-threads` | int | 1 (Blosc codecs) | Blosc threads per shard writer (`blosc` / `blosc-zstd` only). |

Rust hashes barcodes into 4 × shards buckets and packs buckets into shards by site count (memory mode) or temp bytes (stream). A cell stays in one bucket. Python assigns `md5(barcode) % N`. Compare outputs per cell, never per shard file.

`TAPS_STREAM_WRITERS=0` disables feeding writer children during extraction.

---

## Memory

| Option | Type | Default | Description |
|---|---|---|---|
| `--memory-mode` | `auto` \| `stream` \| `memory` | `auto` | `stream`: window results go to compact temp files, then writers read them. `memory`: results stay in RAM; with process writers and more than one shard, each window's payload is piped to the writer children as it finishes. `auto`: memory mode with one streamed writer per shard when the estimated peak (parent + 20% held payload + 16 B/site in the children) fits the budget; otherwise stream. |
| `--no-temp-file` | flag | off | Same as `--memory-mode memory`. Wrapper only; forwarded as `--memory-mode memory`. |
| `--max-memory-gb` | float | unset | Memory budget in GB for the whole process tree (extraction plus writer children). Overrides `--max-memory-fraction`. Minimum 1 GB. |
| `--max-memory-fraction` | float | 0.85 | Share of `MemAvailable` the run may use, clamped to 0.05–1.0. The budget is at least 2 GB. |
| `--temp-dir` | path | `$TMPDIR`, else `/tmp` | Directory for stream-mode temp chunk files, in both engines. Point it at fast local storage, not a network mount. |

The run summary reports the extraction process's peak RSS and the writer children's combined peak. Count both when sizing memory.

**Environment variables** (Rust engine):

| Variable | Default | Effect |
|---|---|---|
| `TMPDIR` | `/tmp` | Where stream-mode temporary chunk files go when `--temp-dir` is not given. |
| `TAPS_STREAM_WRITERS` | on | `0` stops feeding writers during extraction (memory mode hands data over afterwards). |
| `TAPS_STREAM_NICE` | 10 | CPU nice level of streamed writer processes, so extraction gets the CPU first. |
| `TAPS_GPU_SHARE` | on | `0` sends every BAQ batch to the GPU instead of sharing with the CPU. |
| `TAPS_GPU_DEPTH` | 4 | BAQ batches each worker keeps in flight with `--baq-device gpu`. |
| `TAPS_MBIAS_MIN_CALLS` | 500 | Calls a position needs before `--trim auto` uses it. |

---

## Output and compression

| Option | Type | Default | Description |
|---|---|---|---|
| `--compression` | see below | `blosc` (Rust), `gzip` (Python) | HDF5 dataset codec. |

| Codec | What it is | Reading it |
|---|---|---|
| `blosc` | LZ4 level 9 + shuffle (Rust default) | Python: `import hdf5plugin` before `h5py`. R: `rhdf5filters`. |
| `blosc-zstd` | Zstd level 5 + shuffle | Same plugins. Smallest files. |
| `lzf` | HDF5 LZF | R: `rhdf5filters`. Built into h5py. |
| `gzip` | deflate level 1 (Python default) | Plain `rhdf5` / h5py. |
| `gzip-shuffle` | byte shuffle + deflate 1 | Plain `rhdf5`. **Rust engine only**: the Python engine accepts the name but fails when writing. |
| `gzip6` | deflate level 6 | Plain `rhdf5`. **Rust engine only**, as above. |
| `none` | uncompressed | Plain `rhdf5` / h5py. |

The wrapper argparse default is `gzip`; when the user does not pass `--compression`, the wrapper sends `blosc` to Rust. Rates and sizes: [performance.md](performance.md). Schema: [output.md](output.md).

---

## Tuning and profiling

| Option | Type | Default | Description |
|---|---|---|---|
| `--test` | flag | off | Run the machine-tuning probes instead of extracting (needs the Rust engine). `-b` / `-f` optional; together they sample that BAM and predict time and memory. Several `-b` paths or `@FILE` are refused. |
| `--test-quick` | flag | off | Smaller probes (~40 s). |
| `--save-profile` | flag | off | Write `~/.config/taps-sc-extract/machine_profile.toml` (`$XDG_CONFIG_HOME` if set). Later Rust extract runs on this host take default `-t`, writer count, `--decomp-threads` and `--baq-simd` from it (still memory-bounded). |
| `--test-report` | path | none | Also write the profile and report to this file. |
| `--temp-dir-test` | path | `--temp-dir` or `.` | Directory for the synthetic BAM and the temp-filesystem probe. |
| `--profile-json` | path | none | Rust: write per-phase CPU use, per-stage and per-window timings, and per-shard writer reports to JSON. Adds a small amount of processing time. |
| `--log-file` | path | none | Log file. It gets the wrapper's log and, with the Rust engine, everything the engine prints (configuration, progress, summary), which still appears on the terminal too. |
| `-v`, `--verbose` | flag | off | Wrapper debug logging. |

```bash
taps-sc-extract --test --save-profile
taps-sc-extract --test -b sample.bam -f ref.fa -o out/
```

---

## Engine

| Option | Type | Default | Description |
|---|---|---|---|
| `--engine` | `auto` \| `rust` \| `python` | `auto` | `auto`: Rust if the binary is found, else Python. `rust`: require `taps-sc-extract-rs`. `python`: the reference multiprocessing engine (`spawn` workers). |
| `--version` | flag | | Print `taps-sc-extract` version. |

---

## Python engine

`--engine python` (or `auto` with no Rust binary) uses `taps_sc_extract`. Same CLI, these defaults and limits:

| Option | Python engine |
|---|---|
| `-t` / `--workers` | 24 if unset. Pass a positive count (`0` is auto only for Rust). |
| `--decomp-threads` | 1 if unset. |
| `--chunk-size-mb` | 10. |
| `--shards` | 1 if unset. |
| `--compression` | `gzip` (level 1). Codecs: `gzip`, `lzf`, `none`, `blosc`, `blosc-zstd`. `gzip-shuffle` and `gzip6` are Rust-only: the parser accepts them, but the Python writer fails. `blosc` uses hdf5plugin LZ4 at level 1. |
| `--memory-mode` | `memory` (or `--no-temp-file`) keeps results in RAM; any other value, including `auto`, streams to `--temp-dir`. No budget heuristic. |
| `--temp-dir` | Honoured. Stream chunks go under a `taps_extract_*` directory there (system temp if unset). |
| `--max-writer-threads` | `0` (default) = auto: as many as fit available memory, at most 6. A positive count caps it. Writers are threads in the extraction process. |
| `--max-memory-gb`, `--max-memory-fraction`, `--expected-cells`, `--writer-mode`, `--profile-json`, `--caller`, `--baq*`, `--gpu-*` | Parsed, then ignored. |
| `--test` | Always needs the Rust engine. |
| BAQ | Off. `compute_baq` is passed to pysam without a FASTA, so the HMM never runs. Calls match Rust `--no-baq`. |
| `--caller` | Always column pileup. |
| Several inputs | One pileup per file. A finite `--max-depth` is refused. |
| Shard placement | `md5(barcode) % N`. |

Chemistry detection, per-chemistry defaults, trims, site filters, `--context`, cell filters and several-input checks match the Rust engine.

---

## Wrapper vs `taps-sc-extract-rs extract`

The wrapper builds `taps-sc-extract-rs extract ...` and forwards a flag only when you set it, except for a few that it always sends (`-b`, `-f`, `-o`, `--chunk-size-mb`, `--compression`, `--max-writer-threads`, `--writer-mode`, `--chemistry`, `--strand-source`).

| Topic | `taps-sc-extract` | `taps-sc-extract-rs extract` |
|---|---|---|
| Workers flag | `-t` / `--threads` / `--workers` | `-t` / `--workers` only |
| Writer count flag | `--max-writer-threads` / `--writer-threads` | `--max-writer-threads` only |
| Whitelist aliases | `-w` / `--whitelist` / `-a` / `--annot` | `-w` / `--whitelist` |
| `--compression` if unset | Sends `blosc` | clap default `blosc` |
| `--chunk-size-mb` if unset | Sends `0` (auto) | clap default `0` |
| `--temp-dir` | Forwarded when set | `--temp-dir` (default `$TMPDIR`, else `/tmp`) |
| `--log-file`, `-v` | Wrapper logging; `--log-file` also captures the engine's output | No such flags (redirect stderr/stdout) |
| `--no-temp-file` | Becomes `--memory-mode memory` | Use `--memory-mode memory` |
| `--engine`, `--test*` | Wrapper | Use `tune` for `--test` |
| `.xlsx` sheet | Converted to a TSV first | Refused |
| `--start` / `--end` | Not exposed | Restrict to one window on the first contig (0-based, end exclusive) |
| `--keep-temp` | Not exposed | Keep the stream-mode temp directory |
| `--no-ignore-orphans` | Not exposed (orphans always dropped) | Keep paired-but-not-proper reads |

Call the wrapper unless you need a Rust-only extract flag, `stats`, or `baq-selftest`.

---

## Rust binary (`taps-sc-extract-rs`)

```bash
taps-sc-extract-rs extract -b sample.bam -f ref.fa -o out/     # usual extraction
taps-sc-extract-rs stats    -b sample.bam -f ref.fa            # TAPS totals / site dump
taps-sc-extract-rs tune     --save-profile                     # same as taps-sc-extract --test
taps-sc-extract-rs baq-selftest -b sample.bam -f ref.fa        # QUAL memcmp per read
```

`write-shard` is hidden. `extract` spawns it; do not call it yourself.

### `extract`

Same extraction as the wrapper's Rust path. Extra flags: `--start`, `--end`, `--keep-temp`, `--no-ignore-orphans`. Default `--memory-mode auto`, `--compression blosc`, `--shards 0`, `-t 0`, `--chunk-size-mb 0`, `--molecules-per-cell 4`.

### `stats`

TAPS-only totals (and optional per-site TSV). Bisulfite BAMs are refused.

| Option | Default | Notes |
|---|---|---|
| `-t` / `--workers` | 0 (auto from cores and memory) | clap help text that says `min(nproc, 32)` is stale |
| `--chunk-size-mb` | 100 | Not auto |
| `--min-baseq` / `--min-mapq` | 20 / 0 | TAPS hardcoded; no chemistry detection |
| `--max-depth` | 0 | Unlimited |
| `--baq` | `taps` | `--no-baq` turns it off |
| `--caller` | `reads` | |
| `--context` | `both` | |
| `--trim-r1` / `--trim-r2` | 0 / 10 | No `--trim-3p` (always 0) |
| `--pileup` | `htslib` | `noodles` is a slower pure-Rust CIGAR walk |
| `--dump-sites` | none | Per-barcode TSV: barcode, ctx, pos, t, c |
| `--start` / `--end` | none | One window on the first contig |
| `-w` / `--whitelist`, `--complexity`, `--min-complexity`, `--max-complexity` | none | Cell filters, as for extract |
| `--no-overlap-clip`, `--no-ignore-orphans` | off | As for extract |
| `--decomp-threads`, `--baq-simd` | auto | As for extract |

It prints one line to stdout: `cells=N CG=..% (c/total) CH=..% (c/total) CpG=..% CHG=..% CHH=..% elapsed_s=..`. Site filters (`--site-depth-limit`, `--cell-site-cap`, `--blacklist`), GPU BAQ, several inputs and chemistry defaults are not on this subcommand.

### `tune`

Same as `taps-sc-extract --test`. Flags: `--temp-dir` (default `.`), `--out-dir` (default `.`), `--input-dir` (default: the BAM's directory), `-b` / `-f`, `--quick`, `--save-profile`, `--report`.

### `baq-selftest`

Compares SIMD (and optional GPU) BAQ to htslib, read by read. QUAL must be byte-identical. Synthetic reads always; add `-b` / `-f` for real reads.

| Option | Default |
|---|---|
| `-c` / `--contig` | `chr19` |
| `--max-reads` | 2,000,000 |
| `--synthetic-reads` | 200,000 |
| `--gpu` | off |
| `--gpu-backend` | `auto` |
| `--gpu-device` | none |

`install.sh` runs this after the build and stops on a mismatch. With `--gpu`, reads the device does not take still go through the CPU kernel.

### `write-shard` (internal)

Spawned by `extract`. Reads `TAPSCK01` chunks from `--chunk-dir` or `--stdin` and writes one shard `.h5`. clap default `--compression lzf` is overridden by the parent.
