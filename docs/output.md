# Output files

The HDF5 layout, the per-cell dataset format, the run metadata, and how to read the output from R (Amethyst) and Python. [Back to the documentation index](index.md).

## Single file or sharded directory

| `-o` | Shards | What is written |
|---|---|---|
| `out.h5` | 1 | One file, `out.h5`. Naming an `.h5` file also keeps the automatic shard count at 1. |
| `out.h5` | `--shards N` (N > 1) | A directory `out/`, named after the file stem, holding the files in the next row. |
| `out/` (a directory, or any name without `.h5`) | 1 | Rust: `out/output.h5`. The Python engine writes one file at the path given, so give it an `.h5` name. |
| `out/` | N > 1 (the Rust default for large runs) | `out/shard_000.h5` … `shard_NNN.h5` plus `out/master.h5` |

`master.h5` contains no data. It holds HDF5 external links to every cell in the shards, using relative file names, so the directory can be moved or copied as a whole. Point Amethyst at `master.h5`.

- **Placement:** a cell is never split across shards.
- **Rust:** balances shards by size. Barcodes are hashed into 4 buckets per shard, and the buckets are packed into shards by site count (memory mode) or by temp-file bytes (stream mode).
- **Python engine:** uses `md5(barcode) % N`.
- **Comparing outputs:** because of this, which shard holds a cell differs between engines, and can differ between Rust runs (writers fed during extraction deal buckets as data arrives). Always compare outputs per cell, through `master.h5`, never shard by shard.

## Layout

```
out.h5  (or master.h5)
├── metadata/
│   └── version          dataset: b"amethyst2.0.0"   (+ run attributes, below)
├── CG/
│   ├── <barcode_1>/
│   │   └── 1            dataset: one row per CpG position the cell covers
│   └── <barcode_2>/1
└── CH/
    └── <barcode_1>/1    (CHG and CHH positions together)
```

- A group (`/CG` or `/CH`) is created only when some cell has data in it. With `--context CG` or `--context CH`, the other group is absent.
- A cell appears only if it has at least one call in that context.

## Dataset format

Each `/<context>/<barcode>/1` dataset is a structured array sorted by contig (in the run's contig order) and then position:

```python
dtype = [
    ("chr", "S10"),   # contig name, ASCII bytes, e.g. b"chr1"
    ("pos", "<i8"),   # 1-based position of the C (+ strand call) or G (- strand call)
    ("pct", "<f8"),   # 100 * c / (c + t)
    ("t",   "<i8"),   # unmethylated observations
    ("c",   "<i8"),   # methylated observations
]
```

- **Column order:** `t` (unmethylated) comes before `c` (methylated). This is the Amethyst order.
- **Strands are separate positions:** each position is called from one strand only, the C on the + strand or the G on the − strand. A CpG therefore has two positions, `pos` and `pos + 1`.
- **Observation counts:** `t + c` is the number of reads that called the position in that cell, after the read, base and trim filters and the per-cell cap (`--cell-site-cap`). See [settings.md](settings.md).
- **Chunking:** datasets are chunked in up to 65,536 rows.

## Run metadata

Every output file (each shard and `master.h5`) carries the run's settings as string attributes on `/metadata`:

| Attribute | Meaning |
|---|---|
| `chemistry` | `taps` or `bisulfite` |
| `bsbolt_version` | BSBolt version from the BAM header (Rust; empty for TAPS) |
| `trim_r1`, `trim_r2`, `trim_3p` | Trims applied, in bases |
| `trim_mode` | `auto` (M-bias) or `manual` |
| `min_mapq`, `min_baseq` | Read and base filters |
| `baq` | BAQ mode, or `off` (Rust) |
| `max_depth` | htslib pileup cap (`0` = unlimited) |
| `site_depth_limit` | Positions deeper than this were dropped (`0` = off) |
| `cell_site_cap` | Most observations per cell per position (`0` = off) |
| `blacklist` | BED file used, if any |
| `inputs` | Input BAM(s), comma-separated |
| `extractor` | Engine and version, e.g. `taps-sc-extract-rs 0.1.0` or `taps-sc-extract (python)` |

```python
import h5py
with h5py.File("out/master.h5") as f:
    print(dict(f["metadata"].attrs))
```

## Compression and readers

The Rust engine writes `blosc` (LZ4) by default; the Python engine writes `gzip`. Choose with `--compression` (see [performance.md](performance.md#compression)).

| Codec | R | Python (h5py) |
|---|---|---|
| `gzip`, `gzip6`, `gzip-shuffle`, `none` | plain `rhdf5` | plain `h5py` |
| `lzf` | `library(rhdf5filters)` | plain `h5py` |
| `blosc`, `blosc-zstd` | `library(rhdf5filters)` | `import hdf5plugin` before opening |

```python
import hdf5plugin  # registers the Blosc filter; must come before reading
import h5py

with h5py.File("out/master.h5") as f:
    cells = list(f["CG"])
    x = f["CG"][cells[0]]["1"][()]          # numpy structured array
    print(x["chr"][:3], x["pos"][:3], x["t"][:3], x["c"][:3])
```

Read large outputs one cell at a time. A whole run's CH data can be billions of rows.

## Loading into Amethyst (R)

```r
library(amethyst)
library(rhdf5)
library(rhdf5filters)   # for blosc or lzf output (the Rust default)

h5_path <- "out/master.h5"      # or a single out.h5

ls_df    <- h5ls(h5_path)
barcodes <- ls_df$name[ls_df$group == "/CG"]
obj      <- createObject(h5paths = data.frame(path = h5_path, barcode = barcodes))

obj <- indexChr(obj, "CG")
obj <- indexChr(obj, "CH")       # skip for --context CG output
```

- `indexChr` resolves the external links in `master.h5` itself.
- For output written with `--context CG`, index only `"CG"`. Indexing `"CH"` prints a missing-object error for every cell. For CH-only output, take the barcodes from `"/CH"`.

## Other files a run can write

| File | When | Contents |
|---|---|---|
| `<input>.mbias.tsv`, `.png`, `.pdf` | `--trim auto` (the bisulfite default) | Methylation by position from each read end, and the trims chosen. Written beside the output (in the shard directory, or beside a single file), named after the input BAM, or after the list or sheet with several inputs. See [settings.md](settings.md#trims-and-m-bias). |
| `--profile-json PATH` | on request (Rust) | Per-phase, per-stage and per-window timings and writer reports. See [performance.md](performance.md#profiling-a-run). |
| `--log-file PATH` | on request | The run log |
