# Custom Single-Cell TAPS Methylation Extractor
**Detailed Implementation Plan**

## 1. Overview & Goals

Build a Python tool that:

- Takes an Astair-aligned BAM (TAPS / mCtoT chemistry).
- Extracts cell barcodes and calls methylation in **CG** and **CH** (CHG+CHH) contexts.
- Processes **by barcode** (not by genomic coordinate) so that only a limited number of output files are open at once.
- Writes:
  1. Premethyst-style cellCov files (`.CG.cov` / `.CH.cov`)
  2. Amethyst-compatible HDF5
- Scales to **200k+ cells** via barcode batching + multiprocessing.
- Uses Python + `pysam` + `pyfaidx` (reference) + `h5py`.

---

## 2. Key Design Decisions (per user requirements)

| Requirement | Decision |
|-------------|----------|
| Process by barcode, not reference | Primary mode: name-sorted BAM; stream consecutive reads belonging to the same barcode. Optionally build a barcode→offset index for random access. |
| Avoid opening 2×cells files at once | Process barcodes in batches. Only write files for the current batch. |
| Multiprocessing | Split barcode list into batches; each worker processes one batch independently. |
| Output | Both Premethyst cellCov **and** Amethyst HDF5. |
| Language | Python + pysam + pyfaidx + h5py (+ multiprocessing / concurrent.futures). |
| Barcode source | Default: first field of QNAME (`QNAME.split(':')[0]` or `QNAME.split('/')[0]`). Optional: read `CB:Z:` tag if present / requested. |
| Directory hierarchy | Only when writing cellCov files. Default output root = `"."`. User supplies `-o / --outdir` prefix. Hierarchy optional (e.g. `CG/xx/barcode.CG.cov`). |

> **Note on indexing**
> Standard BAM indexes (`.bai`/`.csi`) are coordinate-based. They cannot retrieve "all reads for barcode X".
> `pyfaidx` indexes **FASTA**, not BAM.
> Correct approach:
> 1. Prefer a **name-sorted** BAM (by barcode) and process sequentially.
> 2. Optionally pre-build a simple barcode → list of file offsets index (pickle / sqlite / custom binary) for coordinate-sorted BAMs.

---

## 3. Input / Output Specifications

### 3.1 Inputs

| Input | Description |
|-------|-------------|
| BAM | Astair-aligned, preferably **name-sorted by barcode**. Must be indexed only if using coordinate-sorted fallback. |
| FASTA | Reference genome (indexed with `samtools faidx` / readable by `pyfaidx`). |
| Optional barcode whitelist | Text file, one barcode per line. |
| Optional complexity / cellInfo | For filtering low-coverage barcodes (Premethyst-style). |

### 3.2 Premethyst-style cellCov files

Matches what `premethyst calls2h5` expects (see [premethyst `calls2h5.py`](https://github.com/adeylab/premethyst)):

```
<barcode>.CG.cov
<barcode>.CH.cov
```

Tab-delimited, **no header**:

```
chr    pos    pct    t    c
```

- `pos` — 1-based genomic position
- `t` — unmethylated count
- `c` — methylated count
- `pct` — `100 * c / (c + t)` (or 0 if coverage = 0)
- Files sorted by `(chr, pos)` before writing.

**Directory layout (only when cellCov is requested)**
```
{outdir}/
  CG/
    {barcode[:2]}/          # optional 2-char hierarchy to avoid huge directories
      {barcode}.CG.cov
  CH/
    {barcode[:2]}/
      {barcode}.CH.cov
```
Hierarchy is optional (`--no-hierarchy` writes flat into `{outdir}/`).

### 3.3 Amethyst-compatible HDF5

Two compatible layouts exist:

**A. Classic Premethyst / older Amethyst layout** (produced by `calls2h5.py`):

```
/CG/<barcode>   → structured dataset (chr, pos, pct, t, c)
/CH/<barcode>   → structured dataset (chr, pos, pct, t, c)
```

**B. Current Amethyst / Facet layout** (v1.0+):

```
/CG/<barcode>/1   → base-resolution observations
/CH/<barcode>/1   → base-resolution observations
```

Schema of each dataset (from Facet / Premethyst):

| Field | Type | Meaning |
|-------|------|---------|
| `chr` | bytes / string | Chromosome |
| `pos` | int | 1-based position |
| `pct` | float | methylation % |
| `t` | int | unmethylated count |
| `c` | int | methylated count |

**Recommendation**: Write layout **A** by default (directly consumable by existing `calls2h5` and older Amethyst). Optionally support layout **B** via a flag (`--h5-layout classic|facet`).

Compression: gzip level 9 (or blosc if available) to match Premethyst.

---

## 4. Processing Architecture

```
1. Discover barcodes
   ├─ Scan name-sorted BAM (or load whitelist / pre-built index)
   └─ Produce ordered list of barcodes (optionally filtered)

2. Split barcodes into batches (e.g. 500–2000 barcodes per batch)

3. Multiprocessing pool
   └─ Each worker receives:
        • batch of barcodes
        • path to BAM + FASTA
        • output settings
      Worker:
        a. Open BAM + FASTA (pyfaidx)
        b. For each barcode in batch:
             - Collect all reads for that barcode
             - Call methylation (mCtoT, CG/CH)
             - Accumulate (chrom, pos) → (meth, unmeth)
             - Write .CG.cov / .CH.cov (if requested)
             - Buffer structured arrays for HDF5
        c. Return / write partial HDF5 or records

4. Main process merges partial HDF5 files (or writes final HDF5)
```

### 4.1 Retrieving reads for one barcode

**Primary (name-sorted BAM)**
```python
# Sequential scan; group consecutive identical barcodes
for read in bam.fetch(until_eof=True):   # or iterate without fetch
    bc = extract_barcode(read)
    if bc != current_bc:
        flush_previous_barcode()
        current_bc = bc
    process_read(read)
```

**Optional (coordinate-sorted or random access)**
Pre-build once:

```python
# barcode → list of (file_offset) or (ref, start, end) intervals
# Store as pickle / sqlite / custom index
```

Then for a barcode: seek to each offset or use multiple `bam.fetch()` calls. This is secondary; name-sorted is strongly preferred.

### 4.2 Barcode extraction

```python
def extract_barcode(read, use_cb_tag=False):
    if use_cb_tag and read.has_tag("CB"):
        return read.get_tag("CB")
    # Default: first field of QNAME
    qname = read.query_name
    for sep in (":", "/"):
        if sep in qname:
            return qname.split(sep)[0]
    return qname
```

### 4.3 Methylation calling (mCtoT)

For each aligned base that passes filters:

1. Determine OT / OB from flags + library orientation.
2. Using `pyfaidx` reference:
   - Get reference base and neighboring bases → classify **CG** vs **CH**.
3. Apply positive-readout rules:
   - OT, ref C: read T → methylated; read C → unmethylated
   - OB, ref G: read A → methylated; read G → unmethylated
4. Increment counters for `(chrom, pos, context)`.

Filters (configurable):
- Minimum MAPQ
- Minimum base quality
- Ignore soft-clipped bases (default)
- Optional start/end clip (M-bias)
- Primary alignments only

---

## 5. Command-Line Interface (sketch)

```bash
taps-sc-extract \
  -b sample.name_sorted.bam \
  -f genome.fa \
  -o results/ \
  --write-cellcov \
  --write-h5 \
  --h5-layout classic \
  --batch-size 1000 \
  --threads 16 \
  --min-mapq 10 \
  --min-baseq 20 \
  --use-cb-tag \          # optional
  --barcode-whitelist barcodes.txt \
  --hierarchy             # optional 2-char subdirs for cellCov
```

---

## 6. Module / File Layout

```
taps_sc_extract/
  __init__.py
  cli.py                 # argparse entry point
  barcode.py             # extract_barcode, whitelist loading
  calling.py             # mCtoT logic, context classification
  accumulate.py          # per-barcode counters → cov arrays
  writers/
    cellcov.py           # write .CG.cov / .CH.cov
    h5_writer.py         # classic + facet layouts
  index.py               # optional barcode→offset index builder
  parallel.py            # batch splitting + worker function
  utils.py               # logging, filters, validation
tests/
  test_barcode.py
  test_calling.py
  test_writers.py
  test_end_to_end.py
  data/                  # tiny synthetic BAM + FASTA
```

---

## 7. Unit Tests & Logical Checks

### 7.1 Unit tests

| Test | What it verifies |
|------|------------------|
| `test_extract_barcode_qname` | Correct parsing of `BC:UMI#0`, `BC/1`, plain QNAME |
| `test_extract_barcode_cb_tag` | Prefers CB tag when `--use-cb-tag` |
| `test_context_cpg` | Reference CG → CG context |
| `test_context_ch` | Reference CHG/CHH → CH context |
| `test_mCtoT_OT` | T at ref C = methylated; C = unmethylated |
| `test_mCtoT_OB` | A at ref G = methylated; G = unmethylated |
| `test_softclip_ignored` | Soft-clipped bases do not contribute |
| `test_cov_format` | Written `.cov` matches `chr pos pct t c` and sorts correctly |
| `test_h5_classic_schema` | `/CG/<bc>` and `/CH/<bc>` exist with correct dtype |
| `test_h5_facet_schema` | `/CG/<bc>/1` layout when requested |
| `test_empty_barcode` | Barcode with no informative cytosines produces empty or absent datasets |
| `test_batch_isolation` | Two barcodes in different batches never share open files |

### 7.2 Logical / integration checks

- **Round-trip**: Run extractor on a small known BAM → feed `.cov` folder into `premethyst calls2h5` → confirm Amethyst can open the resulting h5.
- **Consistency**: For the same barcode, sum of `c+t` in `.CG.cov` equals number of CG calls observed while walking reads.
- **No cross-talk**: Methylation counts for barcode A never appear under barcode B.
- **Memory bound**: Peak RSS stays roughly proportional to batch size × average sites per cell, not to total cell number.
- **File-handle limit**: Never open more than `2 × batch_size` cellCov files simultaneously.
- **Sorting**: Every written `.cov` and every HDF5 dataset is sorted by `(chr, pos)`.
- **Chemistry**: On a fully converted spike-in (or synthetic read), observed methylation matches expected mCtoT behavior.

### 7.3 Synthetic test data

Create a tiny FASTA + BAM containing:

- 2–3 barcodes
- Known CpG and CHH sites
- Mix of methylated / unmethylated bases under mCtoT rules
- Soft-clipped and low-quality bases (to test filters)

---

## 8. Performance & Scaling Notes (200k+ cells)

- **Batch size**: Start with 500–2000 barcodes. Tune so each worker stays under a target memory (e.g. 4–8 GB).
- **Name-sorted input is essential** for sequential processing without a huge index.
- **HDF5 writing**: Prefer one final HDF5 written by the main process after workers return structured arrays, or use temporary per-batch HDF5 files that are later merged (h5py can copy datasets).
- **cellCov hierarchy**: Strongly recommended at 200k cells to avoid single directories with hundreds of thousands of entries.
- **Optional future**: Direct streaming into a single HDF5 without ever writing per-cell `.cov` files (still keep `.cov` as an optional output).

---

## 9. References to External Specifications

| Spec | Source |
|------|--------|
| Premethyst cellCov / calls2h5 | [adeylab/premethyst `calls2h5.py`](https://github.com/adeylab/premethyst) — expects `chr pos pct t c` |
| Premethyst bam-extract | [bam_extract.pm](https://github.com/adeylab/premethyst) — barcode from QNAME, CG vs CH, `.meth` / `.cov` intermediates |
| Amethyst HDF5 (classic) | Same as Premethyst `calls2h5` output: `/CG/<bc>`, `/CH/<bc>` |
| Amethyst / Facet HDF5 (v1+) | [amethyst-facet](https://pypi.org/project/amethyst-facet/) — `/[context]/[barcode]/1` with schema `chr, pos, pct, c, t` |
| Astair mCtoT chemistry | Astair documentation — methylated C → T (positive readout) |

---

## 10. Implementation Order

1. Barcode extraction + unit tests
2. Reference context + mCtoT calling logic + unit tests
3. Single-barcode accumulator → `.cov` writer
4. Sequential name-sorted BAM driver (one barcode at a time)
5. Batching + multiprocessing
6. HDF5 writer (classic layout first)
7. Optional hierarchy, whitelist, CB-tag, coordinate-sorted fallback
8. End-to-end test against Premethyst `calls2h5` + Amethyst load

---