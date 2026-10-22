# Settings and their effects

What each setting does to the methylation calls and the output, with measured effects already in the repo. Flag names and defaults: [usage.md](usage.md). Calling tables: [chemistry.md](chemistry.md). [Back to the documentation index](index.md).

Rates below are mapped reads per second on the reference workstation and TAPS test dataset described in [performance.md](performance.md), unless a bisulfite library is named.

---

## Chemistry and defaults

`--chemistry auto` (the default) samples 10,000 mapped primary reads. Any BSBolt `YS` strand tag selects bisulfite; paired-end reads without `YS` select TAPS. An explicit `--chemistry` that contradicts the BAM is an error. Several inputs in one run must all resolve to the same chemistry.

TAPS and bisulfite use separate strand rules and call tables. Calling loops take those tables from the chemistry once per run. Detection keys on `YS` tags (BSBolt 1.2.1 writes `PN:bwa VN:BSB-…`).

Each chemistry fills only the options you did not set:

| Option | TAPS | Bisulfite |
|---|---|---|
| `--min-mapq` | 0 | 10 |
| `--min-baseq` | 20 | 20 |
| `--max-depth` | 0 (unlimited) | 0 (unlimited) |
| `--trim` | `manual` | `auto` |
| `--trim-r1` | 0 | 0 |
| `--trim-r2` | 10 | unused (single-end) |
| `--trim-3p` | 0 | 10 when `--trim manual` |
| BAQ | on (`--baq taps`) | off |
| Supplementary alignments | kept in the pileup, never called | dropped before the pileup |

The FASTA must be the unconverted genome: every used contig must exist at the BAM's length, and about 100 kb from the middle of the first contig must have both C and G at ≥ 5% of A/C/G/T.

Bisulfite is single-end directional only. Paired reads, `*_G2A` / PBAT `YS` tags, or a `YS` that disagrees with the flag (unless `--strand-source flag`) stop the run.

---

## MAPQ and base quality

`--min-mapq` drops whole reads before they enter the pileup (and before BAQ). `--min-baseq` drops individual bases at call time, after BAQ has had a chance to lower qualities.

**TAPS default MAPQ 0.** Properly paired TAPS reads with flags 99 / 147 / 83 / 163 are called even at MAPQ 0.

**Bisulfite default MAPQ 10.** BSBolt assigns MAPQ 0 to alignments that are not bisulfite-unique.

**Base quality 20** is the pileup default in both chemistries. Bases that BAQ lowers below 20 are not called.

On a 304-cell mouse bisulfite library, chr19, calling from BSBolt's `XB` tag with no quality filter read CpG methylation 2.9 percentage points lower than the pileup after `--min-baseq 20`. See [Pileup vs BSBolt XB](#pileup-vs-bsbolt-xb).

---

## Trims and M-bias

Trims drop sequenced bases before calling. They run on BAM SEQ coordinates, soft clips included. Reverse-strand SEQ is reverse-complemented, so the original 5′ end is at the high end of SEQ.

One interval `[lo, hi)` is applied per read:

- 5′: `--trim-r2` bases if the read is read 2 (`BAM_FREAD2`), else `--trim-r1` / `--trim-5p`.
- 3′: `--trim-3p` from the sequenced 3′ end of every read.

With `--trim-3p 0` this is the old 5′-only rule. `--trim-r2` has no effect on single-end reads.

**TAPS manual defaults:** R2 5′ 10 bp (Tn5 / TAPS M-bias on read 2), R1 0, 3′ 0.

**Bisulfite `--trim auto` (the default):** a pass over about 300,000 sampled reads (every read when the selection is that small; otherwise 64 windows dealt by mapped reads) counts CpG and CH methylation by distance from each sequenced end. Only reads at least 40 bp long contribute to the end curves, so the two ends do not mix. The plateau is the mid-read rate on reads at least 60 bp, at least 20 bp from both ends.

The rule, identical in both engines:

1. Walk inward from position 0. While a position's CpG rate is more than 2 percentage points (or 3 standard errors) from the plateau, the trim is that position plus one.
2. Stop at the first 3 consecutive in-band or low-count positions.
3. Cap at 15 bp.
4. A position needs 500 calls (or `TAPS_MBIAS_MIN_CALLS`) before it can set or stop the trim.
5. If CpG lacks data, CH may set the trim (and a warning is logged). If CH wants a longer trim than CpG, CpG is kept and CH is warned.
6. If the end has too few calls, fall back to 0 bp (5′) or the chemistry's 3′ candidate (bisulfite 10, TAPS 0).

An explicit `--trim-5p` / `--trim-r2` / `--trim-3p` wins for that end. The chosen trims are written as HDF5 `/metadata` attributes and as `<bam>.mbias.tsv` / `.png` / `.pdf` next to the output. Plots are drawn by `python -m taps_sc_extract.mbias` (seaborn); a missing plotter logs a note and the run continues.

**Bisulfite `--trim manual`:** 5′ 0 / 3′ 10 bp, for the Tn5 9 bp gap filled with unmethylated C at end repair.

On that 304-cell bisulfite library, the 3′ end accounted for 1.0 percentage point of the gap between unfiltered `XB` calls and the pileup.

---

## BAQ

BAQ (htslib `sam_prob_realn`) realigns each read with a profile HMM and lowers base qualities where the alignment is uncertain. Bases that then fall below `--min-baseq` are not called.

The Python engine leaves BAQ off. Rust `--no-baq` is the setting that matches it.

**Why TAPS needs a masked reference.** On the plain genome a methylated C read as T (A on OB) is a mismatch, so plain BAQ lowers quality at and around methylated sites and drops methylated calls preferentially. The `taps*` modes give the HMM a reference with the convertible base masked to `N`: every C for OT reads, every G for OB reads, in all contexts. Calling still uses the real reference. Details: [chemistry.md](chemistry.md#baq-for-taps).

chr19 of the TAPS test dataset, change in calls against `--no-baq`:

| BAQ | CG methylated | CG unmethylated | CG % | CH methylated | CH unmethylated | CH % |
|---|---:|---:|---:|---:|---:|---:|
| off (`--no-baq`) | 574,228 | 281,919 | 67.07 | 59,592 | 19,984,017 | 0.297 |
| `--baq htslib` | −3.00% | −0.37% | 66.48 | −4.93% | −0.42% | 0.284 |
| `--baq taps-full` | −0.07% | −0.84% | 67.24 | −0.60% | −0.38% | 0.297 |

**Mode trade-off** (full TAPS test dataset, AVX2 kernel). About 72% of reads are all-M and match the masked reference at every A/C/G/T base; for those the HMM mostly lowers end bases that the R2 trim or `--min-baseq` drop anyway.

| BAQ | Reads through the HMM | Processing | End to end | CpG % | CG cell-sites differing from `taps-full` |
|---|---:|---:|---:|---:|---:|
| `taps-full` (HMM on every read) | 99% | 2.19M reads/s | 1.68M reads/s | 71.615 | — |
| `taps` (default) | 28% | 4.48M reads/s | 2.69M reads/s | 71.615 | 172 of 25.9M (0.0007%) |
| `taps-partial` (indel/clip reads only) | 8.5% | 6.41M reads/s | 3.02M reads/s | 71.585 | 15,671 (0.06%) |
| `--no-baq` | 0% | 7.92M reads/s | 3.41M reads/s | 71.543 | 64,056 (0.25%) |
| `htslib` (plain reference) | 99% | 2.19M reads/s | 1.67M reads/s | 70.859 | biased; see the table above |

`--baq taps` is the default for TAPS. Changing its skip rule needs a new frozen chr19 reference.

**Bisulfite.** Default is off. BSBolt, Bismark and MethylDackel do not BAQ. `--baq htslib` is refused (every unmethylated C would look like a mismatch). Passing `--baq taps` (or `taps-full` / `taps-partial`) turns BAQ on.

**SIMD and GPU.** `--baq-simd auto` picks AVX-512 (8 reads/vector), AVX2 (4), or htslib; all give identical qualities (`baq-selftest`). `--baq-device gpu` runs the same HMM on the GPU; results stay identical. Reads longer than 192 bp, or with a band too wide, stay on the CPU. `--test` times GPU vs CPU; the extract default stays `cpu`. GPU rates: [performance.md](performance.md).

---

## Overlap clipping

On by default for TAPS. Where the two mates of a proper pair overlap, htslib 1.21 overlap clipping sets one mate's overlapping base to quality 0 (so it fails `--min-baseq`) and adjusts the other's quality according to whether the two bases agree, so each molecule is called once. Supplementary reads take part. Freeing any read removes its QNAME from the overlap table.

`--no-overlap-clip` keeps both mates' bases. That can double observations at a site, so `--cell-site-cap auto` becomes `--molecules-per-cell` (two calls per molecule) instead of `--molecules-per-cell / 2`.

Bisulfite is single-end; there is nothing to clip.

---

## Depth and site filters

`--max-depth` is the only one of these that changes the pileup. The others remove calls after the window is called, for every cell, and leave other positions untouched.

### `--max-depth` (pileup cap)

All cells together. A finite cap drops a read when it starts at the start of the last accepted read and live reads + 1 (htslib's sentinel) exceed the cap. Dropped reads do not move the position and are removed from the overlap table. The default is **0 = unlimited**.

A cap of 250 (the old default) dropped 805,351 reads on the TAPS test BAM (3.4% of CpG calls). Deep positions are handled by `--site-depth-limit`, which drops a whole position rather than thinning it by read order.

`--max-depth 0` is remapped to unlimited internally. Passing 0 to pysam would leave htslib's 8000; passing `u32::MAX` to htslib would become −1.

The Python engine refuses a finite `--max-depth` with several inputs (it runs one pileup per file). Rust applies the cap to the merged stream.

Frozen chr19 references were made at 250, so tests that compare against them pin `--max-depth 250`.

### `--site-depth-limit`

After calling, drop every call at a position whose kept-cell depth exceeds the limit. Depth = reads of kept cells that passed the read filters (flags, MAPQ, orphans) and whose reference span covers the position, deletions included. Filtered-out cells are not counted.

`auto` (CLI default) = `--molecules-per-cell` × reads per molecule (TAPS 2, bisulfite 1) × cells kept by `--complexity` / `--whitelist`. Without a cell list the limit is off (`0`). Passing `auto` without a cell list is an error.

On a 304-cell mouse bisulfite library the auto limit was 1,216. It removed only collapsed satellite repeats: eight ~100 kb regions.

`0` turns it off. The pileup is unchanged, so calls at other positions equal a run without the filter.

### `--cell-site-cap`

Most observations one cell keeps at one output site. Each site is called on one strand, so `auto` = `--molecules-per-cell / 2`, rounded down and at least 1 (one molecule per haplotype: 2 by default, 4 with `--molecules-per-cell 8`). TAPS with `--no-overlap-clip` multiplies by 2.

A cell over the cap keeps its t:c ratio, scaled to the cap (methylated rounded half up). Cells under the cap are untouched, so unlike `--max-depth` no cell loses data because of read order.

In that 304-cell bisulfite library, the cell-sites with 3 or more observations were PCR duplicates whose read starts differed by 1–3 bp, which exact-position deduplication misses.

**Relax it when real extra copies are expected:** copy-number gains, segmental duplications, aneuploidy (cancer), replicating or polyploid cells, or doublets. Raise `--molecules-per-cell` (8 for one extra genome copy) or set `--cell-site-cap 0`.

The Python API defaults this to off (`0`). The CLI uses `auto`. Frozen chr19 runs pin `--cell-site-cap 0`.

### `--molecules-per-cell`

Default 4 = 2 strands × 2 haplotypes in a diploid, deduplicated library. It is the single knob behind both `auto` limits. Must be at least 1.

### `--blacklist`

BED regions (0-based, half-open; `.gz` accepted) whose calls are dropped for every cell. The pileup is unchanged. Use it for the ENCODE blacklist or other known artefacts.

---

## Cell filters

`--whitelist` and `--complexity` (with optional `--min-complexity` / `--max-complexity`) act at **call time**. Filtered cells' reads stay in the pileup, so depth caps and mate overlaps for kept cells are unchanged. Rust skips only those reads' BAQ. Never dropping them before the pileup is load-bearing: it would free `--max-depth` slots and change kept cells' calls.

Complexity file: barcode in column 2, unique reads in column 4. Bounds are inclusive. With a whitelist, a cell must pass both.

Kept cells' calls are identical to an unfiltered run. Filtered cells are absent from the HDF5.

With several inputs, the kept set is the union of each file's own filter (sheet columns, else the global flags). A barcode in two inputs' complexity files is an error. `--site-depth-limit auto` counts all kept cells.

---

## Context

`--context CG` or `--context CH` skips the other context at the **reference position**, after the trinucleotide is classified and before any per-read work. The pileup, BAQ, depth caps and overlaps stay identical, so the kept context's calls equal a `both` run, cell for cell.

Groups (`/CG`, `/CH`) are created only when a cell has data. A skipped context has no group in shards, single files or `master.h5`. In Amethyst, index only the extracted context (`indexChr(obj, "CG")`).

BAQ still runs on every read, because any read can cover a CpG.

Measured on the TAPS test dataset (CH is ~96% of sites):

| Run | Processing | HDF5 write | End to end | Extraction memory | Writer memory |
|---|---:|---:|---:|---:|---:|
| `both` (default) | 4.48M reads/s | 6.71M reads/s | 2.69M reads/s | 6.3 GB | 10.4 GB |
| `--context CG` | 5.62M reads/s | 139M reads/s | 5.40M reads/s | 1.7 GB | 1.4 GB |
| `--context CH` | 4.60M reads/s | 7.77M reads/s | 2.89M reads/s | 6.6 GB | 9.9 GB |
| `--context CG --no-baq` | 17.0M reads/s | 128M reads/s | 15.0M reads/s | 1.7 GB | 1.3 GB |

CG-only is 101% faster end to end than both contexts (340% with `--no-baq`), with 73% less extraction memory and 87% less writer memory. CH-only is within run-to-run noise of `both`.

---

## Caller

`--caller reads` (Rust default) calls each read once, when htslib would free it from the pileup. It replays the column pileup's depth cap and mate-overlap rules, so per-cell calls match `--caller pileup`.

`--caller pileup` walks htslib's column pileup (the original path, kept for A/B).

On the TAPS test dataset, `reads` processed 36% faster than the column pileup at identical calls. The Python engine is always a column pileup.

---

## Pileup vs BSBolt XB

BSBolt stores its own per-read calls in the `XB` tag (letters `x/X` CG, `y/Y` CHG, `z/Z` CHH; uppercase = methylated). Extracting from `XB` skips the pileup's quality, trim and supplementary filters.

On a 304-cell mouse bisulfite library, chr19:

- Unfiltered `XB` read CpG methylation **3.9 points low** (77.45% against the pileup's 81.34%).
- Low-quality bases explained 2.9 points; the 3′ end explained 1.0 point.
- With the same filters as the extractor, `XB` and the pileup agreed on all but 274 of 1.22 million CpG cell-sites. The remainder is a BSBolt soft-clip bug (an unshifted query passed into `bwa_gen_cigar2`).

`scripts/xb_audit.py` repeats this on any BAM, adding one filter at a time:

```bash
python3 scripts/xb_audit.py -b sample.bam -f ref.fa --h5 out/master.h5 -c chr19 \
    --min-mapq 10 --min-baseq 20 --trim-5p 0 --trim-3p 10 -o prefix
```

Steps: `XB` as is, then supplementary dropped, MAPQ, base quality, trims, then the per-cell cap. It writes `prefix.xb_audit.tsv` and `.png` / `.pdf`.
