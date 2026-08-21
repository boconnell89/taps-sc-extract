# Custom Single-Cell TAPS Methylation Extractor — Phase 1 (HDF5-first, correctness-over-scale)

## Context

The original plan (`basic-plan.md`, this same directory) was written before the real test data was available to check it against. A three-engine review (Grok, Gemini, and a Claude Plan sub-agent) grounded against the real `260731/` BAMs found several of its core assumptions didn't hold: it assumed a name-sorted BAM (the real files are coordinate-sorted with only `.bai`), it assumed astair-specific tags that don't exist, its `.cov`/HDF5 schema tables contradicted each other on field order, and its OT/OB strand rule was too vague to implement correctly (a naive reading would silently misclassify ~half of all calls).

Two follow-up research passes (via Grok against astair's actual source, and Gemini against amethyst's GitHub repo) resolved the remaining open questions with source-level citations, all independently spot-checked in this session: astair's exact OT/OB flag rule, its context-classification and mate-overlap handling, and amethyst's precise HDF5 structural requirements (including a real, already-manifested truncation bug in production data that justifies excluding alt-scaffold contigs). A Claude Plan sub-agent then verified the local environment directly — amethyst R (v1.0.5) and its `.fai`-backed reference are already installed and usable for real round-trip validation, not just aspirational.

The user then set this phase's direction explicitly: treat the reads as astair-aligned (not a red flag — astair's TAPS path literally wraps plain `bwa mem`, confirmed from its source), make annot-file parsing format-flexible, build the Amethyst HDF5 output before the Premethyst `.cov` output, and prioritize a correct pipeline over a scalable one — multiprocessing/200k-cell batching is explicitly deferred. Clarifying answers narrowed scope further: a chr19-restricted run is sufficient to call phase 1 "done," the target file is `taps_ucbTn5_sp2.srt.bam` only, chrM is excluded from output, and — most usefully — the user already has real `astair mbias` output for this exact sample (`260731/taps_ucbTn5_sp2.srt_Mbias.txt`), so no astair install is needed for ground-truth validation. That file was inspected in this session: aggregated across all 100 read-cycle positions it gives CpG≈65.6%, CHG≈0.41%, CHH≈0.41%, pooled≈2.73% — this is the real validation target, not a generic expectation.

The intended outcome of this phase: a correct, single-process Python tool that reads `taps_ucbTn5_sp2.srt.bam` and writes a per-barcode Amethyst-compatible HDF5 file, validated three ways (unit tests, a real amethyst R load, and agreement with the real mbias numbers above) — before any `.cov` writer, batching, or multiprocessing work begins.

## Recommended Approach

### Scope

In scope: BAM → per-barcode `/CG/<bc>/1`, `/CH/<bc>/1` HDF5 datasets, single process, chr19-only as the primary dev/validation target (full-genome run as a checkpoint after chr19 passes, not a blocker). Target file: `/mnt/e/sciMET_TAPS/260731/taps_ucbTn5_sp2.srt.bam`.

Explicitly deferred: Premethyst `.cov` writer and its directory hierarchy, barcode batching, multiprocessing, `--use-cb-tag` (no `CB:Z:` tag exists anywhere in the real data — drop it rather than keep dead code), the `5base_*` sibling files (different chemistry, out of scope), and installing astair locally (the user's already-downloaded `*_Mbias.txt` files supersede that need).

### Methylation-calling spec

Implement as literal SAM-flag-set membership, matching astair's own source exactly (`astair/caller.py` lines 244-245/256-257, verified by direct grep against the cloned repo at `/tmp/astair_research/astair`):

- `flag in {99, 147}` → OT; `flag in {83, 163}` → OB; anything else (secondary/supplementary/orphan/improper-pair) → excluded. This is provably equivalent to `(is_read1 & !is_reverse) | (is_read2 & is_reverse)` = OT (independently confirmed twice this session via from-scratch pysam substitution-rate scripts: OT-flagged reads show ~2.2% C→T against a ~0.02% background floor), but implement the literal flag-membership form so it's auditable line-for-line against astair and implicitly also enforces primary-alignment-only.
- mCtoT rule (astair's default and TAPS-correct): OT + ref `C`, read `T` → methylated; OT + ref `C`, read `C` → unmethylated. OB + ref `G`, read `A` → methylated; OB + ref `G`, read `G` → unmethylated.
- Context from the reference FASTA trinucleotide (`pysam.FastaFile('/mnt/e/refs/mm10/mm10.fa')`, which already has a `.fai`) at the call position — not from the read. The two strands of a CpG dinucleotide are separate, unmerged sites (astair's own convention). Collapse CHG+CHH into `CH` only at the HDF5-write boundary — amethyst has no separate CHG/CHH group.
- Position: pysam is 0-based, Amethyst's `pos` is 1-based. Apply exactly one `+1` at the single point a call is appended to the output buffer — pin this with a unit test.
- Mate overlap and quality filtering: don't hand-roll this — call `bam.pileup(contig, ignore_overlaps=True, min_base_quality=20, stepper='samtools', compute_baq=True, ignore_orphans=True, max_depth=250, min_mapping_quality=0)`, i.e. astair's own default parameters (`astair/caller.py` lines 48-70, 400-402) against the same underlying htslib pileup engine astair uses. This also means soft-clipped bases are excluded for free (pileup only visits aligned, reference-consuming positions) — no separate CIGAR-walk filter needed.

### Processing architecture

One sequential pass, chromosome-by-chromosome, over the canonical contig list only:

```
CANONICAL_CONTIGS = ['chr1'..'chr19', 'chrX', 'chrY']   # chrM excluded per user decision
```
(21 of the BAM's 66 total contigs; the other 44 are underscore-named alt-scaffolds.)

For each contig (in this fixed order): open a fresh per-barcode accumulator dict, run `bam.pileup(contig, ...)` with the parameters above, classify context + OT/OB + mCtoT for each pileup read, accumulate `(barcode, context) → [(pos, meth, unmeth), ...]`, then flush that contig's accumulated rows into the (single, already-open) output HDF5 file before moving to the next contig, discarding the accumulator.

This gets three properties for free, without any explicit sort/dedup step: position-ascending order (pileup visits a contig in order), chromosome-contiguity (one contig fully flushed before the next starts), and no duplicate `(chr,pos)` rows (each pileup column is visited once per contig pass). It also bounds peak memory to roughly one chromosome's worth of calls across all barcodes rather than the whole genome's — chr19 (the dev target) is the smallest canonical chromosome; chr1 (the largest, ~8% of the BAM's 75.1M mapped reads) sets the real ceiling for the eventual full-genome run.

Exclude alt-scaffolds and chrM by never calling `pileup()` on them (read-filtering time), not by writing-then-dropping — cheaper, and the canonical allowlist above already reflects the user's chrM decision.

### HDF5 writer spec

- Structure: `/CG/<barcode>/1`, `/CH/<barcode>/1` only (amethyst 1.0.5's actual expected layout — confirmed against both the installed R package and a real production file; no "classic" flat `/CG/<bc>` layout).
- dtype: `[('chr','S10'), ('pos','<i8'), ('pct','<f8'), ('t','<i8'), ('c','<i8')]` — **`t` before `c`** (this corrects the original plan's §9 table, which had them swapped; confirmed directly against `/home/oconnelb/PBMC2_facet.h5`). `c` = methylated count, `t` = unmethylated count.
- `pct`: amethyst recomputes this from `c`/`t` on load and ignores the stored value, but compute it correctly anyway (`100*c/(c+t)`, or `0.0` if uncovered) — trivial to get right, and real files do store correct values.
- Contigs: canonical-only (21, per above) — this isn't just tidiness. `/home/oconnelb/PBMC2_facet.h5` already contains truncated, silently-colliding alt-scaffold names (`chrUn_GL00`, `chr14_GL00`, etc.) from the `S10` fixed-width field, confirmed by direct inspection in this session — a real bug already in production, not a hypothetical one.
- One open `h5py.File` with resizable datasets (`maxshape=(None,), chunks=True`) created on first touch per barcode/context and extended per contig — this is a single file handle regardless of barcode count, so the original plan's "avoid opening 2×cells files at once" concern doesn't apply to this output path at all (it's specific to separate-file `.cov` writing, which is deferred).

### Annot/whitelist parsing

Generic: `fields = line.rstrip('\n').split('\t'); barcode = fields[0]; extra = fields[1:]`. Never assume a fixed column count or a fixed barcode length (real data mixes 28nt and 30nt barcodes across sub-libraries in the same file). Used in phase 1 only as an optional whitelist filter (`barcode in whitelist_set`), not for any sample-label-driven branching.

### Implementation order

1. `context.py` — reference-trinucleotide classifier, tested against real mm10 trinucleotides pulled directly from `mm10.fa`.
2. `calling.py` — flag classifier (`{99,147}`→OT, `{83,163}`→OB, else None) + mCtoT interpreter; pure functions, no I/O.
3. `barcode.py` — `extract_barcode(qname)` splitting on first `:`; `parse_annot(path)` per the flexible-column spec above.
4. Single-contig pileup driver wired to steps 1-3, run against `--chroms chr19` on the real BAM via the existing `.bai` (no need to materialize a subset file).
5. `h5_writer.py` — resizable per-barcode datasets, dtype/layout per above; wire it to step 4's chr19 output first.
6. Amethyst R round-trip check (see Verification) against the chr19-only file.
7. Cross-check aggregate CG rate against the real `taps_ucbTn5_sp2.srt_Mbias.txt` numbers (see Verification).
8. Only after 1-7 pass: extend the contig loop to all 21 canonical contigs and run the full BAM, watching peak RSS.

### Unit tests

| Test | Verifies |
|---|---|
| `test_context_classifier_real_ref` | Known real mm10 trinucleotides → correct CG/CH |
| `test_flag_classify_ot_ob` | 99/147→OT, 83/163→OB, else→None |
| `test_mCtoT_interpret` | OT+T=meth, OT+C=unmeth, OB+A=meth, OB+G=unmeth |
| `test_extract_barcode_qname` | Real QNAME format `<bc>:<int>#0` parses regardless of barcode length |
| `test_annot_parse_variable_columns` | 1-col, 2-col, N-col annot lines all parse; column 1 always the barcode |
| `test_pos_conversion` | Known 0-based pileup position → correct 1-based output `pos` |
| `test_h5_dtype_and_layout` | `/CG/<bc>/1`, `/CH/<bc>/1` exist with the exact dtype above |
| `test_h5_chr_contiguous_pos_ascending` | Per-barcode rows are chromosome-blocked and position-ascending |
| `test_canonical_contig_allowlist` | No alt-scaffold or chrM rows ever appear in output |
| `test_chr19_real_bam_smoke` | End-to-end run on the real BAM restricted to chr19 produces valid, non-empty output |

## Critical Files

- `/mnt/e/sciMET_TAPS/code/basic-plan.md` — original plan this revision replaces for phase 1
- `/mnt/e/sciMET_TAPS/260731/taps_ucbTn5_sp2.srt.bam` (+ `.bai`) — sole target BAM
- `/mnt/e/sciMET_TAPS/260731/taps_ucbTn5_sp2.srt_Mbias.txt` — real ground-truth validation target (CpG≈65.6%, CHG≈0.41%, CHH≈0.41%, pooled≈2.73%, computed in this session)
- `/mnt/e/sciMET_TAPS/260731/260731_conversion_tn5_sp.annot` — barcode whitelist source (2 columns in practice, parse generically)
- `/mnt/e/refs/mm10/mm10.fa` (+ `.fai`, already present) — reference for context classification
- `/tmp/astair_research/astair/astair/caller.py` — source-of-truth for the exact flag rule and pileup parameters (session-local clone; re-clone from https://github.com/1156054203/astair if this path is gone by implementation time)
- `/home/oconnelb/PBMC2_facet.h5` — real production HDF5 file, used to confirm dtype/field order and the alt-scaffold truncation risk
- `/home/oconnelb/miniforge3/envs/amethyst_r` — R env with amethyst 1.0.5 + rhdf5, used for the round-trip check below

## Verification

1. **Unit tests**: run the table above (pytest) — all pure-function/small-fixture tests, no real BAM needed except `test_chr19_real_bam_smoke`.
2. **chr19 smoke run**: `python -m taps_sc_extract --bam .../taps_ucbTn5_sp2.srt.bam --chroms chr19 --fasta .../mm10.fa -o chr19_test.h5` — confirm it completes, produces a non-empty file, and peak RSS is sane.
3. **Amethyst R round-trip** (env already installed, this is directly runnable):
   ```r
   # /home/oconnelb/miniforge3/envs/amethyst_r/bin/Rscript
   library(amethyst)
   h5paths <- data.frame(barcode = <barcodes from our output>, path = "chr19_test.h5")
   obj <- createObject(h5paths = h5paths)
   obj <- indexChr(obj, type = "CG")   # must not error
   obj <- indexChr(obj, type = "CH")   # must not error
   print(obj@index[["chr_cg"]])        # sanity-check the per-barcode start/count table
   ```
4. **Ground-truth cross-check**: sum our tool's own `c`/`t` counts across all barcodes for chr19, for CG and CH separately, and confirm the resulting rates are in the same ballpark as the real mbias numbers above (CpG≈65.6%, CH≈0.41%) — exact match isn't expected (mbias is genome-wide across all chromosomes and per-read-cycle; ours is chr19-only and per-genomic-position), but a wildly different number (e.g. CG rate near 0% or near 100%, or CH rate above a few percent) indicates a real bug in the flag/context/mCtoT logic, not a scope mismatch.
5. Only after 1-4 pass: repeat step 2 across all 21 canonical contigs on the full BAM as the phase-1.5 checkpoint.
