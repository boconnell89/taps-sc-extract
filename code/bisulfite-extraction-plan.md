# Bisulfite (BSBolt, Ultima SE) extraction in taps-sc-extract — plan

**Status:** final plan, not implemented. Do not test on the pilot BAMs until the other group's `sc-bam-rmdup` has produced merged, deduplicated input. Updated 2026-10-20 with the user's notes (plots as PNG + PDF with a BAM-based prefix, 60 bp plateau floor, pileup-vs-XB report as a headline deliverable, barcode and CHH answers). The TAPS depth cap was raised to 100,000 the same day (§11).
**Date:** 2026-10-20.
**Scope:** extraction only. Merge and deduplication belong to another group (`sc-bam-rmdup --multi`).

**How it was made.** The user's proposal was revised by Claude using measurements on the pilot BAMs, then went through three rounds of adversarial review by Grok 4.6 (xhigh) and Gemini 3.8 Flash (high). Each round reviewed the previous revision and the other reviewer's findings. Claude re-checked every claim adopted here against the BAMs or the source code. Working files (proposal v0, evidence, plan v1–v3, six reviews) are in `scratch/bs-plan/`, which is not in git.

## 1. Facts the plan rests on (verified)

**Input**
- Pilot input: `/mnt/e/sciMETv3/UltimaForLER/*.nsrt.bam`. Single-end, directional (BSBolt strands `W_C2T` and `C_C2T` only), BSBolt 1.2.1, mm10.
- These are name-sorted per-run files. The same library was sequenced in two runs, so the same cells appear in both and deduplication must be joint (upstream).
- Ultima reads are 25–261 bp (mean ~88); 30% are shorter than 60 bp.
- About 0.6 indel operations per read (homopolymer errors).
- Base qualities are binned: Q35 67%, Q25 21%, Q≤19 ~12%.
- Soft clips are rare (~0.6% of reads); there are no hard clips.

**Tags and strands**
- `YS` (`W_C2T`/`C_C2T`) agrees with flag bit 16 on every mapped record checked (183,535 of 183,535). There is no `XG` tag, which only appears from BSBolt 1.6.0.
- BSBolt's `XB` tag holds its per-read calls (`bwa.c` `bwa_gen_cigar2`, `getMethylationContext`):
  - It is a run-length string over aligned query bases (M/=/X/I, soft clips excluded). Each letter is one call, and consecutive letters are common.
  - `x/X` = CG, `y/Y` = CHG, `z/Z` = CHH; uppercase = methylated. This is not Bismark's convention.
  - No base-quality filter and no trim. Contexts come from the unconverted reference.
  - **BSBolt soft-clip bug:** XB reads unshifted query bases in soft-clipped reads. Every site where XB and a reference pileup disagreed in a 20k-read sample was in a soft-clipped read. Excluding those reads, they agree exactly.
- The legacy 2022 script (`/mnt/e/sciMETv2/sciMETv2_call_bsbolt_bam_methylation.py`) can't be the current extractor unmodified:
  - its tokenizer raises `KeyError` on consecutive letters (~99% of these reads);
  - its positioning (`get_reference_positions()` indexed with a match-plus-insertion pointer) would misplace 15.7% of calls after insertions;
  - it never applies its base-quality setting and keeps supplementary alignments.

  The user's current extractor is to be located (§12).
- Supplementary alignments are mostly chimeras on the opposite converted strand (81 of 92 QNAMEs).
- MAPQ 0 = not bisulfite-unique (7.5% of primaries).
- `XC:i:1` (unconverted read) appears on only 16 of 20k reads.

**Existing tool behaviour**
- `--max-depth` caps total live reads across all cells at a position.
- pysam's `max_depth=0` silently means 8000. A Rust cap of 0 drops reads, and `u32::MAX as i32` = −1.
- mm10 chr1 is N for its first 3 Mb.
- `/mnt/e/refs/mm10/BSB_ref.fa` on this host is the unconverted genome despite its name.

**Pilot M-bias** (not deduplicated). The samples are developing mouse brain, so real non-CpG methylation is plausible and CH is treated as signal, not a conversion failure. Reads ≥60 bp, SEQ coordinates, Q≥20:
- The 5′ end is flat.
- At the 3′ end, CG falls to 55–70% over the last ~5 bases and is back at the 82% plateau by 5 bp.
- CH is 4.85–5.6% overall, with a 3′ dip extending to ~11 bp.

## 2. Chemistry: bisulfite kept strictly separate from TAPS

- **Separate modules.** `taps_sc_extract/bisulfite.py` and `taps-sc-extract-rs/src/bisulfite.rs` hold the bisulfite call table and the bisulfite SE strand classifier. The TAPS tables (`MCTOT_LOOKUP`, `call_mctot`) and the TAPS `classify_strand` (exact 99/147/83/163) are not changed or widened.
- **No shared calling function.** Each engine selects a chemistry object once per run (`strand_of(record)`, `call(strand, ref_base, read_base)`). The calling loops never refer to the TAPS globals directly. Python's `extractor.py` and `parallel_extractor.py` are refactored to take this object; today they hard-code `FLAG_STRAND_MAP`/`MCTOT_LOOKUP`, so a half-done port would silently produce zero calls.
- **Bisulfite table.**
  - OT at ref C: read C → methylated (`c`), read T → unmethylated (`t`).
  - OB at ref G: read G → methylated, read A → unmethylated.
  - Read N or any other base → no call. Lowercase read bases are handled as in TAPS.
- **Chemistry is checked against the BAM and recorded in the output.**
  - `--chemistry {taps,bisulfite,auto}`, default `auto`.
  - `auto` inspects the header and a startup sample of 10,000 primary mapped reads. A BSBolt `@PG` or `YS`/`XB` tags → bisulfite. Paired TAPS flags without bisulfite tags → taps. Anything else → error asking for `--chemistry`.
  - An explicit `--chemistry` that contradicts the BAM signature is an error: TAPS on a BSBolt BAM, or bisulfite on a BAM without bisulfite tags and with paired flags.
  - The chemistry, BSBolt version (if any) and the trims used are written as attributes on the root of every output HDF5: each shard, `master.h5`, and the single-file output. Amethyst ignores extra attributes (checked in the load test).

## 3. Library and strand (bisulfite)

- **SE only.**
  - Bisulfite aborts at the first primary mapped read with `BAM_FPAIRED`, after the flag filter, so supplementary reads can't trip it.
  - TAPS aborts only if all 10,000 sampled primary mapped reads lack `BAM_FPAIRED`. Individual orphans in PE data are normal.
  - PE and non-directional bisulfite are non-goals; any `*_G2A` YS is an error at its first occurrence.
- **Strand.** `--strand-source {auto,ys,flag}`:
  - `auto` (default): the strand comes from flag bit 16, and every read's `YS` is verified against it (a cheap aux lookup; YS is the first tag in BSBolt output). The first disagreement aborts the run through `CancelFlag`, so streamed writers are killed and no partial shards are left. This library is 1:1, so any disagreement means something is wrong.
  - `ys`: YS is authoritative; reads without YS are not called and are counted.
  - `flag`: no YS lookup.
  - The BAQ mask (`baq_is_ob`) uses the same source as calling.
  - The check runs identically in Python, in Rust `--caller reads` and in Rust `--caller pileup`.
- **Primary only.** Strand classification happens after the flag filter.

## 4. Read filters and the depth cap (bisulfite)

- **Supplementary (0x800) dropped before push, bisulfite only.**
  - Rust: chemistry-dependent flag filter in `accept_record`.
  - Python: `pileup(flag_filter=0x704 | 0x800, stepper="samtools")` in both `extractor.py` and `parallel_extractor.py`.
  - TAPS is unchanged: there, supplementary reads are uncalled but occupy cap slots, and changing that would move the frozen references.
- **[candidate] `--min-mapq 10`.** Drops MAPQ 0 (bisulfite-ambiguous) and 1–9 (2.3%).
- **[candidate] `--min-baseq 20`.** Ultima bins mean 20 ≡ 21–24: it keeps Q25 and Q35 (88%).
- **`--max-depth 100000`** (the default for both chemistries since 2026-10-20; see §11): a safety cap for satellite repeats. With `--expected-cells N`, the default is `max(100000, 2N)`.
  - `--max-depth 0` = unlimited: remapped to `i32::MAX` before `bam_plp_set_maxcnt`, `ReadPileup::new` and `pysam.pileup(max_depth=…)`. Values above `i32::MAX` are clamped. Never pass 0 to pysam or `u32::MAX` to htslib.
  - Reads dropped by the cap are printed. Rust `--caller reads` counts them already. The column caller and Python gain the same count (the same-start/live test before push) or are documented as not counting.
- **Chemistry-dependent defaults fill in only flags the user didn't set.**
  - The Python wrapper omits unset flags from the Rust command (`option_was_set`, extended to `--max-depth`, `--min-mapq`, `--min-baseq`, the trims and BAQ). Today it always forwards the TAPS argparse defaults.
  - Rust fills bisulfite defaults when clap's value source isn't the command line.
  - The Python engine fills the same defaults on its own namespace.
  - `--chemistry bisulfite --max-depth 250` stays 250.
- Unchanged invariants: MAPQ, flags and supplementary act before push. Barcode, whitelist, complexity and `--context` act at call time.

## 5. Trimming (at extraction)

- **Options** (both engines), counted in sequenced bases of BAM SEQ (`l_qseq`: soft clips included, hard clips excluded):
  - `--trim-5p N`. For SE it is the 5′ trim (`--trim-r1` is accepted as an alias). For TAPS PE the 5′ trim still comes from `--trim-r1`/`--trim-r2` by `BAM_FREAD2`. `--trim-r2` on SE warns only if the user set it.
  - `--trim-3p N` (new).
  - `--trim auto`, the **bisulfite default** (§5.1). Explicitly set ends always win.
  - **Bisulfite SE manual defaults** (used when auto falls back, or when auto is turned off): `--trim-5p 0`, `--trim-3p 10`.
    - The 3′ value matches Tn5: tagmentation leaves a 9 bp gap that end repair fills with unmethylated C, so the sequenced 3′ end reads as unmethylated over ~9–10 bp.
    - On the pilot, the sharp CpG drop ends by ~5 bp, but CG and CH only regain the plateau by ~9–11 bp. The auto rule picked 11 there.
    - Recheck on deduplicated data.
- **One window function** (and `readcall::call_read`'s `[qlo, qhi)`):
  ```
  t5 = trim_r2 if flag & BAM_FREAD2 else trim_r1     # SE: --trim-5p
  t3 = trim_3p
  if t5 + t3 >= len: no calls                          # guard before any subtraction
  reverse (flag & 16): keep t3 <= qpos < len - t5      # sequenced 3′ is low SEQ
  forward:             keep t5 <= qpos < len - t3
  ```
  With `t3 = 0` this is exactly today's 5′ rule (`test_keep_after_5prime_trim` is unchanged). A soft clip at least as long as the trim absorbs it, so no aligned base is removed.
- **Tests:** forward/reverse × leading/trailing clip × 5′/3′ trim; trims meeting; empty SEQ; a 25 bp read; TAPS R1/R2.

### 5.1 `--trim auto`: M-bias trimming with a plot (bisulfite default; chemistry-independent)

**Pre-pass sampler** (one implementation spec, identical in both engines):
- **Windows:** a fixed schedule independent of thread count. 64 windows of 200 kb, evenly spaced over the primary contigs by cumulative length, skipping all-N stretches. All windows are processed; the run doesn't stop early.
- **Reads:** a read counts if its start lies in the window, so each read is counted once. Reads are called **in full**, ignoring window clipping. Production calling clips to windows, which would flatten M-bias at the window edges.
- **Filters:** the run's MAPQ, flag and base-quality filters. No BAQ, depth cap unlimited, no trims. `--context` doesn't skip CH here.
- **Target:** about 300k primary reads (fewer if the BAM is smaller). That's ~2,000+ CG calls per end position. It takes seconds in Rust and well under a minute in Python.

**Count table:** context (CG, CH) × end (5′, 3′) × position p = 0..P−1 (P = 20) × read-length bin (25–39, 40–59, 60–99, 100–149, ≥150), with methylated and unmethylated counts. For TAPS PE, R1 and R2 are tallied separately.

**Variable read lengths (the normalization):**
- End curves count a position only on reads with `len ≥ 2P` (≥ 40 bp). That way a short read's 3′ artifact can't appear in the 5′ curve, or the reverse.
- The plateau comes from positions ≥P from both ends of reads ≥ 60 bp, so short trimmed Illumina libraries still have a plateau. It is kept separate for CG and CH (never pooled). If fewer than `min_calls` plateau calls exist, auto falls back to the chemistry candidate trims with a warning.
- The trim is a fixed base count applied to every read, short reads included, because end-repair artifacts sit at fixed distances from the read end.
- The plot shows the curves by length bin, so a length-dependent bias is visible.

**Rule** (a pure function of the integer count table; the same code is unit-tested in Python and Rust):
1. Positions with fewer than `min_calls = max(500, plateau_n / 100)` calls are "unknown".
2. A position is biased if |rate − plateau| > max(δ, z·SE), with δ_CG = 2 percentage points and z = 3.
3. Walk inward from p = 0. While the position is biased, the trim is p + 1. Stop at the first run of W = 3 consecutive positions that are in band or unknown. Interior blips after recovery don't count.
4. If any of p = 0..4 is unknown while the plateau has enough data, fall back to the chemistry candidate (3′ = 10 for bisulfite SE) and warn.
5. Cap at `max_auto_trim = 15`, with a warning.
6. The trim applied per end comes from **CG**. CH is computed with δ_CH = 1 point. It is used only when no CG position reaches `min_calls`, and otherwise reported. If CH suggests more than CG, warn (e.g. "CH 3′ suggests 11; using CG 5") without raising the trim.

**PE mapping (TAPS):** `trim_r1` and `trim_r2` come from the R1 and R2 5′ curves; `trim_3p` = max of the R1 and R2 3′ results. Auto stays off by default for TAPS until its effect on the frozen references is evaluated (§11).

**Outputs**
- **Log:** the chosen trims, the reason (biased positions, plateau, n), and any warnings.
- **Files**, named `<prefix>.mbias.tsv`, `<prefix>.mbias.png` and `<prefix>.mbias.pdf`, where `<prefix>` is the BAM file name without `.bam`/`.cram`/`.sam` (e.g. `sample.dedup.mbias.png`). They go in the shard directory, or beside a single-file output.
  - **`.tsv`:** the full count table. Python and Rust must write identical tables for the same BAM.
  - **`.png` / `.pdf`** (same figure): CG and CH methylation by position from the 5′ and 3′ ends, one line per read-length bin plus all reads ≥ 40 bp. The plateau and its ±δ band are shown, unknown positions greyed, the chosen trims marked as vertical lines, and the trims and sample size in the title.
  - **Rendering:** `taps_sc_extract/mbias_plot.py` (seaborn on matplotlib, both added as dependencies) draws the figure from the `.tsv`: a 2×2 grid (CG and CH × 5′ and 3′); x axes in whole bases (major ticks every 2, minor every 1); a shared legend. The prototype is `scratch/bs-plan/mbias-prototype/proto.py`. The CLI calls it after either engine. The Rust binary, run directly, writes the `.tsv` and calls `python3 -m taps_sc_extract.mbias_plot <tsv>` if that module is importable; otherwise it prints the command to run.
- **HDF5 root attributes:** `trim_5p`, `trim_3p`, `trim_mode=auto|manual`.

**Tests**
- Synthetic long reads with an injected 3′ step of k bases → k.
- Flat input → 0.
- Short 3′-biased reads mixed with long flat reads → 5′ stays 0. This is the regression that gave 29 without the length floor.
- An interior blip doesn't move the trim.
- Too little data → fallback plus warning.
- Cap warning.
- CH larger than CG → warning, CG used.
- The pure function gives identical output in both languages for hand-built tables.
- Identical `mbias.tsv` from both engines on a small fixture BAM.

## 6. BAQ

- **No bisulfite-specific BAQ is needed.** The TAPS-aware mask (reference C shown as N to OT reads, G to OB reads) neutralises C↔T and G↔A whichever chemistry produced them; only the call table differs. Plain `--baq htslib` is refused for bisulfite.
- **Default `--no-baq` for bisulfite.**
  - Python passes `compute_baq=False`, and is never given the FASTA (that would run unmasked htslib BAQ).
  - BSBolt, Bismark and MethylDackel don't BAQ.
  - Ultima indels defeat the `taps` all-M skip, so cost approaches `taps-full`.
- **Measured in the filter report:** `taps-partial` / `taps` / `taps-full` against off. Homopolymer indels are where BAQ could matter, so the default changes only if indel-adjacent calls measurably improve at an acceptable cost.

## 7. Inputs, reference, contigs, barcodes

- **One BAM:** the merged, deduplicated, coordinate-sorted, indexed output of `sc-bam-rmdup`. No multi-BAM input.
- **Reference check:**
  - every primary contig in the BAM header exists in the FASTA with the same length;
  - ≥100 kb of A/C/G/T sampled from the middle of the first primary contig, skipping N;
  - abort if `min(C, G) / ACGT < 0.05`. Unconverted ≈ 0.19; a C→T or G→A converted FASTA ≈ 0.
  - Filenames are not judged: this host's `BSB_ref.fa` is unconverted and passes.
- **Contigs:** header-based primary detection (done): chr1..19, X, Y; chrM only via `--chroms`.
- **Barcode:** QNAME prefix before the first `:` (`Z0007-ACGAGCCGAA-TTGGACTC`), the same string the complexity file uses (confirmed: column 2 = barcode, column 4 = unique reads).

## 8. XB audit (`scripts/xb_audit.py`, Python; not part of the engines)

- **Decoder:** digits = skipped aligned-query bases (possibly empty); each letter = one call at the current aligned-query index; trailing digits allowed. The walk advances the XB index over M/=/X/I and the reference over M/=/X/D/N; calls happen only on M/=/X; soft clips are skipped. Claude's decoder from the pilot analysis is the starting point.
- **(a) Self-test:** XB length = M+I on every read. On reads without soft clips, each letter must match reference base, read base and context. Mismatches on soft-clipped reads are counted as the BSBolt bug.
- **(b) Oracle:**
  - Slice: reads without soft clips that lie **fully inside** the interval.
  - Pileup settings: `--min-mapq 0 --min-baseq 0`, no trims, `--no-baq`, `--max-depth 0`, supplementary excluded on both sides.
  - Comparison: `(barcode, chr, 1-based pos, t, c)` in `/CG` and `/CH`, which must be identical. CNG with an N in the middle is CHG in XB and CHH here, but both are `/CH`.
  - Whole-contig runs keep one residual class: reads at contig edges that BSBolt couldn't pad.
- **(c) Pileup vs XB calling: a headline deliverable.** The user needs to see how much the reference pileup at production settings changes results compared with XB calling. Compare against Claude's correct XB decoder; the user's own XB extractor can be added later.
  - **Inputs:** a deduplicated chromosome (or the whole BAM) run both ways, cells matched by barcode.
  - **Report:**
    - totals and per-cell CG and CH calls, and CG %/CH % (scatter and correlation per cell);
    - the fraction of sites present in only one method;
    - at shared sites, the fraction with different (t, c);
    - every difference attributed to a cause: trims, base quality, MAPQ, supplementary, soft-clip bug, context classification, other. Attribution works by toggling one setting at a time from the no-filter oracle setting.
  - **Outputs:** a TSV and a PNG/PDF summary with the same naming prefix.

## 9. Filter-impact report (`scripts/bs_filter_report.py`) and freeze gate

- Run on one deduplicated chromosome. Start from a no-filter baseline and add one filter at a time:
  1. MAPQ (0/1/10/20);
  2. base quality (0/10/20/25, plus CH by quality bin to see whether low-Q bases explain part of the CH level);
  3. trims (auto vs 0–15);
  4. depth cap (250 / 8000 / 100k / unlimited);
  5. BAQ (none / partial / taps / full).
- For each step: calls, cells, CG % and CH %, indel-adjacent calls, XC/YC prevalence, runtime.
- **Freeze gate:** the [candidate] defaults become final, and the bisulfite frozen reference is created, only after this table is recorded in `code/rust-implementation-log.md`.

## 10. Validation and implementation order

1. Chemistry modules and separation, bisulfite SE strand with YS verification, bisulfite flag filter, trim window, depth remap, chemistry-dependent defaults (wrapper and engines), chemistry detection and HDF5 attributes. Python first, then Rust, at `--no-baq`, with unit tests (each item in §2–§5).
2. `scripts/xb_audit.py` (a, b) on a deduplicated slice.
3. Reference and library guards; drop counters.
4. `--trim auto`: sampler, count table, pure rule, `<prefix>.mbias.tsv`, PNG and PDF plots, attributes. Python and Rust; table and trim parity tests.
5. Filter-impact report on deduplicated data → finalize defaults → freeze the bisulfite reference.
6. Parity:
   - Python vs Rust per cell (`--no-baq`);
   - `--caller reads` vs `--caller pileup` per cell, including `--max-depth 1/3`, on the BSBolt slice;
   - TAPS frozen references unchanged.
7. Amethyst `createObject` + `indexChr` (with the new attributes), marker loci, and a performance table on deduplicated data. Docs: README chemistry section; AGENTS.md invariants (chemistry separation, bisulfite flag filter, depth remap, trim window, auto-trim rule).

## 11. Non-goals and separate TAPS questions

- **Non-goals (v1):** production XB engine; multi-BAM input; PE or non-directional bisulfite; `XG`; XC/YC read filters; custom bisulfite BAQ.
- **Done for TAPS (2026-10-20):** `--max-depth` default raised from 250 to 100,000 in both engines, and `--max-depth 0` = unlimited, remapped safely. The cap is counted across all cells at a position.
  - On the TAPS test BAM the old cap dropped 805,351 reads. At 100,000 none are dropped, and the output is identical to unlimited.
  - CpG calls +3.5% (31.09M → 32.19M), CH +1%, 2 more cell × context datasets; CpG 71.62% → 71.47%.
  - About 6% slower end to end.
- **Still open for TAPS:**
  - supplementary reads occupying cap slots (irrelevant at 100,000);
  - making `--trim auto` the TAPS default, which would check its fixed R2 5′ trim of 10 but changes TAPS output.

## 12. Open questions for the user
1. The user's current XB extractor: later. The pileup-vs-XB report (§8c) uses Claude's XB decoder meanwhile.

Answered 2026-10-20:
- the barcode is the full QNAME prefix, the same as the complexity file;
- the complexity file is unchanged;
- CHH: developing mouse brain, so non-CpG methylation is plausible;
- the rmdup SE key issue is already handled by that group.
