# Chemistry and calling

How each chemistry's reads are turned into methylation calls, and why BAQ needs a TAPS-aware reference. For what each setting does to the calls, see [settings.md](settings.md). [Back to the documentation index](index.md).

## Choosing the chemistry

`--chemistry auto` (the default) decides from the first 10,000 mapped primary reads of the BAM:

| What the reads look like | Chemistry |
|---|---|
| Any read has a BSBolt `YS` strand tag | `bisulfite` |
| Paired-end reads, no `YS` tags | `taps` |
| Neither | Error: give `--chemistry` |

- **Contradictions are errors:** if you give `--chemistry` explicitly and it contradicts the BAM, the run stops.
- **One chemistry per run:** with several input BAMs, they must all resolve to the same one.
- **Separate rules:** TAPS and bisulfite use their own strand rule and call table, and are never mixed.
- **Defaults:** each chemistry has its own defaults (next section). They fill only the options you did not set.

| Option | TAPS default | Bisulfite default |
|---|---|---|
| `--min-mapq` | 0 | 10 (BSBolt gives 0 to reads that are not bisulfite-unique) |
| `--min-baseq` | 20 | 20 |
| `--trim` | `auto` | `auto` (M-bias) |
| `--trim-r1` / `--trim-5p` | 0 | 0 |
| `--trim-r2` | 10 | (single-end: no read 2) |
| `--trim-3p` | 0 | 10 when `--trim manual` |
| BAQ | on (`--baq taps`) | off |
| Supplementary alignments | kept in the pileup (never called) | dropped before the pileup |

## TAPS

TAPS (TET-assisted pyridine borane sequencing) converts methylated cytosine (5mC and 5hmC) to dihydrouracil, which polymerase reads as **T**. Unmethylated C stays **C**. This is the opposite of bisulfite.

| Strand | Read pairs (SAM flags) | Reference | Read | Call |
|---|---|:---:|:---:|---|
| OT (original top, +) | 99 / 147 | `C` | `T` | methylated (`c += 1`) |
| OT | 99 / 147 | `C` | `C` | unmethylated (`t += 1`) |
| OB (original bottom, −) | 83 / 163 | `G` | `A` | methylated |
| OB | 83 / 163 | `G` | `G` | unmethylated |

- **Which reads are called:** only properly paired reads with exactly these flags. Other reads, such as orphans and supplementary alignments, can sit in the pileup but are never called.
- **Mate overlap:** where two mates of a pair overlap, htslib's overlap clipping keeps one base call. `--no-overlap-clip` turns this off.
- **Trims:** `--trim auto` (the default) measures M-bias and typically trims read 2's first 9 bp (the Tn5 gap, filled with unmethylated C); `--trim manual` drops read 2's first 10 bases. See [settings.md](settings.md#trims-and-m-bias).

## Bisulfite (BSBolt, single-end)

Bisulfite converts **unmethylated** C to U, read as **T**. Methylated C stays **C**. EM-seq reads the same way. The supported input is a single-end, directional, BSBolt-aligned BAM.

| Strand | Flag | BSBolt `YS` | Reference | Read | Call |
|---|---|---|:---:|:---:|---|
| OT | 0 | `W_C2T` | `C` | `C` | methylated |
| OT | 0 | `W_C2T` | `C` | `T` | unmethylated |
| OB | 16 | `C_C2T` | `G` | `G` | methylated |
| OB | 16 | `C_C2T` | `G` | `A` | unmethylated |

**Checks.** The run stops if:
- **The library isn't directional:** reads with `W_G2A` / `C_G2A` strands are an error.
- **Reads are paired:** paired-end bisulfite is not supported.
- **The strand tag disagrees:** a read's `YS` tag contradicts its flag. `--strand-source flag` skips this check, for aligners that don't write `YS`.
- **The reference is converted:** `-f` is a C→T or G→A converted FASTA rather than the normal genome.

**Input:**
- Supplementary alignments are dropped. BSBolt writes them with `bwa -Y`, and most are chimeras on the opposite converted strand.
- The BAM should be merged and deduplicated first. Several deduplicated BAMs, one per barcode group, can be read in one run: see [modes.md](modes.md).

### Pileup calls vs BSBolt's `XB` tag

BSBolt stores its own per-read calls in the `XB` tag. On a 304-cell mouse library (chr19), calling straight from `XB` with no filtering read CpG methylation 3.9 points low (77.45% against 81.34%):
- Low-quality bases explain 2.9 points of that gap.
- The 3′ end explains 1.0 point.

With the same filters, `XB` and the pileup agree on all but 274 of 1.22 million CpG cell-sites. The rest come from a BSBolt soft-clip bug. `scripts/xb_audit.py` reproduces this comparison on any BAM, adding one filter at a time.

## BAQ for TAPS

BAQ (htslib's `sam_prob_realn`) realigns each read to the reference with a profile HMM. It lowers base qualities where the alignment is uncertain, and bases that fall below `--min-baseq` are not called.

**Why plain BAQ is biased for TAPS:** on the plain reference, a methylated C read as T (A on the − strand) counts as a mismatch. So plain BAQ lowers quality at and around methylated sites, and drops methylated calls preferentially.

**What the `taps` modes do instead:** they give the HMM a reference with the convertible base masked to `N`, for every context (CpG and CH):
- every `C` for + strand reads;
- every `G` for − strand reads.

`N` scores as neither a match nor a mismatch, so a converted or unconverted base costs nothing, while indels and other mismatches still lower quality. Only BAQ sees the masked reference; calling uses the real one.

chr19 of the test TAPS dataset, change in calls against `--no-baq`:

| BAQ | CG methylated | CG unmethylated | CG % | CH methylated | CH unmethylated | CH % |
|---|---:|---:|---:|---:|---:|---:|
| off (`--no-baq`) | 574,228 | 281,919 | 67.07 | 59,592 | 19,984,017 | 0.297 |
| `--baq htslib` (plain reference) | −3.00% | −0.37% | 66.48 | −4.93% | −0.42% | 0.284 |
| `--baq taps-full` (masked) | −0.07% | −0.84% | 67.24 | −0.60% | −0.38% | 0.297 |

**What the default skips:** about 72% of reads have an all-M alignment and match the masked reference at every base. For those reads the HMM mostly lowers end bases that the read-2 trim or `--min-baseq` drop anyway. The default `--baq taps` skips them, at a cost of 0.0007% of CG calls against running BAQ on every read. The trade-offs between the modes are in [settings.md](settings.md#baq), and their speed is in [performance.md](performance.md).

The bisulfite default is no BAQ. BSBolt, Bismark and MethylDackel don't apply it either. `--baq htslib` is refused for bisulfite, because it would treat every unmethylated C as a mismatch.

**Exact results whatever the hardware:** the HMM runs on SIMD kernels (8 reads per vector with AVX-512, 4 with AVX2) or optionally on a GPU. Both give base qualities identical to htslib's scalar code for every read. `taps-sc-extract-rs baq-selftest` checks this read by read, and `install.sh` runs it.
