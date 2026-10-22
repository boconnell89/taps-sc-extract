# Bisulfite extraction on the merged LER BAM (2026-10-21)

First measurements of `--chemistry bisulfite` (auto-detected) on real data.

## Input

| | |
|---|---|
| BAM | `AdeyLER1-TTGGACTC.bbrd.q10.bam`: BSBolt 1.2.1, single-end Ultima, two runs merged (one `@RG` each; run 2's QNAMEs end in `_2`), deduplicated, MAPQ ≥ 10, 91.6M reads, 10.2 GB |
| Reference | `/mnt/e/refs/mm10/BSB_ref.fa`, the BSBolt index FASTA the header names (unconverted; contig lengths equal mm10) |
| Cells | `--complexity AdeyLER1-TTGGACTC.complexity.txt --min-complexity 10000`: 304 of 10,073 barcodes, holding 96% of the unique reads |
| Host | Reference workstation (Ryzen 9 5900XT, 16 cores / 32 threads, ~58 GB in WSL2); input, temp and output on the WSL ext4 VHDX. Machine idle (load < 0.6) |

## Results (every run below gave identical totals)

- Settings filled by the chemistry: MAPQ ≥ 10, base quality ≥ 20, no BAQ, supplementary dropped, `--trim auto`.
- `--trim auto`: 5′ 0 bp, 3′ 11 bp, from 316,518 sampled reads in 0.6 s.
- CpG 83.847% (44,408,215 methylated of 52,963,419 calls); CH 5.684% (66,338,843 of 1,167,150,704).
- Output: 1.15 billion cell-sites (49.2M CG, 1.10 billion CH).

The 3′ M-bias depends strongly on read length. In reads under about 100 bp, CpG methylation drops to 58–73% at the last 2 bases and recovers over about 11 bases (plateau 85.5%). Reads of 150 bp or longer barely show it. A likely cause is that only short reads are sequenced up to the Tn5 fill-in end.

## Pileup against BSBolt's XB tags (chr19, 304 cells)

`scripts/xb_audit.py`: 32.9M XB calls from 1.92M reads.

| Calls from XB | CpG % | CH % | CpG cell-sites differing from the pileup |
|---|---:|---:|---:|
| XB as is | 77.45 | 5.01 | 204,279 |
| + supplementary dropped | 77.45 | 5.01 | 204,279 |
| + MAPQ ≥ 10 | 77.45 | 5.01 | 204,279 |
| + base quality ≥ 20 | 80.37 | 5.62 | 133,887 |
| + 3′ trim 11 bp | 81.33 | 5.80 | 274 of 1.22M |
| **taps-sc-extract pileup** | **81.34** | **5.80** | — |

Tag-based extraction without filters underestimates CpG methylation by 3.9 points on chr19, and every cell sits below the pileup. Low base quality accounts for 2.9 points and the 3′ end for 1.0. Supplementary and MAPQ filters change nothing on this BAM. The 274 sites left are BSBolt's soft-clip XB bug (706 misplaced letters) plus a few calls only the pileup makes. Audit run: 9 min, 5.5 GB peak (compared in 5 Mb blocks; see `b6838cb`).

## Performance

Rates use the 91.6M reads. "End to end" is the run's own `elapsed_s` (startup, auto trim, extraction and HDF5 writing). Runs 2 and later had the BAM in the page cache.

| Settings | Runs | Processing | End to end | Peak memory |
|---|---|---:|---:|---|
| Wrapper defaults (`-t` 24, `--shards 1`, stream), cold cache | 1 | 1.69M reads/s | 1.06M reads/s | 15.4 GB |
| Wrapper defaults, warm | 2 | 6.5–6.7M reads/s | 1.88–1.93M reads/s | 15.3 GB |
| Rust binary auto (32 workers, 1 shard, stream) | 2 | 5.4–7.4M reads/s | 1.54–2.03M reads/s | 15.4 GB |
| `-t 32 --shards 16`, stream (auto), 16 writer processes | 2 | 2.0–3.9M reads/s | 1.47–2.03M reads/s | 5.2 GB + writers |
| `-t 32 --shards 16 --memory-mode memory` (5 writer processes) | 2 | 3.9–7.0M reads/s | 2.07–2.51M reads/s | 17.3 GB + 5 × ~1 GB |

Repeats of the same settings differ by up to 2× in processing time, and this was not traced. Treat the ranges as the result.

### Findings

1. **Processing is fast.** At 5.4–7.4M reads/s it beats the TAPS default (4.48M reads/s) because bisulfite has no BAQ, and it is close to TAPS `--no-baq` (7.92M). Reading the 10 GB BAM from a cold cache ran at 188 MB/s and took 54 s.
2. **HDF5 writing limits the run.** The output has 1.15 billion cell-sites (12.5 per read) from 13.3 calls per read; the TAPS test BAM has 16.4 calls per read. With one shard, one writer takes 33 s after extraction. With 16 shards, writing takes 16–23 s.
3. **The wrapper always passes `-t 24 --shards 1`.** The Rust binary therefore treats them as explicit, so `taps-sc-extract` without `--shards` writes a single file with one writer. This is not specific to bisulfite. Rust's own auto setting also picks 1 shard for 304 cells.
4. **The memory planner is conservative in memory mode.** It estimates 28.9 GB for extraction (measured: 17.3 GB), so it allows only 5 writer processes for 16 shards. Streamed writers need one writer per shard, so they stay off and HDF5 writing cannot overlap extraction. In stream mode it estimates a peak of 26.8 of 34 GB but still spills to disk.

### Next

- Shard count: auto by output sites as well as cells (few large cells → several shards).
- The wrapper should forward `-t` and `--shards` only when set.
- Measure memory mode's real per-site cost and correct the planner, so 16 writers stream during extraction as with TAPS.

Each of these needs a TAPS re-check, since the defaults are shared.

## After the fixes (same day)

Changes:
- **Wrapper:** `taps-sc-extract` forwards `-t` and `--shards` only when given.
- **Shard count:** auto now uses estimated output sites / 48M as well as the cell count. It is capped at half the cores, and stays at 1 when `-o` names an `.h5` file.
- **Memory plan:** it now models streamed writers, preferring memory mode with one writer per shard when that fits.
- **`--max-depth auto`** (default): 4 molecules × reads per molecule × cells kept by `--complexity`/`--whitelist`; 100000 when neither is given.

The plain command, with only the complexity cutoff given, now picks 32 workers, 16 shards, memory mode, and 16 writer processes fed during extraction.

| Run (plain command) | Processing | End to end | Peak memory |
|---|---:|---:|---|
| LER, run 1 (BAM partly evicted from cache by the test suite) | 2.08M reads/s | 1.65M reads/s | 6.2 GB extraction |
| LER, runs 2 and 3 | 6.05–6.06M reads/s | 3.41M reads/s | 7.2 GB extraction + 15.9 GB writers (plan: 32.5 GB) |
| TAPS test BAM (53.5M reads, no cell list, depth 100000) | 4.46M reads/s | 2.79M reads/s | 6.4 GB extraction + 10.4 GB writers |

On LER this is 1.8× the previous best warm end-to-end rate with the old wrapper defaults (1.9M reads/s). On the TAPS test BAM the plain command now matches the README's tuned best (`-t 32 --shards 16 --memory-mode memory`, 2.69M reads/s). Before this change it ran with 1 shard and 24 workers. The writers used 13.9 B/site against the 16 assumed.

### Depth limit on LER

With 304 cells the auto limit is 1,216 reads. It removes 0.81% of CpG calls (52,963,419 → 52,532,405), and the CpG rate moves from 83.847% to 83.836%. Almost all of the loss is in four regions:

| Region | CG calls (limit 100000 → 1216) |
|---|---|
| chr2:98.6–98.7 Mb | 400,312 → 50,821 |
| chr9:35.3–35.4 Mb | 59,178 → 8,942 |
| chr9:3.0–3.1 Mb | 68,912 → 56,872 |
| chr6:103.6–103.7 Mb | 16,068 → 4,660 |

The chr2 locus is mm10's known satellite pileup artifact, and the chr9:3 Mb bin sits in pericentromeric satellite. chr19 has no position over 1,216 reads. The limit thins these regions rather than removing them, because htslib drops only reads that start where the last kept read started.
