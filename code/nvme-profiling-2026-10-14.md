# Workstation profiling snapshot — 2026-10-14 9:47pm

As of **9:47pm on 10-14-26**. Rust engine at commit `026c52e` (`main`).
Dataset: `260731/5base_ucbTn5_sp2.srt.bam` (3.5 GB, 53,461,257 mapped records after filters)
on mm10 (`chr1`–`chr19`, `chrX`, `chrY`). 7,351 barcodes. Defaults: `--trim-r1 0 --trim-r2 10`,
`--min-baseq 20`, LZF, 32 workers, `--decomp-threads 0`.

This is a workstation measurement, not a cluster/NFS figure. `--no-baq` is a **test** setting:
Rust BAQ (`sam_prob_realn`) changes which bases pass `min_baseq`, so CG/CH rates from the
no-BAQ run are not comparable to the BAQ-on runs. Production default remains BAQ on.

## Machine

- Host: WSL2 Ubuntu, Linux 6.18.33.1-microsoft-standard-WSL2, 32 cores, ~60 GB RAM
- **E:** Seagate ST8000DM004 8 TB SATA HDD (NTFS “Big Data”), mounted `/mnt/e` via 9p `msize=65536`
- **C:** SPCC M.2 PCIe NVMe. WSL root `/` is `/dev/sdd` ext4 in `ext4.vhdx` on C:
- Inputs for NVMe runs were copied to `/home/oconnelb/taps_io_test/in/` (WSL ext4). That copy is temporary; see `scratch/WSL_IO_TEST_CLEANUP.md` (gitignored).

Printed **Throughput (Reads)** is mapped BAM records divided by **wall** time (process + HDF5 write).
Estimated MB/s is `file_size × mapped/all` over process time (read) or output bytes over write time, not `iostat`.

## Full-genome Rust runs

Shared flags unless noted: `--engine rust -t 32 --decomp-threads 0 --memory-mode memory --compression lzf`.

| Storage | BAQ | Chunk | Shards | Mode | Wall | Process | HDF5 write | Reads/s (wall) | Peak RSS | Notes |
|---|---|---|---|---|---:|---:|---:|---:|---:|---|
| E: 9p HDD | on | 100 Mb (37 win) | 16 | stream | 205.3 s | 116.2 s | 89.1 s | 260,368 | 21.8 GB | disk-backed temps |
| E: 9p HDD | on | 100 Mb | 16 | memory | 171.2 s | 80.0 s | 91.3 s | 312,209 | 31.5 GB | |
| E: 9p HDD | on | 100 Mb | 1 | memory | 247.2 s | 79.1 s | 168.1 s | ~216k | 53.6 GB | one writer |
| WSL NVMe VHDX | on | 100 Mb (37 win) | 16 | memory | 95.0 s | 63.5 s | 31.5 s | 562,565 | 29.7 GB | BAM/FASTA in page cache |
| WSL NVMe VHDX | **off** | **10 Mb (286 win)** | 16 | memory | **45.3 s** | **13.8 s** | 31.6 s | **1,178,919** | 33.7 GB | test: `--no-baq --chunk-size-mb 10` |

BAQ-on methylation (all four BAQ-on rows): CpG **70.859%** (21,499,669 / 30,341,601), CH **1.181%** (9,625,744 / 814,768,999).

`--no-baq` methylation (not comparable): CpG **71.543%** (22,298,196 / 31,167,467), CH **1.239%** (10,173,839 / 821,406,090). Same 53,461,257 mapped records and 7,351 cells.

Logs:

- E: stream — `260731/5base_ucbTn5_sp2_full_rust.extract.log`
- E: memory 16 — `260731/5base_ucbTn5_sp2_full_rust_memory.extract.log`
- E: memory 1 shard — `260731/5base_ucbTn5_sp2_full_rust_memory_noshard.extract.log`
- NVMe BAQ on — `/home/oconnelb/taps_io_test/logs/extract.log`
- NVMe no-BAQ 10 Mb — `/home/oconnelb/taps_io_test/logs/extract_nobaq_10mb.log`

## NVMe pair (same commit, same BAM copy)

Warm page cache (BAM+FASTA prefetched). 16 shards, 16 LZF writers.

| | BAQ on, 100 Mb | `--no-baq`, 10 Mb |
|---|---:|---:|
| Wall | 95.0 s | 45.3 s |
| Process | 63.5 s | 13.8 s |
| HDF5 write | 31.5 s | 31.6 s |
| Est. BAM read | 58.6 MB/s | 270.3 MB/s |
| Est. HDF5 write | 134.1 MB/s | 126.2 MB/s |
| Chunk progress | 20/37 at 52 s; last 17 windows 11 s | ~1,200 chunks/min from 14%–100%; last 6 windows 0.3 s |

`--chunk-size-mb 10` on mm10 canonical contigs is 286 windows (~9 waves of 32). `--chunk-size-mb 100` is 37 windows (~2 waves; last wave 5 windows).

## Disk (WSL `/dev/sdd` during BAQ-on NVMe run)

`iostat` is not installed; 1 Hz `/proc/diskstats` for `sdd` (`logs/sdd_diskstats.tsv`).

| Phase | sdd reads | sdd writes | %util | peak 1 s |
|---|---|---|---:|---|
| Extract (63 s) | 351 reads, 131 MB (almost all in one 349 IOPS / 130 MB/s second) | ~0 | 0.9% | 349 reads/s |
| HDF5 write (34 s) | ~0 | 1.90 GB, 55 MB/s avg | 1.7% | 948 MB/s write, 29.5% util |
| Whole run | 367 reads | 1.90 GB | 1.1% | |

92 of 96 sampled seconds had zero disk reads. The NVMe was not the limiter.

On `/mnt/e` (9p → HDD), sequential BAM read was ~147 MB/s; 8,000 random 4 KiB preads on that BAM were 191 IOPS / 0.8 MB/s. Linux `iostat` does not see E: 9p traffic.

## Limiting factors as of this snapshot

**On WSL ext4 (NVMe VHDX), BAQ on, 100 Mb windows:** process is CPU — `sam_prob_realn` BAQ plus `bam_plp` pileup/calling — plus coarse windows so many of 32 cores idle on the last large chromosomes. HDF5 is ~14.7k LZF datasets (one `/CG|CH/<barcode>/1` per cell), not disk (1.7% util). RAM (29–34 GB peak of 60 GB) is not the cap.

**Same path with `--no-baq` and 10 Mb windows:** process drops to 13.8 s. HDF5 write stays 31.6 s and is **70% of wall**. Occupancy is good (flat chunk rate, 0.3 s tail). Next lever on this box is shard write (HDF5 object creation / LZF / writer count), not BAM I/O.

**If inputs/outputs stay on E:** 9p (`msize=65536`) plus the 5400 RPM HDD dominate wall time (171 s vs 95 s for the same BAQ-on memory/16-shard job).

Rust applies BAQ when FASTA is available. The Python reference sets `compute_baq=True` but does not pass `fastafile=` into `AlignmentFile.pileup()`, so BAQ is a no-op there. `--no-baq` is for isolating CPU cost; it is not a recommendation for production calls.

## Step 0 addendum: window size vs BAQ (2026-10-14 11:30pm)

Same commit, BAM copy, and flags as the NVMe pair above (`-t 32 --decomp-threads 0 --shards 16 --memory-mode memory`, LZF). Two added runs complete a 2×2 table. Logs are `/home/oconnelb/taps_io_test/logs/step0_baq_10mb.log` and `step0_nobaq_100mb.log`.

Process seconds (extraction only):

| | 100 Mb (37 windows) | 10 Mb (286 windows) |
|---|---:|---:|
| BAQ on | 63.5 s | **54.1 s** |
| `--no-baq` | **22.0 s** | 13.8 s |

| Run | Wall | Process | HDF5 write | Reads/s (wall) | Peak RSS |
|---|---:|---:|---:|---:|---:|
| BAQ on, 10 Mb | 86.8 s | 54.1 s | 32.7 s | 615,962 | 32.1 GB |
| `--no-baq`, 100 Mb | 57.0 s | 22.0 s | 35.0 s | 938,590 | 32.5 GB |

- **Window size does not change calls.** BAQ on at 10 Mb gives the same cells and CG/CH totals as at 100 Mb: 21,499,669 / 30,341,601 and 9,625,744 / 814,768,999. `--no-baq` at 100 Mb matches `--no-baq` at 10 Mb: 22,298,196 / 31,167,467 and 10,173,839 / 821,406,090.
- **Smaller windows save 9.4 s with BAQ on and 8.2 s with it off.** That is the 37-window tail at 32 workers.
- **BAQ costs about 40 s of the 54 s process time at 10 Mb (~75%).** It is the main processing cost once windows are balanced.
- **HDF5 write stays at 31–35 s in all four runs.**
