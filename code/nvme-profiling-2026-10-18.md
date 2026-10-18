# NVMe profiling snapshot — 2026-10-18 (read-centric calling)

Reference workstation (Ryzen 9 5900XT, 16C/32T, WSL2 ext4 on NVMe), full test BAM `5base_ucbTn5_sp2.srt.bam` (53,461,257 mapped reads, 7,351 cells), build `72f01cf` plus the uncapped GPU launches. Base command: `taps-sc-extract-rs extract -t 32 --shards 16 --memory-mode memory` with the defaults (`--caller reads`, `--baq taps`, AVX2, `blosc`) unless the run name says otherwise. Raw wall times in seconds; rates are mapped reads / seconds.

| Run | Processing s | Write s | Total s | Processing reads/s | End to end reads/s | Extraction peak GB | Writers peak GB | Output GB | CpG % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| default_1 | 10.53 | 10.36 | 20.89 | 5.08M | 2.56M | 11.1 | 10.7 | 1.49 | 71.615 |
| default_2 | 10.51 | 11.30 | 21.80 | 5.09M | 2.45M | 10.9 | 10.7 | 1.49 | 71.615 |
| default_3 | 9.95 | 11.62 | 21.57 | 5.37M | 2.48M | 11.0 | 10.7 | 1.48 | 71.615 |
| default_4 | 10.56 | 11.27 | 21.83 | 5.06M | 2.45M | 11.0 | 10.7 | 1.48 | 71.615 |
| default_5 | 10.64 | 11.74 | 22.38 | 5.03M | 2.39M | 10.9 | 10.7 | 1.49 | 71.615 |
| stream | 10.76 | 11.73 | 22.49 | 4.97M | 2.38M | 3.7 | 10.7 | 1.48 | 71.615 |
| taps_full | 22.74 | 12.05 | 34.80 | 2.35M | 1.54M | 11.0 | 10.7 | 1.49 | 71.615 |
| taps_partial | 6.86 | 11.58 | 18.43 | 7.80M | 2.90M | 11.1 | 10.7 | 1.49 | 71.585 |
| htslib | 23.03 | 11.74 | 34.78 | 2.32M | 1.54M | 11.1 | 10.7 | 1.47 | 70.859 |
| nobaq | 6.17 | 11.33 | 17.50 | 8.67M | 3.06M | 10.9 | 10.7 | 1.49 | 71.543 |
| simd_off | 16.99 | 11.82 | 28.81 | 3.15M | 1.86M | 10.7 | 10.7 | 1.48 | 71.615 |
| cg | 8.35 | 0.57 | 8.92 | 6.40M | 5.99M | 2.1 | 1.4 | 0.15 | 71.615 |
| ch | 11.26 | 11.53 | 22.79 | 4.75M | 2.35M | 10.6 | 10.3 | 1.34 | — |
| cg_nobaq | 3.21 | 0.58 | 3.78 | 16.68M | 14.13M | 1.9 | 1.4 | 0.15 | 71.543 |
| comp_lzf | 12.22 | 12.72 | 24.94 | 4.37M | 2.14M | 10.9 | 10.6 | 4.01 | 71.615 |
| comp_gzip | 10.33 | 14.83 | 25.16 | 5.17M | 2.13M | 10.9 | 10.7 | 2.06 | 71.615 |
| comp_gzip-shuffle | 10.87 | 17.18 | 28.04 | 4.92M | 1.91M | 11.0 | 10.7 | 1.23 | 71.615 |
| comp_gzip6 | 10.73 | 21.78 | 32.51 | 4.98M | 1.64M | 11.0 | 10.7 | 1.88 | 71.615 |
| comp_blosc-zstd | 11.55 | 23.13 | 34.68 | 4.63M | 1.54M | 10.9 | 10.7 | 1.06 | 71.615 |
| comp_none | 11.25 | 41.54 | 52.79 | 4.75M | 1.01M | 11.0 | 10.7 | 34.40 | 71.615 |
| t16 | 18.47 | 12.47 | 30.94 | 2.89M | 1.73M | 10.2 | 10.4 | 1.49 | 71.615 |
| shards1 | 10.82 | 14.80 | 25.62 | 4.94M | 2.09M | 11.0 | -0.0 | 1.54 | 71.615 |
| cuda | 22.85 | 11.45 | 34.30 | 2.34M | 1.56M | 12.1 | 10.7 | 1.49 | 71.615 |
| cuda_full | 89.40 | 11.30 | 100.70 | 0.60M | 0.53M | 12.6 | 10.7 | 1.48 | 71.615 |
| full_cpu2 | 26.41 | 11.58 | 37.98 | 2.02M | 1.41M | 11.0 | 10.7 | 1.49 | 71.615 |
| pocl | failed (PoCL kernel build error -11) | | | | | | | | |

- `cuda` and `cuda_full` above ran with the 64k launch cap. Uncapped (the new default), measured afterwards: `cuda` 2.65M processing / 1.74M end to end, `cuda_full` 0.76M / 0.66M. With the cap the GPU starved: 2.34–2.41M, and 1.70M at a 16k cap.
- `shards1` writes in-process, so no writer-process peak is reported.
- Limiting factor at the defaults: processing (5.08M reads/s) and HDF5 write (4.73M reads/s) are now about equal. Worker time is 55% BAQ, 23% calling, 13% BAM read and inflate, 7% pileup bookkeeping, 3% setup.
