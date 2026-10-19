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

## Second pass: writer processes fed during extraction (`4e177be`)

Same machine, BAM and commands. "Write s" is now the phase after extraction ends; the transfer to the writers runs during processing.

| Run | Processing s | Write s | Total s | Processing reads/s | End to end reads/s | Extraction peak GB | Writers peak GB | Output GB | CpG % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| default_1 | 11.11 | 7.00 | 18.12 | 4.81M | 2.95M | 6.5 | 10.4 | 1.44 | 71.615 |
| default_2 | 11.95 | 7.07 | 19.02 | 4.47M | 2.81M | 6.0 | 10.4 | 1.44 | 71.615 |
| default_3 | 11.92 | 7.97 | 19.89 | 4.49M | 2.69M | 6.3 | 10.4 | 1.44 | 71.615 |
| default_4 | 11.93 | 8.00 | 19.93 | 4.48M | 2.68M | 5.8 | 10.4 | 1.44 | 71.615 |
| default_5 | 12.26 | 8.28 | 20.54 | 4.36M | 2.60M | 6.3 | 10.4 | 1.44 | 71.615 |
| stream | 10.88 | 10.48 | 21.37 | 4.91M | 2.50M | 3.4 | 10.4 | 1.44 | 71.615 |
| taps_full | 24.50 | 7.43 | 31.93 | 2.18M | 1.67M | 5.7 | 10.4 | 1.44 | 71.615 |
| taps_partial | 8.34 | 9.35 | 17.68 | 6.41M | 3.02M | 9.2 | 10.4 | 1.44 | 71.585 |
| htslib | 24.37 | 7.60 | 31.97 | 2.19M | 1.67M | 5.9 | 10.3 | 1.43 | 70.859 |
| nobaq | 6.75 | 8.92 | 15.67 | 7.92M | 3.41M | 9.1 | 10.4 | 1.45 | 71.543 |
| simd_off | 17.15 | 7.24 | 24.39 | 3.12M | 2.19M | 5.7 | 10.4 | 1.44 | 71.615 |
| cg | 9.51 | 0.39 | 9.89 | 5.62M | 5.40M | 1.7 | 1.4 | 0.15 | 71.615 |
| ch | 11.62 | 6.88 | 18.50 | 4.60M | 2.89M | 6.6 | 9.9 | 1.29 | — |
| cg_nobaq | 3.14 | 0.42 | 3.56 | 17.02M | 15.02M | 1.7 | 1.3 | 0.15 | 71.543 |
| comp_lzf | 11.48 | 8.53 | 20.01 | 4.66M | 2.67M | 6.1 | 10.4 | 3.92 | 71.615 |
| comp_gzip | 12.90 | 9.36 | 22.27 | 4.14M | 2.40M | 6.3 | 10.4 | 2.00 | 71.615 |
| comp_gzip-shuffle | 11.98 | 12.30 | 24.28 | 4.46M | 2.20M | 6.2 | 10.4 | 1.19 | 71.615 |
| comp_gzip6 | 11.68 | 17.28 | 28.96 | 4.58M | 1.85M | 6.1 | 10.4 | 1.82 | 71.615 |
| comp_blosc-zstd | 11.80 | 20.25 | 32.05 | 4.53M | 1.67M | 6.1 | 10.4 | 1.03 | 71.615 |
| comp_none | 12.84 | 33.02 | 45.87 | 4.16M | 1.17M | 5.7 | 10.4 | 34.40 | 71.615 |
| t16 | 18.01 | 7.70 | 25.71 | 2.97M | 2.08M | 4.7 | 10.2 | 1.44 | 71.615 |
| shards1 | 11.33 | 16.13 | 27.46 | 4.72M | 1.95M | 10.9 | -0.0 | 1.44 | 71.615 |
| cuda | 22.84 | 6.97 | 29.82 | 2.34M | 1.79M | 7.1 | 10.4 | 1.44 | 71.615 |
| cuda_full | 71.71 | 6.84 | 78.55 | 0.75M | 0.68M | 9.0 | 10.4 | 1.44 | 71.615 |
| full_cpu2 | 24.27 | 7.47 | 31.74 | 2.20M | 1.68M | 5.8 | 10.3 | 1.44 | 71.615 |

- Defaults, median of five: processing 4.48M, write after extraction 6.71M, end to end 2.69M reads/s (2.60–2.95M). Extraction process peak 6.3 GB (11.0 GB before).
- `cuda` and `cuda_full` run uncapped GPU launches. With the streamed writers, the 3060 processes 2.34M (2.65M without them); end to end 1.79M (1.74M).
- `taps_partial` and `nobaq` peak higher (9.1–9.2 GB): extraction outruns the niced senders and payloads queue in the parent.
- `shards1` and `stream` are not streamed (one shard is written in-process; stream mode uses temp files).
