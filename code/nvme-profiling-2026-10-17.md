# Workstation snapshot 2026-10-17

AMD Ryzen 9 5900XT (16 cores / 32 threads, 3.3 GHz base, up to 4.8 GHz boost, AVX2), ~58 GB RAM in WSL2, ext4 VHDX on NVMe. `5base_ucbTn5_sp2.srt.bam`: 53,461,257 mapped reads, 7,351 cells, mm10. Rust engine (`./install.sh`, native), BAQ kernel AVX2, TAPS mask applied in the kernel. Defaults: `--compression blosc --baq taps --context both`; all runs `-t 32 --shards 16 --memory-mode memory` unless named otherwise. `--profile-json` on.

Rates are mapped reads per second over each phase (processing, HDF5 write) and over the whole run (end to end). Medians where n > 1.

```
default            n=5 proc  3.29M/s write   4.29M/s e2e  1.85M/s ( 110.7M/min) ext  10.7GB wr  10.5GB out  1.49GB  CpG 71.615% (22,263,688 / 31,087,996)
stream             n=1 proc  3.17M/s write   4.94M/s e2e  1.93M/s ( 115.9M/min) ext   3.0GB wr  10.5GB out  1.48GB  CpG 71.615% (22,263,688 / 31,087,996)
t16                n=1 proc  2.55M/s write   3.96M/s e2e  1.55M/s (  93.1M/min) ext  10.0GB wr  10.2GB out  1.49GB  CpG 71.615% (22,263,688 / 31,087,996)
shards1            n=1 proc  3.17M/s write   3.72M/s e2e  1.71M/s ( 102.7M/min) ext  10.8GB wr  -0.0GB out  1.49GB  CpG 71.615% (22,263,688 / 31,087,996)
simd_off           n=1 proc  2.72M/s write   4.40M/s e2e  1.68M/s ( 100.8M/min) ext  10.5GB wr  10.4GB out  1.48GB  CpG 71.615% (22,263,688 / 31,087,996)
taps_full          n=1 proc  1.85M/s write   4.40M/s e2e  1.30M/s (  78.2M/min) ext  10.8GB wr  10.5GB out  1.48GB  CpG 71.615% (22,263,588 / 31,087,846)
taps_partial       n=1 proc  5.11M/s write   4.30M/s e2e  2.33M/s ( 140.0M/min) ext  10.7GB wr  10.5GB out  1.49GB  CpG 71.585% (22,266,336 / 31,104,836)
htslib             n=1 proc  1.83M/s write   4.26M/s e2e  1.28M/s (  76.7M/min) ext  10.8GB wr  10.4GB out  1.47GB  CpG 70.859% (21,499,395 / 30,341,182)
nobaq              n=1 proc  6.17M/s write   4.35M/s e2e  2.55M/s ( 153.1M/min) ext  10.7GB wr  10.5GB out  1.49GB  CpG 71.543% (22,297,279 / 31,166,253)
cg                 n=2 proc  4.29M/s write  89.62M/s e2e  4.09M/s ( 245.4M/min) ext   2.1GB wr   1.3GB out  0.15GB  CpG 71.615% (22,263,688 / 31,087,996)
cg_nobaq           n=1 proc  9.34M/s write  95.47M/s e2e  8.51M/s ( 510.7M/min) ext   2.0GB wr   1.4GB out  0.15GB  CpG 71.543% (22,297,279 / 31,166,253)
ch                 n=1 proc  3.37M/s write   4.82M/s e2e  1.98M/s ( 119.0M/min) ext  10.5GB wr  10.0GB out  1.34GB  CpG not extracted (--context CH)
comp_lzf           n=1 proc  3.05M/s write   4.07M/s e2e  1.74M/s ( 104.6M/min) ext  10.7GB wr  10.4GB out  4.01GB  CpG 71.615% (22,263,688 / 31,087,996)
comp_gzip          n=1 proc  3.32M/s write   3.75M/s e2e  1.76M/s ( 105.6M/min) ext  10.7GB wr  10.4GB out  2.06GB  CpG 71.615% (22,263,688 / 31,087,996)
comp_gzip-shuffle  n=1 proc  3.47M/s write   2.90M/s e2e  1.58M/s (  94.7M/min) ext  10.8GB wr  10.5GB out  1.23GB  CpG 71.615% (22,263,688 / 31,087,996)
comp_gzip6         n=1 proc  3.13M/s write   2.33M/s e2e  1.34M/s (  80.1M/min) ext  10.8GB wr  10.4GB out  1.88GB  CpG 71.615% (22,263,688 / 31,087,996)
comp_blosc-zstd    n=1 proc  3.21M/s write   2.14M/s e2e  1.29M/s (  77.2M/min) ext  10.8GB wr  10.5GB out  1.06GB  CpG 71.615% (22,263,688 / 31,087,996)
comp_none          n=1 proc  3.47M/s write   1.55M/s e2e  1.07M/s (  64.2M/min) ext  10.8GB wr  10.5GB out 34.40GB  CpG 71.615% (22,263,688 / 31,087,996)
```

Per-run spread, default (n=5): end to end 1.78, 1.80, 1.85, 1.88, 1.93M reads/s. `--context CG` (n=2): 4.06, 4.12M reads/s.

- Mask in kernel vs the previous build (separate binary, same day): output identical for every cell in `taps`, `taps-full` and `htslib` modes; extraction peak 13.7 → 11.0 GB (default), and `--context CG` extraction 7.0 → 2.1 GB, stream mode 7.3 → 3.0 GB.
- Worker time at the default: BAQ 37%, calling 36%, read/inflate 14%, pileup 12%; 28.2% of reads run the HMM.
- BAQ stage, AVX2 vs `--baq-simd off`: 180 vs 325 thread-seconds (81% faster).
- One shard vs 16 with `blosc`: write 3.72M vs 4.29M reads/s (with `lzf` on 2026-10-16: 1.8M vs 4.1M).
