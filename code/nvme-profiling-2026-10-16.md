# Workstation snapshot 2026-10-16 (commit a733a66)

AMD Ryzen 9 5900XT (16 cores / 32 threads, 3.3 GHz base, up to 4.8 GHz boost, AVX2), ~58 GB RAM in WSL2, ext4 VHDX on NVMe. `5base_ucbTn5_sp2.srt.bam` (53.5M mapped reads, 7,351 cells), mm10. Rust engine built by `./install.sh`; BAQ kernel AVX2. Default BAQ is `--baq taps` (TAPS-masked, skips reads that already match). One run per row unless numbered; LZF unless noted; `--profile-json` on.

Columns: processing s, HDF5 write s, total wall s, output size, extraction-process peak RSS (MB), writer-processes peak RSS (MB).

```
default_t32_s16_mem_r1 proc   16.1 write   12.2 total   28.3 out   4.01 GB extractMB   13357 writersMB   10670
default_t32_s16_mem_r2 proc   16.5 write   14.0 total   30.5 out   4.01 GB extractMB   13452 writersMB   10625
default_t32_s16_mem_r3 proc   17.5 write   12.9 total   30.4 out   4.01 GB extractMB   13168 writersMB   10664
default_t32_s16_stream proc   17.5 write   11.2 total   28.8 out   3.99 GB extractMB    7318 writersMB   10686
default_t16_s16_mem    proc   23.2 write   13.1 total   36.3 out   4.02 GB extractMB   11796 writersMB   10416
default_t32_s1_mem     proc   15.8 write   29.0 total   44.8 out   4.04 GB extractMB   12732 writersMB      - (one shard is written in-process)
simd_off_t32_s16       proc   20.4 write   12.6 total   33.0 out   4.01 GB extractMB   13255 writersMB   10684
taps_full_t32_s16      proc   28.3 write   13.4 total   41.7 out   4.01 GB extractMB   12742 writersMB   10663
taps_partial_t32_s16   proc   11.1 write   12.5 total   23.6 out   4.01 GB extractMB   13500 writersMB   10698
nobaq_t32_s16          proc    8.1 write   13.7 total   21.8 out   4.02 GB extractMB   10870 writersMB   10681
comp_none              proc   18.1 write   46.3 total   64.4 out  34.40 GB extractMB   13539 writersMB   10642
comp_lzf               proc   22.2 write   12.9 total   35.1 out   4.01 GB extractMB   13002 writersMB   10659
comp_blosc             proc   16.8 write   12.0 total   28.8 out   1.49 GB extractMB   12975 writersMB   10673
comp_gzip              proc   16.8 write   14.6 total   31.4 out   2.06 GB extractMB   13469 writersMB   10694
comp_gzip-shuffle      proc   16.5 write   18.1 total   34.6 out   1.23 GB extractMB   13305 writersMB   10667
comp_gzip6             proc   15.8 write   22.3 total   38.1 out   1.88 GB extractMB   13078 writersMB   10660
comp_blosc-zstd        proc   17.0 write   24.8 total   41.9 out   1.06 GB extractMB   12706 writersMB   10739
```

- Best case: `default_t32_s16_mem` (three runs) 28.3–30.5 s wall. Stream mode is as fast, with 7.3 GB extraction peak.
- 32 vs 16 workers: processing 1.4× faster. 16 vs 1 shard: write 2.2× faster, wall 1.5×.
- SIMD (AVX2) vs `--baq-simd off`: processing 1.2× faster at the default BAQ, 1.5–1.6× with `taps-full`.
- Compression: processing is codec-independent (the `comp_lzf` 22.2 s is run-to-run noise); `blosc` writes as fast as `lzf` with 2.7× smaller output; `none` is disk-bound (34 GB).
- Whole-dataset `baq-selftest`: 53,675,651 reads on 60 contigs identical to htslib for the 4-lane and 8-lane (on AVX2) kernels.
- Worker time split at the default: BAQ 38%, calling 36%, read/inflate 14%, pileup 13%. Calling is the next target.
- Extraction peak is ~0.8 GB above the pre-SIMD build (12.5 → 13.3 GB); not yet investigated.
