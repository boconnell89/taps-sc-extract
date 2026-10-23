# Performance

Measured speed and memory, what each performance option does, and how the defaults were chosen. [Back to the documentation index](index.md).

## How to read these numbers

- **Rates** are mapped reads per second (M = million).
- **Phases:**
  - *Processing* is extraction: BAM read, BAQ, pileup and calling.
  - *HDF5 write* is the write phase that remains after extraction.
  - *End to end* is the whole run.
- **Setup:** unless stated otherwise, rates come from the reference workstation and test dataset below, using the Rust engine.
- **Variation:** single runs vary by about ±4%, and other machines will differ. On a new machine, run `taps-sc-extract --test` first ([modes.md](modes.md#tuning-a-new-machine)).
- **Run summary:** every run prints processing and HDF5 write time, estimated read and write MB/s, and the peak memory of the extraction process and of the writer processes.

## Reference workstation

| | |
|---|---|
| CPU | AMD Ryzen 9 5900XT: 16 cores / 32 threads, AVX2 (no AVX-512) |
| Memory | ~58 GB visible to WSL2 |
| Storage | WSL2 ext4 on NVMe; input, temp and output all on it |
| TAPS test dataset | Paired-end TAPS, 53.5M mapped reads, 7,351 cells, mm10, 3.7 GB BAM |
| Build | `./install.sh` (native); BAQ kernel AVX2 |

Default settings, `-t 32 --shards 16 --memory-mode memory`, median of five runs (end to end ranged 2.60–2.95M reads/s), 2026-10-18:

| Phase | Rate |
|---|---:|
| Processing (writers receive their data alongside) | 4.48M reads/s |
| HDF5 write after extraction (1.4 GB, 16 shards) | 6.71M reads/s |
| **End to end** | **2.69M reads/s (161M reads/min)** |
| Peak memory | 6.3 GB extraction + 10.4 GB writers (stream mode: 3.4 GB + 10.4 GB, 2.50M reads/s) |

The plain command now chooses the same settings by itself: 32 workers, 16 shards by estimated output size, and memory mode with writers fed during extraction. It measured 2.79M reads/s end to end on 2026-10-21.

Faster options on the same machine:

| Options | End to end | Output |
|---|---:|---|
| `--baq taps-partial` | 3.02M reads/s | 0.06% of CG calls differ |
| `--no-baq` | 3.41M reads/s | no BAQ |
| `--context CG` | 5.40M reads/s | CpG only |
| `--context CG --no-baq` | 15.0M reads/s | CpG only, no BAQ |

### Bisulfite library

A merged, deduplicated single-end bisulfite library (91.6M reads, 304 cells), measured 2026-10-21 with the BAM in the page cache:

- **Processing:** about 6M reads/s, faster than TAPS because bisulfite runs no BAQ.
- **End to end:** 3.4M reads/s with the default settings.
- **The limit:** HDF5 writing of its 1.15 billion cell-sites, about 12.5 per read.
- **Cold cache:** reading the 10 GB BAM from disk first ran at about 190 MB/s.

## Options that affect speed

### Compression

`--compression` runs inside the writer processes, so it only changes the HDF5 write phase. One run per codec (`blosc`: median of five), 2026-10-18:

| Codec | End to end | HDF5 write | vs `blosc` | Size vs `none` | R needs |
|---|---:|---:|:---:|:---:|---|
| **`blosc`** (LZ4 level 9; Rust default) | 2.69M reads/s | 6.71M reads/s | — | 96% smaller | `rhdf5filters` |
| `lzf` | 2.67M | 6.27M | 1% slower | 89% smaller | `rhdf5filters` |
| `gzip` (deflate 1; Python default) | 2.40M | 5.71M | 11% slower | 94% smaller | plain `rhdf5` |
| `gzip-shuffle` | 2.20M | 4.35M | 18% slower | 97% smaller | plain `rhdf5` |
| `gzip6` | 1.85M | 3.09M | 31% slower | 95% smaller | plain `rhdf5` |
| `blosc-zstd` | 1.67M | 2.64M | 38% slower | **97% smaller** (smallest) | `rhdf5filters` |
| `none` | 1.17M | 1.62M | 57% slower | 0% (34 GB here) | plain `rhdf5` |

Which codec to pick:
- **`blosc`** writes as fast as `lzf` and its files are 63% smaller.
- **`gzip-shuffle`** is the smallest output that plain `rhdf5` reads without plugins.
- **`none`** is slowest because 34 GB has to reach the disk.

Python readers need `import hdf5plugin` for Blosc ([output.md](output.md#compression-and-readers)).

### Workers (`-t`)

- **What they are:** Rayon threads, each processing genomic windows. Extraction is CPU-bound.
- **Where their time goes,** with the default BAQ: about 54% BAQ, 24% calling, 12% BAM read and inflate, and 8% pileup bookkeeping.
- **Default:** every logical CPU, bounded by the memory budget, or the value saved by `--test --save-profile`.
- **SMT helps:** 32 workers process 4.48M reads/s against 2.97M with 16. End to end that is 2.69M against 2.08M reads/s.

### Decompression threads (`--decomp-threads`)

htslib threads per worker that inflate BGZF blocks ahead of the reader. Leave it unset:
- **One input:** Rust picks 0 at 8 or more workers (inflate on the worker thread), 2 at 5–7 workers and 4 at 1–4. With libdeflate, inflate is only 1–6% of worker time.
- **Several inputs:** the default is 0.
- **Python engine:** uses 1.

### Window size (`--chunk-size-mb`)

- **Rust default, auto:** about 6 windows per worker, clamped to 2–100 Mb (mm10: about 28 Mb at 16 workers, 14 Mb at 32, 4.7 Mb at 96). Enough windows keep every worker busy to the end of the run.
- **No effect on calls:** this was verified per cell at 100, 10 and 2 Mb.
- **Python engine:** 10 Mb.

### Shards (`--shards`)

- **Rust default, auto:** the larger of a cell-count rule and the estimated output sites ÷ 48M, capped at half the CPUs and at the cell count. After extraction, the writer planner may split the output into more files (next section).
  - Cell-count rule: 1 shard up to 10,000 cells, 8 up to 50,000, 16 up to 100,000, 32 above that.
  - The site rule matters for few, deep cells. On the 304-cell bisulfite library above, 16 shards wrote the output several times faster than one, because the shards are written in parallel.
- **Single file:** naming an `.h5` file with `-o` keeps 1 shard.
- **Python engine:** 1 shard.
- **A cell is never split.** So the largest cell sets a floor on write time, and more than about cells ÷ 2 shards adds nothing.
- **Diagnosing writes:** the run summary's `Shard writers:` line shows sites per shard (median and max), the slowest shard, and where writer time went.
- **Network filesystems:** with fewer cells, more shards were slower on one server, so measure there. `--test` probes the output filesystem.

### Writers (`--writer-mode`, `--max-writer-threads`, `--compression-threads`)

- **`--writer-mode process`** (the default with more than one shard): one `write-shard` child process per shard. Each streams its shard one cell at a time in bounded slices.
- **`--writer-mode thread`:** in-process writers. They share one HDF5 lock and so don't write in parallel; auto caps them at 6.
- **`--max-writer-threads`** (0 = auto): planned after extraction, when the shard sizes are known exactly. A writer holds its whole shard (measured: RSS ≈ 1.0 × payload + 0.1 GB; planned at 1.15 × + 0.25 GB). The planner runs as many writers as fit in the memory free at that moment, up to the `--test` profile's count, else the CPU count. With automatic shards it also splits the output into more files, from whole buckets, when that shortens the write (rounds of writers × largest file).
  - On the 40-file, 54.7-billion-site server run, 32 shards of ~21 GB allowed only about 22 writers at once. The planner would write 128 files of ~7 GB with 80 writers.
  - An explicit `--max-writer-threads` or `--shards` is followed as given; if the plan does not fit, the log says `WARNING: ... may run out of memory` with the count that would fit.
  - On network filesystems use 2–4 writers; on spinning disks, 1.
- **`--compression-threads`** (Blosc only, default 1): 1 to 32 threads all wrote equally fast. Most cell datasets are a single chunk too small for Blosc to split, so the parallelism comes from shards.

### Memory mode (`--memory-mode`, `--temp-dir`)

- **`memory`:** window results are piped to the writer processes. With one writer process per shard, the writers start with extraction and take each window's results as it finishes, at low CPU priority, so only the HDF5 writing itself is left afterwards. That was 12% faster end to end than handing the data over after extraction, with about half the extraction memory.
- **`stream`:** window results go to compact temp files under `--temp-dir` (default `$TMPDIR`, else `/tmp`), and each writer reads its own shard's files. Extraction memory is lowest: 3.4 GB against 6.3 GB on the test BAM, at 2.50M against 2.69M reads/s.
- **`auto`** (default): Rust estimates the output sites from the BAM index, at about 14 per mapped read.
  - It prefers memory mode with one streamed writer per shard when the whole process tree fits the memory budget.
  - Otherwise it chooses between memory and stream mode, with as many writers as fit.
  - The config block prints the plan: `Memory plan: ... est. peak ~X GB of Y GB`.
  - The budget is `--max-memory-fraction` (default 0.85 of available memory) or `--max-memory-gb`, and it covers extraction plus writers. The write phase re-reads available memory after extraction.
  - Explicit `--memory-mode memory` is followed even when the estimate exceeds the budget, with a WARNING.
- **Storage:** keep input, temp and output on fast local storage. On the reference machine, the same job ran at 0.31M reads/s from a 9p-mounted HDD against 0.56M reads/s from local NVMe (measured before the 2026-10-15 upgrades).

### Contexts (`--context`)

CH is 96% of sites, so skipping it is the largest single speed-up. Full TAPS test dataset, 2026-10-18:

| Run | Processing | HDF5 write | End to end | Extraction memory | Writer memory |
|---|---:|---:|---:|---:|---:|
| `both` (default) | 4.48M reads/s | 6.71M reads/s | 2.69M reads/s | 6.3 GB | 10.4 GB |
| `--context CG` | 5.62M | 139M | **5.40M** | 1.7 GB | 1.4 GB |
| `--context CH` | 4.60M | 7.77M | 2.89M | 6.6 GB | 9.9 GB |
| `--no-baq` | 7.92M | 5.99M | 3.41M | 9.1 GB | 10.4 GB |
| `--context CG --no-baq` | 17.0M | 128M | **15.0M** | 1.7 GB | 1.3 GB |

BAQ and the pileup still run on every read, since any read can cover a CpG. The kept context's calls are identical to a `both` run.

### BAQ modes and kernels

Speed of each BAQ mode on the full TAPS test dataset (AVX2 kernel). What each mode does to the calls is in [settings.md](settings.md#baq).

| BAQ | Reads through the HMM | Processing | End to end |
|---|---:|---:|---:|
| `taps-full` | 99% | 2.19M reads/s | 1.68M reads/s |
| `taps` (default) | 28% | 4.48M | 2.69M |
| `taps-partial` | 8.5% | 6.41M | 3.02M |
| `--no-baq` | 0% | 7.92M | 3.41M |

- **SIMD kernels** (`--baq-simd`, default `auto`): the AVX2 kernel against htslib's scalar code (`--baq-simd off`) makes processing 44% faster (3.12M → 4.48M reads/s) and end to end 23% faster.
- **AVX-512** runs 8 reads per vector. `--test --save-profile` times each kernel under full load and saves the fastest, because some CPUs lower their clock for AVX-512.

### GPU BAQ (`--baq-device gpu`)

- **How the work is split:** the dense double-precision HMM runs on the GPU. Everything TAPS-specific stays on the CPU: the mask, strand choice, skip check, BAQ_EXTEND and the quality cap. Reads longer than 192 bp, or with very wide bands, also stay on the CPU.
- **Batches:** one service thread pools reads from all workers into large launches. Each worker keeps several batches in flight (`TAPS_GPU_DEPTH`, default 4).
- **Shared with the CPU:** a batch goes to the GPU only if the GPU's queue will clear before the worker needs it back. Otherwise the worker runs it on its own CPU. So `--baq-device gpu` is never much slower than the CPU alone. `TAPS_GPU_SHARE=0` sends every batch to the GPU.
- **Exact:** results are identical to the CPU read by read. `baq-selftest --gpu -b BAM -f FA` checks every read; `install.sh` runs it when a GPU is usable.
- **Backends:**
  - OpenCL is always built. It needs a double-precision GPU and its vendor's OpenCL runtime (Intel Arc / Data Center: Intel's compute runtime; check with `clinfo | grep -i fp64`).
  - CUDA is built when `nvcc` is on the PATH, and HIP when `hipcc` is.
  - `--gpu-backend auto` tries CUDA, then HIP, then OpenCL.

When it helps depends on the GPU's double-precision rate relative to the host CPU:

| Host | BAQ | Device | Processing | End to end |
|---|---|---|---:|---:|
| 16-thread Ryzen 7 1700X + Intel Arc Pro B70 | `taps` | CPU | 1.65M reads/s | 1.18M reads/s |
| same | `taps` | **B70 + CPU (shared)** | **2.57M (+56%)** | **1.58M (+34%)** |
| same | `taps-full` | CPU | 0.72M | 0.61M |
| same | `taps-full` | **B70 + CPU (shared)** | **1.07M (+49%)** | **0.85M (+39%)** |
| Reference workstation + RTX 3060 | `taps` | CPU, AVX2 | 4.86M | 2.97M |
| same | `taps` | 3060 + CPU (shared) | 4.78M | 2.83M |
| same | `taps-full` | CPU, AVX2 | 2.27M | 1.77M |
| same | `taps-full` | **3060 + CPU (shared)** | **2.55M (+12%)** | **1.94M (+10%)** |

Means of two runs each, 2026-10-18. The RTX 3060 has weak double precision and on its own runs at half the CPU's speed. The default stays `cpu`. `--test` times BAQ on the GPU against the CPU and reports which is faster.

### Profiling a run

`--profile-json PATH` (Rust) writes:
- per-phase wall and CPU time, and utilisation;
- the per-stage split of worker time: read/inflate, BAQ, pileup and calling (1 in 16 reads is timed and scaled up);
- per-window timings;
- per-shard writer reports.

It adds about 2% to processing time. Use it to see what limits a run on a new system.

## Performance history

End-to-end and per-phase rates on the reference workstation and TAPS test dataset, `-t 32 --shards 16 --memory-mode memory`, each release's defaults otherwise. "% faster" is relative to the last version before the 2026-10-14 work (`e19e2b4`, 2026-09-25), run with that version's own CLI defaults.

![Throughput relative to the version before 2026-10-14](img/perf_history_throughput.svg)

![Peak memory relative to the version before 2026-10-14](img/perf_history_memory.svg)

Regenerate the plots with `python3 scripts/plot_performance_history.py`.

| Release | End to end | Processing | HDF5 write | Peak memory |
|---|---:|---:|---:|---:|
| Before 2026-10-14 (`e19e2b4`): htslib BAQ, 10 Mb windows, in-process writer threads, gzip | 0.35M reads/s | 0.91M | 0.57M | 33.2 GB |
| 2026-10-14: 100 Mb windows, LZF | 0.56M (59% faster) | 0.84M | 1.70M | 29.7 GB |
| 2026-10-15: balanced windows, libdeflate, shared BAQ references, one writer process per shard, native build | 0.97M (175% faster) | 1.30M | 3.78M | 21.7 GB |
| 2026-10-16: TAPS-aware BAQ (correctness fix) | 1.00M (184% faster) | 1.31M | 4.18M | 23.2 GB |
| 2026-10-16: skip BAQ for reads that match the masked reference | 1.64M (365% faster) | 2.60M | 4.49M | 23.2 GB |
| 2026-10-16: SIMD BAQ (AVX2) | 1.76M (400% faster) | 3.24M | 4.14M | 24.0 GB |
| 2026-10-16: `blosc` default (output 63% smaller) | 1.95M (453% faster) | 3.41M | 4.55M | 23.9 GB |
| 2026-10-17: TAPS mask inside the BAQ kernel | 1.85M (425% faster) | 3.29M | 4.29M | 21.2 GB |
| 2026-10-18: faster BAQ preparation | 2.07M (487% faster) | 3.84M | 4.49M | 21.7 GB |
| 2026-10-18: read-centric calling | 2.45M (595% faster) | 5.08M | 4.73M | 21.7 GB |
| 2026-10-18: writers fed during extraction | **2.69M (663% faster)** | **4.48M** | **6.71M** | **16.7 GB (50% less)** |

Notes on the table:
- **Processing dips in the last row** because processing now includes handing data to the writers, which shortens the write phase after extraction.
- **Peak memory** counts extraction plus writers.
- **BAQ changed in the third row:** the first two rows use htslib BAQ (CpG 70.86%), and later rows use TAPS-aware BAQ (71.62%).
- **The 2026-10-18 rows** are medians of three to five runs.

Dated snapshots with raw timings are in [`code/`](../code/) (`nvme-profiling-*.md`, `bisulfite-profiling-2026-10-21.md`). The change log is [`code/rust-implementation-log.md`](../code/rust-implementation-log.md).
