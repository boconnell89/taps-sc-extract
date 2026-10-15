# taps-sc-extract Rust engine upgrade plan

Combined from the two reviews of `taps-sc-extract-rs`. This is the work order. It is not a claim that any phase is done.

## Baselines

Keep both clocks, and name which clock a number uses.

- Published mm10 job, pre-upgrade: 74.9M mapped reads, 24 workers, about 292 s wall, about 256k reads/s. That clock includes HDF5. Recorded in `code/rust-implementation-log.md`.
- Unpublished high-thread run: about 400k reads/s, 96 threads, memory mode. This run is not in the git history. Four times the threads produced about 1.6× the rate of the 24-worker wall-clock figure. The two clocks are not the same measurement. The 96-thread number is the high-thread baseline.

Every later run records processing seconds and HDF5 seconds separately, at 24 workers and at 96 threads, with BAQ on and with `--no-baq`. The post-run summary already prints those two times plus estimated read MB/s and write MB/s. Read MB/s is BAM file size times the idxstats share of mapped reads on the requested contigs, divided by processing seconds. Write MB/s is the bytes in the `.h5` file, or `master.h5` plus `shard_*.h5`, divided by write seconds. Megabytes are 1,000,000 bytes. Reads/s on the summary still uses total wall time.

## What is already in the tree

- ETA during extraction, on the Python cadence: every 10 chunks, every 15 seconds, or the last chunk.
- Summary lines for processing time, HDF5 write time, read throughput, and write throughput.
- Rust defaults: `--chunk-size-mb 100`, `--compression lzf`. The Python engine stays at 10 Mb and gzip. The Python CLI sends 100 and lzf only when those flags were omitted.
- An `.h5` filename with more than one shard writes `stem/shard_XXX.h5` and `stem/master.h5`. One shard still writes the `.h5` file.
- Public Python-versus-Rust speed figures are removed from `README.md`, the CLI help, and `AGENTS.md` until the phases below are remeasured. `code/rust-implementation-log.md` stays as the dated pre-upgrade record.

## Phase 1 — Defaults that do not change calls

`--decomp-threads 0` becomes the Rust auto value at 8 or more workers. `hts_set_threads(1)` adds a BGZF thread per worker, so `-t 24` is 48 threads. The implementation log already records 0 as the fastest setting at that width. `open_bam` already skips a zero count, which rust-htslib requires. The Python CLI defaults `--decomp-threads` to 1 and always forwards it, so the Python default changes in the same commit.

Remove the second writer clamp. `resolve_writer_threads` already allows up to 32. `assemble_hdf5_from_memory` and `assemble_hdf5_from_temp` then apply `.min(MAX_SHARD_WRITERS)` with that constant at 6, so `--max-writer-threads 16` does nothing. Drop that second `min`.

Retime 24 workers and 96 threads. Calls stay identical.

## Phase 2 — Stop the chr1 stampede

`windows.par_iter()` hands out window 0, 1, 2, so every worker starts on chr1 at a different offset. Each one seeks the BAM on its own, and `refresh_baq_cache` keeps a private full-contig copy for `sam_prob_realn`. chr1 is about 195 MB, so 24 copies are about 4.7 GB, and the HMM walks copies that will not stay in L3 together. At 96 threads that copy set is larger. The 100 Mb window default does not fix it.

Give each worker one contiguous range, balanced by BAM virtual offset from the `.bai`, so the split is equal file bytes. Walk that range in order on the thread-local `IndexedReader`. Keep one shared contig buffer. Calls do not change. Retime both thread counts before touching the inner loop.

## Phase 3 — Profile one window before rewriting it

Run `perf` or `samply` on one 10 Mb window, BAQ on and BAQ off, at the Phase 1–2 settings. The callback rewrite waits on that profile showing `sam_prob_realn` and `bam_plp_next` above the Rust per-base work.

The current split is the debug chr19 run, one thread: 45 s with BAQ, 31 s with `--no-baq`. mCG moved from about 64.37% to about 63.75%.

## Phase 4 — Gate BAQ, then decide whether to keep it

`accept_record` runs `sam_prob_realn` before the MAPQ test and the orphan test. Strand is checked later, per base, and only flags 99, 147, 83, and 163 produce a call. `FLAG_FILTER` drops unmapped, secondary, QC-fail, and duplicate reads. Supplementary reads and anything else that is not OT/OB still pay for the HMM and a pileup copy, then the bases are discarded.

Test strand, MAPQ, and orphans before `sam_prob_realn`. That preserves calls. Rejected reads skip the HMM. Retime, and diff chr19 sites against the current BAQ-on output.

`--no-baq` stays the switch that matches current Python. Python sets `compute_baq=True` and never passes a FASTA into the pileup, so Python does not run the HMM. A further skip for CIGARs with no insertion or deletion is a separate experiment. It is not bit-exact: BAQ also lowers quality beside soft clips and alignment ends. Ship it only if the chr19 site delta is acceptable.

## Phase 5 — Rewrite the callback only if the profile still shows it

Do these in one pass:

- Store the interned barcode id, strand, and trim decision in `bam_pileup_cd` when the read is pushed. htslib copies `cd` onto every column. The QNAME pointer MRU misses whenever the next pileup entry is a different read, which is the normal case above depth 1.
- In the hot loop, push `(cell_id, pos, meth)` into a `Vec`. After the window, sort by `(cell_id, pos)` and compact. That replaces the per-call `FxHashMap` insert. Sites come out sorted, so the writer can stop sorting, and the per-worker map heap shrinks.
- One context-code array per window replaces `classify_trinucleotide` on every column. `call_mctot` can stop uppercasing bases that FASTA fetch and BAM SEQ already uppercase.
- Make `drain_ready` generic instead of `&mut dyn FnMut`, so the column body can inline.

Check cell counts and CG/CH totals against Phase 4. This phase should not move methylation.

The per-window mutexes are not the scaling limit. `window_shard_payloads` sorts and clones before either lock. The memory lock only pushes finished payloads. The stats lock runs once per window. Leave that structure alone unless the profile says otherwise.

## Phase 6 — Compile the C that dominates

Release already uses `lto = true` and `codegen-units = 1`. Cargo’s release profile is already `opt-level = 3`. Add `-C target-cpu=native` for the Rust callback. Build htslib with `-O3 -march=native` as well. `sam_prob_realn` is C, and the Rust flag does not recompile it. `panic = "abort"` is optional and is a binary-size change.

## Phase 7 — A new pileup only if `bam_plp_next` is still the top

The noodles engine is the slower reference. It clones the QNAME into a `HashMap` for every aligned base, then sorts columns. If the profile still shows `bam_plp_next` on top after Phase 4, replace it with a CIGAR walk that records only OT-on-C and OB-on-G, plus a small QNAME-to-mate table for overlap clipping. Re-check overlap parity. The noodles clip is a different algorithm from `bam_plp_init_overlaps`.

## Parity, separate from reads per minute

Equal-start mates are still double-counted. On chr19 5.0–5.1 Mb with `--no-baq`, three sites disagree, all extra Rust `t`. One is a proper pair, flags 99/147, both reads starting at the same position. htslib `overlap_push` only registers a read whose mate has not arrived (`mpos >= pos`). Reproduce that pair against the htslib inside rust-htslib 0.49 and fix the clip. Do this beside Phase 4. It will not show up in the throughput table.

Also record, without blocking the speed work:

- `bam_plp_next` returns a 32-bit position. A contig past 2³¹−1 truncates. mm10 is safe.
- `blosc` uses LZ4 level 9, so that setting is the slow level.
- `.cargo/config.toml` pins `PKG_CONFIG_PATH`, `LIBCLANG_PATH`, and `HDF5_DIR` to this machine’s miniforge prefix. Change that before anyone else builds the crate.

## Leave alone

MD5 per barcode per window, the per-window stats mutex, and mimalloc with `default-features = false`. Smaller chunks will not balance load once ranges are cut by file offset.

## When a speed table can come back

Publish one table only after Phases 1–4 have been rerun on the same BAM, at 24 workers and at 96 threads, in memory mode and in the normal write path. Each row needs processing seconds, HDF5 seconds, reads/s on the processing clock, and reads/s on the wall clock. Put the 96-thread memory-mode run in that table the first time it is published. Until then the README describes the engine and the flags, and it does not quote a speedup.
