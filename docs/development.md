# Development and testing

Building from source, running the tests, and the checks that keep the two engines and every BAQ implementation identical. Contributors (human or AI) should also read [`AGENTS.md`](../AGENTS.md), which lists the invariants. [Back to the documentation index](index.md).

## Layout

- **`taps_sc_extract/`:** the Python package. It holds the CLI wrapper (`cli.py`) and the reference engine (`parallel_extractor.py`, `extractor.py`). The reference engine is the ground truth for calls and must never be removed.
- **`taps-sc-extract-rs/`:** the Rust engine. The wrapper runs it when it is built (`--engine auto`).
- **`tests/`:** pytest suite and the synthetic bisulfite BAM generator (`bsbolt_fixture.py`).
- **`scripts/`:** `xb_audit.py` (pileup vs BSBolt `XB`), `benchmark_ab.py`, `plot_performance_history.py`.
- **`code/`:** design plans, dated profiling snapshots and the Rust change log.

## Building

```bash
./install.sh               # native build + pip install -e .
./install.sh --no-pip      # rebuild the Rust core only
./install.sh --portable    # x86-64-v3 build that runs on other CPUs
```

What `install.sh` does:
- Downloads hts-sys 2.2.1 into `taps-sc-extract-rs/vendor/`, applies the patches in `taps-sc-extract-rs/patches/` (htslib 1.21's read-pair overlap clipping, to match pysam), and re-patches whenever a patch changes.
- Builds with `-ffp-contract=off`, so BAQ arithmetic matches htslib exactly.
- Runs `baq-selftest`, plus `baq-selftest --gpu` when a GPU is usable, and fails on any difference. GPU backends (OpenCL always; CUDA with `nvcc`, HIP with `hipcc`) are described in [performance.md](performance.md#gpu-baq---baq-device-gpu).

Rerun it after every pull.

### Upgrading

```bash
git pull
./install.sh            # same options as the first install (--prefix, --portable, --no-pip)
```

Rerunning is safe. It regenerates `taps-sc-extract-rs/.cargo/config.toml`, rebuilds from clean, reinstalls the CLI and re-patches the vendored hts-sys when its patches changed. The downloaded crate is kept in `vendor/`, so re-patching needs no network. The first install needs `curl`, `patch`, `sha256sum` and access to crates.io.

Checkouts older than 2026-10-15:
- **Config file:** `.cargo/config.toml` is no longer tracked. If you edited it, `git pull` refuses. Discard the edit with `git checkout -- taps-sc-extract-rs/.cargo/config.toml`, pull, and pass a non-default conda location as `./install.sh --prefix /path/to/env`.
- **Re-extract old outputs:** Rust outputs written before that date put about 0.01% of calls under the wrong barcode.

After installing on a new machine, run `taps-sc-extract --test --save-profile`.

## Tests

```bash
pytest tests/ -v                                  # Python and end-to-end tests (both engines)
(cd taps-sc-extract-rs && cargo test --release)   # Rust unit tests
```

- **Synthetic data:** the bisulfite tests build a BSBolt-like BAM at test time. It has known methylation, indels, soft clips, low MAPQ, supplementary alignments and a 3′ bias. Each engine's calls are compared per cell and site against the generator's truth, and against each other.
- **chr19 tests:** these run when the chr19 test data is present (paths in `tests/test_taps_sc_extract.py`). They check:
  - the Rust engine bit for bit against frozen references (one per BAQ mode, made with `--max-depth 250 --cell-site-cap 0`);
  - that Rust `--no-baq` equals the Python engine at every cell and position.
- **Several inputs:** tests split the synthetic BAM by cell into three BAMs and require output byte-identical to the single-BAM run.

## Checks after changing calling code

- **Compare per cell and site**, never by totals or by shard file. Shard placement differs between engines and between runs.
- **Pileup or callers:** `test_read_caller_matches_pileup`, and `stats --dump-sites` diffs of both callers on chr19.
- **BAQ:** `taps-sc-extract-rs baq-selftest -b BAM -f FA`, which compares base qualities read by read.
- **Bisulfite filters or trims:** `scripts/xb_audit.py` on a real BAM. With all filters on, the pileup and BSBolt's `XB` must agree except at BSBolt's soft-clip bug.
- **Large outputs:** stream one cell at a time. A whole run's CH data does not fit in memory as Python objects.
