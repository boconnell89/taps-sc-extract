# taps-sc-extract documentation

Start with the [README](../README.md) for installation and a quickstart. These pages hold the detail.

| Page | What it covers |
|---|---|
| [Running it: examples by use case](modes.md) | Copy-paste commands for TAPS, bisulfite, several input BAMs and sample sheets, cell filters, CG-only output, memory and storage choices, CRAM/SAM input, the Python engine, GPU BAQ, tuning a new machine and profiling |
| [Command-line reference](usage.md) | Every option of `taps-sc-extract` with its default, and the `taps-sc-extract-rs` subcommands |
| [What the settings do](settings.md) | How filters, trims, BAQ, depth limits, per-cell caps, cell filters and contexts change the calls, with measured effects, and when to change the defaults |
| [Chemistry and calling](chemistry.md) | TAPS and bisulfite call rules, chemistry detection and defaults, BAQ for TAPS |
| [Output files](output.md) | HDF5 layout and dtype, run metadata, sharded `master.h5`, reading from Amethyst (R) and Python |
| [Performance](performance.md) | Measured speed and memory, what each performance option does, GPU BAQ, performance history |
| [Troubleshooting](troubleshooting.md) | Common errors and how to fix them |
| [Development and testing](development.md) | Building, tests, and the checks that keep the engines identical |
