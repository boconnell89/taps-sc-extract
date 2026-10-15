# Frozen Rust chr19 reference (barcode fix + htslib 1.21 overlap clipping)

Made 2026-10-15 on branch `rust-upgrade-2026-10-15`, after hts-sys 2.2.1 was
patched with htslib 1.21's `tweak_overlap_quality` (see
`taps-sc-extract-rs/patches/`).

    python3 -m taps_sc_extract -b 260731/taps_ucbTn5_sp2.srt.bam -f /mnt/e/refs/mm10/mm10.fa \
        -o test_chr19_rust_overlapfix_2026-10-15/chr19.h5 -c chr19 -t 24 --engine rust --shards 1

BAQ on (default), trims R1=0 / R2=10, LZF, auto window size.

Against `test_chr19_rust_cellfix_2026-10-15/chr19.h5` (same code, htslib 1.19.1
overlap clipping), 18 CG and 167 CH calls are gone: bases where the mate has a
deletion inside the read-pair overlap now get their quality lowered like other
overlapping bases. With `--no-baq`, Rust now matches the Python engine exactly on
chr19, per cell and position.
