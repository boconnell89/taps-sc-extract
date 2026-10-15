# Frozen Rust chr19 reference (barcode-cache fix)

Made 2026-10-15 from the working tree on top of `026c52e`, after the per-read
barcode cache was keyed on `bam1_t.id` instead of the QNAME pointer.

    python3 -m taps_sc_extract -b 260731/taps_ucbTn5_sp2.srt.bam -f /mnt/e/refs/mm10/mm10.fa \
        -o test_chr19_rust_cellfix_2026-10-15/chr19.h5 -c chr19 -t 24 --engine rust --shards 1

BAQ on (default), trims R1=0 / R2=10, LZF, auto window size.

Against `test_chr19_rust_026c52e/chr19.h5`, every per-position total (t, c summed
over cells) is identical. 79 CG calls (0.009%) and 2,914 CH calls (0.015%) moved
to the correct barcode. The old QNAME-pointer cache could hand a read the
previous read's barcode once htslib recycled that read's pileup node.
