# Frozen Rust chr19 reference (default BAQ `--baq taps`, 2026-10-16)

    python3 -m taps_sc_extract -b 260731/taps_ucbTn5_sp2.srt.bam -f /mnt/e/refs/mm10/mm10.fa \
        -o test_chr19_rust_tapsskip_2026-10-16/chr19.h5 -c chr19 -t 24 --engine rust --shards 1

TAPS-masked BAQ (C masked to N for OT reads, G for OB reads), skipped for reads
whose CIGAR is all M and whose bases equal the masked reference everywhere it is
A/C/G/T. Trims R1=0 / R2=10, LZF, auto window size.

Against `test_chr19_rust_tapsbaq_2026-10-16/` (`--baq taps-full`, HMM on every
read): 4 of 806,646 CG and 964 of 18.8M CH cell-sites differ. Full test BAM
(53.5M reads, `-t 32 --shards 16`): 172 of 25.9M CG cell-sites differ, CpG
71.615% in both, processing 19–23 s vs 41–45 s.
