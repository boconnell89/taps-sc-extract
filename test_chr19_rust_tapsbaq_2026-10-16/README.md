# Frozen Rust chr19 reference (TAPS-aware BAQ, the default from 2026-10-16)

    python3 -m taps_sc_extract -b 260731/taps_ucbTn5_sp2.srt.bam -f /mnt/e/refs/mm10/mm10.fa \
        -o test_chr19_rust_tapsbaq_2026-10-16/chr19.h5 -c chr19 -t 24 --engine rust --shards 1

BAQ on with `--baq taps` (default): OT reads are realigned against the reference
with every C masked to N, OB reads with every G masked. Trims R1=0 / R2=10, LZF,
auto window size. Same code otherwise as `test_chr19_rust_overlapfix_2026-10-15/`,
which `--baq htslib` still reproduces exactly.

Calls against `--no-baq` (chr19):

| BAQ    | CG methylated | CG unmethylated | CG %   | CH methylated | CH unmethylated | CH %   |
|--------|---------------|-----------------|--------|---------------|-----------------|--------|
| off    | 574,228       | 281,919         | 67.071 | 59,592        | 19,984,017      | 0.297  |
| htslib | -3.00%        | -0.37%          | 66.479 | -4.93%        | -0.42%          | 0.284  |
| taps   | -0.07%        | -0.84%          | 67.242 | -0.60%        | -0.38%          | 0.297  |

htslib BAQ treats each methylated C (read as T) as a mismatch and removes
methylated calls preferentially; masking removes that bias.
