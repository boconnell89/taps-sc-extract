"""
Methylation calling logic and SAM flag classification.

Implements astair-compatible flag classification and TAPS (mCtoT) methylation calling.
"""

from typing import Optional

# Exact SAM flag sets matching astair's paired-end directional rules:
# 99:  read paired, read mapped in proper pair, mate reverse strand, first in pair (R1 forward)
# 147: read paired, read mapped in proper pair, read reverse strand, second in pair (R2 reverse)
OT_FLAGS = {99, 147}

# 83:  read paired, read mapped in proper pair, read reverse strand, first in pair (R1 reverse)
# 163: read paired, read mapped in proper pair, mate reverse strand, second in pair (R2 forward)
OB_FLAGS = {83, 163}

# Direct mapping for O(1) strand classification
FLAG_STRAND_MAP = {
    99: 'OT',
    147: 'OT',
    83: 'OB',
    163: 'OB',
}

# SAM flag bits used to locate the sequenced 5' end of R1/R2.
BAM_FREVERSE = 16
BAM_FREAD2 = 128

# Default 5' sequenced-base trims. R2 5' shows strong Tn5/TAPS M-bias;
# R1 is typically clean.
DEFAULT_TRIM_R1 = 0
DEFAULT_TRIM_R2 = 10

# Direct lookup table for mCtoT calling: (strand, ref_base, read_base) -> 0 (unmeth) / 1 (meth)
MCTOT_LOOKUP = {
    ('OT', 'C', 'T'): 1,
    ('OT', 'C', 't'): 1,
    ('OT', 'C', 'C'): 0,
    ('OT', 'C', 'c'): 0,
    ('OB', 'G', 'A'): 1,
    ('OB', 'G', 'a'): 1,
    ('OB', 'G', 'G'): 0,
    ('OB', 'G', 'g'): 0,
}


def classify_strand(flag: int) -> Optional[str]:
    """
    Classify a SAM alignment flag into Original Top ('OT'), Original Bottom ('OB'), or None.

    Non-primary, supplementary, unmapped, orphan, or improper-pair reads
    will not match these flags and return None.
    """
    return FLAG_STRAND_MAP.get(flag)


def keep_after_5prime_trim(
    qpos: int,
    query_length: int,
    flag: int,
    trim_r1: int = DEFAULT_TRIM_R1,
    trim_r2: int = DEFAULT_TRIM_R2,
) -> bool:
    """
    Keep this pileup base if it is past the 5' sequenced-base trim.

    ``qpos`` is the 0-based index into BAM SEQ (reference-oriented). Reverse-
    strand records store SEQ reverse-complemented, so the original sequenced
    5' end is at the high end of SEQ. R2 uses ``trim_r2`` (default 10);
    everything else uses ``trim_r1`` (default 0).
    """
    trim = trim_r2 if (flag & BAM_FREAD2) else trim_r1
    if trim <= 0:
        return True
    if trim >= query_length:
        return False
    if flag & BAM_FREVERSE:
        return qpos < query_length - trim
    return qpos >= trim


def call_mctot(strand: str, ref_base: str, read_base: str) -> Optional[int]:
    """
    Evaluate TAPS (mCtoT) methylation status for a single base call.

    Rules:
    - OT strand (ref 'C'):
        - read 'T' -> 1 (methylated: C was modified to T)
        - read 'C' -> 0 (unmethylated: C was unmodified)
    - OB strand (ref 'G'):
        - read 'A' -> 1 (methylated: G was modified on opposite strand, read as A)
        - read 'G' -> 0 (unmethylated: G was unmodified)

    Returns:
        1 for methylated, 0 for unmethylated, None for non-informative/mismatched bases.
    """
    return MCTOT_LOOKUP.get((strand, ref_base.upper(), read_base))
