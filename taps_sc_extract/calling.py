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
