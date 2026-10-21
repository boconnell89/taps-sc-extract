"""
Bisulfite / EM-seq calling (unmethylated C reads as T).

Kept apart from the TAPS rules in ``calling.py`` (the inverse); ``chemistry.py``
selects one set per run.
"""

from typing import Optional

# Single-end directional libraries (BSBolt W_C2T / C_C2T): primary forward reads
# are the original top strand, primary reverse reads the original bottom strand.
# Paired, secondary and supplementary flags are not classified (supplementary
# reads are filtered out of the pileup for bisulfite).
BISULFITE_SE_STRAND_MAP = {
    0: 'OT',
    16: 'OB',
}

# (strand, ref_base, read_base) -> 1 methylated / 0 unmethylated
BISULFITE_LOOKUP = {
    ('OT', 'C', 'C'): 1,
    ('OT', 'C', 'c'): 1,
    ('OT', 'C', 'T'): 0,
    ('OT', 'C', 't'): 0,
    ('OB', 'G', 'G'): 1,
    ('OB', 'G', 'g'): 1,
    ('OB', 'G', 'A'): 0,
    ('OB', 'G', 'a'): 0,
}


def classify_strand(flag: int) -> Optional[str]:
    return BISULFITE_SE_STRAND_MAP.get(flag)


def call_bisulfite(strand: str, ref_base: str, read_base: str) -> Optional[int]:
    """OT ref C: read C = 1 (methylated), T = 0. OB ref G: read G = 1, A = 0."""
    return BISULFITE_LOOKUP.get((strand, ref_base.upper(), read_base))
