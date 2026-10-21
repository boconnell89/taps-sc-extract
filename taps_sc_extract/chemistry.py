"""
Conversion chemistry: which strand rule and call table a run uses.

TAPS (mC -> T, paired-end) and bisulfite / EM-seq (unmethylated C -> T,
single-end BSBolt) keep separate tables (``calling.py``, ``bisulfite.py``).
Extraction loops take the two dicts from ``tables()`` once per run and never
refer to either module's globals directly.
"""

from typing import Dict, Tuple

from .bisulfite import BISULFITE_LOOKUP, BISULFITE_SE_STRAND_MAP
from .calling import FLAG_STRAND_MAP, MCTOT_LOOKUP

CHEMISTRIES = ("taps", "bisulfite")

# pysam's default pileup flag filter (unmapped, secondary, QC fail, duplicate),
# plus supplementary for bisulfite: BSBolt supplementary records are mostly
# chimeras on the opposite converted strand.
BAM_FSUPPLEMENTARY = 0x800
DEFAULT_FLAG_FILTER = 0x704


def tables(chemistry: str) -> Tuple[Dict[int, str], Dict[tuple, int]]:
    """(flag -> strand map, (strand, ref, read) -> call table) for a chemistry."""
    if chemistry == "taps":
        return FLAG_STRAND_MAP, MCTOT_LOOKUP
    if chemistry == "bisulfite":
        return BISULFITE_SE_STRAND_MAP, BISULFITE_LOOKUP
    raise ValueError(f"unknown chemistry {chemistry!r} (use one of {', '.join(CHEMISTRIES)})")


def flag_filter(chemistry: str) -> int:
    """pysam pileup ``flag_filter``: supplementary is dropped before the pileup
    (so it never takes a --max-depth slot) for bisulfite only. TAPS is unchanged."""
    return DEFAULT_FLAG_FILTER | (BAM_FSUPPLEMENTARY if chemistry == "bisulfite" else 0)
