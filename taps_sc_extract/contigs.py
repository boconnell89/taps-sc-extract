"""Default contigs: the primary chromosomes named in the BAM header.

Same rule as the Rust engine (`window::primary_contigs`): an optional
``chr``/``Chr``/``CHR`` prefix followed by a number (autosomes) or one of X, Y,
Z, W (sex chromosomes). That covers UCSC (``chr1``, ``chrX``) and
Ensembl/NCBI (``1``, ``X``) naming for mammals and birds, and leaves out
mitochondria (chrM, MT), unplaced and unlocalised scaffolds, alt haplotypes,
decoys and spike-ins. Autosomes sort numerically, then X, Y, Z, W, whatever
order the header uses.
"""
import re
from typing import Iterable, List, Optional, Tuple

_PRIMARY = re.compile(r"^(?:chr|Chr|CHR)?(?:([1-9][0-9]{0,2})|([XYZW]))$")
_SEX_ORDER = {"X": 0, "Y": 1, "Z": 2, "W": 3}


def _rank(name: str) -> Optional[Tuple[int, int]]:
    m = _PRIMARY.match(name)
    if not m:
        return None
    if m.group(1):
        return (0, int(m.group(1)))
    return (1, _SEX_ORDER[m.group(2)])


def primary_contigs(names: Iterable[str]) -> List[str]:
    """Primary chromosomes among `names`, in natural order (1, 2, ..., 10, ..., X, Y)."""
    ranked = {}
    for n in names:
        r = _rank(n)
        if r is not None and r not in ranked:
            ranked[r] = n
    return [ranked[r] for r in sorted(ranked)]


def describe_primary(contigs: List[str]) -> str:
    """e.g. '22 autosomes + chrX, chrY'."""
    auto = sum(1 for c in contigs if _rank(c) and _rank(c)[0] == 0)
    sex = [c for c in contigs if _rank(c) and _rank(c)[0] == 1]
    return f"{auto} autosomes" + (f" + {', '.join(sex)}" if sex else "")
