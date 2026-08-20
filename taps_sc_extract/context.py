"""
Reference cytosine context classification.

Determines the sequence context (CG vs CH / CHG / CHH) for cytosines on
both top and bottom strands from the reference genome sequence.
"""

from typing import Optional


def classify_trinucleotide_context(ref_seq: str, pos: int) -> Optional[str]:
    """
    Classify the trinucleotide context of a cytosine at 0-based index `pos` in `ref_seq`.

    Returns:
        'CpG', 'CHG', 'CHH', or None if `pos` is not a C/G or is out of bounds.
    """
    seq_len = len(ref_seq)
    if pos < 0 or pos >= seq_len:
        return None

    base = ref_seq[pos]

    if base in ('C', 'c'):
        # Top strand cytosine: look downstream (pos+1, pos+2)
        n1 = ref_seq[pos + 1].upper() if pos + 1 < seq_len else 'N'
        n2 = ref_seq[pos + 2].upper() if pos + 2 < seq_len else 'N'

        if n1 == 'G':
            return 'CpG'
        elif n2 == 'G' and n1 in ('A', 'C', 'T'):
            return 'CHG'
        else:
            return 'CHH'

    elif base in ('G', 'g'):
        # Bottom strand cytosine: look upstream (pos-1, pos-2)
        # On reverse strand: comp(ref[pos]) = C, comp(ref[pos-1]), comp(ref[pos-2])
        p1 = ref_seq[pos - 1].upper() if pos - 1 >= 0 else 'N'
        p2 = ref_seq[pos - 2].upper() if pos - 2 >= 0 else 'N'

        if p1 == 'C':
            return 'CpG'
        elif p2 == 'C' and p1 in ('A', 'G', 'T'):
            return 'CHG'
        else:
            return 'CHH'

    return None


def classify_context(ref_seq: str, pos: int) -> Optional[str]:
    """
    Classify the context of a cytosine at 0-based index `pos` into 'CG' or 'CH'.

    Returns:
        'CG' for CpG context, 'CH' for non-CpG (CHG or CHH) context, or None.
    """
    seq_len = len(ref_seq)
    if pos < 0 or pos >= seq_len:
        return None

    base = ref_seq[pos]
    if base in ('C', 'c'):
        if pos + 1 < seq_len and ref_seq[pos + 1] in ('G', 'g'):
            return 'CG'
        return 'CH'
    elif base in ('G', 'g'):
        if pos > 0 and ref_seq[pos - 1] in ('C', 'c'):
            return 'CG'
        return 'CH'

    return None
