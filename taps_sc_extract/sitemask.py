"""
Per-position call filters (same rules as the Rust ``sitemask.rs``): positions
covered by more kept-cell reads than ``--site-depth-limit``, and ``--blacklist``
BED regions. The pileup is never changed; a masked position loses its calls in
every cell.

Depth at a position = reads of kept cells that pass the pileup's read filters
(flag, MAPQ, orphans) and whose reference span covers it, deletions included.
"""

import gzip
import heapq
from typing import Dict, Iterable, List, Optional, Set, Tuple

Interval = Tuple[int, int]  # 0-based, half-open

MOLECULES_PER_CELL = 4  # 2 strands x 2 haplotypes (8 for a replicated, G2 genome)


def auto_site_depth_limit(cells: int, molecules_per_cell: int, reads_per_molecule: int) -> int:
    """molecules per cell x reads per molecule x kept cells (Rust ``auto_site_depth_limit``)."""
    return max(1, min(cells * molecules_per_cell * reads_per_molecule, 2**32 - 1))


def cap_counts(t: int, c: int, cap: int) -> Tuple[int, int]:
    """Observations one cell keeps at a site under ``cap``: both counts scaled
    to ``cap`` keeping their ratio, methylated rounded half up (Rust
    ``sitemask::cap_counts``)."""
    n = t + c
    if cap <= 0 or n <= cap:
        return t, c
    c2 = (2 * cap * c + n) // (2 * n)
    return cap - c2, c2


def auto_cell_site_cap(molecules_per_cell: int, reads_per_call: int) -> int:
    """A site is called on one strand: molecules_per_cell / 2 molecules (one per
    haplotype), x calls per molecule (2 only for TAPS with --no-overlap-clip)."""
    return max(1, molecules_per_cell // 2) * max(1, reads_per_call)


def apply_cap(accum, touched, pos_1based: int, cap: int, stats, ctx: str, ctx_tri: str) -> None:
    """Cap each cell touched in this column (pileup columns are complete per position)."""
    for bc in touched:
        counts = accum[bc][pos_1based]
        t, c = cap_counts(counts[0], counts[1], cap)
        if (t, c) != (counts[0], counts[1]):
            for k, key in ((0, "t"), (1, "c")):
                d = counts[k] - (t, c)[k]
                stats[ctx_tri][key] -= d
                stats[ctx][key] -= d
            counts[0], counts[1] = t, c


def merge(intervals: Iterable[Interval]) -> List[Interval]:
    out: List[Interval] = []
    for s, e in sorted(intervals):
        if out and out[-1][1] >= s:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def depth_mask(spans: Iterable[Interval], limit: int) -> List[Interval]:
    """Intervals where more than ``limit`` spans overlap. ``spans`` come in
    start order (as a coordinate-sorted BAM gives them)."""
    ends: List[int] = []
    cur = 0
    masked: List[Interval] = []

    def mark(lo: int, hi: int) -> None:
        if lo >= hi or len(ends) <= limit:
            return
        if masked and masked[-1][1] >= lo:
            masked[-1] = (masked[-1][0], max(masked[-1][1], hi))
        else:
            masked.append((lo, hi))

    def advance(to: int) -> None:
        nonlocal cur
        while ends and ends[0] <= to:
            e = ends[0]
            mark(cur, e)
            cur = max(cur, e)
            heapq.heappop(ends)
        mark(cur, to)
        cur = max(cur, to)

    for s, e in spans:
        if e <= s:
            continue
        advance(s)
        heapq.heappush(ends, e)
    if ends:
        advance(max(ends))
    return masked


def window_depth_mask(bams, contig: str, start: int, end: int, limit: int, flag_filter: int, min_mapq: int,
                      ignore_orphans: bool, whitelist: Optional[Set[str]]) -> List[Interval]:
    """``depth_mask`` over the reads the pileup of ``[start, end)`` would use,
    over every input file together (``bams``: one AlignmentFile or a list)."""
    if not isinstance(bams, (list, tuple)):
        bams = [bams]

    def spans_of(bam):
        for r in bam.fetch(contig, start, end):
            if r.flag & flag_filter or r.mapping_quality < min_mapq:
                continue
            if ignore_orphans and (r.flag & 0x1) and not (r.flag & 0x2):
                continue
            if whitelist is not None:
                q = r.query_name
                i = q.find(":")
                if (q[:i] if i != -1 else q) not in whitelist:
                    continue
            if r.reference_end is not None:
                yield r.reference_start, r.reference_end
    return depth_mask(heapq.merge(*(spans_of(b) for b in bams)), limit)


def load_blacklist(path: str) -> Dict[str, List[Interval]]:
    """BED (contig, start, end; more columns ignored), merged per contig."""
    by: Dict[str, List[Interval]] = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{path} line {i}: expected contig, start, end")
            s, e = int(fields[1]), int(fields[2])
            if e > s:
                by.setdefault(fields[0], []).append((s, e))
    return {c: merge(v) for c, v in by.items()}


def overlapping(regions: List[Interval], start: int, end: int) -> List[Interval]:
    return [(s, e) for s, e in regions if s < end and e > start]


class Cursor:
    """``hit(pos)`` for increasing 0-based positions against merged intervals."""

    def __init__(self, intervals: List[Interval]):
        self.iv = intervals
        self.k = 0

    def hit(self, pos: int) -> bool:
        iv = self.iv
        while self.k < len(iv) and iv[self.k][1] <= pos:
            self.k += 1
        return self.k < len(iv) and iv[self.k][0] <= pos
