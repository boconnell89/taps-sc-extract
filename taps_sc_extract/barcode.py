"""
Cell barcode extraction and whitelist/annot file parsing.
"""

import gzip
import logging
from typing import Optional, Set, Tuple

logger = logging.getLogger(__name__)


def extract_barcode(query_name: str) -> str:
    """
    Extract cell barcode from read QNAME.

    Splits on the first ':' (e.g. '<barcode>:<read_id>#0' -> '<barcode>').
    Falls back to splitting on '/' if ':' is not present.
    """
    if ':' in query_name:
        return query_name.split(':', 1)[0]
    elif '/' in query_name:
        return query_name.split('/', 1)[0]
    return query_name


def parse_annot(path: str) -> Set[str]:
    """
    Parse a barcode whitelist or annotation file with flexible column structure.

    The first column is always taken as the barcode.
    Handles plain text and gzipped files.
    """
    barcodes = set()
    opener = gzip.open if path.endswith('.gz') else open

    with opener(path, 'rt', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split()
            if fields:
                barcodes.add(fields[0])

    return barcodes


def parse_complexity(
    path: str,
    min_complexity: Optional[float] = None,
    max_complexity: Optional[float] = None,
) -> Tuple[Set[str], int]:
    """
    Barcodes from a complexity file whose unique-read count is within
    ``[min_complexity, max_complexity]`` (inclusive; ``None`` = unbounded).

    Format (whitespace-separated): ``rank barcode total_reads unique_reads
    percent_unique``. Column 2 is the barcode, column 4 the unique-read count.
    Empty lines, ``#`` comments and a header line whose 4th column is not a
    number are skipped. Returns (passing barcodes, barcodes in the file).
    Same rules as the Rust ``parse_complexity``.
    """
    passing: Set[str] = set()
    total = 0
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(
                    f"{path} line {i}: expected at least 4 columns "
                    f"(rank, barcode, total, unique), got {len(fields)}"
                )
            try:
                unique = float(fields[3])
            except ValueError:
                if total == 0:
                    continue  # header
                raise ValueError(f"{path} line {i}: unique-read count {fields[3]!r} is not a number")
            total += 1
            if min_complexity is not None and unique < min_complexity:
                continue
            if max_complexity is not None and unique > max_complexity:
                continue
            passing.add(fields[1])
    return passing, total


def load_cell_filter(
    whitelist_path: Optional[str] = None,
    complexity_path: Optional[str] = None,
    min_complexity: Optional[float] = None,
    max_complexity: Optional[float] = None,
) -> Optional[Set[str]]:
    """
    Barcodes to keep from a whitelist and/or a complexity file with min/max
    unique-read cutoffs (intersected when both are given). ``None`` keeps all.
    """
    if complexity_path is None and (min_complexity is not None or max_complexity is not None):
        raise ValueError("--min-complexity / --max-complexity need --complexity FILE")
    whitelist: Optional[Set[str]] = None
    if whitelist_path:
        whitelist = parse_annot(whitelist_path)
        logger.info(f"Whitelist: {len(whitelist):,} barcodes from {whitelist_path}")
    if not complexity_path:
        return whitelist
    passing, total = parse_complexity(complexity_path, min_complexity, max_complexity)
    fmt = lambda v: "-" if v is None else f"{v:g}"
    logger.info(
        f"Complexity filter: {len(passing):,} of {total:,} barcodes with unique reads in "
        f"[{fmt(min_complexity)}, {fmt(max_complexity)}] ({complexity_path})"
    )
    if whitelist is None:
        return passing
    both = passing & whitelist
    logger.info(f"Whitelist and complexity filter together: {len(both):,} barcodes")
    return both
