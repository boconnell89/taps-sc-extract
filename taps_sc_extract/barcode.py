"""
Cell barcode extraction and whitelist/annot file parsing.
"""

import gzip
from typing import Set


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
