"""
Thread-safe, process-safe indexed FASTA reader.

Reads reference genome sequences using .fa and .fai index offsets directly
without C-level library sharing or memory-allocation issues across multiprocessing workers.
"""

from typing import Dict, Tuple


class FastFaiReader:
    """Thread-safe and process-safe indexed FASTA sequence reader."""

    def __init__(self, fa_path: str):
        self.fa_path = fa_path
        self.index: Dict[str, Tuple[int, int, int, int]] = {}
        fai_path = fa_path + ".fai"
        with open(fai_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if parts and len(parts) >= 5:
                    name = parts[0]
                    length = int(parts[1])
                    offset = int(parts[2])
                    line_bases = int(parts[3])
                    line_width = int(parts[4])
                    self.index[name] = (length, offset, line_bases, line_width)

    def get_reference_length(self, contig: str) -> int:
        """Return the sequence length for a given contig."""
        if contig not in self.index:
            raise KeyError(f"Contig '{contig}' not found in FASTA index.")
        return self.index[contig][0]

    def fetch(self, contig: str, start: int = 0, end: int = None) -> str:
        """
        Fetch a sub-sequence from `contig` for 0-based half-open interval [start, end).

        Returns uppercase ASCII DNA string.
        """
        if contig not in self.index:
            raise KeyError(f"Contig '{contig}' not found in FASTA index.")

        length, byte_offset, line_bases, line_width = self.index[contig]
        if end is None or end > length:
            end = length
        if start < 0:
            start = 0
        if start >= end:
            return ""

        start_line = start // line_bases
        start_col = start % line_bases
        start_byte = byte_offset + start_line * line_width + start_col

        end_pos = end - 1
        end_line = end_pos // line_bases
        end_col = end_pos % line_bases
        end_byte = byte_offset + end_line * line_width + end_col + 1

        with open(self.fa_path, "rb") as f:
            f.seek(start_byte)
            raw = f.read(end_byte - start_byte)

        # Strip newlines and carriage returns
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()
