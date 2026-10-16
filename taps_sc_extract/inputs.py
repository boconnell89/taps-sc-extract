"""
Alignment input: BAM or CRAM with an index, built on demand.

Same rules as the Rust engine (``taps-sc-extract-rs/src/input.rs``):

- plain SAM (``.sam``, ``.sam.gz``) is copied to a temporary BAM; it must be
  coordinate-sorted, because windowed extraction needs an index;
- a missing ``.bai`` / ``.csi`` / ``.crai`` is built next to the file, or in a
  temporary directory when that location is not writable;
- CRAM decodes against the ``-f`` FASTA.
"""

import logging
import os
import tempfile
import time
from dataclasses import dataclass
from typing import Optional, Tuple

import pysam

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class AlignmentSource:
    """Picklable description of how to open the alignment file (passed to workers)."""

    path: str
    index: Optional[str] = None
    reference: Optional[str] = None

    def open(self, threads: int = 1) -> pysam.AlignmentFile:
        mode = "rc" if self.path.lower().endswith(".cram") else "rb"
        return pysam.AlignmentFile(
            self.path,
            mode,
            index_filename=self.index,
            reference_filename=self.reference,
            threads=threads,
        )


def _is_sam(path: str) -> bool:
    p = path.lower()
    return p.endswith(".sam") or p.endswith(".sam.gz")


def _sam_to_bam(sam: str, out: str) -> None:
    last = (-1, -1)
    with pysam.AlignmentFile(sam, "r") as src, pysam.AlignmentFile(out, "wb", template=src) as dst:
        for rec in src:
            key = (rec.reference_id if rec.reference_id >= 0 else 2**31, rec.reference_start)
            if key < last:
                raise ValueError(
                    f"{sam} is not coordinate-sorted. Sort it first: samtools sort -o sorted.bam {sam}"
                )
            last = key
            dst.write(rec)


def _has_index(src: AlignmentSource) -> bool:
    # Quiet probe: htslib logs an error for a missing index, which is handled.
    old = pysam.set_verbosity(0)
    try:
        with src.open() as f:
            return f.has_index()
    except (OSError, ValueError):
        return False
    finally:
        pysam.set_verbosity(old)


def prepare_alignment(
    path: str, fasta_path: str
) -> Tuple[AlignmentSource, Optional[tempfile.TemporaryDirectory]]:
    """Return how to open ``path``, and a temp directory the caller must keep alive."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"alignment file {path} does not exist")
    tmp: Optional[tempfile.TemporaryDirectory] = None

    def tmpdir() -> str:
        nonlocal tmp
        if tmp is None:
            tmp = tempfile.TemporaryDirectory(prefix="taps-py-input-")
        return tmp.name

    if _is_sam(path):
        out = os.path.join(tmpdir(), "input.bam")
        logger.info(f"Converting SAM {path} to a temporary BAM in {tmpdir()} ...")
        _sam_to_bam(path, out)
        path = out
    is_cram = path.lower().endswith(".cram")
    src = AlignmentSource(path, None, fasta_path if is_cram else None)
    if _has_index(src):
        return src, tmp

    with src.open() as f:
        longest = max(f.lengths) if f.lengths else 0
    csi = longest >= (1 << 29)
    ext = "crai" if is_cram else ("csi" if csi else "bai")
    flags = ["-c"] if (csi and not is_cram) else []
    logger.info(f"No index for {path}; building .{ext} ...")
    t0 = time.time()
    try:
        pysam.index(*flags, path)
    except pysam.SamtoolsError:
        idx = os.path.join(tmpdir(), f"{os.path.basename(path)}.{ext}")
        pysam.index(*flags, path, idx)
        logger.info(f"  (directory not writable; index written to {idx})")
        src = AlignmentSource(path, idx, src.reference)
    logger.info(f"Index built in {time.time() - t0:.1f}s")
    if not _has_index(src):
        raise OSError(f"could not open {path} with its index")
    return src, tmp
