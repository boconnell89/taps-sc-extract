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


# ---------------------------------------------------------------------------
# Several input BAMs (same rules as taps-sc-extract-rs/src/inputs.rs)

SHEET_COLUMNS = ("bam", "complexity", "min_complexity", "max_complexity", "whitelist")
BARCODE_SAMPLE = 5_000


@dataclass(frozen=True)
class RowFilter:
    """Per-file cell filter from a sample sheet row (None = the global flag)."""

    complexity: Optional[str] = None
    min_complexity: Optional[float] = None
    max_complexity: Optional[float] = None
    whitelist: Optional[str] = None


@dataclass
class InputSpec:
    bams: list
    filters: list
    list_path: Optional[str] = None

    @property
    def name_path(self) -> str:
        """Names outputs after the input: the list or sheet, else the first BAM."""
        return self.list_path or self.bams[0]


def _resolve(base: str, p: str) -> str:
    return p if os.path.isabs(p) else os.path.join(base, p)


def parse_bam_args(args) -> InputSpec:
    """``-b`` values: paths, or one ``@FILE`` naming a list (one path per line)
    or a sheet (TSV whose header has a ``bam`` column)."""
    args = [str(a) for a in args]
    at = [a for a in args if a.startswith("@")]
    if not at:
        if not args:
            raise ValueError("no input BAM (-b)")
        return InputSpec(list(args), [RowFilter()] * len(args))
    if len(args) > 1:
        raise ValueError("-b @FILE takes one list or sheet; do not mix it with other -b paths")
    list_path = at[0][1:]
    base = os.path.dirname(list_path)
    with open(list_path) as f:
        lines = [(i, l.rstrip("\r\n")) for i, l in enumerate(f, 1)]
    lines = [(i, l) for i, l in lines if l.strip() and not l.lstrip().startswith("#")]
    if not lines:
        raise ValueError(f"{list_path}: no input BAMs listed")
    header = [c.strip().lower() for c in lines[0][1].split("\t")]
    if "bam" in header:
        bams, filters = _parse_sheet(list_path, base, header, lines[1:])
    else:
        bams = [_resolve(base, l.strip()) for _, l in lines]
        filters = [RowFilter()] * len(bams)
    if not bams:
        raise ValueError(f"{list_path}: no input BAMs listed")
    seen = set()
    for b in bams:
        if b in seen:
            raise ValueError(f"{list_path}: {b} is listed twice")
        seen.add(b)
    return InputSpec(bams, filters, list_path)


def _parse_sheet(list_path, base, header, rows):
    for c in header:
        if c and c not in SHEET_COLUMNS:
            raise ValueError(f"{list_path}: unknown column {c!r} (columns: {', '.join(SHEET_COLUMNS)})")
    col = {c: i for i, c in enumerate(header) if c}
    bams, filters = [], []
    for line_no, line in rows:
        f = [x.strip() for x in line.split("\t")]

        def get(name):
            i = col.get(name)
            v = f[i] if i is not None and i < len(f) else ""
            return v or None

        def num(name):
            v = get(name)
            if v is None:
                return None
            try:
                return float(v)
            except ValueError:
                raise ValueError(f"{list_path} line {line_no}: {name} {v!r} is not a number")

        bam = get("bam")
        if bam is None:
            raise ValueError(f"{list_path} line {line_no}: empty bam")
        cx, wl = get("complexity"), get("whitelist")
        bams.append(_resolve(base, bam))
        filters.append(RowFilter(
            complexity=_resolve(base, cx) if cx else None,
            min_complexity=num("min_complexity"),
            max_complexity=num("max_complexity"),
            whitelist=_resolve(base, wl) if wl else None,
        ))
    return bams, filters


def xlsx_to_tsv(xlsx_path: str, out_dir: str) -> str:
    """First sheet of an .xlsx sample sheet as a TSV in ``out_dir``, with every
    path made absolute against the workbook's directory (the TSV lives elsewhere)."""
    import openpyxl

    base = os.path.dirname(os.path.abspath(xlsx_path))
    ws = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True).worksheets[0]
    rows = [["" if v is None else str(v).strip() for v in r] for r in ws.iter_rows(values_only=True)]
    rows = [r for r in rows if any(r)]
    if not rows:
        raise ValueError(f"{xlsx_path}: empty sheet")
    header = [c.lower() for c in rows[0]]
    if "bam" not in header:
        raise ValueError(f"{xlsx_path}: the first row needs a 'bam' column")
    path_cols = {i for i, c in enumerate(header) if c in ("bam", "complexity", "whitelist")}
    out = os.path.join(out_dir, os.path.splitext(os.path.basename(xlsx_path))[0] + ".tsv")
    with open(out, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows[1:]:
            r = [(_resolve(base, v) if (i in path_cols and v) else v) for i, v in enumerate(r)]
            # whole numbers from Excel come back as "10000.0": fine for float()
            fh.write("\t".join(r) + "\n")
    return out


def load_cell_filters(spec: InputSpec, whitelist=None, complexity=None, min_c=None, max_c=None):
    """(kept barcodes or None, per-file (kept, listed) counts, per-file listed sets).
    Same rules as Rust ``inputs::load_cell_filters``."""
    from .barcode import load_cell_filter, parse_complexity

    if len(spec.bams) == 1 and spec.filters[0] == RowFilter():
        keep = load_cell_filter(whitelist, complexity, min_c, max_c)
        listed = parse_complexity(complexity)[0] if complexity else None
        return keep, [(None if keep is None else len(keep), None if listed is None else len(listed))], [listed]
    union, per_file, listed_all, owner = set(), [], [], {}
    any_none = any_some = False
    for i, row in enumerate(spec.filters):
        cx = row.complexity or complexity
        keep = load_cell_filter(row.whitelist or whitelist, cx,
                                row.min_complexity if row.min_complexity is not None else min_c,
                                row.max_complexity if row.max_complexity is not None else max_c)
        listed = parse_complexity(cx)[0] if cx else None
        if listed is not None and row.complexity:
            for bc in listed:
                if bc in owner:
                    raise ValueError(f"barcode {bc} is in the complexity files of both {spec.bams[owner[bc]]} "
                                     f"and {spec.bams[i]}; a barcode must belong to one input")
                owner[bc] = i
        per_file.append((None if keep is None else len(keep), None if listed is None else len(listed)))
        if keep is None:
            any_none = True
        else:
            any_some = True
            union |= keep
        listed_all.append(listed if row.complexity else None)
    if any_some and any_none:
        raise ValueError("some inputs have a cell filter (complexity/whitelist) and some do not; give every input one, or none")
    return (union if any_some else None), per_file, listed_all


def check_same_contigs(bams) -> None:
    def sq(p):
        with pysam.AlignmentFile(p, check_sq=False) as f:
            return list(zip(f.references, f.lengths))
    first = sq(bams[0])
    for b in bams[1:]:
        if sq(b) != first:
            raise ValueError(f"{bams[0]} and {b} have different @SQ contigs or lengths; inputs must share one reference")


def check_barcodes_match(bams, listed) -> None:
    """Most of each input's first BARCODE_SAMPLE mapped primary reads must have a
    barcode in its own complexity file (catches a sheet row with the wrong file)."""
    if len(bams) < 2:
        return
    for b, s in zip(bams, listed):
        if s is None:
            continue
        n = hit = 0
        with pysam.AlignmentFile(b) as f:
            for r in f.fetch(until_eof=True):
                if r.flag & 0x904:
                    continue
                n += 1
                q = r.query_name
                i = q.find(":")
                hit += (q[:i] if i != -1 else q) in s
                if n >= BARCODE_SAMPLE:
                    break
        if n and hit * 2 < n:
            raise ValueError(f"only {hit} of {n} sampled reads in {b} have a barcode in its complexity file; "
                             "is the sheet row pointing at the right file?")
