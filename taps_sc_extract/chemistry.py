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


# --- Detection, defaults, guards (same rules as taps-sc-extract-rs chemistry.rs) ---

DETECT_SAMPLE = 10_000
BAM_FPAIRED = 0x1
BAM_FREVERSE = 0x10


def bsbolt_version(header_text: str):
    """BSBolt version from an @PG line (``VN:BSB-x.y.z-...``, or ``bsbolt`` in it)."""
    for line in header_text.splitlines():
        if not line.startswith("@PG"):
            continue
        for field in line.split("\t"):
            if field.startswith("VN:BSB-"):
                return field[len("VN:BSB-"):].split("-")[0]
        if "bsbolt" in line.lower():
            return "unknown"
    return None


def check_ys(read) -> None:
    """A directional BSBolt YS (W_C2T / C_C2T) that agrees with the flag's strand."""
    ys = read.get_tag("YS")
    if "G2A" in ys:
        raise ValueError(
            f"read {read.query_name} has YS:Z:{ys}: non-directional or PBAT bisulfite libraries are not supported"
        )
    reverse = bool(read.flag & BAM_FREVERSE)
    if (ys, reverse) not in (("W_C2T", False), ("C_C2T", True)):
        raise ValueError(
            f"read {read.query_name} has YS:Z:{ys} but flag {read.flag}: the strand tag and the alignment "
            "disagree (use --strand-source flag to trust the flag)"
        )


def detect(bam, verify_ys: bool = True) -> Dict[str, object]:
    """Sample primary mapped reads from the start of an open pysam AlignmentFile.
    With ``verify_ys``, a G2A or flag-disagreeing YS is an error."""
    det = {"sampled": 0, "paired": 0, "with_ys": 0, "bsbolt": bsbolt_version(str(bam.header))}
    for read in bam.fetch(until_eof=True):
        if read.flag & (0x4 | 0x100 | 0x800):
            continue
        det["sampled"] += 1
        if read.flag & BAM_FPAIRED:
            det["paired"] += 1
        if read.has_tag("YS"):
            det["with_ys"] += 1
            if verify_ys:
                check_ys(read)
        if det["sampled"] >= DETECT_SAMPLE:
            break
    return det


def resolve(requested: str, det: Dict[str, object], strand_source: str = "auto") -> str:
    """Resolve --chemistry auto|taps|bisulfite against what the BAM looks like."""
    bs_tags = det["with_ys"] > 0
    all_se = det["sampled"] > 0 and det["paired"] == 0
    if requested == "auto":
        if bs_tags:
            return "bisulfite"
        if det["paired"] > 0:
            return "taps"
        raise ValueError(
            f"cannot tell the chemistry from the BAM ({det['sampled']} single-end reads sampled, no BSBolt YS "
            "tags); pass --chemistry taps or --chemistry bisulfite"
        )
    if requested == "taps":
        if bs_tags:
            raise ValueError("--chemistry taps, but the BAM carries BSBolt bisulfite tags (YS); use --chemistry bisulfite")
        if all_se:
            raise ValueError(f"--chemistry taps needs paired-end reads, but all {det['sampled']} sampled reads are single-end")
        return "taps"
    if requested == "bisulfite":
        if det["paired"] > 0:
            raise ValueError("--chemistry bisulfite supports single-end reads only, but the BAM has paired reads")
        if not bs_tags and det["sampled"] > 0 and strand_source == "auto":
            raise ValueError(
                "--chemistry bisulfite, but no read carries a BSBolt YS strand tag; "
                "pass --strand-source flag to take the strand from the alignment flag"
            )
        return "bisulfite"
    raise ValueError(f"unknown chemistry {requested!r}")


# Reads one molecule can put at a position: both mates for paired-end TAPS
# (--site-depth-limit auto; Rust Chemistry::reads_per_molecule).
READS_PER_MOLECULE = {"taps": 2, "bisulfite": 1}


def defaults(chemistry: str) -> Dict[str, object]:
    """Per-chemistry defaults for options the user did not set (see chemistry.rs)."""
    if chemistry == "taps":
        return {"min_mapq": 0, "min_baseq": 20, "max_depth": 0, "trim_r1": 0, "trim_r2": 10, "trim_3p": 0, "baq": True}
    # Bisulfite SE (BSBolt): drop MAPQ < 10 (0 = not bisulfite-unique), trim
    # 10 bp off the sequenced 3' end (Tn5 gap fill-in), no BAQ.
    return {"min_mapq": 10, "min_baseq": 20, "max_depth": 0, "trim_r1": 0, "trim_r2": 0, "trim_3p": 10, "baq": False}


def check_reference(fasta, bam_lengths: Dict[str, int], contigs) -> None:
    """The FASTA must be the unconverted genome: contigs present with the BAM's
    lengths, and ~100 kb from the middle of the first with C and G >= 5% of A/C/G/T."""
    for c in contigs:
        if c not in bam_lengths:
            continue
        if c not in fasta.index:
            raise ValueError(f"contig {c} is in the BAM but not in the FASTA")
        flen = fasta.get_reference_length(c)
        if flen != bam_lengths[c]:
            raise ValueError(f"contig {c}: FASTA length {flen} differs from the BAM header ({bam_lengths[c]})")
    if not contigs:
        return
    first = contigs[0]
    if first not in fasta.index:
        return
    length = fasta.get_reference_length(first)
    acgt = c_n = g_n = 0
    pos = length // 2
    while acgt < 100_000 and pos < length:
        end = min(pos + 100_000, length)
        seq = fasta.fetch(first, pos, end).upper()
        c, g = seq.count("C"), seq.count("G")
        c_n += c
        g_n += g
        acgt += c + g + seq.count("A") + seq.count("T")
        pos = end
    if acgt >= 10_000 and min(c_n, g_n) < 0.05 * acgt:
        raise ValueError(
            f"the FASTA looks bisulfite-converted on {first} (C {100*c_n/acgt:.1f}%, G {100*g_n/acgt:.1f}% of "
            f"{acgt} bases): pass the unconverted genome the aligner index was built from"
        )


def guard_window(bam, contig: str, start: int, end: int, min_mapq: int, verify_ys: bool) -> None:
    """Bisulfite read guard for one window (Python engine): every read the
    pileup would use must be single-end and, unless --strand-source flag, carry
    a directional YS that agrees with its flag. Same rules as the Rust engine."""
    skip = DEFAULT_FLAG_FILTER | BAM_FSUPPLEMENTARY
    for read in bam.fetch(contig, start, end):
        if read.flag & skip or read.mapping_quality < min_mapq:
            continue
        if read.flag & BAM_FPAIRED:
            raise ValueError(f"read {read.query_name} is paired: --chemistry bisulfite supports single-end reads only")
        if verify_ys:
            if not read.has_tag("YS"):
                raise ValueError(
                    f"read {read.query_name} has no BSBolt YS strand tag "
                    "(use --strand-source flag to take the strand from the flag)"
                )
            check_ys(read)
