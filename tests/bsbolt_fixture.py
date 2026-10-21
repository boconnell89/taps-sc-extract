"""
Synthetic BSBolt-style bisulfite BAM for tests (generated at test time; no data
committed). Single-end, directional (YS W_C2T / C_C2T), Ultima-like reads:
variable length, indels, soft clips, binned qualities, MAPQ 0 reads, opposite-
strand supplementary records, and a 3' end bias (Tn5 gap filled with
unmethylated C). Every read carries a correct BSBolt XB call string.

``expected_calls`` computes the per-cell, per-site (t, c) an extractor should
produce for given filters straight from the generator's truth, independent of
both engines' code.
"""

import os
import random
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import pysam

CONTIG = "chr19"
REGION = (5_000_000, 5_300_000)
CELLS = [f"Z000{k}-{bc}-TTGGACTC" for k, bc in enumerate(
    ["ACGAGCCGAA", "AGAAGAATAC", "TTGGCAATCC", "CAGTTACGGA", "GGATCCTTAA", "TACCGGTTAG"], start=1)]
P_METH = {"CG": 0.8, "CH": 0.05}
BIAS_3P = 10  # bases at the sequenced 3' end forced unmethylated on biased reads


@dataclass
class Call:
    ref_pos: int     # 0-based reference position of the C (OT) / G (OB)
    qpos: int        # index into SEQ
    ctx: str         # 'CG' / 'CH'
    meth: int        # 1 methylated (read C / G), 0 unmethylated (read T / A)
    qual: int


@dataclass
class ReadTruth:
    qname: str
    cell: str
    flag: int
    mapq: int
    seq_len: int
    calls: List[Call] = field(default_factory=list)


def _context_ot(ref: str, i: int) -> Tuple[str, str]:
    """(Amethyst context, BSBolt letter) of a C at ref[i]."""
    n1 = ref[i + 1] if i + 1 < len(ref) else "N"
    n2 = ref[i + 2] if i + 2 < len(ref) else "N"
    if n1 == "G":
        return "CG", "x"
    if n2 == "G":
        return "CH", "y"
    return "CH", "z"


def _context_ob(ref: str, i: int) -> Tuple[str, str]:
    p1 = ref[i - 1] if i >= 1 else "N"
    p2 = ref[i - 2] if i >= 2 else "N"
    if p1 == "C":
        return "CG", "x"
    if p2 == "C":
        return "CH", "y"
    return "CH", "z"


def _qual(rng) -> int:
    r = rng.random()
    return 35 if r < 0.7 else (25 if r < 0.9 else 12)


def make_bsbolt_bam(out_dir: str, fasta_path: str, n_reads: int = 6000, seed: int = 7,
                    flip_ys_of: int = -1) -> Tuple[str, List[ReadTruth]]:
    """Write ``bsbolt_fixture.bam`` (coordinate-sorted, indexed) and return its
    path and the truth for every primary read. ``flip_ys_of`` gives one primary
    read a YS that disagrees with its flag (for the guard test)."""
    rng = random.Random(seed)
    fa = pysam.FastaFile(fasta_path)
    r0, r1 = REGION
    ref = fa.fetch(CONTIG, r0 - 10, r1 + 400).upper()
    off = r0 - 10
    names = [n for n in fa.references if n in {f"chr{i}" for i in range(1, 20)} | {"chrX", "chrY", "chrM"}]
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": n, "LN": fa.get_reference_length(n)} for n in names],
        "PG": [{"ID": "bwa", "PN": "bwa", "VN": "BSB-1.2.1-BWA-fork-0.7.17", "CL": "bwa mem -Y (synthetic)"}],
    }
    tid = names.index(CONTIG)
    unsorted = os.path.join(out_dir, "bsbolt_fixture.unsorted.bam")
    truth: List[ReadTruth] = []
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for i in range(n_reads):
            cell = rng.choice(CELLS)
            reverse = rng.random() < 0.5
            laq = rng.randint(30, 160)          # aligned query bases (M + I)
            indel = rng.random()
            ins = rng.randint(1, 3) if indel < 0.25 else 0
            dele = rng.randint(1, 3) if 0.25 <= indel < 0.5 else 0
            m1 = rng.randint(5, max(5, laq - ins - 5))
            m2 = laq - ins - m1
            if m2 < 1:
                m1, m2 = laq - ins - 1, 1
            s_left = rng.randint(3, 8) if rng.random() < 0.04 else 0
            s_right = rng.randint(3, 8) if rng.random() < 0.04 else 0
            span = m1 + dele + m2
            pos = rng.randint(r0, r1 - span)
            # reference coordinates of M bases
            mref = [pos + k for k in range(m1)] + [pos + m1 + dele + k for k in range(m2)]
            seq = [rng.choice("ACGT") for _ in range(s_left)]
            quals = [_qual(rng) for _ in range(s_left)]
            seg_bases = []  # (qpos, ref_pos or None)
            q = s_left
            for k in range(m1):
                seg_bases.append((q, mref[k])); q += 1
            ins_q = list(range(q, q + ins)); q += ins
            for k in range(m2):
                seg_bases.append((q, mref[m1 + k])); q += 1
            total = s_left + laq + s_right
            seq += ["N"] * (laq)  # placeholders, filled below
            seq += [rng.choice("ACGT") for _ in range(s_right)]
            quals += [_qual(rng) for _ in range(laq + s_right)]
            for qq in ins_q:
                seq[qq] = rng.choice("ACGT")
            biased = rng.random() < 0.3
            calls = []
            for qq, rp in seg_bases:
                b = ref[rp - off]
                # sequenced-end distance: forward reads are sequenced from SEQ start
                d3 = (total - 1 - qq) if not reverse else qq
                if not reverse and b == "C":
                    ctx, _ = _context_ot(ref, rp - off)
                    meth = rng.random() < P_METH[ctx] and not (biased and d3 < BIAS_3P)
                    b = "C" if meth else "T"
                elif reverse and b == "G":
                    ctx, _ = _context_ob(ref, rp - off)
                    meth = rng.random() < P_METH[ctx] and not (biased and d3 < BIAS_3P)
                    b = "G" if meth else "A"
                if rng.random() < 0.002:  # sequencing error
                    b = rng.choice("ACGT")
                seq[qq] = b
            # truth calls and XB
            xb_letters = {}
            for qq, rp in seg_bases:
                b = seq[qq]
                rb = ref[rp - off]
                if not reverse and rb == "C" and b in "CT":
                    ctx, letter = _context_ot(ref, rp - off)
                    meth = int(b == "C")
                elif reverse and rb == "G" and b in "GA":
                    ctx, letter = _context_ob(ref, rp - off)
                    meth = int(b == "G")
                else:
                    continue
                calls.append(Call(rp, qq, ctx, meth, quals[qq]))
                xb_letters[qq - s_left] = letter.upper() if meth else letter
            xb, run = [], 0
            for k in range(laq):
                if k in xb_letters:
                    if run:
                        xb.append(str(run))
                    xb.append(xb_letters[k]); run = 0
                else:
                    run += 1
            if run:
                xb.append(str(run))
            cig = []
            if s_left:
                cig.append((4, s_left))
            cig.append((0, m1))
            if ins:
                cig.append((1, ins))
            if dele:
                cig.append((2, dele))
            cig.append((0, m2))
            if s_right:
                cig.append((4, s_right))
            r = rng.random()
            mapq = 0 if r < 0.05 else (5 if r < 0.08 else 60)
            flag = 16 if reverse else 0
            qname = f"{cell}:{i}"
            ys = "C_C2T" if reverse else "W_C2T"
            if i == flip_ys_of:
                ys = "W_C2T" if reverse else "C_C2T"
            a = pysam.AlignedSegment(out.header)
            a.query_name = qname
            a.flag = flag
            a.reference_id = tid
            a.reference_start = pos
            a.mapping_quality = mapq
            a.cigartuples = cig
            a.query_sequence = "".join(seq)
            a.query_qualities = pysam.qualitystring_to_array("".join(chr(x + 33) for x in quals))
            a.set_tags([("XB", "".join(xb)), ("AS", laq), ("YS", ys)])
            out.write(a)
            truth.append(ReadTruth(qname, cell, flag, mapq, total, calls))
            # Opposite-strand supplementary (chimera); must never be called.
            if rng.random() < 0.03:
                sp = rng.randint(r0, r1 - 60)
                s2 = pysam.AlignedSegment(out.header)
                s2.query_name = qname
                s2.flag = 2048 | (0 if reverse else 16)
                s2.reference_id = tid
                s2.reference_start = sp
                s2.mapping_quality = 60
                s2.cigartuples = [(0, 50)]
                s2.query_sequence = ref[sp - off: sp - off + 50].replace("C", "T") if reverse else ref[sp - off: sp - off + 50].replace("G", "A")
                s2.query_qualities = pysam.qualitystring_to_array("D" * 50)
                s2.set_tags([("XB", "50"), ("YS", "W_C2T" if reverse else "C_C2T")])
                out.write(s2)
    bam = os.path.join(out_dir, "bsbolt_fixture.bam")
    pysam.sort("-o", bam, unsorted)
    pysam.index(bam)
    os.remove(unsorted)
    return bam, truth


def expected_calls(truth: List[ReadTruth], min_mapq: int, min_baseq: int, trim_5p: int, trim_3p: int,
                   contexts=("CG", "CH")) -> Dict[str, Dict[Tuple[str, int], Tuple[int, int]]]:
    """{ctx: {(cell, 1-based pos): (t, c)}} from the generator's truth."""
    acc: Dict[str, Dict[Tuple[str, int], List[int]]] = {"CG": defaultdict(lambda: [0, 0]), "CH": defaultdict(lambda: [0, 0])}
    for rd in truth:
        if rd.mapq < min_mapq:
            continue
        reverse = bool(rd.flag & 16)
        if trim_5p + trim_3p >= rd.seq_len:
            continue
        lo, hi = (trim_3p, rd.seq_len - trim_5p) if reverse else (trim_5p, rd.seq_len - trim_3p)
        for c in rd.calls:
            if c.ctx not in contexts or c.qual < min_baseq or not (lo <= c.qpos < hi):
                continue
            acc[c.ctx][(rd.cell, c.ref_pos + 1)][1 if c.meth else 0] += 1
    return {ctx: {k: tuple(v) for k, v in d.items()} for ctx, d in acc.items()}
