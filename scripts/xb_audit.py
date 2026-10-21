#!/usr/bin/env python3
"""
Compare BSBolt's own per-read calls (the XB tag) with taps-sc-extract's
reference pileup on a bisulfite BAM (code/bisulfite-extraction-plan.md §8).

XB is decoded with a CIGAR walk: digits = aligned query bases skipped (CIGAR
M/=/X/I; soft clips are not in XB), each letter = one call at the current
aligned base (x/X CG, y/Y CHG, z/Z CHH; uppercase = methylated). Calls on
insertion bases cannot occur; calls are placed at the reference position of
their M base.

Reports:
  selftest  Every read's XB length equals its aligned query length; on reads
            without soft clips every letter sits on the strand's reference
            C/G with an informative read base and the right context. Letter
            mismatches on soft-clipped reads are BSBolt's soft-clip bug
            (bwamem.c passes an unshifted query to bwa_gen_cigar2).
  compare   Calls from XB, starting with no filters (what XB-based extraction
            gives) and adding the extractor's filters one at a time, against
            the extractor's HDF5 output. Each step's change is attributed to
            that filter; what remains at the end is attributed to the soft-clip
            bug and anything else (listed).

Usage:
  xb_audit.py -b sample.bam -f ref.fa --h5 out/master.h5 -c chr19 \\
      [--min-mapq 10 --min-baseq 20 --trim-5p 0 --trim-3p 10] -o prefix

Writes <prefix>.xb_audit.tsv (one row per step and context) and
<prefix>.xb_audit.png/.pdf (per-cell methylation, XB vs pileup, and the
differences by cause).
"""

import argparse
import re
import sys
from collections import defaultdict

import h5py
import pysam

try:
    import hdf5plugin  # noqa: F401  (Blosc output)
except ImportError:  # pragma: no cover
    pass

TOKEN = re.compile(r"(\d*)([A-Za-z])")


def decode_xb(read):
    """[(qpos, ref_pos, letter)] from a read's XB tag, or None without XB.
    Raises ValueError if the XB string is longer or shorter than the read's
    aligned query bases (M/=/X/I)."""
    if not read.has_tag("XB"):
        return None
    xb = read.get_tag("XB")
    aligned = []  # (qpos, ref_pos or None) per aligned query base
    q, r = 0, read.reference_start
    for op, n in read.cigartuples:
        if op in (0, 7, 8):
            aligned += [(q + k, r + k) for k in range(n)]
            q += n
            r += n
        elif op == 1:
            aligned += [(q + k, None) for k in range(n)]
            q += n
        elif op == 4:
            q += n
        elif op in (2, 3):
            r += n
    calls, k = [], 0
    for digits, letter in TOKEN.findall(xb):
        k += int(digits) if digits else 0
        if k >= len(aligned):
            raise ValueError(f"{read.query_name}: XB {xb} runs past {len(aligned)} aligned bases")
        qpos, rpos = aligned[k]
        calls.append((qpos, rpos, letter))
        k += 1
    tail = re.search(r"(\d+)$", xb)
    total = k + (int(tail.group(1)) if tail else 0)
    if total != len(aligned):
        raise ValueError(f"{read.query_name}: XB {xb} covers {total} bases, read has {len(aligned)} aligned")
    return calls


def ctx_of(letter):
    return "CG" if letter in "xX" else "CH"


def trim_window(length, reverse, t5, t3):
    if t5 + t3 >= length:
        return 0, 0
    return (t3, length - t5) if reverse else (t5, length - t3)


STEPS = ["XB as is", "+ supplementary dropped", "+ MAPQ", "+ base quality", "+ trims"]


def xb_calls(bam, fa, contig, start, end, args, st, cells=None):
    """{step: {ctx: {(cell, pos1): [t, c]}}} for the cumulative filter steps on
    calls in [start, end) of one contig; adds to the self-test counters `st`.
    Reads that start before `start` were counted by the previous block."""
    out = {s: {"CG": defaultdict(lambda: [0, 0]), "CH": defaultdict(lambda: [0, 0])} for s in STEPS}
    clen = fa.get_reference_length(contig)
    first = start == 0
    if True:
        for read in bam.fetch(contig, start, end):
            if read.flag & (0x4 | 0x100 | 0x200 | 0x400):
                continue
            new = first or read.reference_start >= start
            try:
                calls = decode_xb(read)
            except ValueError:
                st["xb length mismatch"] += new
                continue
            if calls is None:
                st["no XB"] += new
                continue
            st["reads"] += new
            clipped = any(op == 4 for op, _ in read.cigartuples)
            reverse = bool(read.flag & 16)
            cell = read.query_name.split(":", 1)[0]
            if cells is not None and cell not in cells:
                continue
            seq, quals = read.query_sequence, read.query_qualities
            lo, hi = trim_window(read.query_length, reverse, args.trim_5p, args.trim_3p)
            ref = fa.fetch(contig, max(0, read.reference_start - 2), min(clen, read.reference_end + 2)).upper()
            roff = max(0, read.reference_start - 2)
            for qpos, rpos, letter in calls:
                if rpos is None:
                    st["letter on an insertion"] += 1
                    continue
                if not (start <= rpos < end):
                    continue
                # self-test against the reference
                i = rpos - roff
                rb, b = ref[i], seq[qpos].upper()
                if not reverse:
                    ok_base = rb == "C" and (b == "C" if letter.isupper() else b == "T")
                    n1, n2 = ref[i + 1:i + 2], ref[i + 2:i + 3]
                    want = "x" if n1 == "G" else ("y" if n2 == "G" else "z")
                else:
                    ok_base = rb == "G" and (b == "G" if letter.isupper() else b == "A")
                    p1, p2 = ref[i - 1:i] if i >= 1 else "N", ref[i - 2:i - 1] if i >= 2 else "N"
                    want = "x" if p1 == "C" else ("y" if p2 == "C" else "z")
                ok_ctx = letter.lower() == want or (letter.lower() == "y" and want == "z")  # CNG: CHG in BSBolt
                if not (ok_base and ok_ctx):
                    st["letter mismatch (soft-clipped read: BSBolt bug)" if clipped else "letter mismatch (unclipped read)"] += 1
                st["calls"] += 1
                ctx, meth = ctx_of(letter), 1 if letter.isupper() else 0
                key = (cell, rpos + 1)
                keep = [True,
                        not read.flag & 0x800,
                        read.mapping_quality >= args.min_mapq,
                        quals[qpos] >= args.min_baseq,
                        lo <= qpos < hi]
                for s, step in enumerate(STEPS):
                    if all(keep[: s + 1]):
                        out[step][ctx][key][meth] += 1
    return out


def h5_contig(path, contig, contexts=("CG", "CH")):
    """{ctx: {cell: (pos, t, c)}} numpy arrays of one contig's rows. Reads one
    cell at a time, so memory is one cell's dataset plus the contig's rows."""
    import numpy as np
    want = contig.encode()
    got = {}
    with h5py.File(path, "r") as f:
        for ctx in contexts:
            d = {}
            if ctx in f:
                for bc in f[ctx]:
                    x = f[ctx][bc]["1"][()]
                    m = x["chr"] == want
                    if m.any():
                        d[bc] = (np.ascontiguousarray(x["pos"][m]), np.ascontiguousarray(x["t"][m]), np.ascontiguousarray(x["c"][m]))
            got[ctx] = d
    return got


def h5_cells(path):
    with h5py.File(path, "r") as f:
        return {bc for ctx in ("CG", "CH") if ctx in f for bc in f[ctx]}


def h5_block(rows, start, end):
    """{(cell, pos1): (t, c)} for 0-based [start, end) from h5_contig's arrays."""
    import numpy as np
    d = {}
    for bc, (pos, t, c) in rows.items():
        i, j = np.searchsorted(pos, start + 1), np.searchsorted(pos, end + 1)
        for p, tt, cc in zip(pos[i:j].tolist(), t[i:j].tolist(), c[i:j].tolist()):
            d[(bc, p)] = (tt, cc)
    return d


def diff(a, b):
    """Counts of sites only in a, only in b, shared but different, shared equal."""
    only_a = sum(1 for k in a if k not in b)
    only_b = sum(1 for k in b if k not in a)
    shared = [k for k in a if k in b]
    differ = sum(1 for k in shared if tuple(a[k]) != tuple(b[k]))
    return only_a, only_b, differ, len(shared) - differ


def totals(d):
    t = sum(v[0] for v in d.values())
    c = sum(v[1] for v in d.values())
    return t, c


def add_per_cell(agg, d):
    for (cell, _), (t, c) in d.items():
        agg[cell][0] += t
        agg[cell][1] += c


def per_cell_rate(agg):
    return {cell: 100.0 * c / (t + c) for cell, (t, c) in agg.items() if t + c}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-b", "--bam", required=True)
    ap.add_argument("-f", "--fasta", required=True)
    ap.add_argument("--h5", required=True, help="taps-sc-extract output (master.h5 or single file)")
    ap.add_argument("-c", "--chroms", required=True, help="comma-separated contigs the extraction covered")
    ap.add_argument("--start", type=int, default=None)
    ap.add_argument("--end", type=int, default=None)
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--trim-5p", type=int, default=0)
    ap.add_argument("--trim-3p", type=int, default=10)
    ap.add_argument("-o", "--out-prefix", required=True)
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--block", type=int, default=5_000_000, help="bp compared at a time (bounds memory)")
    ap.add_argument("--all-cells", action="store_true", help="count XB for every barcode, not only the HDF5's cells")
    args = ap.parse_args()

    contigs = [c.strip() for c in args.chroms.split(",") if c.strip()]
    cells = None
    if not args.all_cells:
        # Only the cells the extractor kept (e.g. after a complexity cutoff).
        cells = h5_cells(args.h5)
        print(f"Restricting XB to the {len(cells):,} cells in {args.h5} (--all-cells to use every barcode)")

    # Compared block by block: only one block's sites are ever held as dicts.
    st = defaultdict(int)
    acc = {(step, ctx): [0] * 6 for step in STEPS for ctx in ("CG", "CH")}  # t c only_xb only_pileup differ equal
    pt = {ctx: [0, 0] for ctx in ("CG", "CH")}
    cell_xb = {ctx: defaultdict(lambda: [0, 0]) for ctx in ("CG", "CH")}
    cell_pl = {ctx: defaultdict(lambda: [0, 0]) for ctx in ("CG", "CH")}
    fa = pysam.FastaFile(args.fasta)
    with pysam.AlignmentFile(args.bam) as bam:
        for contig in contigs:
            clen = fa.get_reference_length(contig)
            r0 = args.start if args.start is not None else 0
            r1 = min(args.end, clen) if args.end is not None else clen
            rows = h5_contig(args.h5, contig)
            for b0 in range(r0, r1, args.block):
                b1 = min(b0 + args.block, r1)
                steps = xb_calls(bam, fa, contig, b0, b1, args, st, cells)
                for ctx in ("CG", "CH"):
                    pl = h5_block(rows[ctx], b0, b1)
                    t, c = totals(pl)
                    pt[ctx][0] += t
                    pt[ctx][1] += c
                    add_per_cell(cell_pl[ctx], pl)
                    add_per_cell(cell_xb[ctx], steps["XB as is"][ctx])
                    for step in STEPS:
                        d = steps[step][ctx]
                        vals = (*totals(d), *diff(d, pl))
                        acc[(step, ctx)] = [x + y for x, y in zip(acc[(step, ctx)], vals)]
                print(f"  {contig}:{b0:,}-{b1:,} done", file=sys.stderr)

    print("XB self-test:")
    for k in sorted(st):
        print(f"  {k:52s} {st[k]:>12,}")

    rows = []
    print("\nXB calls vs pileup output (sites = cell x position):")
    print(f"  {'step':28s} {'ctx':3s} {'calls t/c':>22s} {'meth %':>7s} {'only XB':>9s} {'only pileup':>11s} {'differ':>8s} {'equal':>9s}")
    for step in STEPS:
        for ctx in ("CG", "CH"):
            t, c, oa, ob, df, eq = acc[(step, ctx)]
            rate = 100.0 * c / (t + c) if t + c else 0.0
            rows.append((step, ctx, t, c, rate, oa, ob, df, eq))
            print(f"  {step:28s} {ctx:3s} {t:>10,}/{c:<11,} {rate:7.3f} {oa:9,} {ob:11,} {df:8,} {eq:9,}")
    for ctx in ("CG", "CH"):
        t, c = pt[ctx]
        print(f"  {'pileup (taps-sc-extract)':28s} {ctx:3s} {t:>10,}/{c:<11,} {100.0*c/(t+c) if t+c else 0:7.3f}")

    with open(f"{args.out_prefix}.xb_audit.tsv", "w") as fh:
        fh.write("step\tcontext\tt\tc\tmeth_pct\tsites_only_xb\tsites_only_pileup\tsites_differ\tsites_equal\n")
        for r in rows:
            fh.write("\t".join(str(x) if not isinstance(x, float) else f"{x:.4f}" for x in r) + "\n")
        for ctx in ("CG", "CH"):
            t, c = pt[ctx]
            fh.write(f"pileup\t{ctx}\t{t}\t{c}\t{100.0*c/(t+c) if t+c else 0:.4f}\t\t\t\t\n")
    print(f"\nwrote {args.out_prefix}.xb_audit.tsv")

    if not args.no_plot:
        plot(cell_xb, cell_pl, rows, args.out_prefix)


def plot(cell_xb, cell_pl, rows, prefix):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid", context="notebook", font_scale=0.9)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
    for ax, ctx in zip(axes[:2], ("CG", "CH")):
        a, b = per_cell_rate(cell_xb[ctx]), per_cell_rate(cell_pl[ctx])
        cells = sorted(set(a) & set(b))
        xs, ys = [a[c] for c in cells], [b[c] for c in cells]
        ax.scatter(xs, ys, s=10, alpha=0.6, color=sns.color_palette("crest")[3], edgecolor="none")
        if xs:
            lo, hi = min(xs + ys), max(xs + ys)
            ax.plot([lo, hi], [lo, hi], color="0.5", lw=0.8, ls="--")
        ax.set_xlabel(f"{ctx} methylation per cell, XB as is (%)")
        ax.set_ylabel(f"{ctx} methylation per cell, pileup (%)")
        ax.set_title(f"{ctx}: {len(cells):,} cells", fontweight="bold")
    ax = axes[2]
    labels, vals_cg, vals_ch = [], [], []
    for step in STEPS:
        r = {row[1]: row for row in rows if row[0] == step}
        labels.append(step)
        vals_cg.append(r["CG"][5] + r["CG"][6] + r["CG"][7])
        vals_ch.append(r["CH"][5] + r["CH"][6] + r["CH"][7])
    y = range(len(labels))
    ax.barh([i + 0.2 for i in y], vals_cg, height=0.4, color=sns.color_palette("crest")[1], label="CG")
    ax.barh([i - 0.2 for i in y], vals_ch, height=0.4, color=sns.color_palette("crest")[4], label="CH")
    ax.set_yticks(list(y), labels)
    ax.invert_yaxis()
    ax.set_xlabel("sites differing from the pileup output")
    ax.set_title("Filters added to XB, one at a time", fontweight="bold")
    ax.legend(frameon=False)
    sns.despine(fig)
    fig.tight_layout()
    fig.savefig(f"{prefix}.xb_audit.png", dpi=150)
    fig.savefig(f"{prefix}.xb_audit.pdf")
    print(f"wrote {prefix}.xb_audit.png / .pdf")


if __name__ == "__main__":
    sys.exit(main())
