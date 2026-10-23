"""
M-bias from a sample of reads, and the automatic trim rule (``--trim auto``).

Design and parameters: code/bisulfite-extraction-plan.md §5.1. The Rust engine
(``mbias.rs``) samples the same windows and writes the same count table; the
trim rule is a pure function of that table (``choose_trims``), identical in
both engines.

Count table rows: (context CG|CH, mate 1|2, end 5|3, length bin, position,
methylated, unmethylated). Positions count from the read's sequenced end in
BAM SEQ coordinates (soft clips included, as the trims are), on reads at least
``2 * P`` long, so a short read's 3' artifact never lands in the 5' curve.
The plateau counts bases at least ``P`` from both ends of reads at least
``PLATEAU_MIN_LEN`` long.
"""

import math
import os
from typing import Dict, List, Optional, Tuple

from .calling import DEFAULT_TRIM_R2

P = 20
PLATEAU_MIN_LEN = 60
LENGTH_BINS = [(40, 59), (60, 99), (100, 149), (150, 1_000_000)]
N_WINDOWS = 64
TARGET_READS = 300_000
MIN_WINDOW = 10_000
DELTA = {"CG": 0.02, "CH": 0.01}
Z = 3.0
STABLE_RUN = 3
MAX_AUTO_TRIM = 15
FALLBACK_HEAD = 5
CANDIDATE_3P = {"taps": 0, "bisulfite": 10}


MIN_CALLS = 500


def min_calls_default(plateau_n: int = 0) -> int:
    """Calls a read-end position needs before it can set or stop the trim.
    Fixed, per position and context: scaling it with the plateau (dominated by
    CH) made it unreachable for CG on real data."""
    env = os.environ.get("TAPS_MBIAS_MIN_CALLS")
    return int(env) if env else MIN_CALLS


def length_bin(n: int) -> Optional[int]:
    for i, (lo, hi) in enumerate(LENGTH_BINS):
        if lo <= n <= hi:
            return i
    return None


def bin_label(i: int) -> str:
    lo, hi = LENGTH_BINS[i]
    return f"{lo}-{hi}" if hi < 1_000_000 else f"{lo}+"


# --------------------------------------------------------------------------
# Sampling windows (deterministic; same formula in mbias.rs)
# --------------------------------------------------------------------------

def sample_windows(contigs: List[Tuple[str, int, int]], start: Optional[int] = None,
                   end: Optional[int] = None) -> List[Tuple[str, int, int]]:
    """Windows to sample from ``contigs`` = [(name, length, mapped reads)].

    If the selection holds at most TARGET_READS mapped reads, every read is
    used (one window per contig, or the --start/--end range). Otherwise
    N_WINDOWS windows are dealt to contigs in proportion to their mapped reads
    (largest remainder; ties to the earlier contig) and spaced evenly within
    each contig; their width aims at TARGET_READS reads in total."""
    if start is not None or end is not None:
        name, length, _ = contigs[0]
        return [(name, start or 0, end if end is not None else length)]
    total_reads = sum(m for _, _, m in contigs)
    if total_reads <= TARGET_READS:
        return [(name, 0, length) for name, length, _ in contigs]
    total_len = sum(l for _, l, m in contigs if m > 0)
    width = max(MIN_WINDOW, int(TARGET_READS / total_reads * total_len / N_WINDOWS))
    quotas = [N_WINDOWS * m / total_reads for _, _, m in contigs]
    counts = [int(q) for q in quotas]
    order = sorted(range(len(contigs)), key=lambda i: (-(quotas[i] - counts[i]), i))
    for i in order[: N_WINDOWS - sum(counts)]:
        counts[i] += 1
    out = []
    for (name, length, _), k in zip(contigs, counts):
        for j in range(k):
            centre = (2 * j + 1) * length // (2 * k)
            ws = max(0, centre - width // 2)
            out.append((name, ws, min(length, ws + width)))
    return out


# --------------------------------------------------------------------------
# Counting (Python engine)
# --------------------------------------------------------------------------

def add_counts(a, b):
    """Add another input's (table, plateau, reads) over the same windows (Rust ``Counts::add``)."""
    table, plateau, n = a
    for k, v in b[0].items():
        e = table.setdefault(k, [0, 0])
        e[0] += v[0]
        e[1] += v[1]
    for ctx in plateau:
        plateau[ctx][0] += b[1][ctx][0]
        plateau[ctx][1] += b[1][ctx][1]
    return table, plateau, n + b[2]


def count_table(bam, fasta, windows, chemistry: str, min_mapq: int, min_baseq: int):
    """Count methylation by position from each sequenced end. Returns
    (table dict, plateau dict, reads used)."""
    from .chemistry import flag_filter, tables

    strand_map, call_table = tables(chemistry)
    skip = flag_filter(chemistry)
    table: Dict[Tuple[str, int, int, int, int], List[int]] = {}
    plateau = {"CG": [0, 0], "CH": [0, 0]}
    nreads = 0
    for contig, ws, we in windows:
        clen = fasta.get_reference_length(contig)
        # One reference fetch per window, padded for reads that run past its end.
        wr0, wr1 = max(0, ws - 2), min(clen, we + 2000)
        wref = fasta.fetch(contig, wr0, wr1).upper()
        for read in bam.fetch(contig, ws, we):
            if not (ws <= read.reference_start < we):
                continue  # count each read once: in the window holding its start
            if read.flag & skip or read.mapping_quality < min_mapq:
                continue
            strand = strand_map.get(read.flag)
            if strand is None:
                continue
            seq, quals = read.query_sequence, read.query_qualities
            if not seq:
                continue
            L = len(seq)
            nreads += 1
            reverse = bool(read.flag & 16)
            mate = 2 if read.flag & 128 else 1
            lb = length_bin(L) if L >= 2 * P else None
            r0 = max(0, read.reference_start - 2)
            r1 = min(clen, read.reference_end + 2)
            if r0 >= wr0 and r1 <= wr1:
                ref = wref[r0 - wr0:r1 - wr0]
            else:
                ref = fasta.fetch(contig, r0, r1).upper()
            q, r = 0, read.reference_start
            for op, n in read.cigartuples:
                if op in (0, 7, 8):
                    for k in range(n):
                        qp, rp = q + k, r + k
                        i = rp - r0
                        rb = ref[i]
                        if strand == "OT" and rb != "C" or strand == "OB" and rb != "G":
                            continue
                        if quals[qp] < min_baseq:
                            continue
                        call = call_table.get((strand, rb, seq[qp]))
                        if call is None:
                            continue
                        if rb == "C":
                            ctx = "CG" if ref[i + 1:i + 2] == "G" else "CH"
                        else:
                            ctx = "CG" if i >= 1 and ref[i - 1] == "C" else "CH"
                        d5, d3 = (L - 1 - qp, qp) if reverse else (qp, L - 1 - qp)
                        idx = 0 if call == 1 else 1
                        if lb is not None:
                            if d5 < P:
                                table.setdefault((ctx, mate, 5, lb, d5), [0, 0])[idx] += 1
                            if d3 < P:
                                table.setdefault((ctx, mate, 3, lb, d3), [0, 0])[idx] += 1
                        if L >= PLATEAU_MIN_LEN and d5 >= P and d3 >= P:
                            plateau[ctx][idx] += 1
                    q += n
                    r += n
                elif op in (1, 4):
                    q += n
                elif op in (2, 3):
                    r += n
    return table, plateau, nreads


# --------------------------------------------------------------------------
# The rule (pure function of the table; identical in mbias.rs)
# --------------------------------------------------------------------------

def curve(table, ctx: str, mate: int, end: int) -> List[Tuple[int, int]]:
    out = []
    for p in range(P):
        m = u = 0
        for lb in range(len(LENGTH_BINS)):
            c = table.get((ctx, mate, end, lb, p))
            if c:
                m += c[0]
                u += c[1]
        out.append((m, u))
    return out


# Each end is compared with its own baseline: positions BASELINE_FROM..P-1 of
# its curve (Rust mbias::baseline). A plateau pooled over mates put TAPS read 1
# and read 2 (about 10 points apart mid-read) both out of band.
BASELINE_FROM = 15
BASELINE_WARN = 0.10  # baseline this far from the plateau: the bias may run past 15 bp


def baseline(cur: List[Tuple[int, int]], min_calls: int) -> Optional[float]:
    """Positions BASELINE_FROM..P-1 pooled after dropping the highest- and
    lowest-rate position (Rust mbias::baseline)."""
    if sum(m + u for m, u in cur[BASELINE_FROM:]) < min_calls:
        return None
    import functools
    tail = [x for x in cur[BASELINE_FROM:] if x[0] + x[1] > 0]
    # Stable sort by rate, compared exactly by cross-multiplication (as in Rust).
    tail.sort(key=functools.cmp_to_key(lambda a, b: (a[0] * (b[0] + b[1]) > b[0] * (a[0] + a[1])) - (a[0] * (b[0] + b[1]) < b[0] * (a[0] + a[1]))))
    mid = tail[1:-1] if len(tail) >= 3 else tail
    m = sum(x[0] for x in mid)
    u = sum(x[1] for x in mid)
    return m / (m + u) if m + u else None


def walk(cur: List[Tuple[int, int]], plateau_rate: float, delta: float, min_calls: int):
    """Trim for one end: walk inward from p = 0; while a position is biased the
    trim is p + 1; stop at the first STABLE_RUN consecutive in-band-or-unknown
    positions. Returns (trim, unknown_in_head)."""
    trim, run, unknown_head = 0, 0, False
    for p, (m, u) in enumerate(cur):
        n = m + u
        if n < min_calls:
            if p < FALLBACK_HEAD:
                unknown_head = True
            run += 1
        else:
            rate = m / n
            se = math.sqrt(rate * (1 - rate) / n)
            if abs(rate - plateau_rate) > max(delta, Z * se):
                trim, run = p + 1, 0
            else:
                run += 1
        if run >= STABLE_RUN:
            break
    return min(trim, MAX_AUTO_TRIM), unknown_head


def choose_trims(table, plateau, chemistry: str, mates=(1,)) -> Dict[str, object]:
    """{'trim_r1', 'trim_r2', 'trim_3p', 'notes': [...]} from the count table."""
    notes = []
    pn = {ctx: sum(plateau[ctx]) for ctx in ("CG", "CH")}
    min_calls = min_calls_default(pn["CG"] + pn["CH"])
    result = {}
    three = []
    for mate in mates:
        for end in (5, 3):
            per_ctx = {}
            for ctx in ("CG", "CH"):
                if pn[ctx] < min_calls:
                    per_ctx[ctx] = None
                    continue
                cur = curve(table, ctx, mate, end)
                pooled = plateau[ctx][0] / pn[ctx]
                rate = baseline(cur, min_calls)
                if rate is None:
                    rate = pooled
                if ctx == "CG" and abs(rate - pooled) > BASELINE_WARN:
                    lab = f"R{mate} " if len(mates) > 1 else ""
                    notes.append(f"{lab}{end}': CG at 15-19 bp is {100 * rate:.1f}% against a {100 * pooled:.1f}% plateau; "
                                 "the bias may extend past 15 bp, check the M-bias plot")
                t, unknown_head = walk(cur, rate, DELTA[ctx], min_calls)
                per_ctx[ctx] = (t, unknown_head)
            label = f"R{mate} " if len(mates) > 1 else ""
            cg = per_ctx["CG"]
            if cg is None or cg[1]:
                fallback = CANDIDATE_3P[chemistry] if end == 3 else (DEFAULT_TRIM_R2 if chemistry == "taps" and mate == 2 else 0)
                ch = per_ctx["CH"]
                if ch is not None and not ch[1]:
                    t = ch[0]
                    notes.append(f"{label}{end}': too few CG calls; CH sets the trim ({t} bp)")
                else:
                    t = fallback
                    notes.append(f"{label}{end}': too few calls near the read end; using the default {t} bp")
            else:
                t = cg[0]
                ch = per_ctx["CH"]
                if ch is not None and not ch[1] and ch[0] > t:
                    notes.append(f"{label}{end}': CH suggests {ch[0]} bp; using CG {t} bp")
            if t >= MAX_AUTO_TRIM:
                notes.append(f"{label}{end}': M-bias reaches the {MAX_AUTO_TRIM} bp cap; check the library")
            if end == 5:
                result["trim_r2" if mate == 2 else "trim_r1"] = t
            else:
                three.append(t)
    result["trim_3p"] = max(three) if three else 0
    result.setdefault("trim_r2", 0 if len(mates) == 1 else result.get("trim_r2", 0))
    result["notes"] = notes
    result["min_calls"] = min_calls
    return result


# --------------------------------------------------------------------------
# Files
# --------------------------------------------------------------------------

TSV_HEADER = "context\tmate\tend\tlength_bin\tposition\tmethylated\tunmethylated"


def write_tsv(path: str, table, plateau, nreads: int, trims: Dict[str, object]) -> None:
    with open(path, "w") as fh:
        fh.write(f"# reads\t{nreads}\n")
        for ctx in ("CG", "CH"):
            fh.write(f"# plateau\t{ctx}\t{plateau[ctx][0]}\t{plateau[ctx][1]}\n")
        fh.write(f"# trims\t5p={trims.get('trim_r1', 0)}\tr2_5p={trims.get('trim_r2', 0)}\t3p={trims.get('trim_3p', 0)}\n")
        fh.write(TSV_HEADER + "\n")
        for key in sorted(table):
            ctx, mate, end, lb, p = key
            m, u = table[key]
            fh.write(f"{ctx}\t{mate}\t{end}\t{bin_label(lb)}\t{p}\t{m}\t{u}\n")


def read_tsv(path: str):
    table, plateau, nreads, trims = {}, {"CG": [0, 0], "CH": [0, 0]}, 0, {}
    labels = {bin_label(i): i for i in range(len(LENGTH_BINS))}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if line.startswith("# reads"):
                nreads = int(f[1])
            elif line.startswith("# plateau"):
                plateau[f[1]] = [int(f[2]), int(f[3])]
            elif line.startswith("# trims"):
                for kv in f[1:]:
                    k, v = kv.split("=")
                    trims[k] = int(v)
            elif line.startswith("context") or line.startswith("#"):
                continue
            else:
                table[(f[0], int(f[1]), int(f[2]), labels[f[3]], int(f[4]))] = [int(f[5]), int(f[6])]
    return table, plateau, nreads, trims


def plot(tsv_path: str, out_prefix: str, title: str) -> None:
    """<prefix>.mbias.png and .pdf from a count table (seaborn)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns
    from matplotlib.ticker import MultipleLocator

    table, plateau, nreads, trims = read_tsv(tsv_path)
    mates = sorted({k[1] for k in table}) or [1]
    sns.set_theme(style="whitegrid", context="notebook", font_scale=0.9)
    rows = []
    for (ctx, mate, end, lb, p), (m, u) in table.items():
        if m + u >= 200:
            rows.append((ctx, mate, end, bin_label(lb) + " bp", p, 100 * m / (m + u)))
    df = pd.DataFrame(rows, columns=["context", "mate", "end", "read length", "position", "methylation"])
    order = [bin_label(i) + " bp" for i in range(len(LENGTH_BINS))]
    palette = dict(zip(order, sns.color_palette("crest", len(order))))
    ncol = 2 * len(mates)
    fig, axes = plt.subplots(2, ncol, figsize=(5.5 * ncol, 7.5), sharex=True, squeeze=False)
    min_calls = min_calls_default(sum(plateau["CG"]) + sum(plateau["CH"]))
    for row, ctx in enumerate(("CG", "CH")):
        pm, pu = plateau[ctx]
        pr = 100 * pm / (pm + pu) if pm + pu else float("nan")
        for mi, mate in enumerate(mates):
            for ei, end in enumerate((5, 3)):
                ax = axes[row][mi * 2 + ei]
                ax.axhspan(pr - 100 * DELTA[ctx], pr + 100 * DELTA[ctx], color="0.9", zorder=0,
                           label=f"plateau ±{100*DELTA[ctx]:.0f} pp")
                ax.axhline(pr, color="0.45", lw=0.9, ls="--", zorder=1)
                sub = df[(df.context == ctx) & (df.mate == mate) & (df.end == end)]
                if len(sub):
                    sns.lineplot(data=sub, x="position", y="methylation", hue="read length", hue_order=order,
                                 palette=palette, lw=1.1, marker="o", markersize=3, markeredgewidth=0, ax=ax,
                                 legend=(row == 0 and mi == len(mates) - 1 and ei == 1))
                cur = curve(table, ctx, mate, end)
                good = [(p, 100 * m / (m + u)) for p, (m, u) in enumerate(cur) if m + u >= min_calls]
                low = [(p, 100 * m / (m + u)) for p, (m, u) in enumerate(cur) if 0 < m + u < min_calls]
                if good:
                    ax.plot(*zip(*good), color="black", lw=2.2, marker="o", ms=4.5, label="all reads ≥40 bp", zorder=5)
                if low:
                    ax.plot(*zip(*low), "o", color="0.7", ms=4.5, label="too few calls", zorder=5)
                t = (trims.get("r2_5p" if mate == 2 else "5p", 0) if end == 5 else trims.get("3p", 0))
                ax.axvline(t - 0.5, color="#c0392b", lw=1.8, label="trim", zorder=6)
                mlabel = f"R{mate} " if len(mates) > 1 else ""
                ax.set_title(f"{ctx} — {mlabel}from sequenced {end}′ end (trim {t} bp)", fontweight="bold")
                ax.set_xlabel("position from read end (bases)" if row == 1 else "")
                ax.set_ylabel(f"{ctx} methylation (%)" if mi == 0 and ei == 0 else "")
                ax.set_xlim(-0.5, P - 0.5)
                ax.xaxis.set_major_locator(MultipleLocator(2))
                ax.xaxis.set_minor_locator(MultipleLocator(1))
                ax.grid(which="minor", axis="x", alpha=0.25)
    seen = {}
    for ax in axes.flat:
        for h, l in zip(*ax.get_legend_handles_labels()):
            seen.setdefault(l, h)
        if ax.get_legend():
            ax.get_legend().remove()
    fig.legend(list(seen.values()), list(seen.keys()), loc="lower center", ncol=5, frameon=False, fontsize=8.5)
    pcg = plateau["CG"]
    pch = plateau["CH"]
    fig.suptitle(
        f"{title}: M-bias — trim 5′ {trims.get('5p', 0)} bp, 3′ {trims.get('3p', 0)} bp   "
        f"({nreads:,} reads; plateau CG {100*pcg[0]/max(1,sum(pcg)):.1f}%, CH {100*pch[0]/max(1,sum(pch)):.2f}%)",
        fontsize=11, fontweight="bold")
    sns.despine(fig)
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    fig.savefig(out_prefix + ".mbias.png", dpi=150)
    fig.savefig(out_prefix + ".mbias.pdf")
    plt.close(fig)


def bam_prefix(bam_path: str) -> str:
    base = os.path.basename(bam_path)
    for ext in (".bam", ".cram", ".sam", ".tsv", ".txt"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


def output_dir(out_path: str, n_shards: int) -> str:
    """Directory for the M-bias files: the shard directory, or beside a single file."""
    if n_shards > 1 or os.path.isdir(out_path) or not out_path.endswith((".h5", ".hdf5")):
        if out_path.endswith((".h5", ".hdf5")) and n_shards > 1:
            return os.path.splitext(out_path)[0]
        return out_path
    return os.path.dirname(out_path) or "."


def main(argv=None) -> int:
    """python -m taps_sc_extract.mbias --plot <tsv> [--title T]  (used by the Rust engine)."""
    import argparse

    ap = argparse.ArgumentParser(description="Plot an M-bias table written by taps-sc-extract.")
    ap.add_argument("--plot", required=True, help="<prefix>.mbias.tsv")
    ap.add_argument("--title", default=None)
    a = ap.parse_args(argv)
    prefix = a.plot[: -len(".mbias.tsv")] if a.plot.endswith(".mbias.tsv") else a.plot
    plot(a.plot, prefix, a.title or os.path.basename(prefix))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
