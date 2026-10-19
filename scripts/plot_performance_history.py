"""Plot the README's Performance History: throughput and peak memory relative to
the last version before the 2026-10-14 work (e19e2b4).

Numbers are the README table's (reference workstation, full test BAM, each
release's defaults). Writes docs/img/perf_history_throughput.svg and
docs/img/perf_history_memory.svg. Run from the repository root:

    python3 scripts/plot_performance_history.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# (label, end to end M reads/s, processing M reads/s, peak memory GB)
RELEASES = [
    ("Before\n10-14\n(e19e2b4)", 0.352, 0.914, 33.2),
    ("100 Mb\nwindows, LZF\n10-14", 0.56, 0.84, 29.7),
    ("Balanced\nwindows\n10-15", 0.97, 1.30, 21.7),
    ("TAPS-aware\nBAQ\n10-16", 1.00, 1.31, 23.2),
    ("Skip matching\nreads' BAQ\n10-16", 1.64, 2.60, 23.2),
    ("SIMD BAQ\n10-16", 1.76, 3.24, 24.0),
    ("blosc\ndefault\n10-16", 1.95, 3.41, 23.9),
    ("Mask in\nBAQ kernel\n10-17", 1.85, 3.29, 21.2),
    ("Faster BAQ\npreparation\n10-18", 2.07, 3.84, 21.7),
    ("Read-centric\ncalling\n10-18", 2.45, 5.08, 21.7),
    ("Writers fed\nduring\nextraction\n10-18", 2.69, 4.48, 16.7),
]

OUT = Path("docs/img")
BLUE, ORANGE, GREEN = "#2f6db5", "#e08a1e", "#3a9a5b"


def style(ax, title, ylabel):
    ax.set_title(title, fontsize=12, loc="left")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color="#dddddd", linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="x", labelsize=7.5)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    labels = [r[0] for r in RELEASES]
    x = range(len(RELEASES))
    e2e0, proc0, mem0 = RELEASES[0][1], RELEASES[0][2], RELEASES[0][3]

    # Throughput relative to the version before 2026-10-14 (e19e2b4).
    e2e = [r[1] / e2e0 for r in RELEASES]
    proc = [r[2] / proc0 for r in RELEASES]
    fig, ax = plt.subplots(figsize=(10, 4.2), facecolor="white")
    ax.plot(x, proc, marker="o", color=ORANGE, label="Processing (read, BAQ, call)")
    ax.plot(x, e2e, marker="o", color=BLUE, linewidth=2.5, label="End to end (incl. HDF5 write)")
    for xi, v in zip(x, e2e):
        ax.annotate(f"{v:.1f}×", (xi, v), textcoords="offset points", xytext=(0, -14), ha="center", fontsize=8, color=BLUE)
    for xi in (len(x) - 2, len(x) - 1):
        ax.annotate(f"{proc[xi]:.1f}×", (xi, proc[xi]), textcoords="offset points", xytext=(0, 7), ha="center", fontsize=8, color=ORANGE)
    ax.axhline(1, color="#999999", linewidth=0.8, linestyle="--")
    ax.set_xticks(list(x), labels)
    ax.set_ylim(0, max(max(proc), max(e2e)) * 1.15)
    style(ax, "Throughput (mapped reads/s) relative to the version before 2026-10-14", "× before 10-14")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "perf_history_throughput.svg")
    plt.close(fig)

    # Peak memory (extraction + writer processes) relative to that version.
    mem = [100 * r[3] / mem0 for r in RELEASES]
    fig, ax = plt.subplots(figsize=(10, 4.2), facecolor="white")
    bars = ax.bar(x, mem, color=[GREEN if i == len(x) - 1 else "#8cc3a0" for i in x], width=0.65)
    for b, r in zip(bars, RELEASES):
        ax.annotate(f"{b.get_height():.0f}%\n{r[3]:.1f} GB", (b.get_x() + b.get_width() / 2, b.get_height()),
                    textcoords="offset points", xytext=(0, 3), ha="center", fontsize=7.5)
    ax.set_xticks(list(x), labels)
    ax.set_ylim(0, 120)
    style(ax, "Peak memory (extraction + writer processes) relative to the version before 2026-10-14", "% of before 10-14")
    fig.tight_layout()
    fig.savefig(OUT / "perf_history_memory.svg")
    plt.close(fig)


if __name__ == "__main__":
    main()
