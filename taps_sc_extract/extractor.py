"""
Main single-cell TAPS methylation extraction driver.

Processes coordinate-sorted BAM files chromosome-by-chromosome,
accumulating methylation calls per barcode and writing to an
Amethyst-compatible HDF5 file.
"""

import gc
import logging
import resource
import time
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pysam

from .contigs import describe_primary, primary_contigs
from .barcode import extract_barcode, load_cell_filter
from .inputs import prepare_alignment
from .calling import (
    pysam_max_depth,
    DEFAULT_TRIM_R1,
    DEFAULT_TRIM_R2,
    trim_window,
)
from .chemistry import flag_filter, guard_window, tables
from .h5_writer import set_run_attrs
from .sitemask import Cursor, apply_cap, load_blacklist, merge, window_depth_mask
from .context import classify_context, classify_trinucleotide_context
from .h5_writer import AmethystH5Writer, METH_DTYPE

logger = logging.getLogger("taps_sc_extract")

# The mm10 primary chromosomes (what `primary_contigs` picks for an mm10 BAM).
# Defaults now come from the BAM header (`contigs.primary_contigs`).
CANONICAL_CONTIGS = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"]


def get_peak_rss_mb() -> float:
    """Return peak resident set size in megabytes."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def extract_methylation(
    bam_path: str,
    fasta_path: str,
    out_h5_path: str,
    chroms: Optional[Iterable[str]] = None,
    whitelist_path: Optional[str] = None,
    complexity_path: Optional[str] = None,
    min_complexity: Optional[float] = None,
    max_complexity: Optional[float] = None,
    min_base_quality: int = 20,
    min_mapq: int = 0,
    max_depth: int = 0,
    compute_baq: bool = True,
    ignore_orphans: bool = True,
    ignore_overlaps: bool = True,
    trim_r1: int = DEFAULT_TRIM_R1,
    trim_r2: int = DEFAULT_TRIM_R2,
    contexts: Tuple[str, ...] = ("CG", "CH"),
    chemistry: str = "taps",
    trim_3p: int = 0,
    verify_ys: bool = True,
    trim_mode: str = "manual",
    site_depth_limit: int = 0,
    blacklist_path: Optional[str] = None,
    cell_site_cap: int = 0,
) -> Dict[str, Any]:
    """
    Extract methylation calls from a TAPS BAM file into an Amethyst HDF5 file.

    Parameters:
        bam_path: Path to coordinate-sorted, indexed BAM file.
        fasta_path: Path to reference FASTA file (indexed with .fai).
        out_h5_path: Output HDF5 file path.
        chroms: Iterable of chromosome names to process (default: the primary chromosomes in the BAM header, `contigs.primary_contigs`).
        whitelist_path: Optional path to barcode whitelist / annotation file.
        complexity_path: Optional complexity file (rank, barcode, total, unique, %).
        min_complexity / max_complexity: Keep barcodes whose unique-read count
            (column 4) is within these bounds (inclusive).
        min_base_quality: Minimum base quality for pileup (default 20).
        min_mapq: Minimum mapping quality (default 0).
        max_depth: Maximum reads per position, all cells together (default 100000; 0 = unlimited).
        compute_baq: Whether to compute Base Alignment Quality (default True).
        ignore_orphans: Whether to ignore unpaired reads (default True).
        ignore_overlaps: Whether to skip overlapping mate bases (default True).
        trim_r1: Drop this many 5' sequenced bases of read 1 (default 0).
        trim_r2: Drop this many 5' sequenced bases of read 2 (default 10).
        contexts: Methylation contexts to extract, any of "CG", "CH" (default
            both). Columns of the other context are skipped and its HDF5 group
            is not written.

    Returns:
        Summary dictionary with counts, methylation percentages, and performance metrics.
    """
    whitelist: Optional[Set[str]] = load_cell_filter(
        whitelist_path, complexity_path, min_complexity, max_complexity
    )

    # Open BAM and FASTA
    # Builds a missing index / converts SAM; keep input_tmp alive for the run.
    source, input_tmp = prepare_alignment(bam_path, fasta_path)
    bam = source.open()
    fasta = pysam.FastaFile(fasta_path)

    # Determine contigs to process
    bam_contigs = set(bam.references)
    if chroms is None:
        target_contigs = primary_contigs(bam.references)
        if not target_contigs:
            raise ValueError(
                "no primary chromosomes (1, 2, ..., X, Y with or without a chr prefix) in the BAM header; pass chroms"
            )
        logger.info(f"Contigs: {describe_primary(target_contigs)} from the BAM header")
    else:
        target_contigs = [c for c in chroms if c in bam_contigs]

    # Index statistics for mapped read counts
    idx_stats = {s.contig: s.mapped for s in bam.get_index_statistics()}
    total_target_reads = sum(idx_stats.get(c, 0) for c in target_contigs)

    logger.info(f"Processing {len(target_contigs)} contig(s) with {total_target_reads:,} total mapped reads.")
    logger.info(f"5' trim: R1={trim_r1} bp, R2={trim_r2} bp; 3' trim: {trim_3p} bp; chemistry: {chemistry}")
    strand_map, call_table = tables(chemistry)
    set_run_attrs({
        "chemistry": chemistry,
        "trim_r1": trim_r1,
        "trim_r2": trim_r2,
        "trim_3p": trim_3p,
        "trim_mode": trim_mode,
        "min_mapq": min_mapq,
        "min_baseq": min_base_quality,
        "max_depth": max_depth,
        "site_depth_limit": site_depth_limit,
        "cell_site_cap": cell_site_cap,
        "blacklist": blacklist_path or "",
        "extractor": "taps-sc-extract (python)",
    })

    total_stats = {
        "CpG": {"c": 0, "t": 0},
        "CHG": {"c": 0, "t": 0},
        "CHH": {"c": 0, "t": 0},
        "CG": {"c": 0, "t": 0},
        "CH": {"c": 0, "t": 0},
    }

    blacklist = load_blacklist(blacklist_path) if blacklist_path else {}
    all_barcodes: Set[str] = set()
    contig_breakdown = []
    start_total_time = time.time()

    with AmethystH5Writer(out_h5_path, mode="w") as writer:
        for contig_idx, contig in enumerate(target_contigs, 1):
            start_contig_time = time.time()
            contig_mapped_reads = idx_stats.get(contig, 0)
            logger.info(f"[{contig_idx}/{len(target_contigs)}] Starting {contig} ({contig_mapped_reads:,} mapped reads)...")

            # Load full contig sequence into memory
            ref_seq = fasta.fetch(contig).upper()
            ref_len = len(ref_seq)

            # Per-contig accumulators: accum[barcode][pos_1based] = [t_count, c_count]
            accum_cg: Dict[str, Dict[int, List[int]]] = defaultdict(lambda: defaultdict(lambda: [0, 0]))
            accum_ch: Dict[str, Dict[int, List[int]]] = defaultdict(lambda: defaultdict(lambda: [0, 0]))

            contig_counts = {
                "CpG": {"c": 0, "t": 0},
                "CHG": {"c": 0, "t": 0},
                "CHH": {"c": 0, "t": 0},
                "CG": {"c": 0, "t": 0},
                "CH": {"c": 0, "t": 0},
            }

            contig_barcodes: Set[str] = set()

            if chemistry == "bisulfite":
                guard_window(bam, contig, 0, ref_len, min_mapq, verify_ys)
            masked = list(blacklist.get(contig, []))
            if site_depth_limit > 0:
                masked += window_depth_mask(bam, contig, 0, ref_len, site_depth_limit, flag_filter(chemistry),
                                            min_mapq, ignore_orphans, whitelist)
            mask = Cursor(merge(masked))
            pileup_iter = bam.pileup(
                contig,
                ignore_overlaps=ignore_overlaps,
                min_base_quality=min_base_quality,
                stepper="samtools",
                flag_filter=flag_filter(chemistry),
                compute_baq=compute_baq,
                ignore_orphans=ignore_orphans,
                max_depth=pysam_max_depth(max_depth),
                min_mapping_quality=min_mapq,
            )

            for col in pileup_iter:
                pos = col.reference_pos  # 0-based
                if pos >= ref_len or mask.hit(pos):
                    continue

                ref_b = ref_seq[pos]
                if ref_b == "C":
                    n1 = ref_seq[pos + 1] if pos + 1 < ref_len else "N"
                    if n1 == "G":
                        ctx_tri = "CpG"
                        ctx = "CG"
                    else:
                        n2 = ref_seq[pos + 2] if pos + 2 < ref_len else "N"
                        if n2 == "G" and n1 in ("A", "C", "T"):
                            ctx_tri = "CHG"
                        else:
                            ctx_tri = "CHH"
                        ctx = "CH"
                elif ref_b == "G":
                    p1 = ref_seq[pos - 1] if pos > 0 else "N"
                    if p1 == "C":
                        ctx_tri = "CpG"
                        ctx = "CG"
                    else:
                        p2 = ref_seq[pos - 2] if pos > 1 else "N"
                        if p2 == "C" and p1 in ("A", "G", "T"):
                            ctx_tri = "CHG"
                        else:
                            ctx_tri = "CHH"
                        ctx = "CH"
                else:
                    continue
                if ctx not in contexts:
                    continue

                target_accum = accum_cg if ctx == "CG" else accum_ch
                pos_1based = pos + 1

                touched = set()
                for pileupread in col.pileups:
                    if pileupread.is_del or pileupread.is_refskip:
                        continue
                    qpos = pileupread.query_position
                    if qpos is None:
                        continue

                    aln = pileupread.alignment
                    strand = strand_map.get(aln.flag)
                    if strand is None:
                        continue
                    lo, hi = trim_window(aln.query_length, aln.flag, trim_r1, trim_r2, trim_3p)
                    if qpos < lo or qpos >= hi:
                        continue

                    read_base = aln.query_sequence[qpos]
                    call = call_table.get((strand, ref_b, read_base))
                    if call is None:
                        continue

                    qname = aln.query_name
                    sep_idx = qname.find(":")
                    bc = qname[:sep_idx] if sep_idx != -1 else qname

                    if whitelist is not None and bc not in whitelist:
                        continue

                    contig_barcodes.add(bc)
                    touched.add(bc)
                    counts = target_accum[bc][pos_1based]
                    if call == 1:
                        counts[1] += 1  # methylated (c)
                        contig_counts[ctx_tri]["c"] += 1
                        contig_counts[ctx]["c"] += 1
                    else:
                        counts[0] += 1  # unmethylated (t)
                        contig_counts[ctx_tri]["t"] += 1
                        contig_counts[ctx]["t"] += 1
                if cell_site_cap and touched:
                    apply_cap(target_accum, touched, pos_1based, cell_site_cap, contig_counts, ctx, ctx_tri)

            # Flush accumulated contig calls to HDF5
            contig_bytes = contig.encode("ascii")
            for ctx, accum in [("CG", accum_cg), ("CH", accum_ch)]:
                for bc, pos_dict in accum.items():
                    if not pos_dict:
                        continue
                    sorted_items = sorted(pos_dict.items())
                    n = len(sorted_items)
                    arr = np.empty(n, dtype=METH_DTYPE)
                    for i, (p, (t, c)) in enumerate(sorted_items):
                        pct = (100.0 * c / (c + t)) if (c + t) > 0 else 0.0
                        arr[i] = (contig_bytes, p, pct, t, c)
                    writer.append_data(ctx, bc, arr)

            all_barcodes.update(contig_barcodes)

            # Update total stats
            for k in ("CpG", "CHG", "CHH", "CG", "CH"):
                total_stats[k]["c"] += contig_counts[k]["c"]
                total_stats[k]["t"] += contig_counts[k]["t"]

            cg_tot = contig_counts["CG"]["c"] + contig_counts["CG"]["t"]
            cg_pct = (100.0 * contig_counts["CG"]["c"] / cg_tot) if cg_tot > 0 else 0.0
            ch_tot = contig_counts["CH"]["c"] + contig_counts["CH"]["t"]
            ch_pct = (100.0 * contig_counts["CH"]["c"] / ch_tot) if ch_tot > 0 else 0.0

            elapsed = time.time() - start_contig_time
            peak_rss = get_peak_rss_mb()
            reads_per_sec = contig_mapped_reads / elapsed if elapsed > 0 else 0.0

            contig_info = {
                "contig": contig,
                "mapped_reads": contig_mapped_reads,
                "elapsed_sec": elapsed,
                "reads_per_sec": reads_per_sec,
                "active_barcodes": len(contig_barcodes),
                "cg_calls": cg_tot,
                "cg_pct": cg_pct,
                "ch_calls": ch_tot,
                "ch_pct": ch_pct,
                "peak_rss_mb": peak_rss,
            }
            contig_breakdown.append(contig_info)

            logger.info(
                f"Finished {contig} in {elapsed:.1f}s ({reads_per_sec:,.0f} reads/s) | "
                f"CG: {contig_counts['CG']['c']:,}/{cg_tot:,} ({cg_pct:.2f}%) | "
                f"CH: {contig_counts['CH']['c']:,}/{ch_tot:,} ({ch_pct:.2f}%) | "
                f"Barcodes: {len(contig_barcodes):,} | Peak RSS: {peak_rss:.1f} MB"
            )

            # Free memory
            del accum_cg, accum_ch, ref_seq
            gc.collect()

    bam.close()
    fasta.close()

    total_elapsed = time.time() - start_total_time
    cg_all = total_stats["CG"]["c"] + total_stats["CG"]["t"]
    cg_all_pct = (100.0 * total_stats["CG"]["c"] / cg_all) if cg_all > 0 else 0.0
    chg_all = total_stats["CHG"]["c"] + total_stats["CHG"]["t"]
    chg_all_pct = (100.0 * total_stats["CHG"]["c"] / chg_all) if chg_all > 0 else 0.0
    chh_all = total_stats["CHH"]["c"] + total_stats["CHH"]["t"]
    chh_all_pct = (100.0 * total_stats["CHH"]["c"] / chh_all) if chh_all > 0 else 0.0
    ch_all = total_stats["CH"]["c"] + total_stats["CH"]["t"]
    ch_all_pct = (100.0 * total_stats["CH"]["c"] / ch_all) if ch_all > 0 else 0.0

    total_calls = cg_all + ch_all
    total_mod = total_stats["CG"]["c"] + total_stats["CH"]["c"]
    pooled_pct = (100.0 * total_mod / total_calls) if total_calls > 0 else 0.0

    reads_per_min = (total_target_reads / (total_elapsed / 60.0)) if total_elapsed > 0 else 0.0
    cells_per_min = (len(all_barcodes) / (total_elapsed / 60.0)) if total_elapsed > 0 else 0.0
    final_peak_rss = get_peak_rss_mb()

    summary = {
        "total_elapsed_sec": total_elapsed,
        "total_elapsed_min": total_elapsed / 60.0,
        "total_mapped_reads": total_target_reads,
        "reads_per_min": reads_per_min,
        "reads_per_sec": total_target_reads / total_elapsed if total_elapsed > 0 else 0.0,
        "total_unique_barcodes": len(all_barcodes),
        "cells_per_min": cells_per_min,
        "peak_rss_mb": final_peak_rss,
        "stats": total_stats,
        "percentages": {
            "CpG": cg_all_pct,
            "CHG": chg_all_pct,
            "CHH": chh_all_pct,
            "CH": ch_all_pct,
            "Pooled": pooled_pct,
        },
        "contig_breakdown": contig_breakdown,
    }

    logger.info("=" * 70)
    logger.info("EXTRACTION SUMMARY & GROUND TRUTH COMPARISON")
    logger.info("=" * 70)
    logger.info(f"Total Wall Time:       {total_elapsed:.1f}s ({total_elapsed/60.0:.2f} min)")
    logger.info(f"Total Mapped Reads:    {total_target_reads:,}")
    logger.info(f"Throughput (Reads):    {reads_per_min:,.0f} reads/min ({total_target_reads/total_elapsed:,.0f} reads/s)")
    logger.info(f"Total Unique Barcodes: {len(all_barcodes):,}")
    if whitelist is not None and not all_barcodes:
        logger.warning(
            "No reads matched the barcode filter (--whitelist / --complexity). "
            "Check that its barcodes use the same format as the QNAME prefix before ':'."
        )
    logger.info(f"Throughput (Cells):    {cells_per_min:,.1f} cells/min")
    logger.info(f"Peak Memory (RSS):     {final_peak_rss:.1f} MB")
    logger.info("-" * 70)
    logger.info(f"CpG Methylation:       {cg_all_pct:.3f}% ({total_stats['CpG']['c']:,} / {cg_all:,})")
    logger.info(f"CHG Methylation:       {chg_all_pct:.3f}% ({total_stats['CHG']['c']:,} / {chg_all:,})")
    logger.info(f"CHH Methylation:       {chh_all_pct:.3f}% ({total_stats['CHH']['c']:,} / {chh_all:,})")
    logger.info(f"CH (Combined):         {ch_all_pct:.3f}% ({total_stats['CH']['c']:,} / {ch_all:,})")
    logger.info(f"Pooled (All C/G sites):{pooled_pct:.3f}% ({total_mod:,} / {total_calls:,})")
    logger.info("=" * 70)

    return summary
