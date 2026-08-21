#!/usr/bin/env python3
"""Python reference stats for one genomic window (PR0 pileup parity)."""
import argparse
import time
from collections import defaultdict

import pysam

from taps_sc_extract.calling import FLAG_STRAND_MAP, MCTOT_LOOKUP
from taps_sc_extract.fasta import FastFaiReader


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-b", "--bam", required=True)
    p.add_argument("-f", "--fasta", required=True)
    p.add_argument("-c", "--chrom", default="chr19")
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--end", type=int, required=True)
    p.add_argument("--min-baseq", type=int, default=20)
    p.add_argument("--min-mapq", type=int, default=0)
    p.add_argument("--max-depth", type=int, default=250)
    p.add_argument(
        "--no-baq",
        action="store_true",
        help="Disable BAQ (still a no-op unless --fasta-baq is set).",
    )
    p.add_argument(
        "--fasta-baq",
        action="store_true",
        help="Pass pysam.FastaFile to pileup so compute_baq actually runs "
        "(current taps-sc-extract production does not do this).",
    )
    p.add_argument(
        "--dump-sites",
        default=None,
        help="Write barcode/ctx/pos/t/c TSV for parity diffs vs Rust.",
    )
    args = p.parse_args()

    fai = FastFaiReader(args.fasta)
    pad = 2
    ref_len = fai.get_reference_length(args.chrom)
    ref_start = max(0, args.start - pad)
    ref_end = min(ref_len, args.end + pad)
    ref_seq = fai.fetch(args.chrom, ref_start, ref_end)

    bam = pysam.AlignmentFile(args.bam, "rb")
    fastafile = pysam.FastaFile(args.fasta) if args.fasta_baq else None
    stats = defaultdict(lambda: {"c": 0, "t": 0})
    barcodes = set()
    sites = defaultdict(lambda: [0, 0]) if args.dump_sites else None
    t0 = time.time()
    pileup_kw = dict(
        ignore_overlaps=True,
        min_base_quality=args.min_baseq,
        stepper="samtools",
        compute_baq=not args.no_baq,
        ignore_orphans=True,
        max_depth=args.max_depth,
        min_mapping_quality=args.min_mapq,
    )
    if fastafile is not None:
        pileup_kw["fastafile"] = fastafile
    pileup = bam.pileup(args.chrom, start=args.start, stop=args.end, **pileup_kw)
    for col in pileup:
        pos = col.reference_pos
        if pos < args.start or pos >= args.end:
            continue
        rel = pos - ref_start
        if rel >= len(ref_seq):
            continue
        ref_b = ref_seq[rel]
        if ref_b == "C":
            n1 = ref_seq[rel + 1] if rel + 1 < len(ref_seq) else "N"
            if n1 == "G":
                ctx_tri, ctx = "CpG", "CG"
            else:
                n2 = ref_seq[rel + 2] if rel + 2 < len(ref_seq) else "N"
                ctx_tri = "CHG" if n2 == "G" and n1 in ("A", "C", "T") else "CHH"
                ctx = "CH"
        elif ref_b == "G":
            p1 = ref_seq[rel - 1] if rel > 0 else "N"
            if p1 == "C":
                ctx_tri, ctx = "CpG", "CG"
            else:
                p2 = ref_seq[rel - 2] if rel > 1 else "N"
                ctx_tri = "CHG" if p2 == "C" and p1 in ("A", "G", "T") else "CHH"
                ctx = "CH"
        else:
            continue
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip:
                continue
            qpos = pr.query_position
            if qpos is None:
                continue
            aln = pr.alignment
            strand = FLAG_STRAND_MAP.get(aln.flag)
            if strand is None:
                continue
            seq = aln.query_sequence
            if seq is None:
                continue
            read_base = seq[qpos]
            call = MCTOT_LOOKUP.get((strand, ref_b, read_base))
            if call is None:
                continue
            qname = aln.query_name
            sep = qname.find(":")
            bc = qname[:sep] if sep != -1 else qname
            barcodes.add(bc)
            key = "c" if call == 1 else "t"
            stats[ctx_tri][key] += 1
            stats[ctx][key] += 1
            if sites is not None:
                slot = sites[(bc, ctx, pos + 1)]
                if call == 1:
                    slot[1] += 1
                else:
                    slot[0] += 1
    dt = time.time() - t0
    cg = stats["CG"]["c"] + stats["CG"]["t"]
    ch = stats["CH"]["c"] + stats["CH"]["t"]
    cpg = stats["CpG"]["c"] + stats["CpG"]["t"]
    def pct(c, t):
        tot = c + t
        return 0.0 if tot == 0 else 100.0 * c / tot
    print(
        f"cells={len(barcodes)} "
        f"CG={pct(stats['CG']['c'], stats['CG']['t']):.4f}% ({stats['CG']['c']}/{cg}) "
        f"CH={pct(stats['CH']['c'], stats['CH']['t']):.4f}% ({stats['CH']['c']}/{ch}) "
        f"CpG={pct(stats['CpG']['c'], stats['CpG']['t']):.4f}% "
        f"elapsed_s={dt:.2f}"
    )
    if sites is not None:
        with open(args.dump_sites, "w", encoding="utf-8") as fh:
            fh.write("barcode\tctx\tpos\tt\tc\n")
            for (bc, ctx, pos) in sorted(sites):
                t, c = sites[(bc, ctx, pos)]
                fh.write(f"{bc}\t{ctx}\t{pos}\t{t}\t{c}\n")
        print(f"wrote site TSV {args.dump_sites}", flush=True)


if __name__ == "__main__":
    main()
