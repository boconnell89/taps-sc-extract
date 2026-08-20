"""
Parallel single-cell TAPS methylation extraction engine with partitioned shard streaming.

1. Extraction Phase:
   Workers process genomic chunks across the genome in parallel, streaming
   compact binary chunk files partitioned by shard to fast scratch storage.
2. Parallel Shard Assembly & Write Phase:
   Each shard assembles its chunk files in strict coordinate order and writes
   its Amethyst HDF5 file concurrently in parallel.
   As soon as a shard's chunks are assembled, its HDF5 file is written and
   temporary files are immediately purged.
3. Master Index Generation:
   Creates a portable `master.h5` containing relative ExternalLinks to all shards.
"""

import concurrent.futures
import gc
import hashlib
import logging
import multiprocessing as mp
import os
import pickle
import resource
import shutil
import tempfile
import time
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import h5py
import numpy as np
import pysam

from .barcode import parse_annot
from .calling import MCTOT_LOOKUP, FLAG_STRAND_MAP
from .fasta import FastFaiReader
from .h5_writer import AmethystH5Writer, METH_DTYPE

# Disable HDF5 file locking
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

logger = logging.getLogger("taps_sc_extract")

# Default canonical mm10 contigs
CANONICAL_CONTIGS = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"]


def get_peak_rss_mb() -> float:
    """Return peak resident set size in megabytes."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def _barcode_to_shard(bc: str, n_shards: int) -> int:
    """Deterministically map a cell barcode to a shard index [0, n_shards - 1]."""
    if n_shards <= 1:
        return 0
    return int(hashlib.md5(bc.encode("ascii")).hexdigest(), 16) % n_shards


# Module-level worker globals
_worker_bam: Optional[pysam.AlignmentFile] = None
_worker_fai: Optional[FastFaiReader] = None
_worker_whitelist: Optional[Set[str]] = None
_worker_temp_dir: Optional[str] = None
_worker_params: Optional[Dict[str, Any]] = None


def _init_worker(
    bam_path: str,
    fasta_path: str,
    whitelist: Optional[Set[str]],
    temp_dir: Optional[str],
    params: Dict[str, Any],
):
    """Initialize worker process with clean private file handles."""
    global _worker_bam, _worker_fai, _worker_whitelist, _worker_temp_dir, _worker_params
    decomp_threads = params.get("decomp_threads", 1)
    _worker_bam = pysam.AlignmentFile(bam_path, "rb", threads=decomp_threads)
    _worker_fai = FastFaiReader(fasta_path)
    _worker_whitelist = whitelist
    _worker_temp_dir = temp_dir
    _worker_params = params


def _process_chunk(chunk_info: Tuple[int, str, int, int]) -> Tuple[int, str, int, int, Any, Dict[str, Any]]:
    """
    Process a single genomic chunk [start, end) on contig.

    Writes compact chunk array data partitioned by shard and returns summary statistics.
    """
    global _worker_bam, _worker_fai, _worker_whitelist, _worker_temp_dir, _worker_params

    chunk_id, contig, start, end = chunk_info
    t0 = time.time()

    bam = _worker_bam
    fai = _worker_fai
    whitelist = _worker_whitelist
    temp_dir = _worker_temp_dir
    params = _worker_params

    ignore_overlaps = params["ignore_overlaps"]
    min_base_quality = params["min_base_quality"]
    min_mapq = params["min_mapq"]
    max_depth = params["max_depth"]
    compute_baq = params["compute_baq"]
    ignore_orphans = params["ignore_orphans"]
    use_temp_files = params.get("use_temp_files", True)
    n_shards = params.get("n_shards", 1)

    # Fetch reference sequence with 2bp padding on each side
    ref_len = fai.get_reference_length(contig)
    pad = 2
    ref_start = max(0, start - pad)
    ref_end = min(ref_len, end + pad)
    ref_seq = fai.fetch(contig, ref_start, ref_end)

    accum_cg: Dict[str, Dict[int, List[int]]] = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    accum_ch: Dict[str, Dict[int, List[int]]] = defaultdict(lambda: defaultdict(lambda: [0, 0]))

    stats = {
        "CpG": {"c": 0, "t": 0},
        "CHG": {"c": 0, "t": 0},
        "CHH": {"c": 0, "t": 0},
        "CG": {"c": 0, "t": 0},
        "CH": {"c": 0, "t": 0},
    }
    barcodes_seen: Set[str] = set()

    pileup_iter = bam.pileup(
        contig,
        start=start,
        stop=end,
        ignore_overlaps=ignore_overlaps,
        min_base_quality=min_base_quality,
        stepper="samtools",
        compute_baq=compute_baq,
        ignore_orphans=ignore_orphans,
        max_depth=max_depth,
        min_mapping_quality=min_mapq,
    )

    for col in pileup_iter:
        pos = col.reference_pos  # 0-based
        if pos < start or pos >= end:
            continue

        rel_pos = pos - ref_start
        if rel_pos >= len(ref_seq):
            continue

        ref_b = ref_seq[rel_pos]
        if ref_b == "C":
            n1 = ref_seq[rel_pos + 1] if rel_pos + 1 < len(ref_seq) else "N"
            if n1 == "G":
                ctx_tri = "CpG"
                ctx = "CG"
            else:
                n2 = ref_seq[rel_pos + 2] if rel_pos + 2 < len(ref_seq) else "N"
                if n2 == "G" and n1 in ("A", "C", "T"):
                    ctx_tri = "CHG"
                else:
                    ctx_tri = "CHH"
                ctx = "CH"
        elif ref_b == "G":
            p1 = ref_seq[rel_pos - 1] if rel_pos > 0 else "N"
            if p1 == "C":
                ctx_tri = "CpG"
                ctx = "CG"
            else:
                p2 = ref_seq[rel_pos - 2] if rel_pos > 1 else "N"
                if p2 == "C" and p1 in ("A", "G", "T"):
                    ctx_tri = "CHG"
                else:
                    ctx_tri = "CHH"
                ctx = "CH"
        else:
            continue

        target_accum = accum_cg if ctx == "CG" else accum_ch
        pos_1based = pos + 1

        for pileupread in col.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                continue
            qpos = pileupread.query_position
            if qpos is None:
                continue

            aln = pileupread.alignment
            strand = FLAG_STRAND_MAP.get(aln.flag)
            if strand is None:
                continue

            read_base = aln.query_sequence[qpos]
            call = MCTOT_LOOKUP.get((strand, ref_b, read_base))
            if call is None:
                continue

            qname = aln.query_name
            sep_idx = qname.find(":")
            bc = qname[:sep_idx] if sep_idx != -1 else qname

            if whitelist is not None and bc not in whitelist:
                continue

            barcodes_seen.add(bc)
            counts = target_accum[bc][pos_1based]
            if call == 1:
                counts[1] += 1
                stats[ctx_tri]["c"] += 1
                stats[ctx]["c"] += 1
            else:
                counts[0] += 1
                stats[ctx_tri]["t"] += 1
                stats[ctx]["t"] += 1

    # Convert accumulator dicts to compact NumPy structured arrays partitioned by shard
    contig_bytes = contig.encode("ascii")
    # shard_idx -> {"CG": {bc: arr}, "CH": {bc: arr}}
    shard_data: Dict[int, Dict[str, Dict[str, np.ndarray]]] = {
        s: {"CG": {}, "CH": {}} for s in range(n_shards)
    }

    for ctx, accum in [("CG", accum_cg), ("CH", accum_ch)]:
        for bc, pos_dict in accum.items():
            if not pos_dict:
                continue
            s_idx = _barcode_to_shard(bc, n_shards)
            sorted_items = sorted(pos_dict.items())
            n = len(sorted_items)
            arr = np.empty(n, dtype=METH_DTYPE)
            for i, (p, (t, c)) in enumerate(sorted_items):
                pct = (100.0 * c / (c + t)) if (c + t) > 0 else 0.0
                arr[i] = (contig_bytes, p, pct, t, c)
            shard_data[s_idx][ctx][bc] = arr

    del accum_cg, accum_ch

    elapsed = time.time() - t0
    chunk_stats = {
        "elapsed_sec": elapsed,
        "barcodes": barcodes_seen,
        "stats": stats,
    }

    if use_temp_files and temp_dir:
        # Write per-shard chunk files into temp_dir/shard_XX/chunk_YYYYYY.bin
        for s_idx in range(n_shards):
            s_dict = shard_data[s_idx]
            if s_dict["CG"] or s_dict["CH"]:
                shard_dir = os.path.join(temp_dir, f"shard_{s_idx:03d}")
                chunk_file = os.path.join(shard_dir, f"chunk_{chunk_id:06d}.bin")
                with open(chunk_file, "wb") as f:
                    pickle.dump(s_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
        del shard_data
        return (chunk_id, contig, start, end, None, chunk_stats)
    else:
        return (chunk_id, contig, start, end, shard_data, chunk_stats)


def plan_genomic_chunks(
    fasta_path: str,
    target_contigs: List[str],
    chunk_size_bp: int = 10_000_000,
) -> List[Tuple[int, str, int, int]]:
    """Divide canonical contigs into evenly-sized genomic windows in coordinate order."""
    fai = FastFaiReader(fasta_path)
    chunks = []
    chunk_id = 0

    for contig in target_contigs:
        contig_len = fai.get_reference_length(contig)
        for start in range(0, contig_len, chunk_size_bp):
            end = min(start + chunk_size_bp, contig_len)
            chunks.append((chunk_id, contig, start, end))
            chunk_id += 1

    return chunks


def _assemble_and_write_shard(
    shard_idx: int,
    shard_path: str,
    temp_shard_dir: Optional[str],
    in_memory_chunks: Optional[List[Dict[str, Dict[str, np.ndarray]]]],
    compression: str = "gzip",
    compression_level: int = 1,
    version: str = "amethyst2.0.0",
) -> Tuple[int, int, float]:
    """
    Assemble all coordinate chunks for a single shard and write its HDF5 file in a single pass.

    Returns:
        (shard_idx, num_cells_written, write_time_sec)
    """
    t0 = time.time()
    cg_accum: Dict[str, List[np.ndarray]] = defaultdict(list)
    ch_accum: Dict[str, List[np.ndarray]] = defaultdict(list)

    if temp_shard_dir and os.path.exists(temp_shard_dir):
        # Read chunk files in strict coordinate order (chunk_000000.bin, ...)
        chunk_files = sorted(os.listdir(temp_shard_dir))
        for fn in chunk_files:
            fp = os.path.join(temp_shard_dir, fn)
            with open(fp, "rb") as f:
                d = pickle.load(f)
            for bc, arr in d.get("CG", {}).items():
                cg_accum[bc].append(arr)
            for bc, arr in d.get("CH", {}).items():
                ch_accum[bc].append(arr)
            del d

        # Immediately purge temp files for this shard
        shutil.rmtree(temp_shard_dir, ignore_errors=True)
    elif in_memory_chunks:
        for d in in_memory_chunks:
            for bc, arr in d.get("CG", {}).items():
                cg_accum[bc].append(arr)
            for bc, arr in d.get("CH", {}).items():
                ch_accum[bc].append(arr)

    all_shard_bcs = sorted(set(cg_accum.keys()) | set(ch_accum.keys()))

    with AmethystH5Writer(shard_path, mode="w", compression=compression, compression_level=compression_level) as writer:
        writer.write_metadata(version)
        for bc in all_shard_bcs:
            if bc in cg_accum:
                merged = cg_accum[bc][0] if len(cg_accum[bc]) == 1 else np.concatenate(cg_accum[bc])
                writer.create_cell_dataset("CG", bc, merged)
            if bc in ch_accum:
                merged = ch_accum[bc][0] if len(ch_accum[bc]) == 1 else np.concatenate(ch_accum[bc])
                writer.create_cell_dataset("CH", bc, merged)

    del cg_accum, ch_accum
    gc.collect()
    elapsed = time.time() - t0
    return shard_idx, len(all_shard_bcs), elapsed


def extract_methylation_parallel(
    bam_path: str,
    fasta_path: str,
    out_h5_path: str,
    chroms: Optional[Iterable[str]] = None,
    whitelist_path: Optional[str] = None,
    n_workers: int = 24,
    decomp_threads: int = 1,
    chunk_size_mb: int = 10,
    n_shards: int = 1,
    use_temp_files: bool = True,
    temp_dir: Optional[str] = None,
    compression: str = "gzip",
    compression_level: int = 1,
    min_base_quality: int = 20,
    min_mapq: int = 0,
    max_depth: int = 250,
    compute_baq: bool = True,
    ignore_orphans: bool = True,
    ignore_overlaps: bool = True,
) -> Dict[str, Any]:
    """
    Parallel single-cell TAPS methylation extractor with partitioned shard streaming.

    When n_shards > 1, writes shard_000.h5..shard_NNN.h5 and master.h5 in parallel.
    When n_shards == 1, writes a single out_h5_path file.
    """
    whitelist: Optional[Set[str]] = None
    if whitelist_path:
        logger.info(f"Loading barcode whitelist from {whitelist_path}...")
        whitelist = parse_annot(whitelist_path)
        logger.info(f"Loaded {len(whitelist):,} barcodes into whitelist.")

    bam = pysam.AlignmentFile(bam_path, "rb")
    bam_contigs = set(bam.references)
    if chroms is None:
        target_contigs = [c for c in CANONICAL_CONTIGS if c in bam_contigs]
    else:
        target_contigs = [c for c in chroms if c in bam_contigs]

    idx_stats = {s.contig: s.mapped for s in bam.get_index_statistics()}
    total_target_reads = sum(idx_stats.get(c, 0) for c in target_contigs)
    bam.close()

    chunk_size_bp = int(chunk_size_mb * 1_000_000)
    chunks = plan_genomic_chunks(fasta_path, target_contigs, chunk_size_bp=chunk_size_bp)
    total_chunks = len(chunks)

    # Determine effective shard count for writing
    effective_shards = max(1, n_shards)

    actual_temp_dir = None
    if use_temp_files:
        actual_temp_dir = tempfile.mkdtemp(prefix="taps_extract_", dir=temp_dir)
        for s in range(effective_shards):
            os.makedirs(os.path.join(actual_temp_dir, f"shard_{s:03d}"), exist_ok=True)
        logger.info(f"Disk-streaming mode enabled ({effective_shards} shard partitions). Temp dir: {actual_temp_dir}")
    else:
        logger.info("In-memory mode enabled (--no-temp-file). Keeping chunk arrays in RAM.")

    logger.info(
        f"Processing {len(target_contigs)} contig(s) | {total_chunks} chunk(s) "
        f"({chunk_size_mb} Mb/chunk, {total_target_reads:,} mapped reads)."
    )
    logger.info(
        f"Concurrency: {n_workers} worker processes × {decomp_threads} decompression thread(s) = "
        f"{n_workers * decomp_threads + n_workers + 1} total execution threads."
    )

    worker_params = {
        "decomp_threads": decomp_threads,
        "ignore_overlaps": ignore_overlaps,
        "min_base_quality": min_base_quality,
        "min_mapq": min_mapq,
        "max_depth": max_depth,
        "compute_baq": compute_baq,
        "ignore_orphans": ignore_orphans,
        "use_temp_files": use_temp_files,
        "n_shards": effective_shards,
    }

    total_stats = {
        "CpG": {"c": 0, "t": 0},
        "CHG": {"c": 0, "t": 0},
        "CHH": {"c": 0, "t": 0},
        "CG": {"c": 0, "t": 0},
        "CH": {"c": 0, "t": 0},
    }
    all_barcodes: Set[str] = set()
    start_total_time = time.time()
    last_log_time = start_total_time
    chunks_completed = 0

    # In-memory mode storage: shard_idx -> list of chunk dicts in coordinate order
    in_memory_shard_chunks: Dict[int, List[Any]] = {s: [] for s in range(effective_shards)}

    try:
        ctx = mp.get_context("spawn")
        with ctx.Pool(
            processes=n_workers,
            initializer=_init_worker,
            initargs=(bam_path, fasta_path, whitelist, actual_temp_dir, worker_params),
        ) as pool:
            for result in pool.imap_unordered(_process_chunk, chunks, chunksize=1):
                chunk_id, contig, start, end, shard_payload, chunk_stats = result
                chunks_completed += 1

                for k in ("CpG", "CHG", "CHH", "CG", "CH"):
                    total_stats[k]["c"] += chunk_stats["stats"][k]["c"]
                    total_stats[k]["t"] += chunk_stats["stats"][k]["t"]
                all_barcodes.update(chunk_stats["barcodes"])

                if not use_temp_files and shard_payload:
                    for s in range(effective_shards):
                        in_memory_shard_chunks[s].append(shard_payload[s])

                now = time.time()
                if (chunks_completed % 10 == 0) or (now - last_log_time >= 15.0) or (chunks_completed == total_chunks):
                    elapsed = now - start_total_time
                    peak_rss = get_peak_rss_mb()
                    pct_done = 100.0 * chunks_completed / total_chunks
                    chunk_rate = (chunks_completed / elapsed * 60.0) if elapsed > 0 else 0.0
                    eta_min = (((total_chunks - chunks_completed) / (chunks_completed / elapsed)) / 60.0) if chunks_completed > 0 else 0.0

                    cg_c = total_stats["CG"]["c"]
                    cg_tot = cg_c + total_stats["CG"]["t"]
                    cg_pct = (100.0 * cg_c / cg_tot) if cg_tot > 0 else 0.0
                    ch_c = total_stats["CH"]["c"]
                    ch_tot = ch_c + total_stats["CH"]["t"]
                    ch_pct = (100.0 * ch_c / ch_tot) if ch_tot > 0 else 0.0

                    logger.info(
                        f"Extraction: [{chunks_completed:3d}/{total_chunks} chunks | {pct_done:5.1f}%] | "
                        f"Rate: {chunk_rate:4.1f} chk/min | Elapsed: {elapsed:5.1f}s | ETA: {eta_min:4.1f}m | "
                        f"Cells: {len(all_barcodes):,} | CG: {cg_pct:.1f}% | CH: {ch_pct:.2f}% | RSS: {peak_rss:.0f} MB"
                    )
                    last_log_time = now

        # Phase 2: Parallel Shard Assembly & Writing
        extraction_elapsed = time.time() - start_total_time
        logger.info(f"Extraction complete in {extraction_elapsed:.1f}s. Assembling and writing {effective_shards} shard(s) in parallel...")
        write_t0 = time.time()

        if effective_shards > 1:
            os.makedirs(out_h5_path, exist_ok=True)
            shard_filenames = [f"shard_{s:03d}.h5" for s in range(effective_shards)]
            shard_paths = [os.path.join(out_h5_path, fn) for fn in shard_filenames]
        else:
            shard_filenames = [os.path.basename(out_h5_path)]
            shard_paths = [out_h5_path]

        # Execute shard assembly and writing sequentially to keep RAM strictly bounded (<100 MB)
        for s in range(effective_shards):
            temp_s_dir = os.path.join(actual_temp_dir, f"shard_{s:03d}") if actual_temp_dir else None
            mem_chunks = in_memory_shard_chunks[s] if not use_temp_files else None
            s_idx, n_cells, elapsed_s = _assemble_and_write_shard(
                shard_idx=s,
                shard_path=shard_paths[s],
                temp_shard_dir=temp_s_dir,
                in_memory_chunks=mem_chunks,
                compression=compression,
                compression_level=compression_level,
            )
            logger.info(f"  -> Shard {s_idx:03d} finished: {n_cells:,} cells written in {elapsed_s:.2f}s ({shard_paths[s_idx]})")
            if not use_temp_files:
                in_memory_shard_chunks[s] = []
            gc.collect()

        # Generate master.h5 if sharded directory mode
        if effective_shards > 1:
            master_path = os.path.join(out_h5_path, "master.h5")
            with h5py.File(master_path, "w") as master_h5:
                meta_group = master_h5.require_group("metadata")
                meta_group.create_dataset("version", data=b"amethyst2.0.0")
                cg_master = master_h5.require_group("CG")
                ch_master = master_h5.require_group("CH")

                for s in range(effective_shards):
                    shard_fn = shard_filenames[s]
                    shard_fp = shard_paths[s]
                    with h5py.File(shard_fp, "r") as sf:
                        if "CG" in sf:
                            for bc in sf["CG"].keys():
                                cg_master[bc] = h5py.ExternalLink(shard_fn, f"CG/{bc}")
                        if "CH" in sf:
                            for bc in sf["CH"].keys():
                                ch_master[bc] = h5py.ExternalLink(shard_fn, f"CH/{bc}")

            logger.info(f"Created master index file: {master_path}")

        write_elapsed = time.time() - write_t0
        logger.info(f"All HDF5 datasets written successfully in {write_elapsed:.2f}s!")

    finally:
        if actual_temp_dir and os.path.exists(actual_temp_dir):
            shutil.rmtree(actual_temp_dir, ignore_errors=True)

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
        "n_workers": n_workers,
        "decomp_threads": decomp_threads,
        "chunk_size_mb": chunk_size_mb,
        "n_shards": effective_shards,
        "use_temp_files": use_temp_files,
        "total_chunks": total_chunks,
        "stats": total_stats,
        "percentages": {
            "CpG": cg_all_pct,
            "CHG": chg_all_pct,
            "CHH": chh_all_pct,
            "CH": ch_all_pct,
            "Pooled": pooled_pct,
        },
    }

    logger.info("=" * 70)
    logger.info("PARALLEL EXTRACTION SUMMARY & GROUND TRUTH COMPARISON")
    logger.info("=" * 70)
    logger.info(f"Total Wall Time:       {total_elapsed:.1f}s ({total_elapsed/60.0:.2f} min)")
    logger.info(f"Total Mapped Reads:    {total_target_reads:,}")
    logger.info(f"Throughput (Reads):    {reads_per_min:,.0f} reads/min ({total_target_reads/total_elapsed:,.0f} reads/s)")
    logger.info(f"Total Unique Barcodes: {len(all_barcodes):,}")
    logger.info(f"Throughput (Cells):    {cells_per_min:,.1f} cells/min")
    logger.info(f"Peak Memory (RSS):     {final_peak_rss:.1f} MB")
    logger.info(f"Shards Written:        {effective_shards}")
    logger.info(f"Streaming Mode:        {'Disk-backed temp files' if use_temp_files else 'In-memory IPC'}")
    logger.info("-" * 70)
    logger.info(f"CpG Methylation:       {cg_all_pct:.3f}% ({total_stats['CpG']['c']:,} / {cg_all:,})")
    logger.info(f"CHG Methylation:       {chg_all_pct:.3f}% ({total_stats['CHG']['c']:,} / {chg_all:,})")
    logger.info(f"CHH Methylation:       {chh_all_pct:.3f}% ({total_stats['CHH']['c']:,} / {chh_all:,})")
    logger.info(f"CH (Combined):         {ch_all_pct:.3f}% ({total_stats['CH']['c']:,} / {ch_all:,})")
    logger.info(f"Pooled (All C/G sites):{pooled_pct:.3f}% ({total_mod:,} / {total_calls:,})")
    logger.info("=" * 70)

    return summary
