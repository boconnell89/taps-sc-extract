"""
Command-line interface for taps_sc_extract.
"""

import argparse
import logging
import sys
from typing import List, Optional

from . import __version__
from .parallel_extractor import CANONICAL_CONTIGS, extract_methylation_parallel


def setup_logging(verbose: bool = False, log_file: Optional[str] = None):
    """Configure detailed timestamped logging to console and optional log file."""
    level = logging.DEBUG if verbose else logging.INFO
    formatter = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    # Clear existing handlers
    root_logger.handlers.clear()

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    # File handler (if specified)
    if log_file:
        file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)
        logging.getLogger("taps_sc_extract").info(f"Detailed logging to: {log_file}")


def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="taps-sc-extract",
        description="Extract single-cell DNA methylation from TAPS BAM into Amethyst HDF5.",
    )
    parser.add_argument(
        "-b", "--bam",
        required=True,
        help="Path to coordinate-sorted, indexed TAPS BAM file (.bam with .bai).",
    )
    parser.add_argument(
        "-f", "--fasta",
        required=True,
        help="Path to reference FASTA file (indexed with .fai).",
    )
    parser.add_argument(
        "-o", "--out",
        required=True,
        help="Output Amethyst-compatible HDF5 file path (.h5).",
    )
    parser.add_argument(
        "-c", "--chroms",
        default=None,
        help=(
            "Comma-separated list of contigs to process (e.g. 'chr19' or 'chr1,chr2'). "
            "Default: canonical mm10 chromosomes (chr1..chr19, chrX, chrY; chrM excluded)."
        ),
    )
    parser.add_argument(
        "-w", "--whitelist", "-a", "--annot",
        dest="whitelist",
        default=None,
        help="Path to optional barcode whitelist or annotation file (first column = barcode).",
    )
    parser.add_argument(
        "-t", "--threads", "--workers",
        dest="workers",
        type=int,
        default=24,
        help="Number of parallel chunk worker processes (default: 24).",
    )
    parser.add_argument(
        "--decomp-threads",
        type=int,
        default=1,
        help="Number of BAM decompression threads per worker process (default: 1).",
    )
    parser.add_argument(
        "--chunk-size-mb",
        type=int,
        default=10,
        help="Genomic chunk size in megabases (default: 10).",
    )
    parser.add_argument(
        "--shards",
        type=int,
        default=1,
        help=(
            "Number of output HDF5 shard files to write in parallel (default: 1). "
            "If >1, -o/--out is treated as a directory where shard_000.h5..shard_NNN.h5 and "
            "a master.h5 (with relative ExternalLinks) will be created."
        ),
    )
    parser.add_argument(
        "--compression",
        choices=["gzip", "lzf", "none"],
        default="gzip",
        help=(
            "HDF5 dataset compression (default: gzip). "
            "gzip is portable to Amethyst/rhdf5 with no extra R packages. "
            "lzf writes much faster but R needs Bioconductor rhdf5filters. "
            "none is fastest and produces the largest files."
        ),
    )
    parser.add_argument(
        "--no-temp-file",
        action="store_true",
        help=(
            "Keep all chunk results directly in memory instead of streaming temporary binary "
            "files to disk. Recommended for high-memory production servers (e.g. >=64 GB RAM)."
        ),
    )
    parser.add_argument(
        "--temp-dir",
        default=None,
        help="Custom directory for temporary chunk streaming (default: system temp /tmp).",
    )
    parser.add_argument(
        "--log-file",
        default=None,
        help="Path to output log file for recording detailed timestamps, progress, and performance metrics.",
    )
    parser.add_argument(
        "--min-baseq",
        type=int,
        default=20,
        help="Minimum base quality for pileup (default: 20).",
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum mapping quality (default: 0).",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=250,
        help="Maximum pileup depth (default: 250).",
    )
    parser.add_argument(
        "--no-baq",
        action="store_true",
        help="Disable BAQ (Base Alignment Quality) computation.",
    )
    parser.add_argument(
        "--no-overlap-clip",
        action="store_true",
        help="Do not ignore overlapping mate read bases.",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser.parse_args(args)


def main(args: Optional[List[str]] = None) -> int:
    parsed = parse_args(args)
    setup_logging(verbose=parsed.verbose, log_file=parsed.log_file)
    logger = logging.getLogger("taps_sc_extract")

    logger.info("=" * 70)
    logger.info(f"taps-sc-extract v{__version__} - Single-Cell TAPS Extractor")
    logger.info(f"Command line: {' '.join(sys.argv)}")
    logger.info(f"Configuration: {vars(parsed)}")
    logger.info("=" * 70)

    chroms_list = None
    if parsed.chroms:
        chroms_list = [c.strip() for c in parsed.chroms.split(",") if c.strip()]

    try:
        extract_methylation_parallel(
            bam_path=parsed.bam,
            fasta_path=parsed.fasta,
            out_h5_path=parsed.out,
            chroms=chroms_list,
            whitelist_path=parsed.whitelist,
            n_workers=parsed.workers,
            decomp_threads=parsed.decomp_threads,
            chunk_size_mb=parsed.chunk_size_mb,
            n_shards=parsed.shards,
            use_temp_files=not parsed.no_temp_file,
            compression=parsed.compression,
            temp_dir=parsed.temp_dir,
            min_base_quality=parsed.min_baseq,
            min_mapq=parsed.min_mapq,
            max_depth=parsed.max_depth,
            compute_baq=not parsed.no_baq,
            ignore_orphans=True,
            ignore_overlaps=not parsed.no_overlap_clip,
        )
        return 0
    except Exception as e:
        logging.getLogger("taps_sc_extract").exception(f"Error during extraction: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
