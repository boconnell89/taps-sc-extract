"""
Command-line interface for taps_sc_extract.
"""

import argparse
import logging
import os
import shutil
import subprocess
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
        choices=["gzip", "lzf", "none", "blosc", "blosc-zstd"],
        default="gzip",
        help=(
            "HDF5 dataset compression (default: gzip). "
            "gzip is portable to Amethyst/rhdf5 with no extra R packages. "
            "blosc uses multithreaded LZ4 (fastest compressed write); "
            "blosc-zstd is nearly as fast with gzip-like size. Both need "
            "Bioconductor rhdf5filters to load in R. "
            "lzf is a fast single-thread filter; none disables compression."
        ),
    )
    parser.add_argument(
        "--compression-threads",
        type=int,
        default=None,
        help=(
            "Blosc worker threads per compress call (default: CPU count / shard "
            "writer count). Ignored unless --compression is blosc or blosc-zstd."
        ),
    )
    parser.add_argument(
        "--max-writer-threads",
        "--writer-threads",
        type=int,
        default=6,
        help=(
            "Maximum parallel shard-writer threads after extraction (default: 6). "
            "Recommended: 6 for local NVMe/SSD, 2-4 for HDD or network filesystems (NFS/Lustre), "
            "1 for strict sequential writing."
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
        "--engine",
        choices=["auto", "rust", "python"],
        default="auto",
        help=(
            "Extraction backend engine (default: auto). "
            "'rust' uses the high-performance taps-sc-extract-rs core (2.6x faster, 56%% less RAM); "
            "'python' runs the reference multiprocessing engine; "
            "'auto' uses rust if available, falling back to python."
        ),
    )
    parser.add_argument(
        "--memory-mode",
        choices=["auto", "stream", "memory"],
        default=None,
        help="Memory mode: stream (disk temporary files), memory (RAM map), or auto (heuristic).",
    )
    parser.add_argument(
        "--max-memory-gb",
        type=float,
        default=None,
        help="Optional memory budget in GB for auto-tuning.",
    )
    parser.add_argument(
        "--expected-cells",
        type=int,
        default=None,
        help="Optional expected cell count (auto-detected from whitelist if provided).",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser.parse_args(args)


def find_rust_binary() -> Optional[str]:
    """Locate the taps-sc-extract-rs executable."""
    # 1. Check PATH
    import shutil
    path = shutil.which("taps-sc-extract-rs")
    if path:
        return path
    # 2. Check local workspace target/release build
    pkg_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    candidate = os.path.join(pkg_dir, "taps-sc-extract-rs", "target", "release", "taps-sc-extract-rs")
    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
        return candidate
    return None


def main(args: Optional[List[str]] = None) -> int:
    parsed = parse_args(args)
    setup_logging(verbose=parsed.verbose, log_file=parsed.log_file)
    logger = logging.getLogger("taps_sc_extract")

    logger.info("=" * 70)
    logger.info(f"taps-sc-extract v{__version__} - Single-Cell TAPS Extractor")
    logger.info(f"Command line: {' '.join(sys.argv)}")
    logger.info(f"Configuration: {vars(parsed)}")
    logger.info("=" * 70)

    # Resolve engine
    rust_bin = find_rust_binary()
    use_rust = False
    if parsed.engine == "rust":
        if not rust_bin:
            logger.error("Requested --engine rust but 'taps-sc-extract-rs' binary was not found.")
            return 1
        use_rust = True
    elif parsed.engine == "auto":
        use_rust = (rust_bin is not None)

    if use_rust and rust_bin:
        logger.info(f"Delegating extraction to Rust acceleration core: {rust_bin}")
        import subprocess

        cmd = [
            rust_bin,
            "extract",
            "-b", parsed.bam,
            "-f", parsed.fasta,
            "-o", parsed.out,
            "-t", str(parsed.workers),
            "--chunk-size-mb", str(parsed.chunk_size_mb),
            "--shards", str(parsed.shards),
            "--compression", parsed.compression,
            "--max-writer-threads", str(parsed.max_writer_threads),
            "--min-baseq", str(parsed.min_baseq),
            "--min-mapq", str(parsed.min_mapq),
            "--max-depth", str(parsed.max_depth),
        ]
        if parsed.chroms:
            cmd.extend(["-c", parsed.chroms])
        if parsed.whitelist:
            cmd.extend(["-w", parsed.whitelist])
        if parsed.decomp_threads is not None:
            cmd.extend(["--decomp-threads", str(parsed.decomp_threads)])
        if parsed.no_baq:
            cmd.append("--no-baq")
        if parsed.no_overlap_clip:
            cmd.append("--no-overlap-clip")
        if parsed.max_memory_gb:
            cmd.extend(["--max-memory-gb", str(parsed.max_memory_gb)])
        if parsed.expected_cells:
            cmd.extend(["--expected-cells", str(parsed.expected_cells)])
        if parsed.memory_mode:
            cmd.extend(["--memory-mode", parsed.memory_mode])
        elif parsed.no_temp_file:
            cmd.extend(["--memory-mode", "memory"])

        try:
            res = subprocess.run(cmd)
            return res.returncode
        except Exception as e:
            logger.exception(f"Error executing Rust backend: {e}")
            return 1

    logger.info("Running Python reference multiprocessing extraction engine.")
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
            decomp_threads=parsed.decomp_threads or 1,
            chunk_size_mb=parsed.chunk_size_mb,
            n_shards=parsed.shards,
            use_temp_files=(parsed.memory_mode != "memory" and not parsed.no_temp_file),
            compression=parsed.compression,
            compression_threads=parsed.compression_threads,
            max_writer_threads=parsed.max_writer_threads,
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
        logger.exception(f"Error during extraction: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
