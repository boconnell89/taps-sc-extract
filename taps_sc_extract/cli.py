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
        default=None,
        help=(
            "Path to coordinate-sorted, indexed TAPS BAM file (.bam with .bai). "
            "Required unless --test runs without a dataset."
        ),
    )
    parser.add_argument(
        "-f", "--fasta",
        default=None,
        help="Path to reference FASTA file (indexed with .fai).",
    )
    parser.add_argument(
        "-o", "--out",
        default=None,
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
        default=None,
        help=(
            "BAM decompression threads per worker (default: Rust picks 4/2/0 by worker "
            "count; the Python engine uses 1)."
        ),
    )
    parser.add_argument(
        "--chunk-size-mb",
        type=int,
        default=10,
        help=(
            "Genomic chunk size in megabases "
            "(default: 10 for the Python engine; Rust sizes windows from the worker count)."
        ),
    )
    parser.add_argument(
        "--shards",
        type=int,
        default=1,
        help=(
            "Number of output HDF5 shard files to write in parallel (default: 1). "
            "If >1, -o/--out is a directory of shard_000.h5..shard_NNN.h5 plus master.h5. "
            "An .h5 file name is written as a directory named after the file stem "
            "(sample.h5 -> sample/)."
        ),
    )
    parser.add_argument(
        "--compression",
        choices=["gzip", "gzip-shuffle", "gzip6", "lzf", "none", "blosc", "blosc-zstd"],
        default="gzip",
        help=(
            "HDF5 dataset compression (default: gzip for Python, lzf for Rust). "
            "gzip (level 1 deflate) is fast and 100%% portable to Amethyst/rhdf5 with no extra packages. "
            "gzip-shuffle applies byte shuffling before level 1 deflate for high ratio and speed. "
            "gzip6 uses standard level 6 deflate. "
            "none disables compression for maximum write speed. "
            "blosc / blosc-zstd / lzf require rhdf5filters in R."
        ),
    )
    parser.add_argument(
        "--compression-threads",
        type=int,
        default=None,
        help=(
            "Rust engine: Blosc threads per shard writer (default 1). Only used with "
            "--compression blosc or blosc-zstd; more threads measured no gain."
        ),
    )
    parser.add_argument(
        "--max-writer-threads",
        "--writer-threads",
        type=int,
        default=0,
        help=(
            "Maximum concurrent shard writers after extraction (default: 0 = auto). "
            "Auto scales up to shard count, available CPU cores, and memory budget. "
            "Suggested: 8-16 for NVMe/SSD, 2-4 for HDD or network filesystems (NFS/Lustre)."
        ),
    )
    parser.add_argument(
        "--profile-json",
        default=None,
        help=(
            "Rust engine: write per-phase CPU use, per-stage and per-window timings, "
            "and shard-writer memory to this JSON file."
        ),
    )
    parser.add_argument(
        "--writer-mode",
        choices=["auto", "process", "thread"],
        default="auto",
        help=(
            "Rust shard writers: process (one child process per shard, each with its own "
            "HDF5 lock), thread (in-process pool sharing one lock), or auto (process when "
            "there is more than one shard)."
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
        "--trim-r2",
        type=int,
        default=10,
        help=(
            "Drop this many bases from the 5' end of read 2 before calling "
            "(default: 10). Removes Tn5/TAPS M-bias on R2."
        ),
    )
    parser.add_argument(
        "--trim-r1",
        type=int,
        default=0,
        help="Drop this many bases from the 5' end of read 1 before calling (default: 0).",
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
            "'rust' uses the taps-sc-extract-rs core; "
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
        help="Memory budget in GB for the whole run (overrides --max-memory-fraction).",
    )
    parser.add_argument(
        "--max-memory-fraction",
        type=float,
        default=None,
        help=(
            "Rust engine: share of available memory the run may use, extraction plus "
            "writer processes (default: 0.6)."
        ),
    )
    parser.add_argument(
        "--expected-cells",
        type=int,
        default=None,
        help="Optional expected cell count (auto-detected from whitelist if provided).",
    )
    test = parser.add_argument_group(
        "tuning (--test)",
        "Measure this machine with a built-in synthetic workload (no data needed) and "
        "recommend -t, --decomp-threads and writer settings. Add -b/-f to also sample "
        "that dataset and predict run time and memory. Needs the Rust engine.",
    )
    test.add_argument("--test", action="store_true", help="Run the tuning probes instead of extracting.")
    test.add_argument("--test-quick", action="store_true", help="Smaller probes (~40 s instead of a few minutes).")
    test.add_argument(
        "--save-profile",
        action="store_true",
        help=(
            "Save the machine profile to ~/.config/taps-sc-extract/machine_profile.toml; "
            "later Rust runs on this host take default -t / writer counts from it."
        ),
    )
    test.add_argument("--test-report", default=None, help="Also write the profile and report to this file.")
    test.add_argument(
        "--temp-dir-test",
        dest="test_temp_dir",
        default=None,
        help="Directory for the synthetic data and temp-filesystem probe (default: --temp-dir or current directory).",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    ns = parser.parse_args(args)
    if not ns.test:
        missing = [f for f, v in (("-b/--bam", ns.bam), ("-f/--fasta", ns.fasta), ("-o/--out", ns.out)) if not v]
        if missing:
            parser.error("the following arguments are required: " + ", ".join(missing))
    elif bool(ns.bam) != bool(ns.fasta):
        parser.error("--test with a dataset needs both -b and -f")
    return ns


RUST_DEFAULT_CHUNK_SIZE_MB = 0  # auto: sized from worker count
RUST_DEFAULT_COMPRESSION = "lzf"


def option_was_set(argv: List[str], option: str) -> bool:
    """True when the user passed `option` or `option=value` on this argv."""
    return any(a == option or a.startswith(option + "=") for a in argv)


def build_rust_command(rust_bin: str, parsed, argv: List[str]) -> List[str]:
    """Argv for taps-sc-extract-rs. Rust defaults apply unless the user set the flag."""
    chunk_mb = (
        parsed.chunk_size_mb
        if option_was_set(argv, "--chunk-size-mb")
        else RUST_DEFAULT_CHUNK_SIZE_MB
    )
    compression = (
        parsed.compression
        if option_was_set(argv, "--compression")
        else RUST_DEFAULT_COMPRESSION
    )
    cmd = [
        rust_bin,
        "extract",
        "-b", parsed.bam,
        "-f", parsed.fasta,
        "-o", parsed.out,
        "-t", str(parsed.workers),
        "--chunk-size-mb", str(chunk_mb),
        "--shards", str(parsed.shards),
        "--compression", compression,
        "--max-writer-threads", str(parsed.max_writer_threads),
        "--writer-mode", parsed.writer_mode,
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
    if parsed.profile_json:
        cmd.extend(["--profile-json", parsed.profile_json])
    if parsed.compression_threads is not None:
        cmd.extend(["--compression-threads", str(parsed.compression_threads)])
    if parsed.max_memory_gb:
        cmd.extend(["--max-memory-gb", str(parsed.max_memory_gb)])
    if parsed.max_memory_fraction:
        cmd.extend(["--max-memory-fraction", str(parsed.max_memory_fraction)])
    if parsed.expected_cells:
        cmd.extend(["--expected-cells", str(parsed.expected_cells)])
    if parsed.memory_mode:
        cmd.extend(["--memory-mode", parsed.memory_mode])
    elif parsed.no_temp_file:
        cmd.extend(["--memory-mode", "memory"])
    cmd.extend(["--trim-r1", str(parsed.trim_r1), "--trim-r2", str(parsed.trim_r2)])
    return cmd


def build_tune_command(rust_bin: str, parsed) -> List[str]:
    """Argv for `taps-sc-extract-rs tune` from `taps-sc-extract --test ...`."""
    temp_dir = parsed.test_temp_dir or parsed.temp_dir or "."
    out_dir = "."
    if parsed.out:
        out_dir = parsed.out if os.path.isdir(parsed.out) else (os.path.dirname(parsed.out) or ".")
    cmd = [rust_bin, "tune", "--temp-dir", temp_dir, "--out-dir", out_dir]
    if parsed.bam:
        cmd.extend(["-b", parsed.bam, "-f", parsed.fasta])
    if parsed.test_quick:
        cmd.append("--quick")
    if parsed.save_profile:
        cmd.append("--save-profile")
    if parsed.test_report:
        cmd.extend(["--report", parsed.test_report])
    return cmd


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
    argv = list(sys.argv[1:] if args is None else args)
    parsed = parse_args(argv)
    setup_logging(verbose=parsed.verbose, log_file=parsed.log_file)
    logger = logging.getLogger("taps_sc_extract")
    if parsed.trim_r1 < 0 or parsed.trim_r2 < 0:
        logger.error("--trim-r1 and --trim-r2 must be >= 0.")
        return 1

    logger.info("=" * 70)
    logger.info(f"taps-sc-extract v{__version__} - Single-Cell TAPS Extractor")
    logger.info(f"Command line: {' '.join(sys.argv)}")
    logger.info(f"Configuration: {vars(parsed)}")
    logger.info("=" * 70)

    # Resolve engine
    rust_bin = find_rust_binary()
    if parsed.test:
        if not rust_bin:
            logger.error("--test needs the Rust engine ('taps-sc-extract-rs'); run ./install.sh first.")
            return 1
        import subprocess

        return subprocess.run(build_tune_command(rust_bin, parsed)).returncode
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

        cmd = build_rust_command(rust_bin, parsed, argv)

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
            decomp_threads=parsed.decomp_threads if parsed.decomp_threads is not None else 1,
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
            trim_r1=parsed.trim_r1,
            trim_r2=parsed.trim_r2,
        )
        return 0
    except Exception as e:
        logger.exception(f"Error during extraction: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
