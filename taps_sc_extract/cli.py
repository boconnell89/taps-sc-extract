"""
Command-line interface for taps_sc_extract.
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
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
        # The log usually sits in the output directory, which may not exist yet.
        log_dir = os.path.dirname(os.path.abspath(log_file))
        os.makedirs(log_dir, exist_ok=True)
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
        nargs="+",
        default=None,
        help=(
            "Coordinate-sorted alignments: BAM, CRAM (decoded with -f) or SAM (converted "
            "to a temporary BAM). A missing .bai/.csi/.crai is built automatically. Several "
            "inputs: list them, or give @FILE naming a list (one path per line) or a sample "
            "sheet (.tsv or .xlsx; header row with a 'bam' column and optional complexity, "
            "min_complexity, max_complexity, whitelist per file). A barcode must belong to "
            "one input; reads are merged per window. Required unless --test runs without a dataset."
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
            "Default: the primary chromosomes in the BAM header (autosomes plus X/Y, or Z/W; chrM, unplaced scaffolds, alts and decoys excluded), e.g. chr1..chr22, chrX, chrY for hg38."
        ),
    )
    parser.add_argument(
        "-w", "--whitelist", "-a", "--annot",
        dest="whitelist",
        default=None,
        help="Path to optional barcode whitelist or annotation file (first column = barcode).",
    )
    parser.add_argument(
        "--complexity",
        default=None,
        help=(
            "Complexity file: rank, barcode, total reads, unique reads, percent (whitespace-separated). "
            "With --min-complexity / --max-complexity, keeps barcodes whose unique-read count "
            "(column 4) is in range. Combined with --whitelist, both must pass."
        ),
    )
    parser.add_argument(
        "--min-complexity", type=float, default=None,
        help="Keep cells with at least this many unique reads (needs --complexity).",
    )
    parser.add_argument(
        "--max-complexity", type=float, default=None,
        help="Keep cells with at most this many unique reads (needs --complexity).",
    )
    parser.add_argument(
        "-t", "--threads", "--workers",
        dest="workers",
        type=int,
        default=None,
        help="Number of parallel chunk workers (default: Rust sizes it from cores and memory; "
        "the Python engine uses 24).",
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
        default=None,
        help=(
            "Number of output HDF5 shard files to write in parallel (default: Rust picks it from "
            "the estimated output sites and cells, 1 when -o names an .h5 file; the Python engine uses 1). "
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
            "HDF5 dataset compression (default: gzip for Python, blosc for Rust). "
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
        default=None,
        help="Minimum base quality for pileup (default: 20).",
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=None,
        help="Minimum mapping quality (default: taps 0, bisulfite 10; BSBolt gives 0 to reads "
        "that are not bisulfite-unique).",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=None,
        help="htslib pileup cap, all cells together (default: 0 = unlimited). A cap thins deep "
        "positions by read order; --site-depth-limit drops them instead.",
    )
    parser.add_argument(
        "--site-depth-limit",
        type=_depth_arg,
        default=None,
        help="Drop every call at positions covered by more kept-cell reads than this: a number "
        "(0 = off) or auto (default) = --molecules-per-cell x reads per molecule (taps 2, "
        "bisulfite 1) x cells kept by --complexity (with its cutoff) / --whitelist; off without "
        "either. The pileup is unchanged.",
    )
    parser.add_argument(
        "--cell-site-cap",
        type=_depth_arg,
        default=None,
        help="Most observations one cell keeps at a site: a number (0 = off) or auto (default) = "
        "--molecules-per-cell / 2 (a site is called on one strand: one molecule per haplotype), "
        "x 2 for TAPS with --no-overlap-clip. Extra observations are scaled away keeping the t:c ratio. "
        "Relax it (raise --molecules-per-cell, or 0) when real extra copies are expected: copy-number "
        "gains, segmental duplications, aneuploidy (e.g. cancer), replicating cells.",
    )
    parser.add_argument(
        "--molecules-per-cell",
        type=int,
        default=None,
        help="Molecules one cell can have at a position, both strands (default 4 = 2 strands x 2 "
        "haplotypes; 8 allows a replicated, G2 genome). Sets --cell-site-cap auto and "
        "--site-depth-limit auto.",
    )
    parser.add_argument(
        "--blacklist",
        default=None,
        help="BED of regions whose calls are dropped for every cell (e.g. the ENCODE blacklist). "
        "The pileup is unchanged.",
    )
    parser.add_argument(
        "--no-baq",
        action="store_true",
        help="Disable BAQ (Base Alignment Quality) computation.",
    )
    parser.add_argument(
        "--baq",
        choices=["taps", "taps-full", "taps-partial", "htslib"],
        default=None,
        help="BAQ mode (Rust engine). The taps modes mask the methylation-convertible "
        "base (C on OT reads, G on OB reads) so C-to-T conversions are not treated as "
        "mismatches. 'taps' (default) skips reads that already match that reference "
        "exactly (~2x faster than 'taps-full', which runs every read); 'taps-partial' "
        "runs only reads with an indel or clip (~4x faster, 0.06%% of CG calls differ). "
        "'htslib' uses the plain reference like samtools/asTair, which lowers "
        "methylated calls.",
    )
    parser.add_argument(
        "--context",
        choices=["both", "CG", "CH"],
        type=lambda v: {"both": "both", "cg": "CG", "ch": "CH"}.get(v.lower(), v),
        default="both",
        help="Methylation contexts to extract (default: both). CG or CH skips the other "
        "context's pileup columns and leaves its HDF5 group (/CG or /CH) out.",
    )
    parser.add_argument(
        "--caller",
        choices=["reads", "pileup"],
        default=None,
        help="Rust engine: how calls are made. 'reads' (default) calls one read at a "
        "time, replaying htslib's pileup depth cap and mate-overlap rules; 'pileup' "
        "uses htslib's column pileup. Same calls; 'reads' is faster.",
    )
    parser.add_argument(
        "--baq-simd",
        choices=["auto", "avx512", "avx2", "off"],
        default=None,
        help="BAQ kernel (Rust engine): auto (default; widest the CPU runs, or the "
        "fastest measured by --test --save-profile), avx512 (8 reads per vector), avx2 "
        "(4), or off (htslib, one read at a time). All give identical results.",
    )
    parser.add_argument(
        "--baq-device",
        choices=["cpu", "gpu"],
        default=None,
        help="Where the BAQ HMM runs (Rust engine): cpu (default) or gpu. Results are "
        "identical; see --gpu-backend.",
    )
    parser.add_argument(
        "--gpu-backend",
        choices=["auto", "opencl", "cuda", "hip"],
        default=None,
        help="GPU backend for --baq-device gpu: auto (CUDA, then HIP, then OpenCL), opencl "
        "(Intel Arc / Data Center GPUs, any device with fp64), cuda (NVIDIA), hip (AMD).",
    )
    parser.add_argument(
        "--gpu-device",
        default=None,
        help="Use the first GPU whose name contains this text.",
    )
    parser.add_argument(
        "--gpu-mem-mb",
        type=int,
        default=None,
        help="GPU memory for BAQ HMM matrices per launch, in MB (default 1024).",
    )
    parser.add_argument(
        "--no-overlap-clip",
        action="store_true",
        help="Do not ignore overlapping mate read bases.",
    )
    parser.add_argument(
        "--chemistry",
        choices=["auto", "taps", "bisulfite"],
        default="auto",
        help="Conversion chemistry. auto (default) reads the BAM: BSBolt YS strand tags mean "
        "bisulfite (single-end, unmethylated C -> T); paired-end reads without them mean taps "
        "(mC -> T). The two use separate strand and call rules and different defaults.",
    )
    parser.add_argument(
        "--strand-source",
        choices=["auto", "flag"],
        default="auto",
        help="Bisulfite strand: auto (flag, checked against every read's BSBolt YS tag; the "
        "run stops at the first disagreement) or flag (no YS check, for aligners without YS).",
    )
    parser.add_argument(
        "--trim-r2",
        type=int,
        default=None,
        help=(
            "Drop this many bases from the 5' end of read 2 before calling "
            "(default: taps 10; no effect on single-end reads). Removes Tn5/TAPS M-bias on R2."
        ),
    )
    parser.add_argument(
        "--trim-r1", "--trim-5p",
        dest="trim_r1",
        type=int,
        default=None,
        help="Drop this many bases from the sequenced 5' end of read 1, or of every "
        "single-end read (default: 0).",
    )
    parser.add_argument(
        "--trim",
        choices=["auto", "manual"],
        default=None,
        help="auto: measure M-bias on a sample of reads and pick the 5' and 3' trims "
        "(writes <bam>.mbias.tsv/.png/.pdf next to the output); manual: use --trim-5p / "
        "--trim-3p / --trim-r2 or the chemistry defaults. Default: auto for bisulfite, "
        "manual for taps. Explicit trim values always win for their end.",
    )
    parser.add_argument(
        "--trim-3p", "--trim-3prime",
        dest="trim_3p",
        type=int,
        default=None,
        help="Drop this many bases from the sequenced 3' end of every read (default: taps 0, "
        "bisulfite 10 for the Tn5 gap filled with unmethylated C at end repair).",
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
            "writer processes (default: 0.85)."
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
RUST_DEFAULT_COMPRESSION = "blosc"


def _depth_arg(value: str):
    """--site-depth-limit: 'auto' or a non-negative read count."""
    if value.lower() == "auto":
        return "auto"
    try:
        n = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a read count or 'auto', got {value!r}")
    if n < 0:
        raise argparse.ArgumentTypeError("--site-depth-limit must be >= 0")
    return n


def run_tee(cmd: List[str], log_path: str) -> int:
    """Run the Rust engine, copying its stdout and stderr (progress, config,
    summary) into the --log-file while still printing them."""
    import threading

    for h in logging.getLogger().handlers:  # the wrapper's own log lines first
        h.flush()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lock = threading.Lock()
    with open(log_path, "ab") as log:
        def pump(src, dst):
            for line in iter(src.readline, b""):
                dst.write(line)
                dst.flush()
                with lock:
                    log.write(line)
                    log.flush()
        threads = [threading.Thread(target=pump, args=(proc.stdout, sys.stdout.buffer)),
                   threading.Thread(target=pump, args=(proc.stderr, sys.stderr.buffer))]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
    return proc.wait()


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
        "-b", *parsed.bam,
        "-f", parsed.fasta,
        "-o", parsed.out,
        "--chunk-size-mb", str(chunk_mb),
        "--compression", compression,
        "--max-writer-threads", str(parsed.max_writer_threads),
        "--writer-mode", parsed.writer_mode,
        "--chemistry", parsed.chemistry,
        "--strand-source", parsed.strand_source,
    ]
    # Forwarded only when set, so the Rust engine applies its own defaults
    # (auto workers and shards, the detected chemistry's filters and trims).
    for flag, value in (
        ("-t", parsed.workers),
        ("--shards", parsed.shards),
        ("--min-baseq", parsed.min_baseq),
        ("--min-mapq", parsed.min_mapq),
        ("--max-depth", parsed.max_depth),
        ("--site-depth-limit", parsed.site_depth_limit),
        ("--cell-site-cap", parsed.cell_site_cap),
        ("--molecules-per-cell", parsed.molecules_per_cell),
        ("--blacklist", parsed.blacklist),
        ("--temp-dir", parsed.temp_dir),
        ("--trim-r1", parsed.trim_r1),
        ("--trim-r2", parsed.trim_r2),
        ("--trim-3p", parsed.trim_3p),
        ("--trim", parsed.trim),
    ):
        if value is not None:
            cmd.extend([flag, str(value)])
    if parsed.chroms:
        cmd.extend(["-c", parsed.chroms])
    if parsed.whitelist:
        cmd.extend(["-w", parsed.whitelist])
    if parsed.complexity:
        cmd.extend(["--complexity", parsed.complexity])
    if parsed.min_complexity is not None:
        cmd.extend(["--min-complexity", f"{parsed.min_complexity:g}"])
    if parsed.max_complexity is not None:
        cmd.extend(["--max-complexity", f"{parsed.max_complexity:g}"])
    if parsed.decomp_threads is not None:
        cmd.extend(["--decomp-threads", str(parsed.decomp_threads)])
    if parsed.no_baq:
        cmd.append("--no-baq")
    if parsed.baq:
        cmd.extend(["--baq", parsed.baq])
    if parsed.baq_simd:
        cmd.extend(["--baq-simd", parsed.baq_simd])
    if parsed.baq_device:
        cmd.extend(["--baq-device", parsed.baq_device])
    if parsed.gpu_backend:
        cmd.extend(["--gpu-backend", parsed.gpu_backend])
    if parsed.gpu_device:
        cmd.extend(["--gpu-device", parsed.gpu_device])
    if parsed.gpu_mem_mb is not None:
        cmd.extend(["--gpu-mem-mb", str(parsed.gpu_mem_mb)])
    if parsed.context != "both":
        cmd.extend(["--context", parsed.context])
    if parsed.caller:
        cmd.extend(["--caller", parsed.caller])
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
    return cmd


def build_tune_command(rust_bin: str, parsed) -> List[str]:
    """Argv for `taps-sc-extract-rs tune` from `taps-sc-extract --test ...`."""
    temp_dir = parsed.test_temp_dir or parsed.temp_dir or "."
    out_dir = "."
    if parsed.out:
        out_dir = parsed.out if os.path.isdir(parsed.out) else (os.path.dirname(parsed.out) or ".")
    cmd = [rust_bin, "tune", "--temp-dir", temp_dir, "--out-dir", out_dir]
    if parsed.bam:
        if len(parsed.bam) > 1 or parsed.bam[0].startswith("@"):
            raise ValueError("--test samples one BAM; give a single -b path")
        cmd.extend(["-b", parsed.bam[0], "-f", parsed.fasta])
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


def auto_trims(parsed, chem, opts, contigs, spec):
    """--trim auto (Python engine): count M-bias on a fixed sample of reads,
    choose trims for the ends the user did not set, and write
    <bam>.mbias.tsv/.png/.pdf beside the output."""
    import pysam
    from . import mbias
    from .fasta import FastFaiReader

    log = logging.getLogger("taps_sc_extract")
    fasta = FastFaiReader(parsed.fasta)
    # Windows from every input's mapped reads together; each input is counted
    # over the same windows and the tables added (same as the Rust engine).
    mapped = defaultdict(int)
    lengths = {}
    for b in spec.bams:
        with pysam.AlignmentFile(b, reference_filename=parsed.fasta) as bam:
            for st in bam.get_index_statistics():
                mapped[st.contig] += st.mapped
            lengths.update(zip(bam.references, bam.lengths))
    sel = [(c, lengths[c], mapped.get(c, 0)) for c in contigs if c in lengths]
    windows = mbias.sample_windows(sel)
    counts = None
    for b in spec.bams:
        with pysam.AlignmentFile(b, reference_filename=parsed.fasta) as bam:
            c = mbias.count_table(bam, fasta, windows, chem, opts["min_mapq"], opts["min_baseq"])
        counts = c if counts is None else mbias.add_counts(counts, c)
    table, plateau, nreads = counts
    mates = (1,) if chem == "bisulfite" else (1, 2)
    chosen = mbias.choose_trims(table, plateau, chem, mates)
    out = {}
    for key in ("trim_r1", "trim_r2", "trim_3p"):
        if getattr(parsed, key) is None:
            out[key] = chosen.get(key, opts[key])
    applied = {**opts, **out}
    shards = parsed.shards if parsed.shards and parsed.shards > 0 else 1
    odir = mbias.output_dir(parsed.out, shards)
    os.makedirs(odir, exist_ok=True)
    prefix = os.path.join(odir, mbias.bam_prefix(spec.name_path))
    trims = {"trim_r1": applied["trim_r1"], "trim_r2": applied["trim_r2"], "trim_3p": applied["trim_3p"]}
    mbias.write_tsv(prefix + ".mbias.tsv", table, plateau, nreads, trims)
    try:
        mbias.plot(prefix + ".mbias.tsv", prefix, mbias.bam_prefix(spec.name_path))
    except Exception as e:  # plotting is a report, never a reason to fail the run
        log.warning(f"M-bias plot not written: {e}")
    log.info(
        f"Auto trim: 5' {applied['trim_r1']} bp, 3' {applied['trim_3p']} bp from {nreads:,} sampled reads "
        f"(M-bias: {prefix}.mbias.tsv/.png/.pdf)"
    )
    for n in chosen["notes"]:
        log.warning(f"Auto trim: {n}")
    return out


def resolve_python_options(parsed, spec):
    """Python engine: detect the chemistry from the BAM and fill unset options
    with its defaults (the Rust engine does the same in chemistry.rs)."""
    import pysam
    from .chemistry import check_reference, defaults, detect, resolve
    from .contigs import primary_contigs
    from .fasta import FastFaiReader

    with pysam.AlignmentFile(spec.bams[0], reference_filename=parsed.fasta) as bam:
        det = detect(bam, parsed.strand_source == "auto")
        lengths = dict(zip(bam.references, bam.lengths))
        contigs = (
            [c.strip() for c in parsed.chroms.split(",") if c.strip()] if parsed.chroms
            else primary_contigs(bam.references)
        )
    chem = resolve(parsed.chemistry, det, parsed.strand_source)
    # Every input must resolve to the same chemistry: runs never mix them.
    for b in spec.bams[1:]:
        with pysam.AlignmentFile(b, reference_filename=parsed.fasta) as bam:
            c = resolve(parsed.chemistry, detect(bam, parsed.strand_source == "auto"), parsed.strand_source)
        if c != chem:
            raise ValueError(f"{spec.bams[0]} is {chem} but {b} is {c}: one run takes one chemistry")
    check_reference(FastFaiReader(parsed.fasta), lengths, contigs)
    if chem == "bisulfite" and parsed.baq == "htslib":
        raise ValueError("--baq htslib treats every unmethylated C as a mismatch in bisulfite data")
    d = defaults(chem)
    opts = {k: (getattr(parsed, k) if getattr(parsed, k) is not None else d[k])
            for k in ("min_mapq", "min_baseq", "max_depth", "trim_r1", "trim_r2", "trim_3p")}
    from .sitemask import MOLECULES_PER_CELL, auto_cell_site_cap
    mpc = parsed.molecules_per_cell or MOLECULES_PER_CELL
    if parsed.cell_site_cap is None or parsed.cell_site_cap == "auto":
        opts["cell_site_cap"] = auto_cell_site_cap(mpc, 2 if chem == "taps" and parsed.no_overlap_clip else 1)
    else:
        opts["cell_site_cap"] = parsed.cell_site_cap
    opts["site_depth_limit"] = parsed.site_depth_limit
    # Cells kept: with several inputs, the union of each input's own filter.
    from .inputs import check_barcodes_match, load_cell_filters
    cells, per_file, listed = load_cell_filters(spec, parsed.whitelist, parsed.complexity,
                                                parsed.min_complexity, parsed.max_complexity)
    check_barcodes_match(spec.bams, listed)
    if len(spec.bams) > 1:
        for b, (kept, n) in zip(spec.bams, per_file):
            logging.getLogger("taps_sc_extract").info(f"  {b}: {'all' if kept is None else kept} cells kept of {'all' if n is None else n} listed")
    opts["keep_cells"] = cells
    if parsed.site_depth_limit is None or parsed.site_depth_limit == "auto":
        from .chemistry import READS_PER_MOLECULE
        from .sitemask import MOLECULES_PER_CELL, auto_site_depth_limit
        mpc = parsed.molecules_per_cell or MOLECULES_PER_CELL
        if cells is not None:
            opts["site_depth_limit"] = auto_site_depth_limit(len(cells), mpc, READS_PER_MOLECULE[chem])
        elif parsed.site_depth_limit == "auto":
            raise ValueError("--site-depth-limit auto needs --complexity or --whitelist to count the cells; give a number instead")
        else:
            opts["site_depth_limit"] = 0
    opts["trim_mode"] = parsed.trim or "auto"  # M-bias trims by default, both chemistries
    if opts["trim_mode"] == "auto":
        opts.update(auto_trims(parsed, chem, opts, contigs, spec))
    logging.getLogger("taps_sc_extract").info(
        f"Chemistry: {chem} ({'detected from the BAM' if parsed.chemistry == 'auto' else 'explicit'}"
        f"{'; BSBolt ' + det['bsbolt'] if det['bsbolt'] else ''})"
    )
    return chem, opts


def main(args: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if args is None else args)
    parsed = parse_args(argv)
    setup_logging(verbose=parsed.verbose, log_file=parsed.log_file)
    logger = logging.getLogger("taps_sc_extract")
    if any(v is not None and v < 0 for v in (parsed.trim_r1, parsed.trim_r2, parsed.trim_3p)):
        logger.error("--trim-r1/--trim-5p, --trim-r2 and --trim-3p must be >= 0.")
        return 1
    # An .xlsx sample sheet becomes a TSV (paths made absolute) for both engines.
    sheet_tmp = None
    if parsed.bam and len(parsed.bam) == 1 and parsed.bam[0].startswith("@") and parsed.bam[0].lower().endswith(".xlsx"):
        from .inputs import xlsx_to_tsv
        sheet_tmp = tempfile.TemporaryDirectory(prefix="taps-sheet-")
        try:
            parsed.bam = ["@" + xlsx_to_tsv(parsed.bam[0][1:], sheet_tmp.name)]
        except Exception as e:
            logger.error(f"sample sheet {parsed.bam[0][1:]}: {e}")
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
            if parsed.log_file:
                return run_tee(cmd, parsed.log_file)
            res = subprocess.run(cmd)
            return res.returncode
        except Exception as e:
            logger.exception(f"Error executing Rust backend: {e}")
            return 1

    logger.info("Running Python reference multiprocessing extraction engine.")
    if not parsed.no_baq:
        logger.info(
            "Note: the Python engine does not apply BAQ (same calls as --no-baq); "
            "TAPS-aware BAQ needs the Rust engine."
        )
    chroms_list = None
    if parsed.chroms:
        chroms_list = [c.strip() for c in parsed.chroms.split(",") if c.strip()]
    try:
        from .inputs import parse_bam_args
        spec = parse_bam_args(parsed.bam)
        chem, opts = resolve_python_options(parsed, spec)
    except Exception as e:
        logger.error(str(e))
        return 1

    try:
        extract_methylation_parallel(
            bam_path=spec.bams,
            keep_cells=opts["keep_cells"],
            fasta_path=parsed.fasta,
            out_h5_path=parsed.out,
            chroms=chroms_list,
            whitelist_path=parsed.whitelist,
            complexity_path=parsed.complexity,
            min_complexity=parsed.min_complexity,
            max_complexity=parsed.max_complexity,
            n_workers=parsed.workers if parsed.workers is not None else 24,
            decomp_threads=parsed.decomp_threads if parsed.decomp_threads is not None else 1,
            chunk_size_mb=parsed.chunk_size_mb,
            n_shards=parsed.shards if parsed.shards is not None else 1,
            use_temp_files=(parsed.memory_mode != "memory" and not parsed.no_temp_file),
            compression=parsed.compression,
            compression_threads=parsed.compression_threads,
            max_writer_threads=parsed.max_writer_threads,
            temp_dir=parsed.temp_dir,
            min_base_quality=opts["min_baseq"],
            min_mapq=opts["min_mapq"],
            max_depth=opts["max_depth"],
            compute_baq=not parsed.no_baq,
            ignore_orphans=True,
            ignore_overlaps=not parsed.no_overlap_clip,
            trim_r1=opts["trim_r1"],
            trim_r2=opts["trim_r2"],
            trim_3p=opts["trim_3p"],
            trim_mode=opts["trim_mode"],
            contexts=("CG", "CH") if parsed.context == "both" else (parsed.context,),
            chemistry=chem,
            verify_ys=parsed.strand_source == "auto",
            site_depth_limit=opts["site_depth_limit"],
            cell_site_cap=opts["cell_site_cap"],
            blacklist_path=parsed.blacklist,
        )
        return 0
    except Exception as e:
        logger.exception(f"Error during extraction: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
