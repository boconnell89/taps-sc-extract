#!/usr/bin/env python3
"""
Automated A/B Performance Benchmark comparing Python and Rust TAPS extraction engines.
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time


def run_timed_command(cmd, name):
    print(f"\n[{name}] Running: {' '.join(cmd)}")
    time_out_file = tempfile.mktemp(prefix="taps_time_")
    time_cmd = ["/usr/bin/time", "-v", "-o", time_out_file] + cmd
    
    t0 = time.time()
    res = subprocess.run(time_cmd, capture_output=True, text=True)
    t1 = time.time()
    
    peak_rss_mb = 0.0
    if os.path.exists(time_out_file):
        with open(time_out_file) as f:
            for line in f:
                if "Maximum resident set size" in line:
                    parts = line.split(":")
                    if len(parts) >= 2:
                        peak_rss_mb = float(parts[1].strip()) / 1024.0
        os.remove(time_out_file)
        
    return {
        "name": name,
        "returncode": res.returncode,
        "elapsed_s": t1 - t0,
        "peak_rss_mb": peak_rss_mb,
        "stdout": res.stdout,
        "stderr": res.stderr,
    }


def main():
    parser = argparse.ArgumentParser(description="TAPS Extract A/B Performance Benchmark (Python vs Rust)")
    parser.add_argument("-b", "--bam", required=True, help="Path to input BAM")
    parser.add_argument("-f", "--fasta", required=True, help="Path to reference FASTA")
    parser.add_argument("-c", "--chroms", default="chr19", help="Contig(s) to benchmark (default: chr19)")
    parser.add_argument("-w", "--whitelist", default=None, help="Optional barcode whitelist")
    parser.add_argument("-t", "--workers", type=int, default=24, help="Worker count (default: 24)")
    parser.add_argument("--shards", type=int, default=8, help="Shard count (default: 8)")
    parser.add_argument("--decomp-threads", type=int, default=1, help="Decompression threads (default: 1)")
    args = parser.parse_args()

    out_py = tempfile.mkdtemp(prefix="taps_ab_py_")
    out_rs = tempfile.mkdtemp(prefix="taps_ab_rs_")

    try:
        # 1. Benchmark Python Engine
        cmd_py = [
            sys.executable, "-m", "taps_sc_extract",
            "-b", args.bam,
            "-f", args.fasta,
            "-o", out_py,
            "-c", args.chroms,
            "-t", str(args.workers),
            "--shards", str(args.shards),
            "--decomp-threads", str(args.decomp_threads),
            "--engine", "python",
            "--no-baq",
        ]
        if args.whitelist:
            cmd_py.extend(["-w", args.whitelist])
        py_res = run_timed_command(cmd_py, "Python Engine")

        # 2. Benchmark Rust Engine
        cmd_rs = [
            sys.executable, "-m", "taps_sc_extract",
            "-b", args.bam,
            "-f", args.fasta,
            "-o", out_rs,
            "-c", args.chroms,
            "-t", str(args.workers),
            "--shards", str(args.shards),
            "--decomp-threads", str(args.decomp_threads),
            "--engine", "rust",
            "--no-baq",
        ]
        if args.whitelist:
            cmd_rs.extend(["-w", args.whitelist])
        rs_res = run_timed_command(cmd_rs, "Rust Engine")

        print("\n" + "=" * 70)
        print("A/B PERFORMANCE BENCHMARK RESULTS")
        print("=" * 70)
        speedup = (py_res["elapsed_s"] / rs_res["elapsed_s"]) if rs_res["elapsed_s"] > 0 else 0.0
        ram_reduction = (1.0 - (rs_res["peak_rss_mb"] / py_res["peak_rss_mb"])) * 100.0 if py_res["peak_rss_mb"] > 0 else 0.0

        print(f"| Metric | Python Engine | Rust Engine | Gain / Delta |")
        print(f"| :--- | :---: | :---: | :---: |")
        print(f"| **Wall Time** | {py_res['elapsed_s']:.2f} s | {rs_res['elapsed_s']:.2f} s | **{speedup:.2f}x faster** |")
        print(f"| **Peak RAM (RSS)** | {py_res['peak_rss_mb']:.1f} MB | {rs_res['peak_rss_mb']:.1f} MB | **{ram_reduction:.1f}% reduction** |")
        print(f"| **Exit Code** | {py_res['returncode']} | {rs_res['returncode']} | {'Pass' if rs_res['returncode'] == 0 else 'Fail'} |")
        print("=" * 70)

    finally:
        shutil.rmtree(out_py, ignore_errors=True)
        shutil.rmtree(out_rs, ignore_errors=True)


if __name__ == "__main__":
    main()
