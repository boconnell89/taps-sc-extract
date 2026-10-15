#!/usr/bin/env bash
# Build taps-sc-extract-rs optimised for this host and install the Python CLI.
#
#   ./install.sh                 # native build for this CPU, then pip install -e .
#   ./install.sh --portable      # x86-64-v3 build that runs on any AVX2 host
#   ./install.sh --prefix DIR    # conda prefix with hdf5, clang/libclang, pkg-config
#   ./install.sh --no-pip        # build the Rust core only
#   ./install.sh --smoke-bam BAM --smoke-fasta FA   # also extract chr19:5.0-5.1 Mb
#
# A native build uses instructions of the CPU it was built on. Rerun this script
# on each host; copying the binary to a different CPU can fail with SIGILL.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CRATE="$ROOT/taps-sc-extract-rs"
PORTABLE=0
PIP=1
PREFIX="${CONDA_PREFIX:-}"
SMOKE_BAM=""
SMOKE_FASTA=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --portable) PORTABLE=1 ;;
        --no-pip) PIP=0 ;;
        --prefix) PREFIX="$2"; shift ;;
        --smoke-bam) SMOKE_BAM="$2"; shift ;;
        --smoke-fasta) SMOKE_FASTA="$2"; shift ;;
        -h|--help) sed -n '2,13p' "$0"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
    shift
done

if [[ -z "$PREFIX" ]]; then
    for p in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3"; do
        if [[ -d "$p" ]]; then PREFIX="$p"; break; fi
    done
fi
missing=()
[[ -n "$PREFIX" && -d "$PREFIX" ]] || missing+=("a conda prefix (activate an env or pass --prefix)")
if [[ ${#missing[@]} -eq 0 ]]; then
    ls "$PREFIX"/lib/libhdf5.so* >/dev/null 2>&1 || missing+=("libhdf5 in $PREFIX/lib")
    ls "$PREFIX"/lib/libclang*.so* >/dev/null 2>&1 || missing+=("libclang in $PREFIX/lib")
fi
command -v cargo >/dev/null || missing+=("cargo (https://rustup.rs)")
for tool in curl patch sha256sum; do
    command -v "$tool" >/dev/null || missing+=("$tool")
done
if [[ ${#missing[@]} -gt 0 ]]; then
    echo "install.sh: missing:" >&2
    printf '  - %s\n' "${missing[@]}" >&2
    echo "Install with: conda install -c conda-forge hdf5 clang libclang pkg-config" >&2
    exit 1
fi

if [[ $PORTABLE -eq 1 ]]; then
    TARGET_CPU="x86-64-v3"
    MARCH="x86-64-v3"
else
    TARGET_CPU="native"
    MARCH="native"
fi
# The generated config carries the target CPU for Rust ([build] rustflags) and
# for the vendored C (forced CFLAGS), so later plain `cargo build --release`
# runs in this checkout produce the same binary.
sed -e "s|@PREFIX@|$PREFIX|g" -e "s|@TARGET_CPU@|$TARGET_CPU|g" -e "s|@MARCH@|$MARCH|g" \
    "$CRATE/.cargo/config.toml.in" > "$CRATE/.cargo/config.toml"
# RUSTFLAGS in the environment would replace [build] rustflags.
unset RUSTFLAGS

CPU_MODEL="$(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2- | sed 's/^ //' || true)"
CPU_FLAGS="$(grep -m1 '^flags' /proc/cpuinfo 2>/dev/null | tr ' ' '\n' | grep -E '^(avx2|avx512f|avx512bw|bmi2|sha_ni)$' | tr '\n' ' ' || true)"
echo "Host:      $(hostname)  ${CPU_MODEL}"
echo "Cores:     $(nproc)  RAM: $(awk '/MemTotal/ {printf "%.0f GB", $2/1048576}' /proc/meminfo)"
echo "CPU flags: ${CPU_FLAGS:-unknown}"
echo "Prefix:    $PREFIX"
echo "Target:    target-cpu=$TARGET_CPU, -march=$MARCH"

# hts-sys 2.2.1 bundles htslib 1.19.1. Patch in htslib 1.21's overlap clipping so
# Rust calls match the Python engine (pysam / htslib 1.21); see
# taps-sc-extract-rs/patches/. Cargo.toml points hts-sys at this copy.
HTS_VER=2.2.1
HTS_DIR="$CRATE/vendor/hts-sys-$HTS_VER"
HTS_PATCH="$CRATE/patches/hts-sys-$HTS_VER-overlap-deletions.patch"
if [[ ! -f "$HTS_DIR/.taps_patched" ]]; then
    echo "Vendoring hts-sys $HTS_VER with $(basename "$HTS_PATCH") ..."
    # crates.io checksum of hts-sys 2.2.1 (Cargo.lock no longer records it once
    # [patch.crates-io] points at the vendored copy).
    want="fc7e68eb880b02c80cfb41e8dc7904062a3ea7e27b7c4556e88d648dd2f038da"
    tmp="$(mktemp -d)"
    curl -sSfL "https://static.crates.io/crates/hts-sys/hts-sys-$HTS_VER.crate" -o "$tmp/hts-sys.crate"
    got="$(sha256sum "$tmp/hts-sys.crate" | cut -d' ' -f1)"
    if [[ "$got" != "$want" ]]; then
        echo "install.sh: hts-sys-$HTS_VER.crate checksum $got does not match crates.io $want" >&2
        exit 1
    fi
    rm -rf "$HTS_DIR"
    mkdir -p "$CRATE/vendor"
    tar -xzf "$tmp/hts-sys.crate" -C "$CRATE/vendor"
    patch -d "$HTS_DIR" -p1 --forward --quiet < "$HTS_PATCH"
    touch "$HTS_DIR/.taps_patched"
    rm -rf "$tmp"
fi

# Flag changes do not always trigger a rebuild of the vendored C, so start clean.
(cd "$CRATE" && cargo clean --release -q && cargo build --release)

BIN="$CRATE/target/release/taps-sc-extract-rs"
cat > "$CRATE/target/release/BUILD_INFO" <<INFO
built=$(date -Is)
host=$(hostname)
cpu_model=$CPU_MODEL
nproc=$(nproc)
target_cpu=$TARGET_CPU
march=$MARCH
prefix=$PREFIX
git=$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)
INFO

"$BIN" --version

if [[ $PIP -eq 1 ]]; then
    python3 -m pip install -e "$ROOT"
fi

if [[ -n "$SMOKE_BAM" && -n "$SMOKE_FASTA" ]]; then
    out="$(mktemp -d)"
    "$BIN" extract -b "$SMOKE_BAM" -f "$SMOKE_FASTA" -o "$out/smoke.h5" \
        -c chr19 --start 5000000 --end 5100000 -t 4 --shards 1 --memory-mode memory 2>&1 | tail -n 4
    rm -rf "$out"
fi
echo "Installed: $BIN ($TARGET_CPU). Rerun install.sh on each new host."
