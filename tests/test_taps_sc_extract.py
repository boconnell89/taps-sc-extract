"""
Unit and integration tests for taps_sc_extract.
"""

import os
import tempfile
import h5py
import numpy as np
import pysam
import pytest

from taps_sc_extract.barcode import extract_barcode, parse_annot
from taps_sc_extract.cli import parse_args
from taps_sc_extract.calling import (
    OB_FLAGS,
    OT_FLAGS,
    call_mctot,
    classify_strand,
    keep_after_5prime_trim,
)
from taps_sc_extract.context import classify_context, classify_trinucleotide_context
from taps_sc_extract.extractor import CANONICAL_CONTIGS, extract_methylation
from taps_sc_extract.h5_writer import AmethystH5Writer, METH_DTYPE


# 1. Context classifier tests
def test_context_classifier_real_ref():
    # Top strand (ref C):
    assert classify_trinucleotide_context("CGA", 0) == "CpG"
    assert classify_context("CGA", 0) == "CG"

    assert classify_trinucleotide_context("CAG", 0) == "CHG"
    assert classify_context("CAG", 0) == "CH"

    assert classify_trinucleotide_context("CTA", 0) == "CHH"
    assert classify_context("CTA", 0) == "CH"

    assert classify_trinucleotide_context("CCC", 0) == "CHH"
    assert classify_context("CCC", 0) == "CH"

    # Bottom strand (ref G):
    # 'TCG' -> G at pos 2. Upstream is C at pos 1 -> CpG
    assert classify_trinucleotide_context("TCG", 2) == "CpG"
    assert classify_context("TCG", 2) == "CG"

    # 'CTG' -> G at pos 2. Upstream is T at 1, C at 0 -> comp is CAG -> CHG
    assert classify_trinucleotide_context("CTG", 2) == "CHG"
    assert classify_context("CTG", 2) == "CH"

    # 'TAG' -> G at pos 2. Upstream is A at 1, T at 0 -> comp is CTA -> CHH
    assert classify_trinucleotide_context("TAG", 2) == "CHH"
    assert classify_context("TAG", 2) == "CH"

    # Edge cases / bounds
    assert classify_context("C", 0) == "CH"
    assert classify_context("G", 0) == "CH"
    assert classify_context("A", 0) is None
    assert classify_context("T", 0) is None
    assert classify_context("CGA", 5) is None
    assert classify_context("CGA", -1) is None

    # Test with real mm10 FASTA if available
    ref_fasta = "/mnt/e/refs/mm10/mm10.fa"
    if os.path.exists(ref_fasta):
        fa = pysam.FastaFile(ref_fasta)
        # Fetch 100bp on chr19
        seq = fa.fetch("chr19", 5000000, 5000100).upper()
        for i, b in enumerate(seq):
            ctx = classify_context(seq, i)
            if b in ("C", "G"):
                assert ctx in ("CG", "CH")
            else:
                assert ctx is None
        fa.close()


# 2. Flag classifier tests
def test_flag_classify_ot_ob():
    # Valid OT flags
    assert classify_strand(99) == "OT"
    assert classify_strand(147) == "OT"

    # Valid OB flags
    assert classify_strand(83) == "OB"
    assert classify_strand(163) == "OB"

    # Excluded flags (unmapped, orphans, secondary, supplementary, etc.)
    for f in [0, 16, 4, 355, 2048, 77, 141, 65, 129, 113, 177]:
        assert classify_strand(f) is None


# 3. mCtoT calling logic tests
def test_mCtoT_interpret():
    # OT strand (ref C): T = meth, C = unmeth
    assert call_mctot("OT", "C", "T") == 1
    assert call_mctot("OT", "C", "C") == 0
    assert call_mctot("OT", "C", "A") is None
    assert call_mctot("OT", "C", "G") is None
    assert call_mctot("OT", "C", "N") is None
    # OT strand with ref G should not make a call
    assert call_mctot("OT", "G", "A") is None
    assert call_mctot("OT", "G", "G") is None

    # OB strand (ref G): A = meth, G = unmeth
    assert call_mctot("OB", "G", "A") == 1
    assert call_mctot("OB", "G", "G") == 0
    assert call_mctot("OB", "G", "T") is None
    assert call_mctot("OB", "G", "C") is None
    assert call_mctot("OB", "G", "N") is None
    # OB strand with ref C should not make a call
    assert call_mctot("OB", "C", "T") is None
    assert call_mctot("OB", "C", "C") is None


def test_keep_after_5prime_trim():
    # R2 forward (flag 163): first 10 SEQ bases are the sequenced 5'.
    flag_r2_fwd = 163
    for q in range(10):
        assert keep_after_5prime_trim(q, 49, flag_r2_fwd, 0, 10) is False
    assert keep_after_5prime_trim(10, 49, flag_r2_fwd, 0, 10) is True
    assert keep_after_5prime_trim(48, 49, flag_r2_fwd, 0, 10) is True

    # R2 reverse (flag 147): sequenced 5' is at the high end of SEQ.
    flag_r2_rev = 147
    assert keep_after_5prime_trim(0, 49, flag_r2_rev, 0, 10) is True
    assert keep_after_5prime_trim(38, 49, flag_r2_rev, 0, 10) is True
    for q in range(39, 49):
        assert keep_after_5prime_trim(q, 49, flag_r2_rev, 0, 10) is False

    # R1 default trim 0 keeps every base.
    for flag in (99, 83):
        for q in range(39):
            assert keep_after_5prime_trim(q, 39, flag, 0, 10) is True

    # R1 reverse with trim_r1=5 drops the last 5 SEQ bases.
    flag_r1_rev = 83
    assert keep_after_5prime_trim(33, 39, flag_r1_rev, 5, 10) is True
    assert keep_after_5prime_trim(34, 39, flag_r1_rev, 5, 10) is False

    # Zero trim keeps everything; trim covering the whole read drops everything.
    assert keep_after_5prime_trim(0, 49, flag_r2_fwd, 0, 0) is True
    assert keep_after_5prime_trim(0, 10, flag_r2_fwd, 0, 10) is False
    assert keep_after_5prime_trim(9, 10, flag_r2_fwd, 0, 10) is False


def test_cli_trim_defaults():
    ns = parse_args(["-b", "x.bam", "-f", "x.fa", "-o", "out.h5"])
    assert ns.trim_r1 == 0
    assert ns.trim_r2 == 10
    ns = parse_args(["-b", "x.bam", "-f", "x.fa", "-o", "out.h5", "--trim-r1", "5", "--trim-r2", "0"])
    assert ns.trim_r1 == 5
    assert ns.trim_r2 == 0


# 4. Barcode extraction tests
def test_extract_barcode_qname():
    # Real QNAME formats
    assert extract_barcode("ATTGGCTCAATGATCCTGGCGCTCCATTCG:77600300#0") == "ATTGGCTCAATGATCCTGGCGCTCCATTCG"
    assert extract_barcode("AACGACGAACGAATGCCTTGCGAATTCGTT:85915554#0") == "AACGACGAACGAATGCCTTGCGAATTCGTT"

    # Variable barcode lengths
    assert extract_barcode("ACGT:123#0") == "ACGT"
    assert extract_barcode("AACGACGAACAGGCCGGCCAAACGTTCC:999#0") == "AACGACGAACAGGCCGGCCAAACGTTCC"

    # Slash separator fallback
    assert extract_barcode("BC123/1") == "BC123"

    # No separator fallback
    assert extract_barcode("BC123") == "BC123"


# 5. Barcode whitelist / annot parsing tests
def test_annot_parse_variable_columns(tmp_path):
    # 1-column file
    p1 = tmp_path / "1col.txt"
    p1.write_text("BC1\nBC2\nBC3\n")
    assert parse_annot(str(p1)) == {"BC1", "BC2", "BC3"}

    # 2-column file
    p2 = tmp_path / "2col.annot"
    p2.write_text("BC1\tsampleA\nBC2\tsampleB\n# comment\n\nBC3\tsampleC\n")
    assert parse_annot(str(p2)) == {"BC1", "BC2", "BC3"}

    # N-column file with mixed whitespace
    pN = tmp_path / "Ncol.txt"
    pN.write_text("BC1   info1   info2   info3\nBC2\tinfo1\tinfo2\n")
    assert parse_annot(str(pN)) == {"BC1", "BC2"}


# 6. Position conversion test
def test_pos_conversion():
    # 0-based pileup position -> 1-based output position
    pysam_pos_0based = 100
    pos_1based = pysam_pos_0based + 1
    assert pos_1based == 101

    rec = np.zeros(1, dtype=METH_DTYPE)
    rec[0]["chr"] = b"chr19"
    rec[0]["pos"] = pos_1based
    rec[0]["pct"] = 100.0
    rec[0]["t"] = 0
    rec[0]["c"] = 1
    assert rec[0]["pos"] == 101


# 7. HDF5 dtype and layout test
def test_h5_dtype_and_layout(tmp_path):
    h5_file = tmp_path / "test.h5"
    bc = "ACGTACGT"

    with AmethystH5Writer(str(h5_file)) as writer:
        cg_data = np.array([
            (b"chr19", 100, 100.0, 0, 1),
            (b"chr19", 200, 0.0, 2, 0)
        ], dtype=METH_DTYPE)
        ch_data = np.array([
            (b"chr19", 150, 0.0, 1, 0)
        ], dtype=METH_DTYPE)

        writer.append_data("CG", bc, cg_data)
        writer.append_data("CH", bc, ch_data)

    with h5py.File(str(h5_file), "r") as f:
        assert "CG" in f
        assert "CH" in f
        assert bc in f["CG"]
        assert bc in f["CH"]
        assert "1" in f["CG"][bc]
        assert "1" in f["CH"][bc]

        ds_cg = f["CG"][bc]["1"]
        assert ds_cg.dtype == METH_DTYPE
        assert ds_cg.dtype.names == ("chr", "pos", "pct", "t", "c")
        # Check field order: t before c
        assert ds_cg.dtype.names[3] == "t"
        assert ds_cg.dtype.names[4] == "c"

        assert len(ds_cg) == 2
        assert ds_cg[0]["chr"] == b"chr19"
        assert ds_cg[0]["pos"] == 100
        assert ds_cg[0]["pct"] == 100.0
        assert ds_cg[0]["t"] == 0
        assert ds_cg[0]["c"] == 1


# 8. HDF5 chromosome contiguity and position ascending test
def test_h5_chr_contiguous_pos_ascending(tmp_path):
    h5_file = tmp_path / "multi_chr.h5"
    bc = "BARCODE_TEST"

    with AmethystH5Writer(str(h5_file)) as writer:
        chr19_data = np.array([
            (b"chr19", 10, 0.0, 1, 0),
            (b"chr19", 50, 100.0, 0, 1),
            (b"chr19", 100, 50.0, 1, 1),
        ], dtype=METH_DTYPE)
        chrX_data = np.array([
            (b"chrX", 5, 0.0, 1, 0),
            (b"chrX", 20, 100.0, 0, 2),
        ], dtype=METH_DTYPE)

        writer.append_data("CG", bc, chr19_data)
        writer.append_data("CG", bc, chrX_data)

    with h5py.File(str(h5_file), "r") as f:
        ds = f["CG"][bc]["1"][:]
        # All chr19 rows first, then all chrX rows
        chrs = [row["chr"].decode("ascii") for row in ds]
        assert chrs == ["chr19", "chr19", "chr19", "chrX", "chrX"]

        # Positions ascending within chr19
        chr19_pos = [row["pos"] for row in ds if row["chr"] == b"chr19"]
        assert chr19_pos == sorted(chr19_pos)

        # Positions ascending within chrX
        chrX_pos = [row["pos"] for row in ds if row["chr"] == b"chrX"]
        assert chrX_pos == sorted(chrX_pos)


# 9. Canonical contig allowlist test
def test_canonical_contig_allowlist():
    expected_21 = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"]
    assert CANONICAL_CONTIGS == expected_21
    assert len(CANONICAL_CONTIGS) == 21

    # Ensure chrM and alt-scaffolds are excluded
    assert "chrM" not in CANONICAL_CONTIGS
    assert "chrUn_GL456239" not in CANONICAL_CONTIGS
    assert "chr1_GL456210_random" not in CANONICAL_CONTIGS


# 10. Real BAM smoke test (chr19)
CHR19_BAM = "/mnt/e/sciMET_TAPS/260731/taps_ucbTn5_sp2.srt.bam"
CHR19_FASTA = "/mnt/e/refs/mm10/mm10.fa"
# Rust chr19 output frozen after the barcode-cache fix (see the README in that
# directory). The older test_chr19_rust_026c52e/ reference has 0.01% of calls
# under the wrong barcode and is kept only as a record.
CHR19_RUST_REFERENCE = "/mnt/e/sciMET_TAPS/test_chr19_rust_cellfix_2026-10-15/chr19.h5"


def _require(*paths):
    for p in paths:
        if not os.path.exists(p):
            pytest.skip(f"{p} not available in environment.")


@pytest.fixture(scope="module")
def chr19_python_run(tmp_path_factory):
    """Python reference engine on chr19 (BAQ is a no-op there), run once per module."""
    _require(CHR19_BAM, CHR19_FASTA)
    out_h5 = tmp_path_factory.mktemp("chr19_py") / "chr19_smoke.h5"
    stats = extract_methylation(
        bam_path=CHR19_BAM,
        fasta_path=CHR19_FASTA,
        out_h5_path=str(out_h5),
        chroms=["chr19"],
    )
    return out_h5, stats


def _run_rust_chr19(out_h5, *extra):
    import subprocess
    from taps_sc_extract.cli import find_rust_binary

    rust_bin = find_rust_binary()
    if not rust_bin:
        pytest.skip("Rust binary not built.")
    subprocess.run(
        [rust_bin, "extract", "-b", CHR19_BAM, "-f", CHR19_FASTA, "-o", str(out_h5),
         "-c", "chr19", "-t", "8", "--shards", "1", "--memory-mode", "memory", *extra],
        check=True,
        capture_output=True,
    )


def _sites_by_cell(path, ctx):
    """{(barcode, pos): (t, c)} and {pos: [t, c]} summed over cells."""
    cells, totals = {}, {}
    with h5py.File(str(path), "r") as f:
        for bc in f[ctx]:
            x = f[ctx][bc]["1"][()]
            for p, t, c in zip(x["pos"].tolist(), x["t"].tolist(), x["c"].tolist()):
                cells[(bc, p)] = (t, c)
                tot = totals.setdefault(p, [0, 0])
                tot[0] += t
                tot[1] += c
    return cells, totals


def test_chr19_real_bam_smoke(chr19_python_run):
    out_h5, stats = chr19_python_run

    assert os.path.exists(out_h5)
    assert os.path.getsize(out_h5) > 0

    with h5py.File(str(out_h5), "r") as f:
        assert "CG" in f
        assert "CH" in f
        cg_barcodes = list(f["CG"].keys())
        assert len(cg_barcodes) > 0

        sample_bc = cg_barcodes[0]
        ds = f["CG"][sample_bc]["1"]
        assert len(ds) > 0
        assert ds.dtype == METH_DTYPE

    # Check that methylation rates are biologically reasonable
    total_stats = stats["stats"]
    cg_tot = total_stats["CG"]["c"] + total_stats["CG"]["t"]
    ch_tot = total_stats["CH"]["c"] + total_stats["CH"]["t"]
    assert cg_tot > 0
    assert ch_tot > 0

    cg_pct = 100.0 * total_stats["CG"]["c"] / cg_tot
    ch_pct = 100.0 * total_stats["CH"]["c"] / ch_tot

    # CpG rate should be high (~50-80%), CH rate should be low (<2%)
    assert 40.0 <= cg_pct <= 85.0
    assert 0.0 <= ch_pct <= 3.0


def test_chr19_real_bam_smoke_rust(tmp_path):
    """Rust engine on chr19 must reproduce the frozen reference sites exactly."""
    _require(CHR19_BAM, CHR19_FASTA, CHR19_RUST_REFERENCE)
    out_h5 = tmp_path / "chr19_rust.h5"
    _run_rust_chr19(out_h5)

    cols = ["chr", "pos", "t", "c"]
    with h5py.File(CHR19_RUST_REFERENCE, "r") as ref, h5py.File(str(out_h5), "r") as new:
        for ctx in ("CG", "CH"):
            assert sorted(ref[ctx].keys()) == sorted(new[ctx].keys())
            for bc in ref[ctx]:
                a = ref[ctx][bc]["1"][()][cols]
                b = new[ctx][bc]["1"][()][cols]
                assert np.array_equal(a, b), f"{ctx}/{bc} differs"


def test_chr19_rust_python_same_cell_per_call(chr19_python_run, tmp_path):
    """Rust --no-baq vs Python on chr19: no call may sit under a different barcode.

    Wherever the two engines agree on a position's total (t, c summed over
    cells), every cell must agree too. Positions whose totals differ are the
    known, separate Rust-has-more-calls gap and are only bounded here.
    """
    py_h5, _ = chr19_python_run
    rs_h5 = tmp_path / "chr19_rust_nobaq.h5"
    _run_rust_chr19(rs_h5, "--no-baq")

    for ctx, max_total_gap in (("CG", 100), ("CH", 1000)):
        py_cells, py_tot = _sites_by_cell(py_h5, ctx)
        rs_cells, rs_tot = _sites_by_cell(rs_h5, ctx)
        gap = {p for p in set(py_tot) | set(rs_tot) if py_tot.get(p) != rs_tot.get(p)}
        moved = [
            k for k in set(py_cells) | set(rs_cells)
            if k[1] not in gap and py_cells.get(k) != rs_cells.get(k)
        ]
        assert not moved, f"{ctx}: {len(moved)} calls under a different barcode, e.g. {moved[:3]}"
        assert len(gap) <= max_total_gap, f"{ctx}: {len(gap)} positions with differing totals"


def test_sharded_hdf5_writing(tmp_path):
    """Test that write_sharded_hdf5_dir creates shards and master.h5 with working external links."""
    from taps_sc_extract.h5_writer import write_sharded_hdf5_dir

    output_dir = tmp_path / "sharded_test"
    cg_data = {
        "cell_1": np.array([(b"chr1", 100, 50.0, 1, 1)], dtype=METH_DTYPE),
        "cell_2": np.array([(b"chr1", 200, 100.0, 0, 2)], dtype=METH_DTYPE),
        "cell_3": np.array([(b"chr2", 300, 0.0, 2, 0)], dtype=METH_DTYPE),
    }
    ch_data = {
        "cell_1": np.array([(b"chr1", 150, 0.0, 3, 0)], dtype=METH_DTYPE),
        "cell_2": np.array([(b"chr1", 250, 0.0, 1, 0)], dtype=METH_DTYPE),
    }

    master_path = write_sharded_hdf5_dir(
        output_dir=str(output_dir),
        cg_data=cg_data,
        ch_data=ch_data,
        n_shards=2,
    )

    assert os.path.exists(master_path)
    assert os.path.exists(output_dir / "shard_000.h5")
    assert os.path.exists(output_dir / "shard_001.h5")

    with h5py.File(master_path, "r") as f:
        assert "CG" in f
        assert "CH" in f
        assert "cell_1" in f["CG"]
        assert "cell_2" in f["CG"]
        assert "cell_3" in f["CG"]
        assert len(f["CG"]["cell_1"]["1"]) == 1
        assert len(f["CH"]["cell_1"]["1"]) == 1
        assert f["CG"]["cell_1"]["1"][0]["pos"] == 100


def _tiny_records(chrom: str, pos: int) -> np.ndarray:
    return np.array([(chrom.encode("ascii"), pos, 50.0, 1, 1)], dtype=METH_DTYPE)


def test_shard_writer_concurrency_after_extract():
    from taps_sc_extract.parallel_extractor import _shard_writer_concurrency

    assert _shard_writer_concurrency(1, use_temp_files=True) == 1
    n_drain = _shard_writer_concurrency(16, use_temp_files=True)
    assert 1 <= n_drain <= 6  # default cap is 6
    assert _shard_writer_concurrency(16, use_temp_files=True, max_writer_threads=3) <= 3
    assert _shard_writer_concurrency(16, use_temp_files=True, max_writer_threads=1) == 1


def test_assemble_and_write_shard_temp_file_order(tmp_path):
    """Temp chunk files must be concatenated in chunk_id / filename order."""
    import pickle
    from taps_sc_extract.parallel_extractor import _assemble_and_write_shard

    shard_dir = tmp_path / "shard_000"
    shard_dir.mkdir()
    # Write chunk 1 first on disk so directory iteration order is not the contract.
    d1 = {"CG": {"cellA": _tiny_records("chr2", 200)}, "CH": {}}
    d0 = {"CG": {"cellA": _tiny_records("chr1", 100)}, "CH": {}}
    with open(shard_dir / "chunk_000001.bin", "wb") as f:
        pickle.dump(d1, f)
    with open(shard_dir / "chunk_000000.bin", "wb") as f:
        pickle.dump(d0, f)

    out_h5 = tmp_path / "shard_000.h5"
    s_idx, n_cells, _elapsed = _assemble_and_write_shard(
        shard_idx=0,
        n_shards=1,
        shard_path=str(out_h5),
        temp_shard_dir=str(shard_dir),
        in_memory_chunks=None,
    )
    assert s_idx == 0
    assert n_cells == 1
    with h5py.File(out_h5, "r") as h5:
        pos = list(h5["CG"]["cellA"]["1"]["pos"])
        chrom = [c.decode() for c in h5["CG"]["cellA"]["1"]["chr"]]
        assert chrom == ["chr1", "chr2"]
        assert pos == [100, 200]


def test_assemble_and_write_shard_in_memory_ignores_completion_order(tmp_path):
    """imap_unordered completion order must not scramble genomic order."""
    from taps_sc_extract.parallel_extractor import _assemble_and_write_shard

    d0 = {"CG": {"cellA": _tiny_records("chr1", 100)}, "CH": {}}
    d1 = {"CG": {"cellA": _tiny_records("chr2", 200)}, "CH": {}}
    # Reverse of genomic order, as imap_unordered can deliver.
    mem = [(1, d1), (0, d0)]
    out_h5 = tmp_path / "mem.h5"
    _s, n_cells, _e = _assemble_and_write_shard(
        shard_idx=0,
        n_shards=1,
        shard_path=str(out_h5),
        temp_shard_dir=None,
        in_memory_chunks=mem,
    )
    assert n_cells == 1
    with h5py.File(out_h5, "r") as h5:
        pos = list(h5["CG"]["cellA"]["1"]["pos"])
        chrom = [c.decode() for c in h5["CG"]["cellA"]["1"]["chr"]]
        assert chrom == ["chr1", "chr2"]
        assert pos == [100, 200]


def test_assemble_lzf_compression(tmp_path):
    """lzf is a supported faster writer path; gzip remains the Amethyst default."""
    import pickle
    from taps_sc_extract.parallel_extractor import _assemble_and_write_shard

    shard_dir = tmp_path / "shard_000"
    shard_dir.mkdir()
    d0 = {"CG": {"cellA": _tiny_records("chr1", 100)}, "CH": {}}
    with open(shard_dir / "chunk_000000.bin", "wb") as f:
        pickle.dump(d0, f)
    out_h5 = tmp_path / "lzf.h5"
    _assemble_and_write_shard(
        shard_idx=0,
        n_shards=1,
        shard_path=str(out_h5),
        temp_shard_dir=str(shard_dir),
        in_memory_chunks=None,
        compression="lzf",
    )
    with h5py.File(out_h5, "r") as h5:
        ds = h5["CG"]["cellA"]["1"]
        assert ds.compression == "lzf"
        assert list(ds["pos"]) == [100]


def test_assemble_blosc_compression(tmp_path):
    """Blosc (hdf5plugin) is the multithreaded compressed writer path."""
    import pickle
    from taps_sc_extract.h5_writer import configure_blosc_threads
    from taps_sc_extract.parallel_extractor import _assemble_and_write_shard

    n = configure_blosc_threads(n_writers=4, compression="blosc")
    assert n >= 1
    assert os.environ["BLOSC_NTHREADS"] == str(n)

    shard_dir = tmp_path / "shard_000"
    shard_dir.mkdir()
    d0 = {"CG": {"cellA": _tiny_records("chr1", 100)}, "CH": {}}
    with open(shard_dir / "chunk_000000.bin", "wb") as f:
        pickle.dump(d0, f)
    out_h5 = tmp_path / "blosc.h5"
    _assemble_and_write_shard(
        shard_idx=0,
        n_shards=1,
        shard_path=str(out_h5),
        temp_shard_dir=str(shard_dir),
        in_memory_chunks=None,
        compression="blosc",
    )
    with h5py.File(out_h5, "r") as h5:
        ds = h5["CG"]["cellA"]["1"]
        filt = ds.id.get_create_plist().get_filter(0)
        assert filt[0] == 32001  # Blosc filter id
        assert filt[3] == b"blosc"
        assert list(ds["pos"]) == [100]


def test_rust_command_uses_rust_defaults_unless_set():
    from taps_sc_extract.cli import build_rust_command, parse_args

    base = ["-b", "in.bam", "-f", "ref.fa", "-o", "out.h5"]
    parsed = parse_args(base)
    cmd = build_rust_command("taps-sc-extract-rs", parsed, base)
    assert cmd[cmd.index("--chunk-size-mb") + 1] == "0"
    assert cmd[cmd.index("--compression") + 1] == "lzf"
    assert "--decomp-threads" not in cmd

    explicit = base + ["--chunk-size-mb", "10", "--compression", "gzip", "--shards", "4"]
    parsed = parse_args(explicit)
    cmd = build_rust_command("taps-sc-extract-rs", parsed, explicit)
    assert cmd[cmd.index("--chunk-size-mb") + 1] == "10"
    assert cmd[cmd.index("--compression") + 1] == "gzip"
    assert cmd[cmd.index("--shards") + 1] == "4"

    equals = base + ["--chunk-size-mb=25", "--compression=none"]
    parsed = parse_args(equals)
    cmd = build_rust_command("taps-sc-extract-rs", parsed, equals)
    assert cmd[cmd.index("--chunk-size-mb") + 1] == "25"
    assert cmd[cmd.index("--compression") + 1] == "none"


def test_sharded_h5_filename_becomes_stem_directory(tmp_path):
    from taps_sc_extract.parallel_extractor import sharded_output_dir

    h5 = tmp_path / "sample.h5"
    h5.write_text("existing file")
    dest = sharded_output_dir(str(h5))
    assert dest == str(tmp_path / "sample")
    assert h5.is_file()

    folder = tmp_path / "already"
    folder.mkdir()
    assert sharded_output_dir(str(folder)) == str(folder)

    blocked = tmp_path / "blocked.h5"
    blocked.write_text("x")
    (tmp_path / "blocked").write_text("not a directory")
    with pytest.raises(FileExistsError):
        sharded_output_dir(str(blocked))

    plain = tmp_path / "plain.txt"
    plain.write_text("x")
    with pytest.raises(FileExistsError):
        sharded_output_dir(str(plain))


def test_h5_compression_kwargs_blosc_zstd():
    from taps_sc_extract.h5_writer import h5_compression_kwargs

    kw = h5_compression_kwargs("blosc-zstd", 1)
    assert kw["compression"] == 32001
    # filter_options: clevel=1, shuffle=1, compressor=zstd(5)
    assert kw["compression_opts"][4] == 1
    assert kw["compression_opts"][6] == 5



def test_test_mode_needs_no_dataset_and_builds_tune_command():
    from taps_sc_extract.cli import build_tune_command, parse_args

    parsed = parse_args(["--test", "--test-quick", "--save-profile"])
    cmd = build_tune_command("taps-sc-extract-rs", parsed)
    assert cmd[:2] == ["taps-sc-extract-rs", "tune"]
    assert "--quick" in cmd and "--save-profile" in cmd
    assert "-b" not in cmd

    parsed = parse_args(["--test", "-b", "in.bam", "-f", "ref.fa", "-o", "out_dir/run.h5"])
    cmd = build_tune_command("taps-sc-extract-rs", parsed)
    assert cmd[cmd.index("-b") + 1] == "in.bam"
    assert cmd[cmd.index("--out-dir") + 1] == "out_dir"

    with pytest.raises(SystemExit):
        parse_args(["--test", "-b", "in.bam"])  # -f missing
    with pytest.raises(SystemExit):
        parse_args(["-b", "in.bam", "-f", "ref.fa"])  # -o missing outside --test
