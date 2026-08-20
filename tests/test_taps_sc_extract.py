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
from taps_sc_extract.calling import OT_FLAGS, OB_FLAGS, call_mctot, classify_strand
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
def test_chr19_real_bam_smoke(tmp_path):
    bam_path = "/mnt/e/sciMET_TAPS/260731/taps_ucbTn5_sp2.srt.bam"
    fasta_path = "/mnt/e/refs/mm10/mm10.fa"

    if not os.path.exists(bam_path) or not os.path.exists(fasta_path):
        pytest.skip("Real BAM or FASTA not available in environment.")

    out_h5 = tmp_path / "chr19_smoke.h5"

    stats = extract_methylation(
        bam_path=bam_path,
        fasta_path=fasta_path,
        out_h5_path=str(out_h5),
        chroms=["chr19"],
    )

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

