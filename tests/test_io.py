"""Tests for sequence I/O parsing."""

import tempfile
from pathlib import Path

from primerviz.io import parse_fasta, read_sequence_input


def test_parse_fasta_single():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(">seq1\n")
        f.write("ATCGATCG\n")
        f.write("ATCGATCG\n")
        f.flush()
        records = parse_fasta(f.name)

    assert len(records) == 1
    assert records[0] == ("seq1", "ATCGATCGATCGATCG")


def test_parse_fasta_multi():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(">first\nATCG\n>second\nGCGC\n")
        f.flush()
        records = parse_fasta(f.name)

    assert len(records) == 2
    assert records[0][0] == "first"
    assert records[1][0] == "second"


def test_read_sequence_input_raw():
    records = read_sequence_input("atcgatcg")
    assert len(records) == 1
    assert records[0] == ("input", "ATCGATCG")


def test_read_sequence_input_fasta_file():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(">myseq\nGCATGCAT\n")
        f.flush()
        records = read_sequence_input(f.name)

    assert len(records) == 1
    assert records[0] == ("myseq", "GCATGCAT")
