"""Tests for primer alignment / binding site finding."""

from primerviz.alignment import find_binding_sites, reverse_complement
from primerviz.models import Direction, Primer


def test_reverse_complement():
    assert reverse_complement("ATCG") == "CGAT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("GCGC") == "GCGC"


def test_exact_match_forward(fwd_primer, template_seq):
    sites = find_binding_sites(fwd_primer, template_seq, max_mismatches=0)
    assert len(sites) > 0
    for s in sites:
        assert s.mismatches == 0
        assert template_seq[s.start : s.end] == fwd_primer.sequence.upper()


def test_no_match():
    primer = Primer(name="nope", sequence="NNNNNNNNNN", direction=Direction.FORWARD)
    template = "ATCGATCGATCGATCGATCG"
    sites = find_binding_sites(primer, template, max_mismatches=0)
    assert len(sites) == 0


def test_mismatch_tolerance():
    primer = Primer(name="mm", sequence="ATCGATCGNN", direction=Direction.FORWARD)
    template = "ATCGATCGAT"
    # With 0 mismatches — no match (NN != AT)
    assert len(find_binding_sites(primer, template, max_mismatches=0)) == 0
    # With 2 mismatches — should match
    sites = find_binding_sites(primer, template, max_mismatches=2)
    assert len(sites) == 1
    assert sites[0].mismatches == 2


def test_reverse_primer_binding():
    # Template: ATCGATCG
    # Reverse primer sequence: CGATCGAT -> RC = ATCGATCG -> matches at pos 0
    primer = Primer(name="rev", sequence="CGATCGAT", direction=Direction.REVERSE)
    template = "ATCGATCG"
    sites = find_binding_sites(primer, template, max_mismatches=0)
    assert len(sites) == 1
    assert sites[0].start == 0
