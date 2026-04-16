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
        if s.strand == Direction.FORWARD:
            assert template_seq[s.start : s.end] == fwd_primer.sequence.upper()
    # At least one hit must be on the forward (sense) strand
    assert any(s.strand == Direction.FORWARD for s in sites)


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
    # With 2 mismatches — should match (forward orientation)
    sites = find_binding_sites(primer, template, max_mismatches=2)
    assert len(sites) == 1
    assert sites[0].mismatches == 2
    assert sites[0].strand == Direction.FORWARD


def test_reverse_primer_binding():
    # Primer sequence: CGATCGAT — its RC (ATCGATCG) matches the sense strand at pos 0
    # → strand=REVERSE (primer anneals to sense strand)
    primer = Primer(name="rev", sequence="CGATCGAT", direction=Direction.REVERSE)
    template = "ATCGATCG"
    sites = find_binding_sites(primer, template, max_mismatches=0)
    assert len(sites) == 1
    assert sites[0].start == 0
    assert sites[0].strand == Direction.REVERSE


def test_both_strands_searched_regardless_of_label():
    """A primer labelled -f that only matches as RC should still be found."""
    primer = Primer(name="mislabelled_fwd", sequence="CGATCGAT", direction=Direction.FORWARD)
    template = "ATCGATCG"
    sites = find_binding_sites(primer, template, max_mismatches=0)
    assert len(sites) == 1
    assert sites[0].strand == Direction.REVERSE
