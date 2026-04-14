"""Tests for primer thermodynamic analysis."""

from primerviz.analysis import (
    analyze_primer,
    calc_3prime_stability,
    calc_gc,
    calc_hairpin,
    calc_mw,
    calc_self_dimer,
    calc_tm,
)


def test_gc_poly_a():
    assert calc_gc("AAAAAAAAAA") == 0.0


def test_gc_poly_gc():
    assert calc_gc("GCGCGCGCGC") == 100.0


def test_gc_mixed():
    gc = calc_gc("ATCGATCG")
    assert gc == 50.0


def test_tm_returns_float():
    tm = calc_tm("ATCGATCGATCGATCGATCG")
    assert isinstance(tm, float)
    assert 40.0 < tm < 80.0  # sanity range for a 20-mer


def test_mw_single_base():
    # Single base: no phosphodiester bonds to subtract
    mw = calc_mw("A")
    assert mw == 331.2


def test_mw_increases_with_length():
    assert calc_mw("ATCG") > calc_mw("AT")


def test_hairpin_returns_float():
    # A self-complementary sequence should form a hairpin
    tm = calc_hairpin("GCATGCATGCATGC")
    assert isinstance(tm, float)


def test_self_dimer_returns_float():
    tm = calc_self_dimer("ATCGATCGATCGATCGATCG")
    assert isinstance(tm, float)


def test_3prime_stability_returns_float():
    dg = calc_3prime_stability("ATCGATCGATCGATCGATCG")
    assert isinstance(dg, float)


def test_analyze_primer_all_fields(fwd_primer):
    props = analyze_primer(fwd_primer)
    assert props.length == 20
    assert 0.0 <= props.gc_percent <= 100.0
    assert props.molecular_weight > 0
    assert isinstance(props.tm, float)
    assert isinstance(props.hairpin_tm, float)
    assert isinstance(props.self_dimer_tm, float)
    assert isinstance(props.three_prime_stability, float)
