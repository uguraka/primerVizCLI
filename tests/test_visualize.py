"""Smoke tests for visualization rendering."""

from primerviz.alignment import find_binding_sites
from primerviz.analysis import analyze_primer
from primerviz.models import AnalysisResult, Direction, Primer
from primerviz.visualize import render


def _make_result():
    template = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    fwd = Primer(name="fwd_1", sequence="ATCGATCGATCGATCGATCG", direction=Direction.FORWARD)
    rev = Primer(name="rev_1", sequence="CGATCGATCGATCGATCGAT", direction=Direction.REVERSE)

    bindings = find_binding_sites(fwd, template) + find_binding_sites(rev, template)
    props = {
        fwd.name: analyze_primer(fwd),
        rev.name: analyze_primer(rev),
    }
    return AnalysisResult(
        template_name="test_template",
        template_seq=template,
        bindings=bindings,
        primer_props=props,
    )


def test_render_rich_no_crash():
    result = _make_result()
    render(result, use_rich=True)


def test_render_ascii_no_crash():
    result = _make_result()
    render(result, use_rich=False)
