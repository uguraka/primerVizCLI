"""Shared fixtures for primerviz tests."""

import pytest

from primerviz.models import Direction, Primer


@pytest.fixture
def fwd_primer():
    return Primer(name="fwd_1", sequence="ATCGATCGATCGATCGATCG", direction=Direction.FORWARD)


@pytest.fixture
def rev_primer():
    return Primer(name="rev_1", sequence="CGATCGATCGATCGATCGAT", direction=Direction.REVERSE)


@pytest.fixture
def template_seq():
    return "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
