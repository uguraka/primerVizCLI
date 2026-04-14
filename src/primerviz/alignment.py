"""Find primer binding sites on a template sequence."""

from __future__ import annotations

from .models import BindingSite, Direction, Primer

_COMPLEMENT = str.maketrans("ATCGatcg", "TAGCtagc")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def _count_mismatches(query: str, target: str) -> tuple[int, list[int]]:
    """Compare two equal-length sequences. Return (count, positions)."""
    positions = [i for i, (a, b) in enumerate(zip(query, target)) if a != b]
    return len(positions), positions


def find_binding_sites(
    primer: Primer,
    template_seq: str,
    max_mismatches: int = 0,
) -> list[BindingSite]:
    """Slide primer along the template (both strands) and find binding sites.

    For forward primers: match directly against the sense strand (5'→3').
    For reverse primers: reverse-complement the primer, then match against
    the sense strand — the binding position is where the RC'd primer aligns.
    """
    template = template_seq.upper()
    seq = primer.sequence.upper()

    if primer.direction == Direction.REVERSE:
        seq = reverse_complement(seq)

    plen = len(seq)
    sites: list[BindingSite] = []

    for i in range(len(template) - plen + 1):
        window = template[i : i + plen]
        n_mis, mis_pos = _count_mismatches(seq, window)
        if n_mis <= max_mismatches:
            sites.append(
                BindingSite(
                    primer=primer,
                    start=i,
                    end=i + plen,
                    mismatches=n_mis,
                    mismatch_positions=mis_pos,
                )
            )

    return sites
