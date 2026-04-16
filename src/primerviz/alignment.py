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
    """Slide primer along both strands of the template to find binding sites.

    The primer's direction label (-f / -r) is user metadata for naming only.
    Actual binding orientation is determined by which form matches:

      strand=FORWARD  — primer sequence matches the sense strand directly
                        (primer anneals to antisense strand, extends →)
      strand=REVERSE  — primer's RC matches the sense strand
                        (primer anneals to sense strand, extends ←)

    Both orientations are always searched so that mislabelled or ambiguous
    primers are still found.
    """
    template = template_seq.upper()
    seq = primer.sequence.upper()
    seq_rc = reverse_complement(seq)

    plen = len(seq)
    sites: list[BindingSite] = []

    for i in range(len(template) - plen + 1):
        window = template[i : i + plen]

        # Check sense-strand match (primer binds antisense strand)
        n_mis, mis_pos = _count_mismatches(seq, window)
        if n_mis <= max_mismatches:
            sites.append(
                BindingSite(
                    primer=primer,
                    start=i,
                    end=i + plen,
                    mismatches=n_mis,
                    strand=Direction.FORWARD,
                    mismatch_positions=mis_pos,
                )
            )
            # If both orientations match at the same position (palindrome),
            # only record the forward match to avoid a duplicate.
            continue

        # Check antisense-strand match (primer RC matches sense strand)
        n_mis_rc, mis_pos_rc = _count_mismatches(seq_rc, window)
        if n_mis_rc <= max_mismatches:
            sites.append(
                BindingSite(
                    primer=primer,
                    start=i,
                    end=i + plen,
                    mismatches=n_mis_rc,
                    strand=Direction.REVERSE,
                    mismatch_positions=mis_pos_rc,
                )
            )

    return sites
