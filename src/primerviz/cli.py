"""CLI entry point for primerVizCLI."""

from __future__ import annotations

from itertools import combinations

import click

from .alignment import find_binding_sites
from .analysis import analyze_primer, calc_primer_dimer
from .io import read_sequence_input
from .models import AnalysisResult, Direction, Primer
from .visualize import render


@click.command()
@click.option(
    "--template", "-t", required=True,
    help="Template sequence or path to a FASTA/text file.",
)
@click.option(
    "--forward", "-f", multiple=True,
    help="Forward primer sequence or file (repeatable).",
)
@click.option(
    "--reverse", "-r", multiple=True,
    help="Reverse primer sequence or file (repeatable).",
)
@click.option(
    "--mismatches", "-m", default=0, show_default=True,
    help="Maximum allowed mismatches when finding binding sites.",
)
@click.option(
    "--plain", is_flag=True, default=False,
    help="Force plain ASCII output (no colors).",
)
def main(
    template: str,
    forward: tuple[str, ...],
    reverse: tuple[str, ...],
    mismatches: int,
    plain: bool,
) -> None:
    """Visualize primer binding sites and thermodynamic properties."""
    if not forward and not reverse:
        raise click.UsageError("At least one --forward or --reverse primer is required.")

    # Parse template
    templates = read_sequence_input(template)

    # Build primer objects
    primers: list[Primer] = []
    for i, fwd_val in enumerate(forward):
        for name, seq in read_sequence_input(fwd_val):
            label = name if name != "input" else f"fwd_{i + 1}"
            primers.append(Primer(name=label, sequence=seq, direction=Direction.FORWARD))
    for i, rev_val in enumerate(reverse):
        for name, seq in read_sequence_input(rev_val):
            label = name if name != "input" else f"rev_{i + 1}"
            primers.append(Primer(name=label, sequence=seq, direction=Direction.REVERSE))

    for tpl_name, tpl_seq in templates:
        # Find binding sites
        bindings = []
        for primer in primers:
            bindings.extend(find_binding_sites(primer, tpl_seq, max_mismatches=mismatches))

        # Compute primer properties
        primer_props = {p.name: analyze_primer(p) for p in primers}

        # Compute primer-dimer Tms for all pairs
        dimer_tms: dict[tuple[str, str], float] = {}
        for p1, p2 in combinations(primers, 2):
            tm = calc_primer_dimer(p1.sequence, p2.sequence)
            dimer_tms[(p1.name, p2.name)] = tm

        result = AnalysisResult(
            template_name=tpl_name,
            template_seq=tpl_seq,
            bindings=bindings,
            primer_props=primer_props,
            primer_dimer_tms=dimer_tms,
        )

        render(result, use_rich=not plain)
