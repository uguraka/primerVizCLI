"""CLI entry point for primerVizCLI."""

from __future__ import annotations

from itertools import combinations

import click

from .alignment import find_binding_sites
from .analysis import analyze_primer, calc_primer_dimer, run_qpcr_checks
from .io import invalid_bases, read_sequence_input
from .models import AnalysisResult, Direction, Primer
from .visualize import render


def _check_seq(label: str, seq: str, param: str) -> None:
    bad = invalid_bases(seq)
    if not bad:
        return
    pos, base = bad[0]
    extra = f" (and {len(bad) - 1} more)" if len(bad) > 1 else ""
    raise click.BadParameter(
        f"invalid base '{base}' at position {pos}{extra} — only A, C, G, T are allowed.",
        param_hint=f"--{param} '{label}'",
    )


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
    "--probe", "-p", default=None,
    help="Probe sequence or file (optional; enables qPCR report automatically).",
)
@click.option(
    "--mismatches", "-m", default=0, show_default=True,
    help="Maximum allowed mismatches when finding binding sites.",
)
@click.option(
    "--qpcr", is_flag=True, default=False,
    help="Run qPCR compatibility checks and show a pass/warn/fail report.",
)
@click.option(
    "--plain", is_flag=True, default=False,
    help="Force plain ASCII output (no colors).",
)
def main(
    template: str,
    forward: tuple[str, ...],
    reverse: tuple[str, ...],
    probe: str | None,
    mismatches: int,
    qpcr: bool,
    plain: bool,
) -> None:
    """Visualize primer binding sites and thermodynamic properties."""
    if not forward and not reverse:
        raise click.UsageError("At least one --forward or --reverse primer is required.")

    templates = read_sequence_input(template)
    for tpl_name, tpl_seq in templates:
        _check_seq(tpl_name, tpl_seq, "template")

    # Build primer objects
    primers: list[Primer] = []
    for i, fwd_val in enumerate(forward):
        for name, seq in read_sequence_input(fwd_val):
            label = name if name != "input" else f"fwd_{i + 1}"
            _check_seq(label, seq, "forward")
            primers.append(Primer(name=label, sequence=seq, direction=Direction.FORWARD))
    for i, rev_val in enumerate(reverse):
        for name, seq in read_sequence_input(rev_val):
            label = name if name != "input" else f"rev_{i + 1}"
            _check_seq(label, seq, "reverse")
            primers.append(Primer(name=label, sequence=seq, direction=Direction.REVERSE))

    # Build probe object (if supplied)
    probe_primer: Primer | None = None
    if probe:
        records = read_sequence_input(probe)
        name, seq = records[0]
        label = name if name != "input" else "probe"
        _check_seq(label, seq, "probe")
        probe_primer = Primer(name=label, sequence=seq, direction=Direction.FORWARD)

    for tpl_name, tpl_seq in templates:
        # Primer binding sites
        bindings = []
        for primer in primers:
            bindings.extend(find_binding_sites(primer, tpl_seq, max_mismatches=mismatches))

        # Primer properties + dimer Tms
        primer_props = {p.name: analyze_primer(p) for p in primers}
        dimer_tms: dict[tuple[str, str], float] = {}
        for p1, p2 in combinations(primers, 2):
            dimer_tms[(p1.name, p2.name)] = calc_primer_dimer(p1.sequence, p2.sequence)

        # Probe binding sites + properties
        probe_bindings = []
        probe_props = None
        if probe_primer:
            probe_bindings = find_binding_sites(probe_primer, tpl_seq, max_mismatches=mismatches)
            probe_props = analyze_primer(probe_primer)

        result = AnalysisResult(
            template_name=tpl_name,
            template_seq=tpl_seq,
            bindings=bindings,
            primer_props=primer_props,
            primer_dimer_tms=dimer_tms,
            probe_bindings=probe_bindings,
            probe_props=probe_props,
        )

        # qPCR report — enabled explicitly or implicitly when a probe is given
        if qpcr or probe_primer:
            result.qpcr_report = run_qpcr_checks(result)

        render(result, use_rich=not plain)
