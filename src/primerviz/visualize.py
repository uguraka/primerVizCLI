"""Terminal visualization of primer binding — Rich and ASCII modes."""

from __future__ import annotations

import sys

from .models import AnalysisResult, BindingSite, Direction


def render(result: AnalysisResult, use_rich: bool = True) -> None:
    """Dispatch to the appropriate renderer."""
    if use_rich:
        try:
            _render_rich(result)
            return
        except ImportError:
            pass
    _render_ascii(result)


# ---------------------------------------------------------------------------
# Rich renderer
# ---------------------------------------------------------------------------

def _render_rich(result: AnalysisResult) -> None:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table
    from rich.text import Text

    console = Console()
    tpl = result.template_seq
    tpl_len = len(tpl)

    # Header
    console.print()
    console.print(
        Panel(f"[bold]{result.template_name}[/bold]  ({tpl_len} bp)", expand=False)
    )

    # Ruler
    ruler = _build_ruler(tpl_len)
    console.print(Text(ruler, style="dim"))

    # Template strand
    tpl_line = Text()
    tpl_line.append("5'- ", style="dim")
    tpl_line.append(tpl)
    tpl_line.append(" -3'", style="dim")
    console.print(tpl_line)

    # Binding lines — forward above, reverse below
    fwd_sites = [b for b in result.bindings if b.primer.direction == Direction.FORWARD]
    rev_sites = [b for b in result.bindings if b.primer.direction == Direction.REVERSE]

    for site in fwd_sites:
        console.print(_binding_line_rich(site, tpl_len, ">"))
    for site in rev_sites:
        console.print(_binding_line_rich(site, tpl_len, "<"))

    console.print()

    # Properties table
    table = Table(title="Primer Properties", show_header=True, header_style="bold cyan")
    table.add_column("Name")
    table.add_column("Dir", justify="center")
    table.add_column("Len", justify="right")
    table.add_column("Tm (C)", justify="right")
    table.add_column("GC%", justify="right")
    table.add_column("MW", justify="right")
    table.add_column("Hairpin Tm", justify="right")
    table.add_column("Self-dimer Tm", justify="right")
    table.add_column("3' dG", justify="right")

    for name, props in result.primer_props.items():
        # Find direction from bindings or primers
        direction = _direction_for(name, result)
        table.add_row(
            name,
            direction,
            str(props.length),
            f"{props.tm:.1f}",
            f"{props.gc_percent:.1f}",
            f"{props.molecular_weight:.1f}",
            f"{props.hairpin_tm:.1f}",
            f"{props.self_dimer_tm:.1f}",
            f"{props.three_prime_stability:.2f}",
        )

    console.print(table)

    # Primer-dimer table (if any pairs)
    if result.primer_dimer_tms:
        dimer_table = Table(title="Primer-Dimer Tm (C)", show_header=True, header_style="bold yellow")
        dimer_table.add_column("Primer 1")
        dimer_table.add_column("Primer 2")
        dimer_table.add_column("Tm (C)", justify="right")
        for (p1, p2), tm in result.primer_dimer_tms.items():
            dimer_table.add_row(p1, p2, f"{tm:.1f}")
        console.print(dimer_table)

    console.print()


def _binding_line_rich(site: BindingSite, tpl_len: int, char: str) -> "Text":
    from rich.text import Text

    padding = "    "  # matches "5'- " prefix
    line = Text()
    line.append(padding)
    line.append(" " * site.start)

    primer_len = site.end - site.start
    for i in range(primer_len):
        if i in site.mismatch_positions:
            line.append("x", style="bold red")
        else:
            line.append(char, style="bold green")

    label = f"  {site.primer.name} (Tm: {_get_tm(site)}C)"
    if site.mismatches > 0:
        label += f" [{site.mismatches} mm]"
    line.append(label, style="dim")
    return line


def _get_tm(site: BindingSite) -> str:
    """Helper to format Tm — we don't store it on BindingSite, so re-compute."""
    from .analysis import calc_tm
    return f"{calc_tm(site.primer.sequence):.1f}"


# ---------------------------------------------------------------------------
# ASCII fallback
# ---------------------------------------------------------------------------

def _render_ascii(result: AnalysisResult) -> None:
    tpl = result.template_seq
    tpl_len = len(tpl)
    out = sys.stdout

    out.write(f"\nTemplate: {result.template_name} ({tpl_len} bp)\n")
    out.write(_build_ruler(tpl_len) + "\n")
    out.write(f"5'- {tpl} -3'\n")

    fwd_sites = [b for b in result.bindings if b.primer.direction == Direction.FORWARD]
    rev_sites = [b for b in result.bindings if b.primer.direction == Direction.REVERSE]

    for site in fwd_sites:
        out.write(_binding_line_ascii(site, ">") + "\n")
    for site in rev_sites:
        out.write(_binding_line_ascii(site, "<") + "\n")

    out.write("\nPrimer Properties:\n")
    out.write(f"{'Name':<12} {'Dir':<5} {'Len':>4} {'Tm':>7} {'GC%':>6} {'MW':>9} {'HP Tm':>7} {'SD Tm':>7} {'3dG':>7}\n")
    out.write("-" * 72 + "\n")

    for name, props in result.primer_props.items():
        direction = _direction_for(name, result)
        out.write(
            f"{name:<12} {direction:<5} {props.length:>4} "
            f"{props.tm:>7.1f} {props.gc_percent:>6.1f} {props.molecular_weight:>9.1f} "
            f"{props.hairpin_tm:>7.1f} {props.self_dimer_tm:>7.1f} "
            f"{props.three_prime_stability:>7.2f}\n"
        )

    if result.primer_dimer_tms:
        out.write("\nPrimer-Dimer Tm (°C):\n")
        for (p1, p2), tm in result.primer_dimer_tms.items():
            out.write(f"  {p1} / {p2}: {tm:.1f}\n")

    out.write("\n")


def _binding_line_ascii(site: BindingSite, char: str) -> str:
    padding = "    "  # matches "5'- "
    primer_len = site.end - site.start
    arrows = ""
    for i in range(primer_len):
        arrows += "x" if i in site.mismatch_positions else char
    label = f"  {site.primer.name}"
    return padding + " " * site.start + arrows + label


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _build_ruler(length: int) -> str:
    """Build a position ruler line: tick marks every 10 bp."""
    padding = "    "  # matches "5'- "
    ruler_chars = []
    for i in range(1, length + 1):
        if i % 10 == 0:
            label = str(i)
            # Place label so the last digit lands on the tick
            ruler_chars.append(label[-1])
        elif (i + 1) % 10 == 0 and i + 1 <= length:
            # Peek: next position is a tick — skip if label will overlap
            ruler_chars.append(" ")
        else:
            ruler_chars.append(" ")

    # Build a proper ruler with numbers positioned above ticks
    nums = []
    ticks = []
    i = 0
    while i < length:
        if (i + 1) % 10 == 0:
            label = str(i + 1)
            # Right-align the number above the tick
            start = max(0, len(nums) - (len(label) - 1))
            # Overwrite previous spaces with the number
            nums = nums[:start]
            nums.extend(list(label))
            ticks.append("|")
        else:
            nums.append(" ")
            ticks.append(" ")
        i += 1

    return padding + "".join(nums) + "\n" + padding + "".join(ticks)


def _direction_for(name: str, result: AnalysisResult) -> str:
    """Look up direction string for a primer name from the bindings list."""
    for b in result.bindings:
        if b.primer.name == name:
            return "fwd" if b.primer.direction == Direction.FORWARD else "rev"
    # Fallback: infer from name
    if "fwd" in name.lower() or "forward" in name.lower():
        return "fwd"
    if "rev" in name.lower() or "reverse" in name.lower():
        return "rev"
    return "?"
