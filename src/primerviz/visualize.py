"""Terminal visualization of primer binding — Rich and ASCII modes."""

from __future__ import annotations

import sys
from dataclasses import dataclass

from .models import AnalysisResult, BindingSite, Direction, QpcrReport


# ---------------------------------------------------------------------------
# Amplicon detection
# ---------------------------------------------------------------------------

@dataclass
class _AmpliconView:
    seq: str
    offset: int    # 0-based start in the original template
    tpl_len: int
    is_full: bool  # True when no valid amplicon was found


def _get_amplicon(result: AnalysisResult) -> _AmpliconView:
    """Trim view to the outermost fwd/rev binding pair (the amplicon).

    Falls back to the full template when only one primer type is present.
    """
    tpl = result.template_seq
    tpl_len = len(tpl)

    fwd_sites = [b for b in result.bindings if b.strand == Direction.FORWARD]
    rev_sites = [b for b in result.bindings if b.strand == Direction.REVERSE]

    if not fwd_sites or not rev_sites:
        return _AmpliconView(seq=tpl, offset=0, tpl_len=tpl_len, is_full=True)

    amp_start = min(s.start for s in fwd_sites)
    amp_end   = max(s.end   for s in rev_sites)

    if amp_start >= amp_end:
        return _AmpliconView(seq=tpl, offset=0, tpl_len=tpl_len, is_full=True)

    return _AmpliconView(
        seq=tpl[amp_start:amp_end],
        offset=amp_start,
        tpl_len=tpl_len,
        is_full=False,
    )


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

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
    amp = _get_amplicon(result)
    amp_len = len(amp.seq)

    # ── Header ──────────────────────────────────────────────────────────────
    if amp.is_full:
        title = f"[bold]{result.template_name}[/bold]  ({amp.tpl_len} bp)"
    else:
        title = (
            f"[bold]{result.template_name}[/bold]  "
            f"amplicon: {amp_len} bp  "
            f"[dim](pos {amp.offset + 1}..{amp.offset + amp_len} of {amp.tpl_len} bp)[/dim]"
        )
    console.print()
    console.print(Panel(title, expand=False))

    # ── Sequence view ────────────────────────────────────────────────────────
    ruler = _build_ruler(amp_len, base=amp.offset + 1)
    console.print(Text(ruler, style="dim"))

    tpl_line = Text()
    tpl_line.append("5'- ", style="dim")
    tpl_line.append(amp.seq)
    tpl_line.append(" -3'", style="dim")
    console.print(tpl_line)

    # Forward primers (>>>)
    for site in (b for b in result.bindings if b.strand == Direction.FORWARD):
        console.print(_primer_line_rich(site, amp.offset, ">", "bold green"))

    # Probe (===)
    for site in result.probe_bindings:
        console.print(_probe_line_rich(site, amp.offset))

    # Reverse primers (<<<)
    for site in (b for b in result.bindings if b.strand == Direction.REVERSE):
        console.print(_primer_line_rich(site, amp.offset, "<", "bold green"))

    console.print()

    # ── Primer properties table ──────────────────────────────────────────────
    table = Table(title="Primer Properties", show_header=True, header_style="bold cyan")
    table.add_column("Name")
    table.add_column("Role", justify="center")
    table.add_column("Len", justify="right")
    table.add_column("Tm (C)", justify="right")
    table.add_column("GC%", justify="right")
    table.add_column("MW", justify="right")
    table.add_column("Hairpin Tm", justify="right")
    table.add_column("Self-dimer Tm", justify="right")
    table.add_column("3' dG", justify="right")

    for name, props in result.primer_props.items():
        role = _direction_for(name, result)
        table.add_row(
            name, role, str(props.length),
            f"{props.tm:.1f}", f"{props.gc_percent:.1f}",
            f"{props.molecular_weight:.1f}", f"{props.hairpin_tm:.1f}",
            f"{props.self_dimer_tm:.1f}", f"{props.three_prime_stability:.2f}",
        )

    if result.probe_props and result.probe_bindings:
        p = result.probe_props
        table.add_row(
            result.probe_bindings[0].primer.name, "probe", str(p.length),
            f"{p.tm:.1f}", f"{p.gc_percent:.1f}",
            f"{p.molecular_weight:.1f}", f"{p.hairpin_tm:.1f}",
            f"{p.self_dimer_tm:.1f}", f"{p.three_prime_stability:.2f}",
        )

    console.print(table)

    # ── Primer-dimer table ───────────────────────────────────────────────────
    if result.primer_dimer_tms:
        dt = Table(title="Primer-Dimer Tm (C)", show_header=True, header_style="bold yellow")
        dt.add_column("Primer 1")
        dt.add_column("Primer 2")
        dt.add_column("Tm (C)", justify="right")
        for (p1, p2), tm in result.primer_dimer_tms.items():
            dt.add_row(p1, p2, f"{tm:.1f}")
        console.print(dt)

    # ── qPCR report card ─────────────────────────────────────────────────────
    if result.qpcr_report:
        _render_qpcr_rich(console, result.qpcr_report)

    console.print()


def _primer_line_rich(site: BindingSite, offset: int, char: str, style: str):
    from rich.text import Text

    padding = "    "
    line = Text()
    line.append(padding)
    line.append(" " * (site.start - offset))

    for i in range(site.end - site.start):
        line.append(
            "x" if i in site.mismatch_positions else char,
            style="bold red" if i in site.mismatch_positions else style,
        )

    label = f"  {site.primer.name} (Tm: {_get_tm(site)}C)"
    if site.mismatches > 0:
        label += f" [{site.mismatches} mm]"
    line.append(label, style="dim")
    return line


def _probe_line_rich(site: BindingSite, offset: int):
    from rich.text import Text

    padding = "    "
    line = Text()
    line.append(padding)
    line.append(" " * (site.start - offset))

    for i in range(site.end - site.start):
        line.append(
            "x" if i in site.mismatch_positions else "=",
            style="bold red" if i in site.mismatch_positions else "bold magenta",
        )

    line.append(f"  {site.primer.name} [probe] (Tm: {_get_tm(site)}C)", style="dim magenta")
    return line


def _render_qpcr_rich(console, report: QpcrReport) -> None:
    from rich.table import Table

    _STYLE = {"pass": "bold green", "warn": "bold yellow", "fail": "bold red"}
    _ICON  = {"pass": "[PASS]", "warn": "[WARN]", "fail": "[FAIL]"}

    table = Table(title="qPCR Compatibility Report", show_header=True, header_style="bold white")
    table.add_column("", width=6, no_wrap=True)
    table.add_column("Criterion")
    table.add_column("Value")
    table.add_column("Threshold", style="dim")

    for chk in report.checks:
        s = _STYLE[chk.status]
        table.add_row(
            f"[{s}]{_ICON[chk.status]}[/{s}]",
            chk.name,
            f"[{s}]{chk.value}[/{s}]",
            chk.detail,
        )
    console.print(table)


def _get_tm(site: BindingSite) -> str:
    from .analysis import calc_tm
    return f"{calc_tm(site.primer.sequence):.1f}"


# ---------------------------------------------------------------------------
# ASCII fallback
# ---------------------------------------------------------------------------

def _render_ascii(result: AnalysisResult) -> None:
    amp = _get_amplicon(result)
    amp_len = len(amp.seq)
    out = sys.stdout

    if amp.is_full:
        out.write(f"\nTemplate: {result.template_name} ({amp.tpl_len} bp)\n")
    else:
        out.write(
            f"\nTemplate: {result.template_name}  "
            f"amplicon: {amp_len} bp  "
            f"(pos {amp.offset + 1}..{amp.offset + amp_len} of {amp.tpl_len} bp)\n"
        )

    out.write(_build_ruler(amp_len, base=amp.offset + 1) + "\n")
    out.write(f"5'- {amp.seq} -3'\n")

    for site in (b for b in result.bindings if b.strand == Direction.FORWARD):
        out.write(_primer_line_ascii(site, amp.offset, ">") + "\n")
    for site in result.probe_bindings:
        out.write(_probe_line_ascii(site, amp.offset) + "\n")
    for site in (b for b in result.bindings if b.strand == Direction.REVERSE):
        out.write(_primer_line_ascii(site, amp.offset, "<") + "\n")

    out.write("\nPrimer Properties:\n")
    out.write(
        f"{'Name':<12} {'Role':<6} {'Len':>4} {'Tm':>7} {'GC%':>6} "
        f"{'MW':>9} {'HP Tm':>7} {'SD Tm':>7} {'3dG':>7}\n"
    )
    out.write("-" * 76 + "\n")

    for name, props in result.primer_props.items():
        role = _direction_for(name, result)
        out.write(
            f"{name:<12} {role:<6} {props.length:>4} "
            f"{props.tm:>7.1f} {props.gc_percent:>6.1f} {props.molecular_weight:>9.1f} "
            f"{props.hairpin_tm:>7.1f} {props.self_dimer_tm:>7.1f} "
            f"{props.three_prime_stability:>7.2f}\n"
        )

    if result.probe_props and result.probe_bindings:
        p = result.probe_props
        pname = result.probe_bindings[0].primer.name
        out.write(
            f"{pname:<12} {'probe':<6} {p.length:>4} "
            f"{p.tm:>7.1f} {p.gc_percent:>6.1f} {p.molecular_weight:>9.1f} "
            f"{p.hairpin_tm:>7.1f} {p.self_dimer_tm:>7.1f} "
            f"{p.three_prime_stability:>7.2f}\n"
        )

    if result.primer_dimer_tms:
        out.write("\nPrimer-Dimer Tm (C):\n")
        for (p1, p2), tm in result.primer_dimer_tms.items():
            out.write(f"  {p1} / {p2}: {tm:.1f}\n")

    if result.qpcr_report:
        _render_qpcr_ascii(out, result.qpcr_report)

    out.write("\n")


def _primer_line_ascii(site: BindingSite, offset: int, char: str) -> str:
    padding = "    "
    arrows = "".join("x" if i in site.mismatch_positions else char for i in range(site.end - site.start))
    return padding + " " * (site.start - offset) + arrows + f"  {site.primer.name}"


def _probe_line_ascii(site: BindingSite, offset: int) -> str:
    padding = "    "
    marks = "".join("x" if i in site.mismatch_positions else "=" for i in range(site.end - site.start))
    return padding + " " * (site.start - offset) + marks + f"  {site.primer.name} [probe]"


def _render_qpcr_ascii(out, report: QpcrReport) -> None:
    _ICON = {"pass": "[PASS]", "warn": "[WARN]", "fail": "[FAIL]"}
    out.write("\nqPCR Compatibility Report:\n")
    out.write("-" * 76 + "\n")
    for chk in report.checks:
        out.write(f"  {_ICON[chk.status]:<7} {chk.name:<30} {chk.value:<20} {chk.detail}\n")


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _build_ruler(length: int, base: int = 1) -> str:
    """Two-line ruler with ticks at multiples of 10 in genomic coordinates."""
    padding = "    "
    nums:  list[str] = [" "] * length
    ticks: list[str] = [" "] * length

    for i in range(length):
        pos = base + i
        if pos % 10 == 0:
            ticks[i] = "|"
            label = str(pos)
            label_start = i - (len(label) - 1)
            for j, ch in enumerate(label):
                idx = label_start + j
                if 0 <= idx < length:
                    nums[idx] = ch

    return padding + "".join(nums) + "\n" + padding + "".join(ticks)


def _direction_for(name: str, result: AnalysisResult) -> str:
    for b in result.bindings:
        if b.primer.name == name:
            return "fwd" if b.strand == Direction.FORWARD else "rev"
    if "fwd" in name.lower() or "forward" in name.lower():
        return "fwd"
    if "rev" in name.lower() or "reverse" in name.lower():
        return "rev"
    return "?"
