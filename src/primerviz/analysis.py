"""Primer thermodynamic analysis using primer3-py."""

from __future__ import annotations

import primer3

from .models import AnalysisResult, Direction, Primer, PrimerProperties, QpcrCheck, QpcrReport

# Average molecular weights of deoxyribonucleotides (g/mol, anhydrous)
_NUC_MW = {"A": 331.2, "T": 322.2, "G": 347.2, "C": 307.2}
_WATER_MW = 18.02


# ---------------------------------------------------------------------------
# Core thermodynamic calculations
# ---------------------------------------------------------------------------

def calc_tm(seq: str) -> float:
    """Nearest-neighbor Tm (C) using SantaLucia parameters."""
    return primer3.calc_tm(seq)


def calc_gc(seq: str) -> float:
    """GC content as a percentage."""
    seq = seq.upper()
    gc = sum(1 for b in seq if b in "GC")
    return (gc / len(seq)) * 100.0 if seq else 0.0


def calc_mw(seq: str) -> float:
    """Molecular weight of a single-stranded DNA oligo (g/mol)."""
    seq = seq.upper()
    total = sum(_NUC_MW.get(b, 0.0) for b in seq)
    if len(seq) > 1:
        total -= _WATER_MW * (len(seq) - 1)
    return round(total, 1)


def calc_hairpin(seq: str) -> float:
    """Hairpin Tm (C)."""
    return round(primer3.calc_hairpin(seq).tm, 1)


def calc_self_dimer(seq: str) -> float:
    """Homodimer Tm (C)."""
    return round(primer3.calc_homodimer(seq).tm, 1)


def calc_primer_dimer(seq1: str, seq2: str) -> float:
    """Heterodimer Tm (C) between two sequences."""
    return round(primer3.calc_heterodimer(seq1, seq2).tm, 1)


def calc_3prime_stability(seq: str) -> float:
    """Delta-G (kcal/mol) of the last 5 bases."""
    result = primer3.calc_end_stability(seq[-5:], seq)
    return round(result.dg / 1000.0, 2)


def analyze_primer(primer: Primer) -> PrimerProperties:
    """Compute all thermodynamic properties for a primer."""
    seq = primer.sequence.upper()
    return PrimerProperties(
        tm=round(calc_tm(seq), 1),
        gc_percent=round(calc_gc(seq), 1),
        length=len(seq),
        molecular_weight=calc_mw(seq),
        hairpin_tm=calc_hairpin(seq),
        self_dimer_tm=calc_self_dimer(seq),
        three_prime_stability=calc_3prime_stability(seq),
    )


# ---------------------------------------------------------------------------
# qPCR helper checks
# ---------------------------------------------------------------------------

def gc_clamp(seq: str, n: int = 3) -> int:
    """Count G/C bases in the last *n* bases of a sequence (3' clamp)."""
    return sum(1 for b in seq.upper()[-n:] if b in "GC")


def mono_run(seq: str) -> int:
    """Return the length of the longest mono-nucleotide run."""
    if not seq:
        return 0
    max_run = current = 1
    for i in range(1, len(seq)):
        if seq[i].upper() == seq[i - 1].upper():
            current += 1
            max_run = max(max_run, current)
        else:
            current = 1
    return max_run


def _check(name: str, value: str, status: str, detail: str = "") -> QpcrCheck:
    return QpcrCheck(name=name, value=value, status=status, detail=detail)


# ---------------------------------------------------------------------------
# qPCR report
# ---------------------------------------------------------------------------

def run_qpcr_checks(result: AnalysisResult) -> QpcrReport:
    """Evaluate all qPCR design criteria and return a structured report."""
    checks: list[QpcrCheck] = []

    fwd_sites = [b for b in result.bindings if b.strand == Direction.FORWARD]
    rev_sites = [b for b in result.bindings if b.strand == Direction.REVERSE]

    # ── Amplicon ────────────────────────────────────────────────────────────
    amp_len = 0
    amp_gc = 0.0
    amp_seq = ""

    if fwd_sites and rev_sites:
        amp_start = min(s.start for s in fwd_sites)
        amp_end   = max(s.end   for s in rev_sites)
        amp_seq   = result.template_seq[amp_start:amp_end]
        amp_len   = len(amp_seq)
        amp_gc    = round(calc_gc(amp_seq), 1)

        if 80 <= amp_len <= 200:
            status = "pass"
        elif 60 <= amp_len <= 300:
            status = "warn"
        else:
            status = "fail"
        checks.append(_check(
            "Amplicon size", f"{amp_len} bp", status,
            "ideal 80-200 bp, warn 60-300 bp",
        ))
        checks.append(_check(
            "Amplicon GC%", f"{amp_gc}%",
            "pass" if 40 <= amp_gc <= 60 else ("warn" if 35 <= amp_gc <= 65 else "fail"),
            "ideal 40-60%",
        ))
    else:
        checks.append(_check(
            "Amplicon size", "N/A", "warn",
            "need both a forward and reverse primer to determine amplicon",
        ))
        checks.append(_check("Amplicon GC%", "N/A", "warn", "amplicon not determined"))

    # ── Per-primer checks ───────────────────────────────────────────────────
    for name, props in result.primer_props.items():
        seq = _primer_seq_for(name, result)
        prefix = name

        # Tm range
        tm_status = (
            "pass" if 58 <= props.tm <= 65 else
            "warn" if 55 <= props.tm <= 68 else "fail"
        )
        checks.append(_check(f"{prefix} Tm", f"{props.tm}C", tm_status, "ideal 58-65C, warn 55-68C"))

        # GC%
        gc_status = (
            "pass" if 40 <= props.gc_percent <= 60 else
            "warn" if 35 <= props.gc_percent <= 65 else "fail"
        )
        checks.append(_check(f"{prefix} GC%", f"{props.gc_percent}%", gc_status, "ideal 40-60%"))

        if seq:
            # 3' GC clamp (last 3 bases)
            clamp = gc_clamp(seq, 3)
            clamp_status = "pass" if 1 <= clamp <= 2 else "warn"
            checks.append(_check(
                f"{prefix} 3' clamp", f"{clamp}/3 G/C", clamp_status,
                "1-2 G/C in last 3 bases is ideal",
            ))

            # Mono-nucleotide runs
            run = mono_run(seq)
            run_status = "pass" if run <= 3 else ("warn" if run == 4 else "fail")
            checks.append(_check(
                f"{prefix} mono-run", f"{run} bp", run_status,
                "max 3 consecutive identical bases",
            ))

        # Self-dimer
        sd_status = "pass" if props.self_dimer_tm <= 40 else ("warn" if props.self_dimer_tm <= 45 else "fail")
        checks.append(_check(f"{prefix} self-dimer Tm", f"{props.self_dimer_tm}C", sd_status, "ideal <=40C, warn <=45C"))

        # Hairpin
        hp_status = "pass" if props.hairpin_tm <= 40 else ("warn" if props.hairpin_tm <= 45 else "fail")
        checks.append(_check(f"{prefix} hairpin Tm", f"{props.hairpin_tm}C", hp_status, "ideal <=40C, warn <=45C"))

    # ── Tm delta between fwd and rev ────────────────────────────────────────
    fwd_names = {b.primer.name for b in fwd_sites}
    rev_names = {b.primer.name for b in rev_sites}
    fwd_tms = [result.primer_props[n].tm for n in fwd_names if n in result.primer_props]
    rev_tms = [result.primer_props[n].tm for n in rev_names if n in result.primer_props]

    if fwd_tms and rev_tms:
        delta = round(abs(sum(fwd_tms) / len(fwd_tms) - sum(rev_tms) / len(rev_tms)), 1)
        delta_status = "pass" if delta <= 1 else ("warn" if delta <= 3 else "fail")
        checks.append(_check("Tm delta (fwd-rev)", f"{delta}C", delta_status, "ideal <=1C, warn <=3C"))

    # ── Primer-dimer pairs ──────────────────────────────────────────────────
    for (p1, p2), tm in result.primer_dimer_tms.items():
        pd_status = "pass" if tm <= 40 else ("warn" if tm <= 45 else "fail")
        checks.append(_check(f"Primer-dimer {p1}/{p2}", f"{tm}C", pd_status, "ideal <=40C, warn <=45C"))

    # ── Probe checks ────────────────────────────────────────────────────────
    if result.probe_bindings and result.probe_props:
        probe_props = result.probe_props
        probe_seq = result.probe_bindings[0].primer.sequence

        # Probe Tm vs average primer Tm
        all_primer_tms = [p.tm for p in result.primer_props.values()]
        avg_primer_tm = sum(all_primer_tms) / len(all_primer_tms) if all_primer_tms else 0
        probe_delta = round(probe_props.tm - avg_primer_tm, 1)
        pd_tm_status = "pass" if probe_delta >= 8 else ("warn" if probe_delta >= 5 else "fail")
        checks.append(_check(
            "Probe Tm above primers",
            f"+{probe_delta}C ({probe_props.tm}C vs avg {avg_primer_tm:.1f}C)",
            pd_tm_status,
            "probe should be >=8C above primer Tm",
        ))

        # Probe must not start with G
        starts_with_g = probe_seq.upper().startswith("G")
        checks.append(_check(
            "Probe 5' base", probe_seq[0].upper(),
            "fail" if starts_with_g else "pass",
            "probe must not start with G (quenching interference)",
        ))

        # Probe within amplicon
        if amp_len > 0:
            amp_start = min(s.start for s in fwd_sites)
            amp_end   = max(s.end   for s in rev_sites)
            for pb in result.probe_bindings:
                within = amp_start <= pb.start and pb.end <= amp_end
                checks.append(_check(
                    "Probe within amplicon",
                    f"pos {pb.start + 1}..{pb.end}",
                    "pass" if within else "fail",
                    f"amplicon is pos {amp_start + 1}..{amp_end}",
                ))

        # Probe overlap with primers
        overlapping = [
            ps for ps in result.bindings
            for pb in result.probe_bindings
            if _overlaps(pb, ps)
        ]
        if overlapping:
            names = ", ".join({ps.primer.name for ps in overlapping})
            checks.append(_check(
                "Probe-primer overlap", f"overlaps with: {names}", "fail",
                "probe and primers cannot bind at the same position",
            ))
        else:
            checks.append(_check(
                "Probe-primer overlap", "none", "pass",
                "probe does not overlap any primer",
            ))

    return QpcrReport(checks=checks, amplicon_gc=amp_gc, amplicon_len=amp_len)


def _primer_seq_for(name: str, result: AnalysisResult) -> str:
    for b in result.bindings:
        if b.primer.name == name:
            return b.primer.sequence
    return ""


def _overlaps(a, b) -> bool:
    """Return True if two BindingSites overlap on the template coordinate system."""
    return a.start < b.end and b.start < a.end
