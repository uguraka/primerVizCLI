"""Primer thermodynamic analysis using primer3-py."""

from __future__ import annotations

import primer3

from .models import Primer, PrimerProperties

# Average molecular weights of deoxyribonucleotides (g/mol, anhydrous)
_NUC_MW = {"A": 331.2, "T": 322.2, "G": 347.2, "C": 307.2}
_WATER_MW = 18.02


def calc_tm(seq: str) -> float:
    """Nearest-neighbor Tm (°C) using SantaLucia parameters."""
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
    # Subtract water for each phosphodiester bond (len - 1)
    if len(seq) > 1:
        total -= _WATER_MW * (len(seq) - 1)
    return round(total, 1)


def calc_hairpin(seq: str) -> float:
    """Hairpin Tm (°C). Returns 0.0 if no significant structure."""
    result = primer3.calc_hairpin(seq)
    return round(result.tm, 1)


def calc_self_dimer(seq: str) -> float:
    """Homodimer Tm (°C)."""
    result = primer3.calc_homodimer(seq)
    return round(result.tm, 1)


def calc_primer_dimer(seq1: str, seq2: str) -> float:
    """Heterodimer Tm (°C) between two sequences."""
    result = primer3.calc_heterodimer(seq1, seq2)
    return round(result.tm, 1)


def calc_3prime_stability(seq: str) -> float:
    """Delta-G (kcal/mol) of the last 5 bases — measures 3' end stability."""
    result = primer3.calc_end_stability(seq[-5:], seq)
    return round(result.dg / 1000.0, 2)  # primer3 returns cal/mol


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
