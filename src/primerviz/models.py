"""Data models for primer analysis and visualization."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class Direction(Enum):
    FORWARD = "forward"
    REVERSE = "reverse"


@dataclass
class Primer:
    name: str
    sequence: str
    direction: Direction


@dataclass
class PrimerProperties:
    tm: float
    gc_percent: float
    length: int
    molecular_weight: float
    hairpin_tm: float
    self_dimer_tm: float
    three_prime_stability: float  # delta-G (kcal/mol)


@dataclass
class BindingSite:
    primer: Primer
    start: int  # 0-based position on template
    end: int
    mismatches: int
    mismatch_positions: list[int] = field(default_factory=list)


@dataclass
class AnalysisResult:
    template_name: str
    template_seq: str
    bindings: list[BindingSite]
    primer_props: dict[str, PrimerProperties]  # keyed by primer name
    primer_dimer_tms: dict[tuple[str, str], float] = field(default_factory=dict)
