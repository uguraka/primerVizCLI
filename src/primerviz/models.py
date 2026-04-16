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
    strand: Direction = Direction.FORWARD   # FORWARD = primer matches sense strand (binds antisense, >>>)
                                            # REVERSE = primer RC matches sense strand (binds sense, <<<)
    mismatch_positions: list[int] = field(default_factory=list)


# ---------------------------------------------------------------------------
# qPCR
# ---------------------------------------------------------------------------

@dataclass
class QpcrCheck:
    name: str        # human-readable criterion name
    value: str       # formatted measured value
    status: str      # "pass" | "warn" | "fail"
    detail: str = "" # optional explanation / threshold info


@dataclass
class QpcrReport:
    checks: list[QpcrCheck]
    amplicon_gc: float   # GC% of the amplicon sequence
    amplicon_len: int    # bp


@dataclass
class AnalysisResult:
    template_name: str
    template_seq: str
    bindings: list[BindingSite]
    primer_props: dict[str, PrimerProperties]  # keyed by primer name
    primer_dimer_tms: dict[tuple[str, str], float] = field(default_factory=dict)
    # Optional probe (hydrolysis) — kept separate from primer bindings
    probe_bindings: list[BindingSite] = field(default_factory=list)
    probe_props: PrimerProperties | None = None
    # qPCR report — populated when --qpcr flag is used or probe is given
    qpcr_report: QpcrReport | None = None
