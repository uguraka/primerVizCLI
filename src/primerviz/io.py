"""Sequence input parsing — FASTA files and raw sequences."""

from __future__ import annotations

import os
from pathlib import Path


def parse_fasta(filepath: str | Path) -> list[tuple[str, str]]:
    """Parse a FASTA file into a list of (name, sequence) tuples."""
    records: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []

    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts)))
                name = line[1:].split()[0]  # first word after >
                parts = []
            else:
                parts.append(line.upper())

    if name is not None:
        records.append((name, "".join(parts)))

    return records


def read_sequence_input(value: str) -> list[tuple[str, str]]:
    """Auto-detect input: if a file path, parse it; otherwise treat as raw sequence."""
    if os.path.isfile(value):
        ext = Path(value).suffix.lower()
        if ext in (".fa", ".fasta", ".fna", ".fas"):
            return parse_fasta(value)
        # Plain text file: one sequence per line (no headers)
        with open(value) as fh:
            seq = "".join(line.strip().upper() for line in fh if line.strip())
            name = Path(value).stem
            return [(name, seq)]

    # Raw sequence string
    return [("input", value.upper())]


def invalid_bases(seq: str) -> list[tuple[int, str]]:
    """Return (1-based position, base) for every non-ACGT character in seq."""
    return [(i + 1, b) for i, b in enumerate(seq) if b not in "ACGT"]
