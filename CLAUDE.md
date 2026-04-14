# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

primerVizCLI is a CLI tool for visualizing DNA primer binding to template sequences. It shows where primers anneal and displays thermodynamic properties (Tm via primer3-py, GC%, MW, hairpin, self-dimer, primer-dimer, 3' stability).

## Development Environment

- Python 3.12, venv at `.venv/`
- Activate: `.venv/Scripts/activate` (Windows) or `source .venv/bin/activate` (Unix)
- Windows-primary (PyCharm, cp1254 console encoding — avoid Unicode symbols like °, Δ in terminal output)

## Build & Run

```bash
pip install -e ".[dev]"       # install in dev mode with test deps
primerviz --help              # show CLI usage
primerviz -t SEQ -f FWD -r REV              # run with Rich output
primerviz -t SEQ -f FWD --plain             # force ASCII fallback
pytest                        # run all tests
pytest tests/test_analysis.py::test_gc_mixed  # run a single test
```

## Architecture

Uses `src/` layout. Entry point: `primerviz = primerviz.cli:main` (Click command).

- **`cli.py`** — Click CLI. Parses `--template/-t`, `--forward/-f` (repeatable), `--reverse/-r` (repeatable), `--mismatches/-m`, `--plain`. Orchestrates I/O → alignment → analysis → rendering.
- **`models.py`** — Dataclasses: `Primer`, `PrimerProperties`, `BindingSite`, `AnalysisResult`. `Direction` enum (FORWARD/REVERSE).
- **`analysis.py`** — All thermodynamic calculations via `primer3-py` (`calc_tm`, `calc_hairpin`, `calc_homodimer`, `calc_heterodimer`, `calc_end_stability`). MW computed from nucleotide weights. `analyze_primer()` is the convenience wrapper.
- **`alignment.py`** — Sliding-window binding site finder. Reverse primers are reverse-complemented before matching against the sense strand. Configurable mismatch threshold.
- **`io.py`** — Auto-detects input: FASTA file (by extension), plain text file, or raw sequence string.
- **`visualize.py`** — Dual renderer: Rich (colored tables/panels) with ASCII fallback. `render()` dispatches based on `use_rich` flag, falls back to ASCII if Rich import fails.

## Key Dependencies

- `click` — CLI framework
- `primer3-py` — Tm and thermodynamic calculations (wraps primer3 C library)
- `rich` — Terminal formatting (optional at runtime, ASCII fallback exists)
