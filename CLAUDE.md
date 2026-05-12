# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project: primerVizCLI

A command-line tool for visualizing DNA primer binding sites on template sequences with full thermodynamic analysis and qPCR compatibility reporting.

**Key features:**
- Primer binding visualization (both strands searched automatically)
- Amplicon zoom when both forward and reverse primers are detected
- Thermodynamic properties (Tm, GC%, MW, hairpin, self-dimer, 3' stability) via primer3-py
- Probe (hydrolysis) support with overlap and Tm delta checks
- qPCR report card (15+ design criteria evaluated as pass/warn/fail)
- Rich terminal output with ASCII fallback
- Flexible input: raw sequences, FASTA files, or plain text files

## Quick Start

### Installation & Setup
```bash
pip install -e ".[dev]"   # Install with test dependencies
```

### Running the CLI
```bash
primerviz -t TEMPLATE -f FORWARD -r REVERSE [OPTIONS]
primerviz -t template.fasta -f fwd_primer.fasta -r rev_primer.fasta
primerviz -t template.fasta -f FWD -r REV -p PROBE  # With probe (auto-enables qPCR)
primerviz -t template.fasta -f FWD -r REV --qpcr     # Explicit qPCR report
primerviz -t template.fasta -f FWD -r REV --plain    # ASCII output (no Rich)
primerviz -t template.fasta -f FWD -r REV -m 2       # Allow up to 2 mismatches
```

Requires Python 3.10+.

### Testing
```bash
pytest                                    # Run all tests
pytest tests/test_analysis.py             # Run a single test module
pytest tests/test_analysis.py::test_gc_mixed  # Run a single test
pytest -v                                 # Verbose output
```

## Architecture

The tool follows a **linear processing pipeline**:

```
[CLI Input] → [I/O Parsing] → [Alignment] → [Analysis] → [Rendering]
  cli.py        io.py            alignment.py   analysis.py   visualize.py
```

### Module Responsibilities

**cli.py** — Click command entry point
- Parses 7 CLI options (template, forward/reverse primers, probe, mismatches, qpcr flag, plain flag)
- Validates at least one primer is provided
- Orchestrates the entire pipeline: I/O → alignment → analysis → rendering
- Handles multiple templates (if FASTA with multiple records) and multiple primers (-f/-r repeatable)

**io.py** — Sequence input parsing
- Detects if input is a file path or raw sequence string
- Parses FASTA files (recognizes .fa, .fasta, .fna, .fas extensions)
- Reads plain text files (one sequence per line, no headers)
- Returns list of (name, sequence) tuples

**alignment.py** — Primer binding site discovery
- Core function: `find_binding_sites(primer, template, max_mismatches=0)`
- Slides primer along both strands of template (length of template - primer length + 1 positions)
- Matches both primer sequence AND its reverse complement
- Determines actual binding strand based on which form matches (not user's -f/-r label)
  - `strand=FORWARD`: primer sequence matches sense strand (binds antisense, displays as `>>>`)
  - `strand=REVERSE`: primer's RC matches sense strand (binds sense, displays as `<<<`)
- Counts mismatches and records mismatch positions for visualization
- Returns list of `BindingSite` objects with position, strand, and mismatch info

**analysis.py** — Thermodynamic calculations & qPCR evaluation
- **Thermodynamics** (via primer3-py):
  - `calc_tm()` — Nearest-neighbor Tm (SantaLucia parameters)
  - `calc_gc()` — GC content percentage
  - `calc_mw()` — Molecular weight accounting for nucleotide masses and phosphodiester bonds
  - `calc_hairpin()` — Hairpin Tm (self-complementarity)
  - `calc_self_dimer()` — Homodimer Tm
  - `calc_primer_dimer(seq1, seq2)` — Heterodimer Tm between two primers
  - `calc_3prime_stability()` — Delta-G (kcal/mol) of 3' end
- **qPCR checks** via `run_qpcr_checks(AnalysisResult)`:
  - Amplicon size (80-200 bp ideal) and GC% (40-60% ideal)
  - Per-primer Tm (58-65°C ideal), GC%, 3' GC clamp, mono-nucleotide runs
  - Self-dimer and hairpin Tm (≤40°C ideal)
  - Tm delta between forward and reverse (≤1°C ideal)
  - Primer-dimer pairs (≤40°C ideal)
  - Probe Tm above primers (≥+8°C ideal), 5' base (not G), amplicon containment, overlap with primers
- Returns `QpcrReport` with list of `QpcrCheck` objects (name, value, status, detail)

**models.py** — Data structures
- `Primer`: name, sequence, direction (FORWARD/REVERSE)
- `BindingSite`: primer, start/end positions, mismatches, strand, mismatch_positions
- `PrimerProperties`: tm, gc_percent, length, molecular_weight, hairpin_tm, self_dimer_tm, three_prime_stability
- `QpcrCheck`: name, value, status (pass/warn/fail), detail
- `QpcrReport`: list of checks, amplicon_gc, amplicon_len
- `AnalysisResult`: template_name, template_seq, bindings, primer_props, primer_dimer_tms, probe_bindings, probe_props, qpcr_report

**visualize.py** — Terminal output rendering
- Public function: `render(result, use_rich=True)`
  - Attempts Rich renderer; falls back to ASCII if Rich import fails
- **Rich renderer** (`_render_rich()`):
  - Panel header with template name and amplicon info
  - Ruler showing genomic coordinates (multiples of 10)
  - Template sequence with 5'/3' markers
  - Primer lines with position, Tm, mismatch count (using `>` for FORWARD, `<` for REVERSE, `x` for mismatches)
  - Probe lines (using `=` for matches, `x` for mismatches)
  - Properties table (Tm, GC%, MW, hairpin Tm, self-dimer Tm, 3' dG)
  - Primer-dimer Tm table (if multiple primers)
  - qPCR report table (if enabled) with pass/warn/fail icons and thresholds
- **ASCII renderer** (`_render_ascii()`): Text-only tables and formatting, no colors
- **Amplicon zoom** (`_get_amplicon()`):
  - When both FORWARD and REVERSE binding sites exist, view is trimmed to the amplicon region
  - Amplicon start = min(forward binding starts), Amplicon end = max(reverse binding ends)
  - Falls back to full template if only one primer type found
  - Coordinate offsets are tracked to display correct genomic positions in ruler

### Data Flow Example

For input: `primerviz -t template.fasta -f fwd.fasta -r rev.fasta -p probe.fasta --qpcr`

1. **io.py**: Parse all 4 FASTA files → list of (name, sequence) tuples
2. **cli.py**: Construct `Primer` objects for 2 primers + 1 probe
3. **alignment.py**: Search template for fwd/rev/probe binding sites (3 calls to `find_binding_sites()`)
4. **analysis.py**: 
   - Call `analyze_primer()` for each primer/probe → `PrimerProperties`
   - Call `calc_primer_dimer()` for each primer pair (e.g., fwd-rev)
   - Call `run_qpcr_checks()` → `QpcrReport` with 15+ checks
5. **models.py**: Assemble `AnalysisResult` with all results
6. **visualize.py**: Render AnalysisResult with Rich formatter (or ASCII fallback)

## Key Design Patterns

**Strand Detection**: The `-f` (forward) and `-r` (reverse) labels are user metadata only. Actual binding orientation is discovered:
- If primer sequence matches the sense strand at position i, it binds the antisense strand (FORWARD, `>>>`)
- If primer's reverse complement matches the sense strand at position i, it binds the sense strand (REVERSE, `<<<`)
- Both forms are always searched; a "mislabeled" forward primer that only matches as RC will still be found as REVERSE

**Amplicon**: A valid amplicon requires both forward and reverse binding sites. The amplicon spans from the leftmost forward binding to the rightmost reverse binding. When rendering, the view is automatically zoomed to the amplicon region.

**Mismatch Handling**: Mismatches are tallied during alignment and positions are stored. In visualization, mismatched bases are marked with `x` instead of `>/<` or `=`.

**qPCR Criteria**: Each check has three status levels (pass/warn/fail) with inclusive ranges. For example:
- Amplicon size: PASS 80-200 bp, WARN 60-300 bp, FAIL outside WARN range
- Primer Tm: PASS 58-65°C, WARN 55-68°C, FAIL outside WARN range

## Testing

**Test organization**:
- `conftest.py`: Pytest fixtures (fwd_primer, rev_primer, template_seq)
- `test_alignment.py`: Reverse complement, exact matches, mismatch tolerance, strand detection
- `test_analysis.py`: GC%, Tm, MW, hairpin, dimer, 3' stability calculations
- `test_io.py`: FASTA parsing (single/multiple records), raw sequence input, file detection
- `test_visualize.py`: Smoke tests for Rich and ASCII renderers (no crash)

**Running tests**:
```bash
pytest                          # All tests
pytest -v                       # Verbose
pytest tests/test_alignment.py  # Single module
pytest -k test_gc_mixed         # By keyword
```

## Dependencies

- `click>=8.0` — CLI argument parsing and command structure
- `primer3-py>=2.0` — Thermodynamic calculations (Tm, hairpin, dimer via primer3 C library)
- `rich>=13.0` — Colored terminal output (optional; ASCII fallback if unavailable)
- `pytest>=7.0` — Test framework (dev only)

## Notes

- **Primer labeling**: Always remember that `-f` and `-r` are hints, not truth. The actual binding strand is determined by sequence matching during alignment.
- **Coordinate systems**: All positions in `BindingSite` are 0-based on the template. Visualization adds 1 for 1-based genomic coordinate display.
- **Rich optional**: The tool gracefully falls back to ASCII if Rich is not available or fails to import. Both renderers should produce equivalent information.
- **qPCR thresholds**: The pass/warn/fail ranges are intentionally inclusive to allow users to understand margin of safety. A check can only have one status (not multiple).
