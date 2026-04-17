# primerVizCLI

A command-line tool for visualizing DNA primer binding sites on a template sequence, with full thermodynamic analysis and qPCR compatibility reporting.

## Features

- **Binding visualization** — shows where primers anneal on the template (both strands searched automatically)
- **Amplicon zoom** — when both a forward and reverse primer are found, the view is trimmed to the amplicon region with genomic coordinates
- **Thermodynamic properties** — Tm (nearest-neighbor via primer3), GC%, molecular weight, hairpin, self-dimer, 3' stability
- **Probe support** — visualize a hydrolysis probe (`-p`) alongside primers; checks overlap, Tm delta, 5' G, and amplicon containment
- **qPCR report card** — pass/warn/fail evaluation of all design criteria (`--qpcr` or automatic with `-p`)
- **Rich terminal output** with plain ASCII fallback (`--plain`)
- **Flexible input** — raw sequences, FASTA files, or plain text files; repeatable primer options

## Installation

```bash
git clone https://github.com/uguraka/primerVizCLI.git
cd primerVizCLI
pip install -e ".[dev]"
```

Requires Python 3.10+.

## Usage

```bash
primerviz --help
```

```
Options:
  -t, --template TEXT       Template sequence or FASTA file  [required]
  -f, --forward TEXT        Forward primer sequence or file (repeatable)
  -r, --reverse TEXT        Reverse primer sequence or file (repeatable)
  -p, --probe TEXT          Probe sequence or file (enables qPCR report)
  -m, --mismatches INTEGER  Max allowed mismatches  [default: 0]
  --qpcr                    Run qPCR compatibility report
  --plain                   Force plain ASCII output
```

### Basic primer visualization

```bash
primerviz \
  -t "AAAGTATGATTGGCGTCGGTACAATCGTGGATGTCTCCAATAGCAAAAGTAGTTTGCGAGGGTGAGTTGTCTCTACCATCCAGCAAACTTATCGGATTCATCGTCCTAAATCTTTACTCACGCTAATATCAGTGGGAGAA" \
  -f "TTCTCCCACTGATATTAGCGTGA" \
  -r "AAAGTATGATTGGCGTCGGTAC"
```

```
+------------------+
| input  (140 bp)  |
+------------------+
            10        20   ...
             |         |
5'- AAAGTATGATTGGCGTCGGTAC...GGGAGAA -3'
    >>>>>>>>>>>>>>>>>>>>>>  rev_1 (Tm: 58.8C)
                                   <<<<<<<<<<<<<<<<<<<<<<<  fwd_1 (Tm: 59.0C)
```

### With FASTA file input

```bash
primerviz -t template.fasta -f primer_fwd.fasta -r primer_rev.fasta
```

### qPCR report

```bash
primerviz -t template.fasta -f FWD -r REV --qpcr
```

### TaqMan probe + automatic qPCR report

```bash
primerviz -t template.fasta -f FWD -r REV -p PROBE
```

```
5'- ...amplicon... -3'
    >>>>>>>>>>>>>>>>>>>>>>>    fwd (Tm: 58.8C)
              ===========      probe [probe] (Tm: 70.2C)
                          <<<  rev (Tm: 59.0C)

qPCR Compatibility Report:
  [PASS]  Amplicon size         140 bp       ideal 80-200 bp
  [PASS]  Amplicon GC%          42.9%        ideal 40-60%
  [PASS]  fwd_1 Tm              59.0C        ideal 58-65C
  [PASS]  Tm delta (fwd-rev)    0.2C         ideal <=1C
  [FAIL]  Probe Tm above primers +1.2C       probe should be >=8C above primer Tm
  ...
```

### Mismatch tolerance

```bash
primerviz -t template.fasta -f FWD -r REV -m 2
```

Mismatched positions are shown as `x` in the binding line.

## qPCR Criteria

| Criterion | Pass | Warn | Fail |
|---|---|---|---|
| Amplicon size | 80–200 bp | 60–300 bp | outside warn |
| Amplicon GC% | 40–60% | 35–65% | outside warn |
| Primer Tm | 58–65°C | 55–68°C | outside warn |
| Tm delta (fwd/rev) | ≤1°C | ≤3°C | >3°C |
| GC% | 40–60% | 35–65% | outside warn |
| 3' GC clamp (last 3 bp) | 1–2 G/C | 0 or 3 G/C | — |
| Mono-nucleotide run | ≤3 bp | 4 bp | ≥5 bp |
| Self-dimer / hairpin Tm | ≤40°C | ≤45°C | >45°C |
| Primer-dimer Tm | ≤40°C | ≤45°C | >45°C |
| Probe Tm above primers | ≥+8°C | ≥+5°C | <+5°C |
| Probe 5' base | not G | — | G |
| Probe within amplicon | yes | — | no |
| Probe-primer overlap | none | — | any |

## Development

```bash
pip install -e ".[dev]"   # install with test dependencies
pytest                    # run all tests
pytest tests/test_analysis.py::test_gc_mixed  # run a single test
```

## Architecture

```
src/primerviz/
  cli.py        Click entry point; wires I/O → alignment → analysis → render
  models.py     Dataclasses: Primer, BindingSite, AnalysisResult, QpcrReport
  analysis.py   Thermodynamics via primer3-py; run_qpcr_checks()
  alignment.py  Sliding-window search on both strands; strand auto-detected
  io.py         FASTA / plain-text / raw-sequence input parser
  visualize.py  Rich renderer + ASCII fallback; amplicon zoom; probe lines
```

## Dependencies

| Package | Purpose |
|---|---|
| `click` | CLI argument parsing |
| `primer3-py` | Nearest-neighbor Tm and thermodynamic calculations |
| `rich` | Colored terminal output (optional; ASCII fallback if unavailable) |
