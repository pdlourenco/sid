# Cross-Language Reference Data

This directory contains canonical test vectors used to verify numerical
equivalence across the MATLAB, Python, and Julia implementations of sid.

## Generating reference data

Run from the repository root in MATLAB:

```matlab
run('testdata/generate_reference.m')
```

This produces JSON files with full double-precision outputs for a curated
set of test cases. Each JSON file contains:

- `function` — the sid function name
- `params` — algorithm parameters used
- `input` — input data arrays (not seeds, to avoid RNG differences)
- `output` — expected outputs at 16 significant digits
- `tolerance` — per-field relative tolerances for cross-language comparison

## Validating against reference data

Each language validates its own outputs against the committed JSON files:

```matlab
% MATLAB / Octave
run('testdata/validate_reference.m')
```

The validator reads each `reference_*.json`, calls the corresponding sid
function with the stored input data, and checks that outputs match within
the specified tolerances.  This runs automatically in CI via
`cross-validate.yml` (currently Octave; Python and Julia will be added
in future phases).

## Format

All numeric arrays are stored as JSON arrays of numbers with full double
precision.  Complex arrays are split into `_real` and `_imag` components.
