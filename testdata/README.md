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
- `output` — expected outputs at full double precision (shortest round-trip)
- `tolerance` — per-field relative tolerances for cross-language comparison

### Canonical environment

`reference_*.json` is regenerated and committed **only by CI** — the MATLAB job in
[`../.github/workflows/tests.yml`](../.github/workflows/tests.yml), pinned to
**MATLAB R2025a** (`matlab-actions/setup-matlab@v2.7.0`), which commits any change
as `Update cross-language reference data`. That job is the single source of truth.

Regenerating locally is fine for inspection, but **do not commit the result.** The
stored *inputs* are RNG-stable and reproduce exactly, but the computed *outputs*
differ by ~1 ULP across MATLAB versions (and more for ill-conditioned paths — the
LTV-IO solve drifts ~1e-4, still within its 1% `A_rel`/`B_rel`/`Cost_rel`
tolerance). Because `jsonencode` emits shortest-round-trip decimals, a 1-ULP output
change rewrites the entire file, so a local regen under any engine other than the
pinned R2025a churns most `reference_*.json` for no real numerical change. Let CI
regenerate.

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
