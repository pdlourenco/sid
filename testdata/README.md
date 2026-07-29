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

CI is the canonical generator: the MATLAB job in
[`../.github/workflows/tests.yml`](../.github/workflows/tests.yml), pinned to
**MATLAB R2025a** (`matlab-actions/setup-matlab@v2.7.0`), regenerates the vectors on
push to `main` and commits any change as `Update cross-language reference data`.
Those committed bytes are the reference.

**Never commit regen *churn*** — files whose values changed only by engine ULPs. The
stored *inputs* are RNG-stable and reproduce exactly, but computed *outputs* differ
by ~1 ULP across MATLAB versions (more for ill-conditioned paths — the LTV-IO solve
drifts ~1e-4, still within its 1% `A_rel`/`B_rel`/`Cost_rel` tolerance), and because
`jsonencode` emits shortest-round-trip decimals a 1-ULP change rewrites the whole
file. So a local regen under any engine other than the pinned R2025a churns most
`reference_*.json` for no real numerical change — leave those untouched.

A PR that *changes* reference semantics or *adds* a vector is different: commit
exactly the files the change affects (flag them in the PR, and note the canonical
R2025a round-trip CI will apply post-merge), and leave the rest untouched. This is
[ADR-0002](../docs/decisions/ADR-0002-contract-artifact-hardening.md) rules 3–4 —
the vector and its consumer land together, and payloads change only by regeneration
— see #147, #174, and #175 for the pattern in practice.

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

## Serialization conventions (avoid these round-trip traps)

`jsonencode`/`jsondecode` and `numpy` do not round-trip every array shape the
way you might expect. When adding a vector:

- **Store matrices as genuine 2-D arrays, never a `1×n` row.** MATLAB
  `jsondecode` collapses a stored row vector to a **column**, silently
  transposing a `(p_y × n)` observation matrix — the Octave validator then
  errors while the Python consumer's `H[newaxis, :]` masks it. Use a `2×n`
  (or larger) `H`, or reshape defensively in *both* consumers.
- **Store multi-trajectory input as true 3-D `(N × ny × L)`.** A squeezed
  `(N × L)` matrix is ambiguous — the estimators read it as MIMO (`ny = L`),
  not multi-trajectory. MATLAB `(d1 × d2 × d3)` maps to numpy `(d1, d2, d3)`
  directly, so 3-D storage round-trips cleanly and dispatches by `ndims == 3`.
- **Variable-length trajectories** are stored as cells of unequal-size arrays;
  `jsondecode` returns a cell and `json.load` a list, so build a list of
  per-trajectory arrays in the Python consumer.
- **Floor every field that can hold a zero or near-zero element with an
  `<field>_atol`** (ADR-0002 rule 1). `atol = 0` makes a near-zero element's
  pass condition `|diff| <= rtol·|expected|`, which any cross-engine ULP
  difference fails — latent until the canonical R2025a regen meets the Octave
  validator. Tolerance blocks may be hand-edited (generator **and** committed
  block) in the same PR that adds the vector, per #176's amended policy.
