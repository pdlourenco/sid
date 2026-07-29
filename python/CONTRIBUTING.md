# Contributing to sid — Python

This guide covers Python-specific contribution standards for the sid toolbox.
For general project guidelines, see the root [CONTRIBUTING.md](../CONTRIBUTING.md).

> **⚠ Read this first.** Before writing or modifying any algorithmic code,
> read the [Specification as Source of Truth](../CONTRIBUTING.md#specification-as-source-of-truth)
> section in the root contributing guide. `spec/SPEC.md` is the binding
> contract for this implementation — Python conforms to the spec, not to
> the MATLAB port. If MATLAB and the spec disagree, the spec wins.

Please ensure that `pytest python/tests/ -v` passes before submitting — the
CI pipeline checks automatically.

---

## Docstring Standard

Every Python module (public and private) **must** follow the NumPy-style
docstring template below. This ensures consistency with the MATLAB function
header format, enables IDE documentation popups, and keeps a clear link
between code and the algorithm specification.

### Canonical Template

```python
# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Brief one-line description of the module."""

from __future__ import annotations

import numpy as np

from sid._results import FreqResult


def function_name(
    in1: np.ndarray,
    in2: np.ndarray | None = None,
    *,
    option_name: int | None = None,
    sample_time: float = 1.0,
) -> FreqResult:
    """Brief one-line description.

    Extended description paragraph(s). What the function does, context,
    and any important notes.

    Parameters
    ----------
    in1 : ndarray, shape (N, ny) or (N, ny, L)
        Description. Constraints.
    in2 : ndarray or None, optional
        Description. Use ``None`` for time-series mode. Default is ``None``.
    option_name : int, optional
        Description. Default is ``min(N // 10, 30)``.
    sample_time : float, optional
        Sample time in seconds. Default is ``1.0``.

    Returns
    -------
    FreqResult
        Frozen dataclass with fields:

        - **frequency** (*ndarray, shape (nf,)*) -- Frequency vector, rad/sample.
        - **response** (*ndarray or None*) -- Complex frequency response.
        - ...

    Raises
    ------
    SidError
        If data contains NaN or Inf values (code: ``'non_finite'``).
    ValueError
        If ``in1`` and ``in2`` have incompatible shapes.

    Examples
    --------
    Basic SISO usage:

    >>> result = function_name(y, u)  # doctest: +SKIP
    >>> result.response.shape
    (128,)

    With options:

    >>> result = function_name(y, u, option_name=50)  # doctest: +SKIP

    Notes
    -----
    **Algorithm:**

    1. Step description.
    2. Step description.

    **Specification:** SPEC.md §X.Y -- Section Title

    References
    ----------
    .. [1] Author, "Title", Publisher, Year. Sections X.Y.

    See Also
    --------
    related_function_1 : Brief description.
    related_function_2 : Brief description.

    Changelog
    ---------
    YYYY-MM-DD : First version by Author Name.
    """
    ...
```

### Section Order (fixed)

Sections must appear in this exact order within the docstring:

1. **One-line summary** — first line of docstring
2. **Extended description** — one or more paragraphs
3. **`Parameters`** — all positional and keyword arguments
4. **`Returns`** — return value(s)
5. **`Raises`** — exceptions that may be raised
6. **`Examples`** — runnable code snippets
7. **`Notes`** — algorithm description and SPEC.md cross-reference
8. **`References`** — academic citations
9. **`See Also`** — related functions
10. **`Changelog`** — entries in `YYYY-MM-DD : Description by Author.` format

### Required vs Optional Sections

| Section | Required? |
|---------|-----------|
| One-line summary | Always |
| Extended description | Always (can be brief for simple helpers) |
| Parameters | Always (unless the function takes no arguments) |
| Returns | Always (unless the function returns nothing) |
| Raises | Only if the function raises exceptions |
| Examples | Always |
| Notes | Only for non-trivial algorithms or SPEC.md cross-references |
| References | Only when citing academic papers |
| See Also | Always |
| Changelog | Always |

### Key Rules

- **Section headings are always plural**: `Parameters`, `Returns`, `Examples`,
  `References`. Exception: `See Also` and `Changelog` keep their
  traditional casing.
- **SPEC.md cross-reference**: if the function implements or relates to a
  section of `SPEC.md`, add a line in the `Notes` section:
  `**Specification:** SPEC.md §2 -- Blackman-Tukey Spectral Analysis`
- **Copyright block** is a comment block at the **top of every file**, before
  the module docstring. It uses `#` comment lines.
- **Changelog dates** use ISO 8601 format (`YYYY-MM-DD`).
- **`See Also`** appears exactly once per docstring.
- **Keyword-only arguments**: all optional parameters use keyword-only syntax
  (after `*` in the signature). Positional arguments come first.
- **Type hints**: all function signatures include type annotations.

### Mapping from MATLAB Header Sections

| MATLAB Section | Python Section | Notes |
|----------------|---------------|-------|
| `% FUNCTIONNAME Brief description.` | One-line summary | |
| Usage signatures | *(implicit from signature + type hints)* | Not a separate section |
| Extended description | Extended description | |
| `INPUTS:` | `Parameters` | Positional arguments |
| `NAME-VALUE OPTIONS:` | `Parameters` | Keyword-only arguments (after `*`) |
| `OUTPUTS:` | `Returns` | |
| `EXAMPLES:` | `Examples` | Use `>>> ` doctest format |
| `ALGORITHM:` | `Notes` | Under `**Algorithm:**` sub-heading |
| `SPECIFICATION:` | `Notes` | Under `**Specification:**` line |
| `REFERENCES:` | `References` | Use `.. [1]` format |
| `See also:` | `See Also` | |
| `Changelog:` | `Changelog` | |
| Copyright block | File-level `#` comment | Before module docstring |

---

## Naming Conventions

Python functions follow the same `sid` + `Domain` + `Method` pattern as MATLAB,
translated to snake_case:

```
sid.freq_bt          # sidFreqBT
sid.ltv_disc         # sidLTVdisc
sid.bode_plot        # sidBodePlot
sid.model_order      # sidModelOrder
```

### Module and function naming

| MATLAB | Python module | Python function |
|--------|--------------|-----------------|
| `sidFreqBT.m` | `sid/freq_bt.py` | `freq_bt()` |
| `sidLTVdisc.m` | `sid/ltv_disc.py` | `ltv_disc()` |
| `sidBodePlot.m` | `sid/bode_plot.py` | `bode_plot()` |

### Private helpers

Internal helper functions live in `sid/_internal/` and use the same snake_case
convention (e.g., `validate_data`, `hann_win`, `ltv_cosmic_solve`). The leading
underscore on `_internal` makes these private by Python convention.

Many existing helpers happen to correspond one-to-one with a file in
`matlab/sid/private/`, because that is how the port was originally built — but
**that correspondence is not a rule**. Helper decomposition is a private
implementation detail each port chooses for itself; only the *public* catalogue
and the spec's behaviour are shared (`docs/roadmap.md`,
[ADR-0001](../docs/decisions/ADR-0001-spec-is-the-contract.md)). Do not add a
Python helper merely to mirror a MATLAB one, and do not treat the existing
pairing as a constraint on new work.

### Result field naming

All result struct fields use snake_case:

| MATLAB field | Python field |
|-------------|-------------|
| `Response` | `response` |
| `FrequencyHz` | `frequency_hz` |
| `NoiseSpectrum` | `noise_spectrum` |
| `NumTrajectories` | `num_trajectories` |
| `SampleTime` | `sample_time` |

### Reserved word handling

MATLAB's `Lambda` parameter becomes `lambda_` in Python (PEP 8 trailing
underscore convention for reserved words).

---

## Code Style

| Rule | Value |
|------|-------|
| Indentation | 4 spaces (no tabs) |
| Line length | 100 characters max |
| Line endings | LF |
| Charset | UTF-8 |
| Formatter | ruff format |
| Linter | ruff check |
| Type hints | Required on all public function signatures |
| Imports | No star imports; group: stdlib, third-party, local |

### Translating spec notation to NumPy

[`spec/SPEC.md`](../spec/SPEC.md) states formulas in mathematical notation with
1-based indices. These are the recurring translation gotchas:

| Concern | Rule |
|---|---|
| Indexing | Spec formulas are **1-indexed**; NumPy is 0-indexed. Adjust bin/step extraction rather than transcribing indices literally. |
| FFT convention | `np.fft.fft` matches the spec's unscaled forward transform (no `1/N` on the forward direction). |
| Linear solves | A spec `A⁻¹b` (MATLAB `A\b`) is `np.linalg.solve(A, b)` — never `inv(A) @ b`. |
| Data layout | `(N, ny)` for signals and `(N, ny, L)` for multi-trajectory, matching the spec's stated shapes. |

Reserved-word escaping (`Lambda` -> `lambda_`) is covered under
[Reserved word handling](#reserved-word-handling); RNG incompatibility is
covered under [Test categories](#test-categories).

### Plotting conventions

Plotting functions must not make matplotlib a hard dependency of the package:

- **Lazy-import matplotlib** inside the function body, not at module scope —
  `import sid` must work without matplotlib installed.
- **Accept an optional `ax=` keyword** so a caller can embed the plot in a
  figure they created.
- **Return a `dict` of handles** rather than bare figure/axes objects.

### Inline Comments

Code comments within function bodies should make the mathematical intent
clear and link back to the specification. Follow these guidelines:

**Section separators.** Use `# ---- Name ----` to mark major computational
phases. These should correspond to steps in the function's `Notes` section:

```python
# ---- Build data matrices (SPEC.md §8.3.2) ----
D, Xp = ltv_build_data_matrices(X, U)
```

**SPEC.md cross-references.** When a code block implements a specific
equation or algorithm step from SPEC.md, cite the section number:

```python
# Schur complement forward pass (SPEC.md §8.3.4, Eq. 8.3):
#   Lambda(k) = S(k) - lambda(k-1)^2 * Lambda(k-1)^{-1}
Lbd[:, :, k] = S[:, :, k] - lam[k - 1] ** 2 * np.linalg.solve(Lbd[:, :, k - 1], I)
```

**Mathematical steps.** Annotate non-obvious operations — matrix
inversions, Schur complements, spectral transformations, and
regularization terms. Write the formula in comment notation before
the code that implements it:

```python
# G(w) = Phi_yu(w) / Phi_u(w)  -- transfer function estimate
G = Phi_yu / Phi_u

# Phi_v(w) = Phi_y(w) - |Phi_yu(w)|^2 / Phi_u(w)  -- noise spectrum
Phi_v = Phi_y - np.abs(Phi_yu) ** 2 / Phi_u
```

**Variable-to-notation mapping.** When a variable name differs from the
mathematical notation in SPEC.md, state the correspondence on first use:

```python
# Lbd corresponds to Lambda_k in SPEC.md §8.3 (forward Schur complement)
Lbd = np.zeros((d, d, N))
```

**Dimensions.** Annotate array dimensions on the line that creates or
returns them, using trailing comments:

```python
Ryy = sid_cov(y, y, M)  # (M+1, ny, ny) biased auto-covariance
```

**What not to comment.** Do not comment self-explanatory operations
(loop counters, standard imports, trivial assignments). Focus
comments on *why*, not *what*:

```python
# Bad: increment k by 1
k = k + 1

# Good: skip the first segment (it has incomplete overlap)
k = k + 1
```

---

## Testing

### Running tests

```bash
# All tests
pytest python/tests/ -v

# Single test file
pytest python/tests/test_freq_bt.py -v
```

### Test structure

Tests use **pytest** and follow auto-discovery conventions. To add a new
test, create a file matching the `test_*.py` pattern in `python/tests/`.

Each test file tests one public function or one private helper. Test
functions use descriptive names:

```python
def test_hann_win_symmetry():
    """Hann window is symmetric: w(tau) == w(-tau)."""
    ...

def test_freq_bt_siso_known_system():
    """SISO BT on AR(1) recovers known frequency response."""
    ...
```

### Test categories

Every public function has **two kinds of tests**:

1. **Unit tests** (in `test_<function>.py`) — assert the rules of the function's
   [`spec/SPEC.md`](../spec/SPEC.md) section, so a failure reads "this violates
   the spec" rather than "this differs from MATLAB". Reusing the *scenarios* from
   the corresponding MATLAB test file (`matlab/tests/test_sid<Function>.m`) is
   encouraged — the same plants and edge cases are worth exercising in both ports
   — but derive the **assertions and tolerances from the spec**, not from what
   MATLAB happens to produce. Random data will differ anyway (MATLAB and NumPy
   RNGs are incompatible).

2. **Equivalence tests** (in `test_cross_validation.py`) — load JSON
   reference data from `testdata/` and verify that the Python function produces
   the same output to the vector's stored tolerance. These use stored input data
   (not seeds), so they are RNG-independent. They check that the ports *agree*;
   they do **not** prove either satisfies the spec — two ports can drift
   identically and still match ([ADR-0001](../docs/decisions/ADR-0001-spec-is-the-contract.md)).
   That is why category 1 is not optional.

### Shared fixtures

`python/tests/conftest.py` provides shared fixtures:

```python
@pytest.fixture
def rng():
    """Reproducible random number generator."""
    return np.random.default_rng(42)

def load_reference(name: str) -> dict:
    """Load a JSON reference file from testdata/."""
    ...
```

---

## Examples (Jupyter Notebooks)

> **⚠ Read this first.** The binding contract for every example is
> [`spec/EXAMPLES.md`](../spec/EXAMPLES.md). It defines the physical plant
> catalog, the `util_msd*` helper API, and — for every example — the required
> pedagogical sections, the `sid.*` call graph, and the required plots and
> prints. Python notebooks and MATLAB scripts are parallel implementations of
> the same spec. If you are porting, adding, or modifying an example, start
> there; the conventions in this file cover only the Python-specific notebook
> mechanics.

Examples live in `python/examples/` as **Jupyter notebooks** (`.ipynb`). Each
MATLAB example script (`matlab/examples/example*.m`) maps to one notebook.
Notebooks are the Python equivalent of MATLAB's `%%`-sectioned scripts — they
combine narrative, code, and inline plots in a single runnable document.

### Naming convention

| MATLAB example | Python notebook |
|---|---|
| `exampleSISO.m` | `example_siso.ipynb` |
| `exampleFreqMap.m` | `example_freq_map.ipynb` |
| `exampleLTVdisc.m` | `example_ltv_disc.ipynb` |
| `exampleOutputCOSMIC.m` | `example_output_cosmic.ipynb` |

Pattern: drop the `example` prefix camelCase, convert to `example_snake_case.ipynb`.

### Auto-discovery

The example runner and CI discover notebooks by globbing `example_*.ipynb`
in `python/examples/`. **Do not add entries to a hardcoded list** — it defeats
auto-discovery and creates merge conflicts. This mirrors the MATLAB convention
where `runAllExamples.m` uses `dir('example*.m')`.

### Notebook structure

Each notebook should:

1. **Title cell** (Markdown) — `# Example: <title>` and a one-paragraph
   description of the topic the example's [`spec/EXAMPLES.md`](../spec/EXAMPLES.md)
   §3 entry assigns it. Exact prose is yours (§0.3: header prose is SHOULD).

2. **Setup cell** — imports and `%matplotlib inline`:

   ```python
   import numpy as np
   import sid

   %matplotlib inline
   ```

3. **One code cell per spec section** — the ordered list of pedagogical sections
   and each section's topic are **MUST**-match items in
   [`spec/EXAMPLES.md`](../spec/EXAMPLES.md) §0.3, so follow the example's §3
   entry section by section. Use Markdown cells between code cells for narrative.

4. **Self-contained** — every notebook must run top-to-bottom without external
   data files. All data is generated inline, using the plant parameters and
   topology fixed in `spec/EXAMPLES.md` §1 (also MUST-match items).

5. **No output committed** — clear all cell outputs before committing. CI
   validates that notebooks execute without error, but stored outputs bloat
   the repository and cause noisy diffs.

### Checklist (per notebook)

Examples are a contract too: [`spec/EXAMPLES.md`](../spec/EXAMPLES.md) binds the
plant, the ordered sections, the `sid.*` call graph and its options, the plot
kinds, and the print semantics — so a notebook is written **from its §3 entry**,
not by translating `matlab/examples/*.m`. Per §0.3 the prose, styling, RNG seed,
and cosmetics are deliberately yours.

1. Read the example's `spec/EXAMPLES.md` §3 entry, plus §1 for the plant
   parameters and §2 for the helper definitions it uses.
2. Create the notebook with the spec's ordered sections and section topics.
3. Write each section's code so the `sid.*` calls and their options match the
   spec's call graph; use `sid.*_plot()` or matplotlib for the specified plot
   kinds, and print the specified quantities.
4. Add Markdown narrative in your own words for each section.
5. Run the notebook top-to-bottom to verify it executes cleanly.
6. Clear outputs, then commit.

The MATLAB example is a useful cross-check for "did I read the spec entry the
same way?" — but where the two disagree, the spec entry governs and the MATLAB
example is the one that needs fixing.

### CI validation

Notebooks are validated in CI using `pytest --nbmake python/examples/` to
ensure they execute without errors. Only execution is checked — output cell
content is not asserted.

### `python/examples/README.md`

Maintain an `examples/README.md` with an index table and per-notebook
descriptions, mirroring `matlab/examples/README.md`. The table should list
notebook name, description, and which `sid` functions are demonstrated.

---

## Implementation Workflow

**Implement from [`spec/SPEC.md`](../spec/SPEC.md), not from the MATLAB source.**
The spec is the sole contract; the MATLAB port is another consumer of it, not a
reference implementation ([ADR-0001](../docs/decisions/ADR-0001-spec-is-the-contract.md),
`CLAUDE.md` §3). Transcribing MATLAB logic into Python is the **copy-by-copy
porting** red flag in [`docs/REVIEW_CONTEXT.md`](../docs/REVIEW_CONTEXT.md): it
propagates MATLAB's bugs into Python, and cross-validation will *not* catch it —
both ports drift from the spec identically and the vectors still pass.

When adding a function to the Python port, follow this order:

1. **Read the spec section.** Find the normative rules for the function:
   formulas, defaults, input handling, edge cases, error conditions, output
   fields, and numerical diagnostics. All of it is contract, not just the
   formulas. If the spec is ambiguous or silent on something you need, **stop and
   ask** — do not infer the answer from MATLAB's behaviour and do not invent one
   (`CLAUDE.md` §3). Resolving it means a spec change first.

2. **Implement the private dependencies, then the public function**, each from
   the spec rules of the behaviour it carries, with a docstring following the
   template above. Helper decomposition is yours to choose — it need not mirror
   `matlab/sid/private/`. **Each new private helper gets its own test file**
   (`test_<helper>.py`), so a building block is tested in isolation before it is
   composed.

3. **Write unit tests from the spec.** Derive the acceptance criteria from the
   spec rules, so a test failure means "this violates the spec", not "this
   differs from MATLAB". Sharing *scenarios* with `matlab/tests/test_sid*.m` is
   fine and encouraged — the same plants and the same edge cases are worth
   exercising in both ports — but the assertions come from the spec.

4. **Add a cross-validation case** in `test_cross_validation.py` against the
   stored `testdata/` vectors. This checks agreement; it does **not** prove
   either port satisfies the spec (`CLAUDE.md` §3), so it complements step 3
   rather than replacing it. Vectors are produced **only** by regenerating
   `testdata/generate_reference.m` under MATLAB — never hand-write or hand-edit a
   payload, and never leave a vector no consumer reads
   ([ADR-0002](../docs/decisions/ADR-0002-contract-artifact-hardening.md)). If a
   new function needs a vector, add its block to the generator.

5. **Run the full suite.** `pytest python/tests/ -v` must be green.

6. **Update `__init__.py`.** Export the new function from `sid/__init__.py`.

**Using the MATLAB source.** Reading it is legitimate — as a cross-check once
you have your own spec-derived implementation, or to understand *why* a spec rule
exists. What is not legitimate is treating it as the source of behaviour. If your
implementation and MATLAB's disagree, that is a finding, not something to
paper over: re-read the spec, and if the spec is on your side, **MATLAB is wrong**
— file an issue against it rather than reproducing the discrepancy. If you
genuinely cannot tell which is right, **stop and ask** — do not pick a side
unilaterally (`CLAUDE.md` §3).

---

## Dependencies

| Package | Minimum | Role |
|---------|---------|------|
| numpy | >= 1.22 | Array operations, FFT |
| scipy | >= 1.8 | Linear algebra, signal processing |
| matplotlib | >= 3.5 | Plotting (optional; required for examples) |
| pytest | >= 7.0 | Testing (dev only) |
| nbmake | >= 1.0 | Notebook execution testing (dev only) |
| ruff | latest | Linting and formatting (dev only) |
