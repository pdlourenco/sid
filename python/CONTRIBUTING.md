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

| MATLAB private function | Python module |
|------------------------|---------------|
| `sidCov.m` | `sid/_internal/cov.py` |
| `sidValidateData.m` | `sid/_internal/validate_data.py` |
| `sidLTVcosmicSolve.m` | `sid/_internal/ltv_cosmic_solve.py` |

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

Every ported function has **two kinds of tests**:

1. **Unit tests** (in `test_<function>.py`) — port the logic from the
   corresponding MATLAB test file (`matlab/tests/test_sid<Function>.m`).
   Use the same test scenarios, tolerances, and assertions. The random
   data will differ (MATLAB and NumPy RNGs are incompatible), but the
   test *structure* and *acceptance criteria* should match.

2. **Equivalence tests** (in `test_cross_validation.py`) — load JSON
   reference data from `testdata/` (generated by MATLAB) and verify
   that the Python function produces the same output to `rtol=1e-10`.
   These use stored input data (not seeds), so they are RNG-independent.

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
   description matching the MATLAB example's `%%` header comment.

2. **Setup cell** — imports and `%matplotlib inline`:

   ```python
   import numpy as np
   import sid

   %matplotlib inline
   ```

3. **One code cell per MATLAB `%%` section** — mirror the MATLAB structure
   section by section. Use Markdown cells between code cells for narrative
   (replacing MATLAB `%%` section comments and `%` inline explanations).

4. **Self-contained** — every notebook must run top-to-bottom without external
   data files. All data is generated inline (matching the MATLAB examples).

5. **No output committed** — clear all cell outputs before committing. CI
   validates that notebooks execute without error, but stored outputs bloat
   the repository and cause noisy diffs.

### Porting checklist (per notebook)

When porting a MATLAB example to a notebook:

1. Read the MATLAB example from top to bottom.
2. Create the notebook with matching sections.
3. Translate MATLAB code to Python, using the ported `sid.*` functions.
4. Replace MATLAB plotting idioms (`figure; hold on; semilogx; ...`) with
   matplotlib equivalents or `sid.*_plot()` functions.
5. Add Markdown narrative explaining what each section does.
6. Run the full notebook to verify it executes cleanly.
7. Clear outputs, then commit.

### CI validation

Notebooks are validated in CI using `pytest --nbmake python/examples/` to
ensure they execute without errors. Only execution is checked — output cell
content is not asserted.

### `python/examples/README.md`

Maintain an `examples/README.md` with an index table and per-notebook
descriptions, mirroring `matlab/examples/README.md`. The table should list
notebook name, description, and which `sid` functions are demonstrated.

---

## Porting Workflow

When porting a MATLAB function to Python, follow this order:

1. **Port private dependencies first.** For each private helper used by the
   target function, read the MATLAB source line by line, write the Python
   equivalent, and write its test file.

2. **Port the public function.** Read the MATLAB source, write the Python
   function with a full docstring following the template above.

3. **Port the MATLAB tests.** Translate the test logic from the corresponding
   `matlab/tests/test_sid*.m` to pytest.

4. **Write equivalence tests.** Add a test case to `test_cross_validation.py`
   that loads the JSON reference and compares outputs.

5. **Run the full suite.** `pytest python/tests/ -v` must be green.

6. **Update `__init__.py`.** Export the new function from `sid/__init__.py`.

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
