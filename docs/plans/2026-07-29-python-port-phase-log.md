# Python port — phase log (historical record)

**Date archived:** 2026-07-29
**Status:** executed (historical record)
**Origin:** the "Porting Workflow", "Phased Roadmap", "Dependency Graph" and
"Timeline" sections of `docs/roadmap_python.md`, reproduced **verbatim** below,
extracted when the roadmaps were restructured along the contract/implementation
layers (issue #195, ADR-0006). This is the phase log only — the former file was
491 lines, of which these are ~266. Its other sections were relocated (naming
convention and catalogue merged into [`docs/roadmap.md`](../roadmap.md);
"Technical Notes" merged into
[`python/CONTRIBUTING.md`](../../python/CONTRIBUTING.md); out-of-scope items
merged into the roadmap's list) or deliberately deleted ("Result Types" and the
`FreqResult` field list, which duplicated the spec and `sid/_results.py`; the
MATLAB↔Python private-helper mapping, which was already drifting and is the
port-to-port correspondence ADR-0001 rejects). That file no longer exists.

---

> **Historical document — do not treat as current guidance.**
>
> This is the execution log of the Python port, preserved for the decision
> trail. Two aspects of it are **contradicted by current project doctrine** and
> are reproduced here only because the record is archived rather than rewritten:
>
> - **"The MATLAB code in `matlab/sid/` is the ground truth for numerical
>   behaviour"** (in the preamble below) is **wrong under current doctrine.** It
>   predates [ADR-0001](../decisions/ADR-0001-spec-is-the-contract.md).
>   [`spec/SPEC.md`](../../spec/SPEC.md) is the sole contract; implementations
>   conform to the spec, **not to each other**. If MATLAB and the spec disagree,
>   MATLAB is wrong. Cross-language numerical agreement is a *consequence* of
>   each port independently satisfying the spec, never a goal pursued by copying
>   one port into another. See also `CLAUDE.md` §3.
> - **The "Porting Workflow"** below ("read MATLAB source, write Python") is the
>   copy-by-copy-porting pattern that
>   [`docs/REVIEW_CONTEXT.md`](../REVIEW_CONTEXT.md) now flags as a red flag: it
>   propagates MATLAB's bugs into Python and makes cross-validation agree for the
>   wrong reason. New work implements from the spec.
>
> Preserved verbatim below, disclaimed rather than edited.

---

## Original preamble (verbatim)

> This document tracks all phases required to achieve feature parity with the
> MATLAB/Octave v1.0 implementation. The authoritative algorithm specification
> is `spec/SPEC.md`. The MATLAB code in `matlab/sid/` is the ground truth for
> numerical behaviour.

---

## Porting Workflow

For each public function, the porting order is:

1. **Port private dependencies** — read MATLAB source, write Python, write tests
2. **Port the public function** — read MATLAB source, write Python with docstring
3. **Port MATLAB tests** — translate test logic to pytest
4. **Write equivalence tests** — load JSON reference data, compare to `rtol=1e-10`
5. **Run full suite** — verify green

This ensures every building block is tested in isolation before composition.

---

## Phased Roadmap

### Phase P1 — Scaffolding ✅

- `python/pyproject.toml` — packaging (numpy, scipy, optional matplotlib)
- `python/sid/__init__.py` — public API exports
- `python/sid/_results.py` — frozen dataclasses for all result types
- `python/sid/_exceptions.py` — `SidError` exception class
- `python/sid/_internal/__init__.py`
- `python/CONTRIBUTING.md` — Python coding standards and docstring template
- `.github/scripts/check_python_headers.py` — docstring validation
- `.github/workflows/python-lint.yml` — ruff + header check CI
- `python/tests/conftest.py` — shared fixtures, tolerance helpers

### Phase P2 — `freq_bt` SISO Core ✅

Port the Blackman-Tukey estimator and all its private dependencies.

**Private helpers (in order):**

| Step | MATLAB source | Python target | Tests |
|------|--------------|---------------|-------|
| 2.1 | `sidValidateData.m` | `_internal/validate_data.py` | `test_validate.py` |
| 2.2 | `sidHannWin.m` | `_internal/hann_win.py` | `test_hann_win.py` |
| 2.3 | `sidCov.m` | `_internal/cov.py` | `test_cov.py` |
| 2.4 | `sidDFT.m` | `_internal/dft.py` | `test_dft.py` |
| 2.5 | `sidIsDefaultFreqs.m` | `_internal/is_default_freqs.py` | (tested via windowed_dft) |
| 2.6 | `sidWindowedDFT.m` | `_internal/windowed_dft.py` | `test_windowed_dft.py` |
| 2.7 | `sidUncertainty.m` | `_internal/uncertainty.py` | `test_uncertainty.py` |

**Public function:**

| Step | MATLAB source | Python target | Tests |
|------|--------------|---------------|-------|
| 2.8 | `sidFreqBT.m` | `freq_bt.py` | `test_freq_bt.py` |

### Phase P3 — Time Series + MIMO ✅

Extend `freq_bt` to handle `u=None` (time series mode) and multi-channel
input/output (MIMO). No new files — extends Phase P2 code.

- Time series: `freq_bt(y, None)` returns output spectrum only
- MIMO: matrix-valued covariances, spectral matrix inversion, per-element uncertainty

Tests: additional cases in `test_freq_bt.py`.

### Phase P4 — Multi-Trajectory ✅

Extend `validate_data`, `cov`, `freq_bt` to accept 3D arrays `(N, nch, L)` and
`list[ndarray]` for variable-length trajectories.

- Ensemble-averaged covariances: `R_ens(tau) = (1/L) sum_l R_l(tau)`
- Variance reduction: `1/(N*L)` in uncertainty formulas
- `num_trajectories` field in output

Tests: `test_multi_trajectory.py`, additional cases in `test_freq_bt.py`.

### Phase P5 — `freq_etfe`, `freq_btfdr`, `detrend` ✅

| Step | MATLAB source | Python target | Tests |
|------|--------------|---------------|-------|
| 5.1 | `sidFreqETFE.m` | `freq_etfe.py` | `test_freq_etfe.py` |
| 5.2 | `sidFreqBTFDR.m` | `freq_btfdr.py` | `test_freq_btfdr.py` |
| 5.3 | `sidDetrend.m` | `detrend.py` | `test_detrend.py` |

### Phase P6 — `spectrogram`, `freq_map` ✅

| Step | MATLAB source | Python target | Tests |
|------|--------------|---------------|-------|
| 6.1 | `sidSpectrogram.m` | `spectrogram.py` | `test_spectrogram.py` |
| 6.2 | `sidFreqMap.m` | `freq_map.py` | `test_freq_map.py` |

`freq_map` supports two algorithms:
- `algorithm='bt'` (default): calls `freq_bt` per segment
- `algorithm='welch'`: Welch's method per segment

### Phase P7 — `ltv_disc` Core (COSMIC) ✅

**Private helpers:**

| Step | MATLAB source | Python target |
|------|--------------|---------------|
| 7.1 | `matlab/examples/util_msd.m` | `python/examples/util_msd.py` |
| 7.2 | `sidLTVbuildDataMatrices.m` | `_internal/ltv_build_data_matrices.py` |
| 7.3 | `sidLTVbuildBlockTerms.m` | `_internal/ltv_build_block_terms.py` |
| 7.4 | `sidLTVcosmicSolve.m` | `_internal/ltv_cosmic_solve.py` |
| 7.5 | `sidLTVevaluateCost.m` | `_internal/ltv_evaluate_cost.py` |

**Public function:**

| Step | MATLAB source | Python target | Tests |
|------|--------------|---------------|-------|
| 7.6 | `sidLTVdisc.m` | `ltv_disc.py` | `test_ltv_disc.py` |

### Phase P7a — Variable-Length Trajectories ✅

| Step | MATLAB source | Python target | Tests |
|------|--------------|---------------|-------|
| 7a.1 | `sidLTVbuildDataMatricesVarLen.m` | extend `ltv_build_data_matrices.py` | `test_ltv_disc_var_len.py` |

### Phase P7b — Bayesian Uncertainty ✅

**Private helpers:**

| Step | MATLAB source | Python target |
|------|--------------|---------------|
| 7b.1 | `sidLTVuncertaintyBackwardPass.m` | `_internal/ltv_uncertainty_backward_pass.py` |
| 7b.2 | `sidEstimateNoiseCov.m` | `_internal/estimate_noise_cov.py` |
| 7b.3 | `sidExtractStd.m` | `_internal/extract_std.py` |

Tests: `test_ltv_disc_uncertainty.py`

### Phase P7c — `ltv_disc_tune`, `ltv_disc_frozen` ✅

| Step | MATLAB source | Python target | Tests |
|------|--------------|---------------|-------|
| 7c.1 | `sidLTVdiscTune.m` | `ltv_disc_tune.py` | `test_ltv_disc_tune.py` |
| 7c.2 | `sidLTVdiscFrozen.m` | `ltv_disc_frozen.py` | `test_ltv_disc_frozen.py` |

### Phase P8 — Output-COSMIC ✅

Depends on P2 (frequency-domain) and P7 (COSMIC core).

**Private helpers:**

| Step | MATLAB source | Python target |
|------|--------------|---------------|
| 8.1 | `sidLTVblkTriSolve.m` | `_internal/ltv_blk_tri_solve.py` |

**Public functions:**

| Step | MATLAB source | Python target | Tests |
|------|--------------|---------------|-------|
| 8.2 | `sidLTIfreqIO.m` | `lti_freq_io.py` | `test_lti_freq_io.py` |
| 8.3 | `sidLTVStateEst.m` | `ltv_state_est.py` | `test_ltv_state_est.py` |
| 8.4 | `sidModelOrder.m` | `model_order.py` | `test_model_order.py` |
| 8.5 | `sidLTVdiscIO.m` | `ltv_disc_io.py` | `test_ltv_disc_io.py` |

### Phase P9 — Workflow Utilities ✅

**Private helpers:**

| Step | MATLAB source | Python target |
|------|--------------|---------------|
| 9.1 | `sidFreqDomainSim.m` | `_internal/freq_domain_sim.py` |

**Public functions:**

| Step | MATLAB source | Python target | Tests |
|------|--------------|---------------|-------|
| 9.2 | `sidResidual.m` | `residual.py` | `test_residual.py` |
| 9.3 | `sidCompare.m` | `compare.py` | `test_compare.py` |

### Phase P10 — Plotting ✅

All plotting functions lazy-import matplotlib and accept an optional `ax=`
keyword argument for embedding in user-created figures.

| Step | MATLAB source | Python target | Tests |
|------|--------------|---------------|-------|
| 10.1 | `sidBodePlot.m` | `bode_plot.py` | `test_plotting.py` |
| 10.2 | `sidSpectrumPlot.m` | `spectrum_plot.py` | `test_plotting.py` |
| 10.3 | `sidMapPlot.m` | `map_plot.py` | `test_plotting.py` |
| 10.4 | `sidSpectrogramPlot.m` | `spectrogram_plot.py` | `test_plotting.py` |

### Phase P10a — Examples (Jupyter Notebooks) ✅

Each MATLAB example script is ported to a Jupyter notebook in `python/examples/`.
Notebooks are the natural Python equivalent of MATLAB's `%%`-sectioned scripts:
they combine narrative, code, and inline plots in a single runnable document.

**Example catalog:**

| MATLAB example | Python notebook | Functions demonstrated | Depends on |
|---|---|---|---|
| `exampleSISO.m` | `example_siso.ipynb` | `freq_bt`, `bode_plot`, `spectrum_plot` | P2, P10 |
| `exampleETFE.m` | `example_etfe.ipynb` | `freq_etfe`, `bode_plot`, `spectrum_plot` | P5, P10 |
| `exampleFreqDepRes.m` | `example_freq_dep_res.ipynb` | `freq_btfdr`, `freq_bt`, `bode_plot` | P5, P10 |
| `exampleCoherence.m` | `example_coherence.ipynb` | `freq_bt`, `bode_plot` | P2, P10 |
| `exampleMethodComparison.m` | `example_method_comparison.ipynb` | `freq_bt`, `freq_btfdr`, `freq_etfe` | P5, P10 |
| `exampleMIMO.m` | `example_mimo.ipynb` | `freq_bt` (MIMO mode) | P3, P10 |
| `exampleFreqMap.m` | `example_freq_map.ipynb` | `freq_map`, `map_plot` | P6, P10 |
| `exampleSpectrogram.m` | `example_spectrogram.ipynb` | `spectrogram`, `spectrogram_plot` | P6, P10 |
| `exampleLTVdisc.m` | `example_ltv_disc.ipynb` | `ltv_disc`, `ltv_disc_tune`, `ltv_disc_frozen` | P7c, P10 |
| `exampleMultiTrajectory.m` | `example_multi_trajectory.ipynb` | `freq_bt`, `freq_map`, `spectrogram`, `ltv_disc` | P4, P6, P7 |
| `exampleOutputCOSMIC.m` | `example_output_cosmic.ipynb` | `freq_bt`, `model_order`, `ltv_disc_io` | P8, P10 |

**Notebook conventions:**
- One notebook per MATLAB example, mirroring its structure section by section
- Markdown cells for narrative (replacing MATLAB `%%` section comments)
- Each code cell corresponds to one MATLAB `%%` section
- Inline plots via `%matplotlib inline` (no `figure;` / `hold on` boilerplate)
- All notebooks can be run top-to-bottom without external data files
- `python/examples/README.md` provides the index and descriptions

**CI validation:**
- Notebooks are validated using `nbval` or `pytest --nbmake` in CI to ensure
  they execute without errors (output cells are not checked, only execution)

### Phase P11 — Cross-Validation + CI ✅

- `python/tests/test_cross_validation.py` — load JSON from `testdata/`, call
  Python functions, assert `rtol=1e-10`
- `.github/workflows/python-tests.yml` — pytest on Python 3.10–3.13
- Update `.github/workflows/cross-validate.yml` — add `validate-python` job

---

## Dependency Graph

```
P1 (scaffolding)
├── P2 (freq_bt core)
│   ├── P3 (time series + MIMO)
│   │   └── P4 (multi-trajectory)
│   ├── P5 (freq_etfe, freq_btfdr, detrend)
│   │   └── P6 (spectrogram, freq_map)
│   │       └── P10 (plotting)
│   │           └── P10a (examples / notebooks)
│   ├── P8 (Output-COSMIC) ←── also needs P7
│   └── P9 (residual, compare) ←── also needs P7
├── P7 (ltv_disc core)
│   ├── P7a (variable-length)
│   ├── P7b (uncertainty)
│   │   └── P7c (tune, frozen)
│   ├── P8 (Output-COSMIC)
│   └── P9 (residual, compare)
└── P11 (cross-validation + CI) ←── after all phases
```

---

## Timeline

| Phase | Description | Dependencies | Status |
|-------|-------------|-------------|--------|
| P1 | Scaffolding | — | ✅ |
| P2 | `freq_bt` SISO core | P1 | ✅ |
| P3 | Time series + MIMO | P2 | ✅ |
| P4 | Multi-trajectory | P3 | ✅ |
| P5 | `freq_etfe`, `freq_btfdr`, `detrend` | P2 | ✅ |
| P6 | `spectrogram`, `freq_map` | P2, P5 | ✅ |
| P7 | `ltv_disc` core (COSMIC) | P1 | ✅ |
| P7a | Variable-length trajectories | P7 | ✅ |
| P7b | Bayesian uncertainty | P7 | ✅ |
| P7c | `ltv_disc_tune`, `ltv_disc_frozen` | P7, P7b | ✅ |
| P8 | Output-COSMIC | P2, P7 | ✅ |
| P9 | Workflow utilities | P2, P7 | ✅ |
| P10 | Plotting | P2, P6 | ✅ |
| P10a | Examples (Jupyter notebooks) | P10, varies per notebook | ✅ |
| P11 | Cross-validation + CI | all | ✅ |

