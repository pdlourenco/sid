# sid v0.1.0 — python — Release Notes

**Release date:** 2026-04-11
**License:** MIT

## Overview

sid is a free, open-source system identification toolbox. This release is
the first public release of the **Python** port. It is numerically
equivalent to the MATLAB/Octave reference implementation to within
floating-point tolerance, and every public function is cross-validated in
CI against JSON reference vectors produced by the MATLAB port.

The toolbox covers two complementary paths: a frequency-domain path built on
Blackman-Tukey spectral analysis (Ljung, 1999), and a state-space path built
on the COSMIC algorithm — a closed-form, O(N)-complexity solver for discrete
linear time-varying system identification (Carvalho et al., 2022). Both paths
support SISO, MIMO, time series, and multi-trajectory data through a unified
API.

The package depends only on **NumPy** and **SciPy**; **Matplotlib** is an
optional extra needed for the plotting functions and for running the example
notebooks.

## Release Contents

```
sid-v0.1.0-python/
├── README.md                              # Overview, features, references
├── LICENSE                                # MIT License
├── spec/
│   ├── SPEC.md                            # Full algorithm specification
│   ├── EXAMPLES.md                        # Example-suite specification
│   └── cosmic/                            # COSMIC theory documents
├── testdata/                              # Cross-language reference vectors (JSON)
├── binder/                                # MyBinder environment definition
│   ├── runtime.txt                        # python-3.11
│   ├── requirements.txt                   # numpy, scipy, matplotlib
│   └── postBuild                          # pip install -e ./python[plot]
└── python/
    ├── README.md                          # Installation, quick start, function reference
    ├── CONTRIBUTING.md                    # Python-specific contribution guide
    ├── pyproject.toml                     # Packaging metadata (sid-toolbox 0.1.0)
    ├── sid/                               # Public API
    │   ├── __init__.py                    # Public re-exports, __version__ = "0.1.0"
    │   ├── _exceptions.py                 # SidError exception class
    │   ├── _results.py                    # Frozen dataclasses for result types
    │   ├── freq_bt.py                     # Blackman-Tukey spectral estimation
    │   ├── freq_btfdr.py                  # BT with frequency-dependent resolution
    │   ├── freq_etfe.py                   # Empirical transfer function estimate
    │   ├── freq_map.py                    # Time-varying frequency response map
    │   ├── spectrogram.py                 # Short-time FFT spectrogram
    │   ├── ltv_disc.py                    # COSMIC: LTV state-space identification
    │   ├── ltv_disc_tune.py               # Regularisation tuning (validation, frequency)
    │   ├── ltv_disc_frozen.py             # Frozen transfer function G(ω,k)
    │   ├── ltv_disc_io.py                 # Output-COSMIC: partial-observation LTV ID
    │   ├── ltv_state_est.py               # Batch LTV state estimation (RTS smoother)
    │   ├── lti_freq_io.py                 # LTI realisation from I/O freq. response
    │   ├── model_order.py                 # Model order estimation (Hankel SVD)
    │   ├── detrend.py                     # Polynomial detrending
    │   ├── residual.py                    # Residual analysis (whiteness + independence)
    │   ├── compare.py                     # Model output comparison (NRMSE fit)
    │   ├── bode_plot.py                   # Bode diagram with confidence bands
    │   ├── spectrum_plot.py               # Power spectrum with confidence bands
    │   ├── map_plot.py                    # Time-frequency colour map
    │   ├── spectrogram_plot.py            # Spectrogram colour map
    │   └── _internal/                     # Internal helper modules
    │       ├── cov.py                     #   Biased cross-covariance estimation
    │       ├── dft.py                     #   Discrete Fourier transform
    │       ├── hann_win.py                #   Hann lag window
    │       ├── windowed_dft.py            #   Windowed FFT (fast + direct paths)
    │       ├── uncertainty.py             #   Asymptotic variance formulas
    │       ├── validate_data.py           #   Input validation and parsing
    │       ├── is_default_freqs.py        #   Default frequency grid detection
    │       ├── freq_domain_sim.py         #   Frequency-domain model simulation
    │       ├── estimate_noise_cov.py      #   Noise covariance estimation
    │       ├── extract_std.py             #   Standard deviation extraction
    │       ├── test_msd.py                #   Mass-spring-damper test system
    │       ├── ltv_build_data_matrices.py #   COSMIC data matrix construction (+ var-len)
    │       ├── ltv_build_block_terms.py   #   Block tridiagonal term assembly
    │       ├── ltv_cosmic_solve.py        #   COSMIC forward-backward solve
    │       ├── ltv_blk_tri_solve.py       #   Generic block tridiagonal solver
    │       ├── ltv_evaluate_cost.py       #   COSMIC cost function evaluation
    │       └── ltv_uncertainty_backward_pass.py # Posterior covariance recursion
    ├── examples/                          # Jupyter notebooks
    │   ├── README.md                      # Examples guide
    │   ├── util_msd.py                    # SMD helpers (ZOH, LTV stack, Duffing RK4)
    │   ├── example_siso.ipynb             # Basic SISO frequency response
    │   ├── example_etfe.ipynb             # Empirical transfer function
    │   ├── example_freq_dep_res.ipynb     # Frequency-dependent resolution
    │   ├── example_coherence.ipynb        # Coherence analysis
    │   ├── example_method_comparison.ipynb # BT vs. BTFDR vs. ETFE
    │   ├── example_mimo.ipynb             # MIMO estimation
    │   ├── example_freq_map.ipynb         # Time-varying frequency maps
    │   ├── example_spectrogram.ipynb      # Spectrogram
    │   ├── example_ltv_disc.ipynb         # LTV identification with COSMIC
    │   ├── example_multi_trajectory.ipynb # Multi-trajectory ensemble averaging
    │   └── example_output_cosmic.ipynb    # Partial-observation identification
    └── tests/                             # pytest suite (Python 3.10–3.13 in CI)
        ├── conftest.py
        ├── test_cross_validation.py       # MATLAB ↔ Python numerical equivalence
        └── test_*.py
```

## Installation

From PyPI-style install against the source tree:

```bash
pip install -e ./python
```

With plotting support (adds Matplotlib):

```bash
pip install -e "./python[plot]"
```

For development (tests + linting + notebook validation):

```bash
pip install -e "./python[dev]"
```

Or clone the full repository:

```bash
git clone https://github.com/pdlourenco/sid.git
cd sid
pip install -e "./python[plot]"
```

### Try it in the browser

Every example notebook can be launched on MyBinder without any local
installation — click the **Binder** badge in
[`python/README.md`](README.md) to open a ready-to-run Jupyter environment.

## Quick Start

```python
import numpy as np
from scipy.signal import lfilter
import sid

# SISO frequency response estimation
N = 1000
Ts = 0.01
rng = np.random.default_rng(42)
u = rng.standard_normal(N)
y = lfilter([1], [1, -0.9], u) + 0.1 * rng.standard_normal(N)

result = sid.freq_bt(y, u, sample_time=Ts)
sid.bode_plot(result)

# Time series (output spectrum only)
result_ts = sid.freq_bt(y, None)
sid.spectrum_plot(result_ts)
```

See the [examples guide](examples/README.md) for ETFE, frequency-dependent
resolution, MIMO, time-varying maps, spectrograms, LTV state-space
identification, multi-trajectory ensemble averaging, and partial-observation
Output-COSMIC.

## Function Reference

All public functions are snake_case. MATLAB `Name,Value` pairs map to
Python keyword arguments (e.g. `'SampleTime', 0.01` → `sample_time=0.01`).
Every public function returns an immutable `@dataclass(frozen=True)` with
tab-completable attributes.

### Frequency-Domain Estimation

| Function | MATLAB equivalent | Description |
|---|---|---|
| `sid.freq_bt` | `sidFreqBT` / `spa` | Blackman-Tukey frequency response and noise spectrum estimation. Supports SISO, MIMO, time series, and multi-trajectory inputs (3D arrays `(N, n_ch, L)` and lists of ndarrays). |
| `sid.freq_btfdr` | `sidFreqBTFDR` / `spafdr` | Blackman-Tukey with frequency-dependent resolution. Adapts the correlation window length across frequency. |
| `sid.freq_etfe` | `sidFreqETFE` / `etfe` | Empirical transfer function estimate via FFT ratio, with optional smoothing. |

### Time-Frequency Analysis

| Function | MATLAB equivalent | Description |
|---|---|---|
| `sid.freq_map` | `sidFreqMap` / `tfestimate` | Sliding-window time-varying frequency response map using Blackman-Tukey or Welch windowing. Produces time-varying transfer function, coherence, and noise spectrum. |
| `sid.spectrogram` | `sidSpectrogram` / `spectrogram` | Short-time FFT spectrogram with configurable window, overlap, and FFT length. |

### State-Space Identification (COSMIC)

| Function | MATLAB equivalent | Description |
|---|---|---|
| `sid.ltv_disc` | `sidLTVdisc` | COSMIC algorithm for discrete LTV system identification from full-state measurements. Identifies A(k), B(k) from single or multiple trajectories (including variable-length via lists of ndarrays) with time-varying or constant regularisation `lambda_`. Returns system matrices and optional Bayesian posterior uncertainty. |
| `sid.ltv_disc_io` | `sidLTVdiscIO` | Partial-observation LTV identification (Output-COSMIC). Alternates between state estimation (RTS smoother) and dynamics identification (COSMIC) to jointly recover A(k), B(k), and hidden states. Includes fast path when H is full-rank. |
| `sid.ltv_disc_tune` | `sidLTVdiscTune` | Regularisation tuning via validation-based grid search, L-curve, or frequency-response consistency scoring (frozen transfer function vs. `sid.freq_map` with Mahalanobis distance). |
| `sid.ltv_disc_frozen` | `sidLTVdiscFrozen` | Frozen transfer function G(ω, k) from identified A(k), B(k) with uncertainty propagated from the COSMIC posterior, enabling direct comparison with non-parametric frequency estimates. |
| `sid.lti_freq_io` | `sidLTIfreqIO` | LTI realisation from I/O frequency response via Ho-Kalman algorithm with H-basis transform. Used by `sid.ltv_disc_io` for initialisation; also usable standalone. |
| `sid.ltv_state_est` | `sidLTVStateEst` | Batch LTV state estimation via RTS smoother given A(k), B(k), H, R, and optional Q. Usable standalone or as a building block for `sid.ltv_disc_io`. |
| `sid.model_order` | `sidModelOrder` | Model order estimation from any `sid.freq_*` result via Hankel SVD. Supports automatic gap detection and user-specified threshold. |

### Analysis and Validation

| Function | MATLAB equivalent | Description |
|---|---|---|
| `sid.detrend` | `sidDetrend` | Polynomial detrending. Removes constant, linear, or higher-order trends. Multi-channel support. |
| `sid.residual` | `sidResidual` | Residual analysis with whiteness (autocorrelation) and independence (cross-correlation with input) tests, with confidence bounds. Works with both frequency-domain and state-space models. |
| `sid.compare` | `sidCompare` | Model output simulation and comparison. Computes NRMSE fit metric. Supports frequency-domain and state-space models, multi-channel, and multi-trajectory data. |

### Plotting

Plotting functions lazy-import Matplotlib and accept an optional `ax=`
keyword for embedding in user-created figures.

| Function | MATLAB equivalent | Description |
|---|---|---|
| `sid.bode_plot` | `sidBodePlot` | Bode diagram with shaded confidence bands. |
| `sid.spectrum_plot` | `sidSpectrumPlot` | Power spectrum plot with shaded confidence bands. |
| `sid.map_plot` | `sidMapPlot` | Time-frequency colour map for `sid.freq_map` results. |
| `sid.spectrogram_plot` | `sidSpectrogramPlot` | Spectrogram colour map for `sid.spectrogram` results. |

### Result Types

Every public function returns a frozen dataclass from `sid._results`, with
tab-completable fields accessed via dot notation (`result.response`,
`result.fit`, etc.).

| Result type | Returned by | Key fields |
|---|---|---|
| `FreqResult` | `freq_bt`, `freq_btfdr`, `freq_etfe` | `response`, `noise_spectrum`, `coherence`, `response_std`, `frequency` |
| `FreqMapResult` | `freq_map` | `time`, `response`, `noise_spectrum`, `coherence`, `frequency` |
| `SpectrogramResult` | `spectrogram` | `time`, `frequency`, `power`, `power_db`, `complex_stft` |
| `LTVResult` | `ltv_disc`, `ltv_disc_tune` | `a`, `b`, `lambda_`, `cost`, `a_std`, `b_std` |
| `LTVIOResult` | `ltv_disc_io` | `a`, `b`, `x`, `cost`, `iterations` |
| `FrozenResult` | `ltv_disc_frozen` | `response`, `response_std`, `time_steps`, `frequency` |
| `CompareResult` | `compare` | `predicted`, `measured`, `fit`, `residual` |
| `ResidualResult` | `residual` | `residual`, `auto_corr`, `cross_corr`, `whiteness_pass`, `independence_pass` |

## Key Features

- **Numerically equivalent to MATLAB/Octave.** Every public function is
  cross-validated in CI against JSON reference vectors produced by the
  MATLAB implementation (`python/tests/test_cross_validation.py`). The
  binding contract is [`spec/SPEC.md`](../spec/SPEC.md); both language
  ports derive from the spec independently, not by copying one
  implementation into the other.

- **Closed-form, O(N) LTV identification.** COSMIC solves the regularised
  least-squares problem in a single forward-backward pass — no iterations,
  no convergence tuning. Complexity is linear in the number of time steps
  and cubic in the state dimension, independent of the number of
  trajectories.

- **Partial-observation identification.** Output-COSMIC extends the
  algorithm to systems where only a subset of states is measured. An
  alternating minimisation between dynamics identification and state
  estimation converges to a joint optimum, with a fast path that avoids
  iteration entirely when the observation matrix is full-rank.

- **Multi-trajectory support.** Frequency-domain functions accept 3D
  arrays `(N, n_ch, L)` and lists of ndarrays for ensemble averaging with
  1/L variance reduction. COSMIC pools multiple trajectories sharing the
  same time-varying dynamics, including variable-length trajectories.

- **Bayesian uncertainty quantification.** Per-timestep posterior
  covariance of the identified system matrices via a backward recursion
  at O(N) cost. The posterior has matrix-normal structure and the MAP
  estimate is independent of the noise covariance. Uncertainty propagates
  through to frozen transfer functions for comparison with non-parametric
  frequency estimates.

- **Immutable, tab-completable results.** Every public function returns a
  `@dataclass(frozen=True)` with type-annotated attributes, so result
  fields autocomplete in IPython/Jupyter and cannot be mutated after
  construction.

- **Minimal runtime dependencies.** NumPy ≥ 1.22 and SciPy ≥ 1.8 are the
  only required runtime dependencies. Matplotlib is an optional extra
  needed only for the plotting functions and the example notebooks.

## Compatibility

| Dependency | Minimum version | Notes |
|---|---|---|
| Python | 3.10 | Tested in CI on 3.10, 3.11, 3.12, 3.13 |
| NumPy | 1.22 | Required |
| SciPy | 1.8 | Required |
| Matplotlib | 3.5 | Optional (`sid-toolbox[plot]`) — needed for plotting functions and example notebooks |

## Testing

```bash
# Unit tests
pytest python/tests/ -v

# Notebook execution (requires nbmake)
pytest --nbmake python/examples/ -v

# Cross-language validation (testdata/*.json reference vectors)
pytest python/tests/test_cross_validation.py -v
```

## Out of Scope for this version

- Online/recursive COSMIC (deferred to a future release)
- Parametric identification: ARX, ARMAX, state-space subspace methods
- LPV identification
- Unknown or time-varying observation matrix H
- Alternative regularisation norms
- MATLAB/Octave and Julia implementations — MATLAB/Octave ships as a
  separate v0.1.0 release (`v0.1-matlab`); Julia is planned for a later
  version.

## References

- Ljung, L. (1999). *System Identification: Theory for the User*, 2nd ed.
  Prentice Hall.
- Blackman, R. B. & Tukey, J. W. (1959). *The Measurement of Power Spectra*.
  Dover.
- Carvalho, M., Soares, C., Lourenco, P., and Ventura, R. (2022). "COSMIC:
  fast closed-form identification from large-scale data for LTV systems."
  [arXiv:2112.04355v2](https://arxiv.org/abs/2112.04355v2)
- Laszkiewicz, P., Carvalho, M., Soares, C., and Lourenco, P. (2025). "The
  impact of modeling approaches on controlling safety-critical, highly
  perturbed systems: the case for data-driven models."
  [arXiv:2509.13531](https://arxiv.org/abs/2509.13531)
