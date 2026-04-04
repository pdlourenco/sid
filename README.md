# sid — Open-Source System Identification for MATLAB/Octave

![Tests](https://github.com/pdlourenco/sid-matlab/actions/workflows/tests.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![MATLAB R2016b+](https://img.shields.io/badge/MATLAB-R2016b%2B-orange.svg)
![GNU Octave 8+](https://img.shields.io/badge/GNU_Octave-8%2B-blue.svg)

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=pdlourenco/sid-matlab&file=examples/exampleSISO.m)

**sid** is a free, open-source toolbox for **system identification** — covering both **non-parametric frequency response estimation** and **time-varying state-space identification** — implemented in pure MATLAB/Octave code with **zero toolbox dependencies**. It provides drop-in replacements for `spa`, `spafdr`, and `etfe` that you can use without a licence, and implements the **COSMIC** (**C**losed-form **O**ptimal data-driven linear time-varying **S**yste**M** **I**dentifi**C**ation) algorithm for identifying linear time-varying (LTV) state-space models from data.

All estimation functions natively support **multiple trajectories** — pass L independent experiments as a 3D array `(N x n x L)` or, for the time-domain functions, as a cell array of variable-length sequences. Frequency-domain estimates are ensemble-averaged for lower variance; COSMIC pools trajectories into a single regularized least-squares problem. sid also exposes **asymptotic uncertainty estimates** and built-in **confidence-band plotting** out of the box. It supports SISO, MIMO, and time-series (output-only) modes, all with a unified result struct and consistent API.

sid is designed from the ground up to run on **GNU Octave** as a first-class target — not just MATLAB with Octave as an afterthought. Every function is tested on both platforms in CI, so Octave users get the same reliability as MATLAB users. If you work in an environment where MATLAB licences are limited or unavailable, sid and Octave give you a fully open-source path to system identification.

## Features

### Frequency-domain estimation
- **Blackman-Tukey spectral analysis** (`sidFreqBT`) — the workhorse estimator, with configurable window size and automatic defaults
- **Frequency-dependent resolution** (`sidFreqBTFDR`) — vary the smoothing bandwidth across the frequency axis
- **Empirical transfer function estimate** (`sidFreqETFE`) — maximum resolution via FFT ratio, with optional smoothing
- **Time-varying frequency maps** (`sidFreqMap`) — sliding-window frequency response maps using Blackman-Tukey or Welch inner estimators; `sidSpectrogram` for short-time FFT spectrograms
- **Asymptotic uncertainty** — standard deviations and squared coherence for every estimate (Ljung, 1999), rendered as shaded confidence bands by `sidBodePlot`, `sidSpectrumPlot`, `sidMapPlot`, and `sidSpectrogramPlot`

### Time-varying state-space identification (COSMIC)
- **LTV identification** (`sidLTVdisc`) — the COSMIC algorithm: identify time-varying A(k), B(k) from state measurements with O(N) complexity, automatic or manual regularization, and optional Bayesian uncertainty quantification
- **Output-COSMIC** (`sidLTVdiscIO`) — identify A(k), B(k) from output-only observations when the full state is not measured, with alternating state-dynamics optimization and automatic LTI initialization
- **Frozen transfer function** (`sidLTVdiscFrozen`) — compute instantaneous G(w,k) = (e^{jw}I - A(k))^{-1} B(k) with propagated uncertainty bands
- **Regularization tuning** (`sidLTVdiscTune`) — select lambda via held-out validation loss, L-curve, or frequency-response consistency scoring
- **LTI realization** (`sidLTIfreqIO`) — Ho-Kalman realization from I/O frequency response data
- **State estimation** (`sidLTVStateEst`) — batch LTV state estimation via RTS smoother
- **Model order estimation** (`sidModelOrder`) — Hankel singular value analysis for selecting state dimension

### Multi-trajectory support
- **Frequency-domain**: pass `(N x n x L)` arrays or cell arrays of variable-length sequences — spectral estimates are ensemble-averaged across L trajectories, reducing variance by a factor of L
- **COSMIC**: pool L trajectories (uniform or variable-length) into a single identification problem — data matrices are stacked across trajectories while the regularization structure is shared

### Analysis, validation, and preprocessing
- `sidDetrend` — polynomial detrending (replaces `detrend`)
- `sidResidual` — residual whiteness and independence diagnostics (replaces `resid`)
- `sidCompare` — model output comparison with NRMSE fit metric (replaces `compare`)

### Cross-cutting
- **SISO, MIMO, and time-series modes** — unified API across all estimation functions
- **Validated against MATLAB's System Identification Toolbox** — comparison tests for `spa`, `spafdr`, and `etfe` run in CI
- **34 test suites** with continuous integration on both MATLAB and GNU Octave

## How It Works

**Frequency-domain path.** The core spectral estimators use the **Blackman-Tukey method**: compute biased cross-covariances between input and output, apply a Hann lag window, then transform via FFT. The transfer function is the cross-spectrum / input auto-spectrum ratio; asymptotic variance formulas (Ljung, 1999) provide per-frequency uncertainty. When multiple trajectories are provided, covariances are ensemble-averaged before forming the ratio, reducing variance by a factor of L without sacrificing frequency resolution.

**State-space path.** The **COSMIC algorithm** (Carvalho et al., 2022) identifies discrete-time LTV models x(k+1) = A(k)x(k) + B(k)u(k) by solving a block-tridiagonal regularized least-squares problem in O(N) time. Multiple trajectories — including variable-length sequences — are pooled into the data matrices. When only outputs are observed, **Output-COSMIC** alternates between state estimation (RTS smoother) and dynamics identification, converging to a joint optimum. Bayesian uncertainty quantification propagates through to frozen transfer functions G(w,k) for direct comparison with non-parametric frequency estimates.

See [SPEC.md](SPEC.md) for the full mathematical derivation.

## Function Comparison

**Frequency-domain estimation:**

| sid function | Replaces | Description |
|---|---|---|
| `sidFreqBT` | `spa` | Blackman-Tukey frequency response and noise spectrum estimation |
| `sidFreqBTFDR` | `spafdr` | Blackman-Tukey with frequency-dependent resolution |
| `sidFreqETFE` | `etfe` | Empirical transfer function estimate (FFT ratio) |

**Time-frequency analysis:**

| sid function | Replaces | Description |
|---|---|---|
| `sidFreqMap` | `tfestimate`, `mscohere` | Time-varying frequency response map (BT or Welch) |
| `sidSpectrogram` | `spectrogram` | Short-time FFT spectrogram |

**State-space identification (COSMIC):**

| sid function | Replaces | Description |
|---|---|---|
| `sidLTVdisc` | — | Discrete LTV state-space identification (COSMIC algorithm) |
| `sidLTVdiscIO` | — | LTV identification from partial observations (Output-COSMIC) |
| `sidLTVdiscTune` | — | Regularization tuning via validation loss, L-curve, or frequency consistency |
| `sidLTVdiscFrozen` | — | Frozen transfer function G(w,k) with uncertainty propagation |
| `sidLTIfreqIO` | — | LTI realization from I/O frequency response (Ho-Kalman) |
| `sidLTVStateEst` | — | Batch LTV state estimation (RTS smoother) |
| `sidModelOrder` | — | Model order estimation via Hankel SVD |

**Plotting:**

| sid function | Replaces | Description |
|---|---|---|
| `sidBodePlot` | — | Bode diagram with shaded confidence bands |
| `sidSpectrumPlot` | — | Power spectrum plot with shaded confidence bands |
| `sidMapPlot` | — | Time-frequency color map for `sidFreqMap` results |
| `sidSpectrogramPlot` | — | Spectrogram color map |

**Analysis and validation:**

| sid function | Replaces | Description |
|---|---|---|
| `sidDetrend` | `detrend` | Remove mean, linear, or polynomial trends from signals |
| `sidResidual` | `resid` | Residual analysis with auto/cross-correlation diagnostics |
| `sidCompare` | `compare` | Simulate model output and compute fit percentage |

All estimation functions support both positional and name-value calling conventions.

## Installation

Clone the repository:

```bash
git clone https://github.com/pdlourenco/sid-matlab.git
```

Then add it to your MATLAB or Octave path:

```matlab
run('/path/to/sid-matlab/sidInstall.m')
```

To make the path persistent across sessions, add the line above to your [`startup.m`](https://www.mathworks.com/help/matlab/ref/startup.html) file (or [`.octaverc`](https://docs.octave.org/latest/Startup-Files.html) for Octave). No `pkg install` is needed — sid is a plain directory of `.m` files that works on both platforms.

## Quick Start

```matlab
% Generate example data: first-order system with noise
N = 1000; Ts = 0.01;
u = randn(N, 1);
y = filter([1], [1 -0.9], u) + 0.1 * randn(N, 1);

% Estimate frequency response (Blackman-Tukey)
result = sidFreqBT(y, u, 'SampleTime', Ts);

% Plot Bode diagram and noise spectrum with confidence bands
sidBodePlot(result);
sidSpectrumPlot(result);
```

Time-series mode (output spectrum only, no input signal):

```matlab
result_ts = sidFreqBT(y, []);
sidSpectrumPlot(result_ts);
```

See the [examples guide](examples/README.md) for more usage patterns, including ETFE, frequency-dependent resolution, MIMO, time-varying maps, spectrograms, and LTV state-space identification.

## Compatibility

| Platform | Version | Status |
|---|---|---|
| **MATLAB** | R2016b or later | Tested in CI |
| **GNU Octave** | 8.0 or later | Tested in CI |

No toolboxes are required. The entire codebase uses only core MATLAB/Octave functions (`fft`, `filter`, `conv`, etc.), so it runs anywhere the base language does — including MATLAB Online.

## Documentation

- [**SPEC.md**](SPEC.md) — Full algorithm specification with mathematical derivations
- [**Roadmap**](docs/sid_matlab_roadmap.md) — Development phases and planned features
- [**COSMIC uncertainty derivation**](docs/cosmic_uncertainty_derivation.md) — Bayesian posterior covariance for LTV identification
- [**COSMIC online recursion**](docs/cosmic_online_recursion.md) — Recursive/streaming formulation of the COSMIC algorithm
- [**COSMIC automatic tuning**](docs/cosmic_automatic_tuning.md) — Regularization parameter selection via validation and L-curve
- [**Output-COSMIC**](docs/cosmic_output.md) — LTV identification from partial (output-only) observations
- [**Multi-trajectory spectral theory**](docs/multi_trajectory_spectral_theory.md) — Averaging frequency estimates across repeated experiments
- [**Examples**](examples/README.md) — Runnable scripts demonstrating typical workflows

## References

- Ljung, L. (1999). *System Identification: Theory for the User*, 2nd ed. Prentice Hall.
- Blackman, R. B. & Tukey, J. W. (1959). *The Measurement of Power Spectra*. Dover.
- Carvalho, M., Soares, C., Lourenço, P., and Ventura, R. (2022). "COSMIC: fast closed-form identification from large-scale data for LTV systems." [arXiv:2112.04355v2](https://arxiv.org/abs/2112.04355v2)
- Łaszkiewicz, P., Carvalho, M., Soares, C., and Lourenço, P. (2025). "The impact of modeling approaches on controlling safety-critical, highly perturbed systems: the case for data-driven models." [arXiv:2509.13531](https://arxiv.org/abs/2509.13531)

## Contributing

Contributions are welcome via issues and pull requests. See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on function headers, code style, and testing. Please ensure that `tests/runAllTests.m` passes on both MATLAB and Octave before submitting.

## License

MIT License. Copyright (c) 2026 Pedro Lourenço. See [LICENSE](LICENSE).
