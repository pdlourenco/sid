# sid — Open-Source System Identification for MATLAB/Octave

![Tests](https://github.com/pdlourenco/sid-matlab/actions/workflows/tests.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![MATLAB R2016b+](https://img.shields.io/badge/MATLAB-R2016b%2B-orange.svg)
![GNU Octave 8+](https://img.shields.io/badge/GNU_Octave-8%2B-blue.svg)

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=pdlourenco/sid-matlab&file=examples/exampleSISO.m)

**sid** is a free, open-source toolbox for **non-parametric frequency response estimation** — the same class of algorithms found in MATLAB's commercial [System Identification Toolbox](https://www.mathworks.com/products/sysid.html), but implemented in pure MATLAB/Octave code with **zero toolbox dependencies**. It provides drop-in replacements for `spa`, `spafdr`, and `etfe` that you can use without a licence.

Unlike the commercial toolbox, sid exposes **asymptotic uncertainty estimates** and built-in **confidence-band plotting** out of the box, making it straightforward to assess estimation quality. It supports SISO, MIMO, and time-series (output-only) modes, all with a unified result struct and consistent API.

sid is designed from the ground up to run on **GNU Octave** as a first-class target — not just MATLAB with Octave as an afterthought. Every function is tested on both platforms in CI, so Octave users get the same reliability as MATLAB users. If you work in an environment where MATLAB licences are limited or unavailable, sid and Octave give you a fully open-source path to frequency-domain system identification.

## Features

- **Blackman-Tukey spectral analysis** (`sidFreqBT`) — the workhorse estimator, with configurable window size and automatic defaults
- **Frequency-dependent resolution** (`sidFreqBTFDR`) — vary the smoothing bandwidth across the frequency axis
- **Empirical transfer function estimate** (`sidFreqETFE`) — maximum resolution via FFT ratio, with optional smoothing
- **Multi-trajectory averaging** — pool frequency estimates across repeated experiments for lower variance
- **Time-varying analysis** — `sidFreqMap` for sliding-window frequency response maps (Blackman-Tukey or Welch), `sidSpectrogram` for short-time FFT spectrograms
- **LTV state-space identification** (`sidLTVdisc`) — the COSMIC algorithm for identifying time-varying A(k), B(k) matrices from state measurements, with automatic or manual regularization tuning and optional Bayesian uncertainty quantification
- **Output-COSMIC** (`sidLTVdiscIO`) — LTV identification from partial (output-only) observations, with variable-length trajectory support
- **LTI realization** (`sidLTIfreqIO`) — Ho-Kalman realization from input-output frequency response data
- **State estimation** (`sidLTVStateEst`) — batch LTV state estimation via RTS smoother
- **Model order estimation** (`sidModelOrder`) — Hankel singular value analysis for selecting state dimension
- **Frozen transfer function** (`sidLTVdiscFrozen`) — compute instantaneous G(w,k) = (e^{jw}I - A(k))^{-1} B(k) with propagated uncertainty bands
- **Analysis and validation** — `sidDetrend` for signal preprocessing, `sidResidual` for residual diagnostics, `sidCompare` for model-vs-data comparison
- **Asymptotic uncertainty** — standard deviations and squared coherence returned for every frequency-domain estimate (Ljung, 1999)
- **Confidence-band plotting** — `sidBodePlot`, `sidSpectrumPlot`, `sidMapPlot`, and `sidSpectrogramPlot` render shaded confidence bands out of the box
- **SISO, MIMO, and time-series modes** — unified API across all frequency-domain estimation functions
- **Validated against MATLAB's System Identification Toolbox** — comparison tests for `spa`, `spafdr`, and `etfe` run in CI
- **34 test suites** with continuous integration on both MATLAB and GNU Octave

## How It Works

The core frequency-domain estimators use the **Blackman-Tukey method**: compute biased cross-covariances between input and output, apply a Hann lag window for spectral smoothing, then transform to the frequency domain via FFT (or direct DFT for custom frequency grids). The transfer function estimate is the ratio of cross-spectrum to input auto-spectrum, and the noise spectrum is obtained by subtraction. Asymptotic variance formulas from Ljung (1999) provide uncertainty estimates at each frequency, which are rendered as shaded confidence bands in the plotting functions. For time-varying state-space identification, the toolbox implements the **COSMIC algorithm** (Carvalho et al., 2022): regularized least-squares estimation of A(k), B(k) matrices with optional Bayesian uncertainty, supporting both full-state and output-only observations. See [SPEC.md](SPEC.md) for the full mathematical derivation.

## Function Comparison

**Frequency-domain estimation:**

| sid function | Replaces | Description |
|---|---|---|
| `sidFreqBT` | `spa` | Blackman-Tukey frequency response and noise spectrum estimation |
| `sidFreqBTFDR` | `spafdr` | Blackman-Tukey with frequency-dependent resolution |
| `sidFreqETFE` | `etfe` | Empirical transfer function estimate (FFT ratio) |

**Time-varying analysis:**

| sid function | Replaces | Description |
|---|---|---|
| `sidFreqMap` | `tfestimate`, `mscohere` | Time-varying frequency response map (BT or Welch) |
| `sidSpectrogram` | `spectrogram` | Short-time FFT spectrogram |
| `sidLTVdisc` | — | Discrete LTV state-space identification (COSMIC algorithm) |
| `sidLTVdiscIO` | — | LTV identification from partial observations (Output-COSMIC) |
| `sidLTIfreqIO` | — | LTI realization from I/O frequency response (Ho-Kalman) |
| `sidLTVStateEst` | — | Batch LTV state estimation (RTS smoother) |
| `sidModelOrder` | — | Model order estimation via Hankel SVD |
| `sidLTVdiscTune` | — | Regularization tuning via validation loss or L-curve |
| `sidLTVdiscFrozen` | — | Frozen transfer function G(w,k) with uncertainty propagation |

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
