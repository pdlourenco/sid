# sid — MATLAB/Octave Implementation

![Tests](https://github.com/pdlourenco/sid/actions/workflows/tests.yml/badge.svg)
![Lint](https://github.com/pdlourenco/sid/actions/workflows/lint.yml/badge.svg)
![MATLAB R2016b+](https://img.shields.io/badge/MATLAB-R2016b%2B-orange.svg)
![GNU Octave 8+](https://img.shields.io/badge/GNU_Octave-8%2B-blue.svg)

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=pdlourenco/sid&file=matlab/examples/exampleSISO.m)

MATLAB/Octave implementation of the sid system identification toolbox.
Pure MATLAB/Octave code with **zero toolbox dependencies**. Provides
drop-in replacements for `spa`, `spafdr`, and `etfe` that you can use
without a System Identification Toolbox licence, and implements the
**COSMIC** algorithm for identifying linear time-varying (LTV) state-space
models from data.

sid is designed from the ground up to run on **GNU Octave** as a first-class
target — not just MATLAB with Octave as an afterthought. Every function is
tested on both platforms in CI.

## Installation

Clone the repository:

```bash
git clone https://github.com/pdlourenco/sid.git
```

Then add it to your MATLAB or Octave path:

```matlab
run('/path/to/sid/matlab/sidInstall.m')
```

To make the path persistent across sessions, add the line above to your
[`startup.m`](https://www.mathworks.com/help/matlab/ref/startup.html) file
(or [`.octaverc`](https://docs.octave.org/latest/Startup-Files.html) for
Octave). No `pkg install` is needed — sid is a plain directory of `.m` files
that works on both platforms.

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

See the [examples guide](examples/README.md) for more usage patterns, including
ETFE, frequency-dependent resolution, MIMO, time-varying maps, spectrograms,
LTV state-space identification, multi-trajectory ensemble averaging, and
partial-observation Output-COSMIC.

## Function Reference

### Frequency-domain estimation

| sid function | Replaces | Description |
|---|---|---|
| `sidFreqBT` | `spa` | Blackman-Tukey frequency response and noise spectrum estimation |
| `sidFreqBTFDR` | `spafdr` | Blackman-Tukey with frequency-dependent resolution |
| `sidFreqETFE` | `etfe` | Empirical transfer function estimate (FFT ratio) |

### Time-frequency analysis

| sid function | Replaces | Description |
|---|---|---|
| `sidFreqMap` | `tfestimate`, `mscohere` | Time-varying frequency response map (BT or Welch) |
| `sidSpectrogram` | `spectrogram` | Short-time FFT spectrogram |

### State-space identification (COSMIC)

| sid function | Replaces | Description |
|---|---|---|
| `sidLTVdisc` | — | Discrete LTV state-space identification (COSMIC algorithm) |
| `sidLTVdiscIO` | — | LTV identification from partial observations (Output-COSMIC) |
| `sidLTVdiscTune` | — | Regularization tuning via validation loss, L-curve, or frequency consistency |
| `sidLTVdiscFrozen` | — | Frozen transfer function G(w,k) with uncertainty propagation |
| `sidLTIfreqIO` | — | LTI realization from I/O frequency response (Ho-Kalman) |
| `sidLTVStateEst` | — | Batch LTV state estimation (RTS smoother) |
| `sidModelOrder` | — | Model order estimation via Hankel SVD |

### Plotting

| sid function | Replaces | Description |
|---|---|---|
| `sidBodePlot` | — | Bode diagram with shaded confidence bands |
| `sidSpectrumPlot` | — | Power spectrum plot with shaded confidence bands |
| `sidMapPlot` | — | Time-frequency color map for `sidFreqMap` results |
| `sidSpectrogramPlot` | — | Spectrogram color map |

### Analysis and validation

| sid function | Replaces | Description |
|---|---|---|
| `sidDetrend` | `detrend` | Remove mean, linear, or polynomial trends from signals |
| `sidResidual` | `resid` | Residual analysis with auto/cross-correlation diagnostics |
| `sidCompare` | `compare` | Simulate model output and compute fit percentage |

All estimation functions support both positional and name-value calling conventions.

## Result Types

Every estimation function returns a struct with documented fields. Run
`help sidResultTypes` for the complete field-by-field reference, or see
the table below for a quick overview.

| Result type | Produced by | Key fields |
|---|---|---|
| FreqResult | `sidFreqBT`, `sidFreqBTFDR`, `sidFreqETFE` | `.Response`, `.NoiseSpectrum`, `.Coherence`, `.ResponseStd` |
| FreqMapResult | `sidFreqMap` | `.Time`, `.Response`, `.NoiseSpectrum`, `.Coherence` |
| SpectrogramResult | `sidSpectrogram` | `.Time`, `.Frequency`, `.Power`, `.PowerDB`, `.Complex` |
| LTVResult | `sidLTVdisc`, `sidLTVdiscTune` | `.A`, `.B`, `.Lambda`, `.Cost` (+ `.AStd`, `.BStd`, `.P` with uncertainty) |
| LTVIOResult | `sidLTVdiscIO` | `.A`, `.B`, `.X`, `.H`, `.R`, `.Cost`, `.Iterations` |
| FrozenResult | `sidLTVdiscFrozen` | `.Response`, `.ResponseStd`, `.TimeSteps` |
| CompareResult | `sidCompare` | `.Predicted`, `.Measured`, `.Fit`, `.Residual` |
| ResidualResult | `sidResidual` | `.Residual`, `.AutoCorr`, `.CrossCorr`, `.WhitenessPass` |

All result structs include `.Method` (identifier string) and metadata fields
such as `.SampleTime`, `.DataLength`, and `.NumTrajectories`.

## Compatibility

| Platform | Version | Status |
|---|---|---|
| **MATLAB** | R2016b or later | Tested in CI |
| **GNU Octave** | 8.0 or later | Tested in CI |

No toolboxes are required. The entire codebase uses only core MATLAB/Octave
functions (`fft`, `filter`, `conv`, etc.), so it runs anywhere the base
language does — including MATLAB Online.

## Testing

```matlab
run('matlab/tests/runAllTests.m')
run('matlab/examples/runAllExamples.m')
```

## Documentation

- [Algorithm specification](../spec/SPEC.md)
- [Examples guide](examples/README.md)
- [Development roadmap](../docs/roadmap.md)
- [Contributing guide](CONTRIBUTING.md)
