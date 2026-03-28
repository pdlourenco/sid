# sid — Open-Source System Identification for MATLAB/Octave

![Tests](https://github.com/pdlourenco/sid-matlab/actions/workflows/tests.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![MATLAB R2016b+](https://img.shields.io/badge/MATLAB-R2016b%2B-orange.svg)
![GNU Octave 8+](https://img.shields.io/badge/GNU_Octave-8%2B-blue.svg)

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=pdlourenco/sid-matlab&file=examples/exampleSISO.m)

**sid** is a free, open-source toolbox for **non-parametric frequency response estimation** — the same class of algorithms found in MATLAB's commercial [System Identification Toolbox](https://www.mathworks.com/products/sysid.html), but implemented in pure MATLAB/Octave code with **zero toolbox dependencies**. It provides drop-in replacements for `spa`, `spafdr`, and `etfe` that you can use without a licence.

Unlike the commercial toolbox, sid exposes **asymptotic uncertainty estimates** and built-in **confidence-band plotting** out of the box, making it straightforward to assess estimation quality. It supports SISO, MIMO, and time-series (output-only) modes, all with a unified result struct and consistent API.

sid is designed from the ground up to run on **GNU Octave** as a first-class target — not just MATLAB with Octave as an afterthought. Every function is tested on both platforms in CI, so Octave users get the same reliability as MATLAB users. If you work in an environment where MATLAB licences are limited or unavailable, sid and Octave give you a fully open-source path to frequency-domain system identification.

## Function Comparison

| sid function | Replaces (System Identification Toolbox) | Description |
|---|---|---|
| `sidFreqBT` | `spa` | Blackman-Tukey frequency response and noise spectrum estimation |
| `sidFreqBTFDR` | `spafdr` | Blackman-Tukey with frequency-dependent resolution |
| `sidFreqETFE` | `etfe` | Empirical transfer function estimate (FFT ratio) |
| `sidBodePlot` | — | Bode diagram with shaded confidence bands |
| `sidSpectrumPlot` | — | Power spectrum plot with shaded confidence bands |
| `sidFreqBTMap` | — | Time-varying frequency response map |
| `sidSpectrogram` | `spectrogram` | Short-time FFT spectrogram |
| `sidMapPlot` | — | Time-frequency color map for `sidFreqBTMap` results |
| `sidSpectrogramPlot` | — | Spectrogram color map |

All estimation functions return uncertainty estimates and squared coherence (SISO), and support both positional and name-value calling conventions.

## Installation

Clone the repository:

```bash
git clone https://github.com/pdlourenco/sid-matlab.git
```

Then add it to your MATLAB or Octave path:

```matlab
run('/path/to/sid-matlab/sidInstall.m')
```

To make the path persistent across sessions, add the line above to your [`startup.m`](https://www.mathworks.com/help/matlab/ref/startup.html) file.

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

See [`examples/`](examples/) for more usage patterns, including MIMO and time-series estimation.

## Compatibility

| Platform | Version | Status |
|---|---|---|
| **MATLAB** | R2016b or later | Tested in CI |
| **GNU Octave** | 8.0 or later | Tested in CI |

No toolboxes are required. The entire codebase uses only core MATLAB/Octave functions (`fft`, `filter`, `conv`, etc.), so it runs anywhere the base language does — including MATLAB Online.

## Documentation

- [**SPEC.md**](SPEC.md) — Full algorithm specification with mathematical derivations (reference: Ljung, *System Identification: Theory for the User*, 1999)
- [**Roadmap**](docs/sid_matlab_roadmap.md) — Development phases and planned features
- [**Examples**](examples/) — Runnable scripts demonstrating typical workflows

## Contributing

Contributions are welcome via issues and pull requests. Please ensure that `tests/runAllTests.m` passes on both MATLAB and Octave before submitting — the CI pipeline checks both platforms automatically.

## License

MIT License. Copyright (c) 2026 Pedro Lourenço. See [LICENSE](LICENSE).
