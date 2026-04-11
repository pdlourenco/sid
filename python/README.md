# sid — Python Implementation

![Python Tests](https://github.com/pdlourenco/sid/actions/workflows/python-tests.yml/badge.svg)
![Python Lint](https://github.com/pdlourenco/sid/actions/workflows/python-lint.yml/badge.svg)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pdlourenco/sid/main?labpath=python%2Fexamples)

Python implementation of the sid system identification toolbox. Numerically
identical to the MATLAB/Octave reference implementation — every public function
produces the same output to within floating-point tolerance.

Requires **Python 3.10+**, **NumPy**, and **SciPy**. Matplotlib is optional
(needed only for plotting functions and examples).

## Installation

From the repository root:

```bash
pip install -e ./python
```

With plotting support:

```bash
pip install -e "./python[plot]"
```

For development (tests + linting + notebook validation):

```bash
pip install -e "./python[dev]"
```

## Quick Start

```python
import numpy as np
import sid

# Generate example data: first-order system with noise
N = 1000
Ts = 0.01
rng = np.random.default_rng(42)
u = rng.standard_normal(N)
y = np.convolve(u, [1.0], mode="full")[:N]  # placeholder; use scipy.signal.lfilter
from scipy.signal import lfilter
y = lfilter([1], [1, -0.9], u) + 0.1 * rng.standard_normal(N)

# Estimate frequency response (Blackman-Tukey)
result = sid.freq_bt(y, u, sample_time=Ts)

# Plot Bode diagram and noise spectrum with confidence bands
sid.bode_plot(result)
sid.spectrum_plot(result)
```

Time-series mode (output spectrum only, no input signal):

```python
result_ts = sid.freq_bt(y, None)
sid.spectrum_plot(result_ts)
```

See the [examples](examples/README.md) for more usage patterns, including
ETFE, frequency-dependent resolution, MIMO, time-varying maps, spectrograms,
LTV state-space identification, multi-trajectory ensemble averaging, and
partial-observation Output-COSMIC.

## Function Reference

### Frequency-domain estimation

| Python function | MATLAB equivalent | Description |
|---|---|---|
| `sid.freq_bt` | `sidFreqBT` / `spa` | Blackman-Tukey frequency response and noise spectrum |
| `sid.freq_btfdr` | `sidFreqBTFDR` / `spafdr` | Blackman-Tukey with frequency-dependent resolution |
| `sid.freq_etfe` | `sidFreqETFE` / `etfe` | Empirical transfer function estimate (FFT ratio) |

### Time-frequency analysis

| Python function | MATLAB equivalent | Description |
|---|---|---|
| `sid.freq_map` | `sidFreqMap` / `tfestimate` | Time-varying frequency response map (BT or Welch) |
| `sid.spectrogram` | `sidSpectrogram` / `spectrogram` | Short-time FFT spectrogram |

### State-space identification (COSMIC)

| Python function | MATLAB equivalent | Description |
|---|---|---|
| `sid.ltv_disc` | `sidLTVdisc` | Discrete LTV state-space identification (COSMIC algorithm) |
| `sid.ltv_disc_io` | `sidLTVdiscIO` | LTV identification from partial observations (Output-COSMIC) |
| `sid.ltv_disc_tune` | `sidLTVdiscTune` | Regularization tuning via validation loss or frequency consistency |
| `sid.ltv_disc_frozen` | `sidLTVdiscFrozen` | Frozen transfer function G(w,k) with uncertainty propagation |
| `sid.lti_freq_io` | `sidLTIfreqIO` | LTI realization from I/O frequency response (Ho-Kalman) |
| `sid.ltv_state_est` | `sidLTVStateEst` | Batch LTV state estimation (RTS smoother) |
| `sid.model_order` | `sidModelOrder` | Model order estimation via Hankel SVD |

### Plotting

| Python function | MATLAB equivalent | Description |
|---|---|---|
| `sid.bode_plot` | `sidBodePlot` | Bode diagram with shaded confidence bands |
| `sid.spectrum_plot` | `sidSpectrumPlot` | Power spectrum plot with confidence bands |
| `sid.map_plot` | `sidMapPlot` | Time-frequency color map for `freq_map` results |
| `sid.spectrogram_plot` | `sidSpectrogramPlot` | Spectrogram color map |

### Analysis and preprocessing

| Python function | MATLAB equivalent | Description |
|---|---|---|
| `sid.detrend` | `sidDetrend` | Remove mean, linear, or polynomial trends from signals |
| `sid.residual` | `sidResidual` | Residual analysis with auto/cross-correlation diagnostics |
| `sid.compare` | `sidCompare` | Simulate model output and compute fit percentage |

All estimation functions return frozen dataclasses with tab-completable fields.
Options are passed as keyword arguments (MATLAB name-value pairs become
`snake_case` kwargs).

### Result types

Every public function returns a typed, frozen dataclass from `sid._results`:

| Result type | Returned by | Key fields |
|---|---|---|
| `FreqResult` | `freq_bt`, `freq_etfe`, `freq_btfdr` | `response`, `noise_spectrum`, `coherence`, `frequency` |
| `FreqMapResult` | `freq_map` | `response`, `noise_spectrum`, `coherence`, `time`, `frequency` |
| `SpectrogramResult` | `spectrogram` | `power`, `power_db`, `complex_stft`, `time`, `frequency` |
| `LTVResult` | `ltv_disc` | `a`, `b`, `a_std`, `b_std`, `lambda_`, `cost` |
| `LTVIOResult` | `ltv_disc_io` | `a`, `b`, `x`, `cost`, `iterations` |
| `FrozenResult` | `ltv_disc_frozen` | `response`, `response_std`, `time_steps`, `frequency` |
| `ResidualResult` | `residual` | `residual`, `auto_corr`, `whiteness_pass`, `independence_pass` |
| `CompareResult` | `compare` | `predicted`, `measured`, `fit`, `residual` |

All result objects are immutable (frozen) and support tab-completion in
IPython/Jupyter. Access fields with dot notation: `result.response`,
`result.fit`, etc.

## Examples

Examples are provided as **Jupyter notebooks** in [`python/examples/`](examples/README.md).
Each notebook mirrors a MATLAB example script, combining narrative, code, and
inline plots:

| Notebook | Description |
|---|---|
| `example_siso.ipynb` | Basic SISO frequency response estimation |
| `example_etfe.ipynb` | Empirical transfer function estimate |
| `example_freq_dep_res.ipynb` | Frequency-dependent resolution |
| `example_coherence.ipynb` | Coherence analysis and signal quality |
| `example_method_comparison.ipynb` | Comparing BT, BTFDR, and ETFE |
| `example_mimo.ipynb` | Multi-input multi-output systems |
| `example_freq_map.ipynb` | Time-varying frequency response maps |
| `example_spectrogram.ipynb` | Short-time FFT spectrogram |
| `example_ltv_disc.ipynb` | LTV state-space identification (COSMIC) |
| `example_multi_trajectory.ipynb` | Multi-trajectory ensemble averaging |
| `example_output_cosmic.ipynb` | Partial-observation Output-COSMIC |

Launch any notebook locally:

```bash
jupyter notebook python/examples/example_siso.ipynb
```

Or try them online — click the **Binder** badge above to launch a ready-to-run
Jupyter environment in the browser (no installation required).

## Compatibility

| Dependency | Minimum version |
|---|---|
| Python | 3.10 |
| NumPy | 1.22 |
| SciPy | 1.8 |
| Matplotlib | 3.5 (optional) |

## Testing

```bash
# Unit tests
pytest python/tests/ -v

# Notebook execution (requires nbmake)
pytest --nbmake python/examples/

# Cross-language validation (requires testdata/*.json from MATLAB CI)
pytest python/tests/test_cross_validation.py -v
```

## Documentation

- [Algorithm specification](../spec/SPEC.md)
- [Examples](examples/README.md)
- [Python development roadmap](../docs/roadmap_python.md)
- [Contributing guide](CONTRIBUTING.md)
- [MATLAB reference implementation](../matlab/README.md)
