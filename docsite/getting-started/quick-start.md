# Quick start

Two equivalent identifications of a first-order system with output noise —
one in each language. Both produce numerically identical results to within
floating-point tolerance.

## Python

```python
import numpy as np
import sid
from scipy.signal import lfilter

N, Ts = 1000, 0.01
rng = np.random.default_rng(42)
u = rng.standard_normal(N)
y = lfilter([1], [1, -0.9], u) + 0.1 * rng.standard_normal(N)

result = sid.freq_bt(y, u, sample_time=Ts)
sid.bode_plot(result)
sid.spectrum_plot(result)
```

Time-series mode (output spectrum only):

```python
result_ts = sid.freq_bt(y, None)
sid.spectrum_plot(result_ts)
```

## MATLAB / Octave

```matlab
N = 1000; Ts = 0.01;
u = randn(N, 1);
y = filter([1], [1 -0.9], u) + 0.1 * randn(N, 1);

result = sidFreqBT(y, u, 'SampleTime', Ts);
sidBodePlot(result);
sidSpectrumPlot(result);
```

Time-series mode:

```matlab
result_ts = sidFreqBT(y, []);
sidSpectrumPlot(result_ts);
```

## Where next

- [Function index](../api/index.md) — every public function in both languages
- [Examples](../examples/index.md) — full notebooks and scripts
- [Specification](../spec/index.md) — math behind every estimator
