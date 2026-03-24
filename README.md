# sid — Open-Source System Identification for MATLAB/Octave

Non-parametric frequency response estimation, with no toolbox dependencies.

## Quick Start

```matlab
% Add to path
run('/path/to/sid/sidInstall.m')

% Generate example data: first-order system with noise
N = 1000;
u = randn(N, 1);
y = filter([1], [1 -0.9], u) + 0.1 * randn(N, 1);

% Estimate frequency response
result = sidFreqBT(y, u);

% Plot
sidBodePlot(result);
sidSpectrumPlot(result);
```

## Functions

| Function | Description | Status |
|----------|-------------|--------|
| `sidFreqBT` | Frequency response via Blackman-Tukey (replaces `spa`) | Implemented |
| `sidFreqBTFDR` | Blackman-Tukey with frequency-dependent resolution (replaces `spafdr`) | Stub |
| `sidFreqETFE` | Empirical transfer function estimate (replaces `etfe`) | Stub |
| `sidBodePlot` | Bode diagram with confidence bands | Implemented |
| `sidSpectrumPlot` | Power spectrum plot with confidence bands | Implemented |

## Requirements

- MATLAB R2016b or later, **or** GNU Octave 8+
- No toolbox dependencies

## Algorithm

See [SPEC.md](SPEC.md) for the full algorithm specification.

## License

MIT. See [LICENSE](LICENSE).
