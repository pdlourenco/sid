# sid — Open-Source System Identification Toolbox

![MATLAB/Octave Tests](https://github.com/pdlourenco/sid-matlab/actions/workflows/tests.yml/badge.svg)
![MATLAB/Octave Lint](https://github.com/pdlourenco/sid-matlab/actions/workflows/lint.yml/badge.svg)
![Python Tests](https://github.com/pdlourenco/sid-matlab/actions/workflows/python-tests.yml/badge.svg)
![Python Lint](https://github.com/pdlourenco/sid-matlab/actions/workflows/python-lint.yml/badge.svg)
![Cross-Language Validation](https://github.com/pdlourenco/sid-matlab/actions/workflows/cross-validate.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**sid** is a free, open-source toolbox for **system identification** — covering
both **non-parametric frequency response estimation** and **time-varying
state-space identification**. All implementations share a single mathematical
specification and cross-language reference test vectors to ensure numerical
consistency.

## Features

- **Blackman-Tukey spectral analysis** — frequency response and noise spectrum estimation with configurable window size
- **Frequency-dependent resolution** — vary the smoothing bandwidth across the frequency axis
- **Empirical transfer function estimate** — maximum resolution via FFT ratio, with optional smoothing
- **Time-varying frequency maps and spectrograms** — sliding-window analysis for non-stationary signals
- **LTV state-space identification (COSMIC)** — identify time-varying A(k), B(k) with O(N) complexity, automatic regularization, and Bayesian uncertainty
- **Partial-observation identification (Output-COSMIC)** — identify dynamics from output-only measurements
- **Multi-trajectory support** — ensemble averaging for frequency estimates; pooled least-squares for state-space
- **Asymptotic uncertainty estimates** — confidence bands for all estimation functions
- **Analysis and validation** — detrending, residual diagnostics, and model comparison
- **SISO, MIMO, and time-series modes** — unified API across all estimation functions

## Language Implementations

| Language | Status | Directory | README | Requirements |
|----------|--------|-----------|--------|--------------|
| MATLAB/Octave | Stable | [`matlab/`](matlab/) | [README](matlab/README.md) | MATLAB R2016b+ or GNU Octave 8.0+ |
| Python | Planned | [`python/`](python/) | [README](python/README.md) | TBD |
| Julia | Planned | [`julia/`](julia/) | [README](julia/README.md) | TBD |

See each language's README for installation, quick start, function reference, and compatibility details.

<details>
<summary><strong>Tip:</strong> Only need one language? Use a sparse checkout.</summary>

```bash
git clone --no-checkout https://github.com/pdlourenco/sid-matlab.git sid
cd sid
git sparse-checkout init --cone
git sparse-checkout set spec testdata matlab   # replace 'matlab' with your language
git checkout
```

This downloads only the shared specification, test data, and your chosen implementation.
</details>

## How It Works

**Frequency-domain path.** The core spectral estimators use the **Blackman-Tukey method**: compute biased cross-covariances between input and output, apply a Hann lag window, then transform via FFT. The transfer function is the cross-spectrum / input auto-spectrum ratio; asymptotic variance formulas (Ljung, 1999) provide per-frequency uncertainty. When multiple trajectories are provided, covariances are ensemble-averaged before forming the ratio, reducing variance by a factor of L without sacrificing frequency resolution.

**State-space path.** The **COSMIC algorithm** (Carvalho et al., 2022) identifies discrete-time LTV models x(k+1) = A(k)x(k) + B(k)u(k) by solving a block-tridiagonal regularized least-squares problem in O(N) time. Multiple trajectories — including variable-length sequences — are pooled into the data matrices. When only outputs are observed, **Output-COSMIC** alternates between state estimation (RTS smoother) and dynamics identification, converging to a joint optimum. Bayesian uncertainty quantification propagates through to frozen transfer functions G(w,k) for direct comparison with non-parametric frequency estimates.

See [SPEC.md](spec/SPEC.md) for the full mathematical derivation.

## Documentation

- [**SPEC.md**](spec/SPEC.md) — Full algorithm specification with mathematical derivations
- [**Roadmap**](docs/roadmap.md) — Development phases and planned features
- [**COSMIC uncertainty derivation**](spec/cosmic/uncertainty_derivation.md) — Bayesian posterior covariance for LTV identification
- [**COSMIC online recursion**](spec/cosmic/online_recursion.md) — Recursive/streaming formulation of the COSMIC algorithm
- [**COSMIC automatic tuning**](spec/cosmic/automatic_tuning.md) — Regularization parameter selection via validation and L-curve
- [**Output-COSMIC**](spec/cosmic/output.md) — LTV identification from partial (output-only) observations

## References

- Ljung, L. (1999). *System Identification: Theory for the User*, 2nd ed. Prentice Hall.
- Blackman, R. B. & Tukey, J. W. (1959). *The Measurement of Power Spectra*. Dover.
- Carvalho, M., Soares, C., Lourenco, P., and Ventura, R. (2022). "COSMIC: fast closed-form identification from large-scale data for LTV systems." [arXiv:2112.04355v2](https://arxiv.org/abs/2112.04355v2)
- Laszkiewicz, P., Carvalho, M., Soares, C., and Lourenco, P. (2025). "The impact of modeling approaches on controlling safety-critical, highly perturbed systems: the case for data-driven models." [arXiv:2509.13531](https://arxiv.org/abs/2509.13531)

## Contributing

Contributions are welcome via issues and pull requests. See [CONTRIBUTING.md](CONTRIBUTING.md) for general guidelines and links to language-specific contributing guides.

## License

MIT License. Copyright (c) 2026 Pedro Lourenço. See [LICENSE](LICENSE).
