# Uncertainty quantification

Every estimator in `sid` returns confidence information alongside its
point estimates. The non-parametric estimators use asymptotic variance
formulas (Ljung, 1999) per frequency; the COSMIC state-space identifier
propagates a Bayesian posterior through to *A(k), B(k)* standard
deviations and, optionally, to frozen transfer functions G(ω, k).

## Frequency-domain

Per-frequency variance for the Blackman-Tukey estimator follows from the
Hann window's effective number of independent averages and the
input-output coherence. See [Specification](../spec/index.md) §3.

## State-space (COSMIC)

The COSMIC posterior is derived in
[Uncertainty derivation](../spec/cosmic/uncertainty-derivation.md).
Activate it by passing `Uncertainty=True` (Python) /
`'Uncertainty', true` (MATLAB) to
[`ltv_disc`](../api/python/ltv_disc.md) /
[`sidLTVdisc`](../api/matlab/sidLTVdisc.md), which adds `AStd`, `BStd`,
and a per-step row covariance `P` to the result.
