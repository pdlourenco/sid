# How sid works

`sid` follows two complementary paths: a **frequency-domain** path for
non-parametric estimation, and a **time-domain** path for parametric
state-space identification. Both paths share a single mathematical
specification and cross-language reference test vectors.

## Frequency-domain path

The core spectral estimators use the **Blackman-Tukey method**: compute
biased cross-covariances between input and output, apply a Hann lag
window, then transform via FFT. The transfer function is the
cross-spectrum / input auto-spectrum ratio; asymptotic variance formulas
(Ljung, 1999) provide per-frequency uncertainty. When multiple
trajectories are provided, covariances are ensemble-averaged before
forming the ratio, reducing variance by a factor of *L* without
sacrificing frequency resolution.

See [`freq_bt`](../api/python/freq_bt.md) /
[`sidFreqBT`](../api/matlab/sidFreqBT.md) for the entry point, and the
[Specification](../spec/index.md) for the full mathematical derivation.

## State-space path

The **COSMIC algorithm** (Carvalho et al., 2022) identifies discrete-time
LTV models *x(k+1) = A(k) x(k) + B(k) u(k)* by solving a block-tridiagonal
regularized least-squares problem in O(N) time. Multiple trajectories —
including variable-length sequences — are pooled into the data matrices.
When only outputs are observed, **Output-COSMIC** alternates between
state estimation (RTS smoother) and dynamics identification, converging
to a joint optimum. Bayesian uncertainty quantification propagates
through to frozen transfer functions G(ω, k) for direct comparison with
non-parametric frequency estimates.

See [`ltv_disc`](../api/python/ltv_disc.md) /
[`sidLTVdisc`](../api/matlab/sidLTVdisc.md) for the entry point, and the
[COSMIC notes](../spec/cosmic/automatic-tuning.md) for derivations.
