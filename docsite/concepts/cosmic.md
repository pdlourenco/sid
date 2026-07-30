# COSMIC — LTV state-space identification

COSMIC (Carvalho et al., 2022) is the regularized least-squares algorithm
behind every `ltv_disc*` function. It identifies discrete-time LTV state
space matrices *A(k), B(k)* from state and input trajectories — or, in
its output-only variant, from output measurements alone — in O(N) time
via a block-tridiagonal solver.

## Entry points

- [`ltv_disc`](../api/python/ltv_disc.md) / [`sidLTVdisc`](../api/matlab/sidLTVdisc.md) — full-state observation
- [`ltv_disc_io`](../api/python/ltv_disc_io.md) / [`sidLTVdiscIO`](../api/matlab/sidLTVdiscIO.md) — partial-observation (Output-COSMIC)
- [`ltv_disc_tune`](../api/python/ltv_disc_tune.md) / [`sidLTVdiscTune`](../api/matlab/sidLTVdiscTune.md) — automatic regularization tuning
- [`ltv_disc_frozen`](../api/python/ltv_disc_frozen.md) / [`sidLTVdiscFrozen`](../api/matlab/sidLTVdiscFrozen.md) — frozen transfer functions G(ω, k)

## Derivations

- [Automatic tuning](../spec/cosmic/automatic-tuning.md)
- [Online recursion](../spec/cosmic/online-recursion.md)
- [Output-COSMIC](../spec/cosmic/output.md)
- [Uncertainty derivation](../spec/cosmic/uncertainty-derivation.md)
