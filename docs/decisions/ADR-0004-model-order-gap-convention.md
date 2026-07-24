# ADR-0004: Model-order gap search — machine-eps floor, cliff exclusion, retained cap

- **Status:** Proposed — _2026-07-24_
- **Deciders:** Pedro Lourenço (maintainer), implementer session
- **Related:** issue #160 (tracks this cluster), #139 / PR #157 (lag-0 Hankel fix that precedes it), [SPEC §8.12.12](../../spec/SPEC.md), [ADR-0001](ADR-0001-spec-is-the-contract.md)

## Context

`sidModelOrder` picks the state dimension as the largest gap in the Hankel singular
values. SPEC §8.12.12 nominally specified `n = argmax_k σ_k/σ_{k+1}` over
`k = 1..r-1`, but that rule is un-runnable on an *estimated* frequency response:
the noisy singular-value tail has no exact zeros, so an unguarded argmax wanders
into the tail. Both ports therefore grew two undocumented guards — a noise floor
and a `floor(m/2)` cap — and they drifted apart in the process. Issue #160 split
this cluster out of #139 and proposed a direction: *add a data-aware relative
noise floor, then lift the r/2 cap and reconcile the spec.*

Two facts, established numerically on the post-#157 tree (see
[analysis 2026-07-24](../analyses/2026-07-24-model-order-gap-search.md)), reshaped
that direction:

1. **A real cross-language divergence.** The floor marks the last resolvable
   singular value `σ_L`. MATLAB used the **1-based** index `lastSig = L` as the
   ratio bound and so *included* the resolvable→floor cliff `σ_L/σ_{L+1}`
   (a resolvable value over a numerical zero — an enormous, meaningless ratio);
   Python used the **0-based** `last_sig = L-1` and *excluded* it. On the exact
   finite-zero plant `G(z)=(z+0.5)/(z²−1.2728z+0.81)` this yields Python `n=2`
   (correct) vs MATLAB `n=3` on an identical `Σ`.
2. **The reconstruction fix removed the need for a data-aware floor.** The lag-1
   IFFT reconstruction merged in #139/#157 suppresses finite-grid/DC artifacts to
   near-`ε` magnitudes. On that tree a machine-eps floor already recovers the true
   order on exact plants (finite-zero → 2, third-order → 3), and on noisy BT
   estimates a machine-eps, a `1e-8`, and a `1e-3` floor all agree. A `1e-3`
   relative floor is therefore *never* helpful and, on near-exact data, actively
   harmful: it cuts the small-but-real `σ` and, combined with cliff exclusion,
   collapses the finite-zero plant to `n=1`.

The `floor(m/2)` cap is separately load-bearing: noise-dominated data (e.g. output
unrelated to input) produces a slowly decaying spectrum with *no* cliff, every
value stays above any sane floor, and an uncapped search drifts to `n ≈ m`
(white-noise BT estimate → `n=40` of 42).

## Decision

The default gap search is defined on singular-value **magnitudes**, not array
indices: `L =` number of resolvable values above the machine-eps floor
`σ_1·√m·ε`; search `k = 1 … K` with `K = min(L−1, floor(m/2))`; `n = argmax σ_k/σ_{k+1}`.
The `L−1` bound **excludes** the resolvable→floor cliff ratio, the machine-eps
floor is kept (no data-aware floor, no new option), and the `floor(m/2)` cap is
**retained, not lifted**. MATLAB and Python compute from the same count `L`, so
they return identical `n` for identical `Σ`.

## Consequences

- **MATLAB output changes** on low-order plants: the finite-zero reference goes
  `3 → 2`, matching Python and the true order. This is a deliberate, recorded
  behavior change, pinned by a cross-language reference vector on a strictly-proper
  finite-zero plant (the plant #160 recommends over the current biproper
  test-case 7).
- **Python is essentially unchanged** — its 0-based `last_sig` already equals
  `L-1`; only the comment/framing updates. The behavioral fix is a one-line MATLAB
  correction (`lastSig → lastSig-1`, i.e. count-based `L-1`).
- **No new public surface.** No `'NoiseFloor'` option; `'Threshold'` remains the
  escape hatch for anyone needing an explicit magnitude cutoff or rank-count
  semantics.
- **The `floor(m/2)` cap becomes contract**, documented in SPEC with its
  rationale, rather than an undocumented implementation detail. Detectable order is
  bounded to `floor(m/2)` — acceptable for the COSMIC use case and the price of
  robustness on no-cliff data.
- **Issue #160's "then lift the cap" direction is closed as superseded**, with the
  reason recorded here rather than lost in a closed thread.
- The in-code `issue #160` NOTE markers (both ports) are removed by the fixing PR,
  per ADR-0003.

## Alternatives considered

- **Data-aware relative floor (`ρ = 1e-3`, exposed as `'NoiseFloor'`).** The
  original #160 direction, and the choice initially approved for this ADR *before*
  re-measuring on the post-#157 tree. Rejected because the reconstruction fix
  shrank artifacts below any useful relative floor: `1e-3` changes nothing on noisy
  estimates and, on exact low-order plants, removes a genuine mode and collapses
  the order (finite-zero → `n=1`, third-order → `n=2`). A larger floor cannot
  distinguish a weak real mode from an artifact once artifacts are `O(ε)`.
- **Lift the `floor(m/2)` cap.** The second half of #160's direction. Rejected
  because the floor and the cap guard *orthogonal* regimes: the floor handles the
  sharp-cliff (low-order) case, the cap handles the no-cliff (noise-dominated)
  case. With the cap lifted, a white-noise BT estimate returns `n=40` of 42 — no
  floor engages because there is no resolvable/floor boundary to find.
- **Adopt MATLAB's include-cliff behavior as the cross-language reference.** Would
  align the ports by making Python *include* the cliff. Rejected because counting
  resolvable values (`n ≈ L`) is rank-count semantics, which `'Threshold'` already
  provides; it also returns the wrong order (`3`, not `2`) on the finite-zero
  plant, since the cliff `σ_L/σ_{L+1}` divides a real value by a numerical zero and
  always dominates the argmax.
- **Leave the cap and floor undocumented (spec keeps `k = 1..r-1`).** Rejected:
  the spec would keep describing an algorithm neither port runs, and the
  cross-language divergence would have no authoritative resolution — exactly the
  drift ADR-0001 exists to prevent.
