# Changelog

All notable changes to **sid** are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the project uses
per-language release tags (`vX.Y.Z-matlab`, `vX.Y.Z-python`) that move together.

Numerical equivalence across the MATLAB/Octave and Python ports is a consequence
of each conforming to [`spec/SPEC.md`](spec/SPEC.md); changes below apply to both
ports unless noted.

## [0.2.0] — 2026-07-26

The correctness-and-verification release. It remediates every finding of the
2026-07-12 repo-wide review across both ports, closes the cross-language
verification gaps, and hardens the contract-artifact machinery. Several fixes
change numerical outputs — see **Changed** for the ones to be aware of when
upgrading.

### Changed (behavioural — review before upgrading)

- **COSMIC posterior covariance rescaled** (`sidLTVdisc`, `sidLTVdiscIO`): the
  reported posterior `P(k)` is now `P_scaled / N`, and the smoothness term uses
  the effective weight `Nλ` consistent with the `1/√N` data scaling. The
  previous reconstruction **overstated** `P` (hence `AStd`/`BStd`) by up to a
  factor `N`, worst in the mid-`λ` regime; the reported `Cost` magnitude also
  changes (it now uses the objective the alternation provably decreases). Point
  estimates `A(k)`/`B(k)` are unaffected. (#137, SPEC §8.9/§8.12.2)
- **Output-COSMIC trust region** (`sidLTVdiscIO`, `'TrustRegion'`): rewritten
  from a broken fused loop into a correct two-level μ-homotopy. Previously
  `TrustRegion = 1` left the cost ~two decades *worse* than off; it now helps on
  hard partial-observation cases. When enabled it runs materially more inner
  iterations (per the normative budget `MaxIter × (⌈log₂(1/ε_μ)⌉ + 2)`); `off`
  remains the default and is unchanged. (#138, SPEC §8.12.4)
- **`sidLTVdiscFrozen` of an `sidLTVdiscIO` result** now returns the **output**
  transfer function `H(e^{jω}I − A)⁻¹B` (`p_y × q`) with uncertainty folded
  through `H`, instead of the `n × q` state response (dimensionally wrong for
  `H ≠ I`). Full-state results (`H = I`) are unchanged. (#144, SPEC §8.11.1)
- **Welch one-sided PSD** (`sidFreqMap` Welch path) no longer doubles the Nyquist
  bin. (#142, SPEC §6.5)

### Fixed

- **`sidLTIfreqIO` stabilization** no longer blows up (entries ~1e14) on
  defective / integrator-chain dynamics: eigenvalues are reflected/clamped in
  real Schur form rather than by eigenvector inversion. Requesting an order
  above the Hankel numerical rank now raises `sid:orderExceedsRank` instead of
  propagating `inf`/`NaN`. (#144, SPEC §8.13.1)
- **Degenerate excitation** is handled uniformly across `sidFreqBT`,
  `sidFreqBTFDR`, `sidFreqETFE`, and the Welch path: the MIMO singular-solve
  crash is fixed (rcond-guarded → `NaN` + `σ_G = Inf` sentinel + one warning),
  with an absolute input-excitation floor and NaN-row propagation. (#141/#143,
  SPEC §3.3/§5.2/§10.3)
- Variable-length smoother terminal block (#134); crashes on documented MIMO /
  multi-output / 3-D / plotting input shapes (#135); model-order Hankel lag-0
  off-by-one (#139); multi-trajectory fit/residual (#140); complex-input
  validation (#136).
- `sidLTVdiscTune` validation mode accepts 1-D / 2-D single-trajectory data
  rather than requiring 3-D. (#189, SPEC §8.11)
- Cross-validation infrastructure: mis-generated reference `P`, four orphaned
  reference vectors, and divergent tolerance handling. (#145a–c)

### Added

- **Numerical diagnostics** (stable warning identifiers, part of the contract):
  `sid:notConverged` (base alternation hit the iteration cap), `sid:stabilized`
  (stabilization moved eigenvalues), `sid:orderExceedsRank`.
- **Generator-provenance metadata**: every `testdata/reference_*.json` carries a
  `provenance` block (generator, source-commit SHA/date), gated by both
  validators. (#172, ADR-0002)

### Testing & verification

- Cross-language reference vectors grew 22 → 32; **every `Verified by: none`
  cross-vector gap is discharged**. New coverage in both ports for the L-curve
  auto-λ corner selection and ill-conditioning warning (#120), the Welch
  inner-path (#122), the plotting confidence-band formulas (#123), the
  frequency-domain simulation helper and `sidLTIfreqIO` port symmetry (#124),
  and multi-trajectory / var-len / tuning / frozen-of-IO / `H ≠ I` paths (#145d).
- `Verified by:` verification annotations across `spec/SPEC.md` §1–§15 (#113).
- Contract-artifact hardening (ADR-0002): absolute-tolerance floors, authoritative
  stored tolerances, no-orphan-artifacts, payload-only-by-regeneration, and
  structural gates — all now enforced.

### Project / governance

- Adopted the disciplined-project-seed governance layer: `CLAUDE.md`, PR
  template, ADR infrastructure (`docs/decisions/`), `REVIEW_CONTEXT.md` +
  pre-push self-review, `docs/DESIGN.md`, and CI hardening (least-privilege
  workflows, LF `.gitattributes`, local-CI runner). (#117)

## [0.1.0] — 2026-04-11

Initial release: frequency-domain (Blackman-Tukey) and COSMIC LTV state-space
identification for MATLAB/Octave and Python, built on a shared specification
with cross-language reference vectors. No MathWorks toolbox dependencies.
