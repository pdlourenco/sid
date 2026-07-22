# Remediation plan for the repo-wide review findings (issues #134–#145)

**Date:** 2026-07-22
**Scope:** the twelve issues filed from the repo-wide review
([`docs/analyses/2026-07-12-repo-wide-review.md`](../analyses/2026-07-12-repo-wide-review.md),
PR #133): #134–#145.
**Status:** proposed

This plan sequences the fixes so that (1) the validation infrastructure that
proves later fixes is repaired first, (2) spec-convention decisions are made
before the code that depends on them, and (3) no two concurrently open PRs
touch the same file.

---

## 1. Issue inventory

| Issue | Review finding | Severity | Languages | Needs spec edit | Blocks / blocked by |
|---|---|---|---|---|---|
| #134 | C1 | critical — silent wrong states | both | no | — |
| #135 | C2/C3/C5 | critical — crashes on documented shapes | both | no | — |
| #136 | C4 | critical — silent complex→real cast | Python | no | — |
| #137 | S1 | major — posterior covariance mis-scaled | both | yes (§8.9.2, §8.12.2/3) | blocked by #145a |
| #138 | S2 | major — trust region harmful | both | yes (§8.12.4 decision) | after #137 (same file) |
| #139 | S3 | major — model order +1 | both | minor (DC handling) | — |
| #140 | S4 | major — multi-traj fit optimistic | both | no (§15.3 already clear) | — |
| #141 | S5 | major — BTFDR missing guards | both | via #143 | with #143 |
| #142 | S6 | major — Welch scaling conventions | both | yes (§6.5/§6.6/§7.3) | — |
| #143 | S7 | major — σ conventions inconsistent | both | yes (§10.2/§10.3) | pairs with #141 |
| #144 | S8/S9 | major — stabilize blowup; frozen H omitted | both | yes (IO-frozen contract) | — |
| #145 | S10 | major — cross-validation blind spots | infra | no | #145a blocks #137 |

---

## 2. Phases

### Phase 1 — Repair the validation infrastructure (#145a–c)

The cross-language machinery must be trustworthy before it can adjudicate the
numerical fixes in later phases.

1. **#145a — fix the reference generator's double scaling.**
   `testdata/generate_reference.m` divides `S` by `N` before calling
   `sidLTVuncertaintyBackwardPass`, while production passes `S` directly; the
   stored `P` in `reference_cosmic_internals.json` matches neither
   implementation. Remove the extra division, regenerate the JSON, and extend
   `python/tests/test_cross_validation.py` to actually read `D, Xl, S, T, C, P`
   from that file (today only `cost/fidelity/regularization` are read, so the
   whole uncertainty path has zero effective cross-language validation).
   *Note:* regenerate against the **current** (pre-#137) backward pass — this
   phase pins existing behavior; #137 later changes it deliberately.
2. **#145b — wire the four orphaned reference JSONs** (`reference_detrend`,
   `reference_freq_domain_sim`, `reference_ltv_io`, `reference_uncertainty`)
   into `test_cross_validation.py` so the Python ports of `detrend`,
   `freq_domain_sim`, `ltv_disc_io`, and the frequency-domain uncertainty are
   checked against MATLAB vectors.
3. **#145c — make the stored `tolerance` blocks authoritative.** Both consumers
   (`testdata/validate_reference.m`, `test_cross_validation.py`) read the JSON
   tolerances, with one agreed default (rtol 1e-6, atol 0) for fields that lack
   an entry. Remove the per-language hardcoded overrides.

**Exit criteria:** cross-validation reads every stored field of every JSON;
both consumers agree on tolerances; full local + CI pipeline green.

### Phase 2 — Self-contained correctness fixes

No spec decisions needed; independent files; one PR each, any order. Planned
order below is by severity and effort (smallest first).

1. **#136 — complex `u` validation (Python).** Move the `iscomplexobj(u)`
   check in `validate_data.py` ahead of the `asarray` cast (mirroring `y`),
   delete the two dead post-cast checks, add a complex-`u` unit test.
2. **#134 — var-len smoother terminal block (both).** For each trajectory with
   `Nl < N`, substitute the terminal-form diagonal block `HᵀR⁻¹H + Q⁻¹` at
   position `Kl−1` (`python/sid/ltv_state_est.py`, `matlab/sid/sidLTVStateEst.m`).
   Add the var-len smoother-vs-dense-LSQ reference test in both languages and a
   var-len cross-language vector (also listed in #145d).
3. **#135 — crash fixes on documented shapes (both).**
   - `freq_map` storage branch: condition on `ny == 1` only (Python + MATLAB).
   - `cov.py` `L == 1` branch: handle the still-3-D array.
   - `spectrogram.py`: branch on `x.ndim == 3`, not `n_traj > 1`.
   - `map_plot.py` / `spectrum_plot.py`: reduce MIMO responses to the plotted
     channel (`[:, :, 0, 0]` / `[:, 0, 0]`).
   - `sidMapPlot.m`: accept `sidSpectrogram` results per SPEC §7.5; the
     `'spectrum'` branch must read the field spectrogram results actually have.
   - Edge crashes: validate `Lsub` before the NFFT default; error (or warn +
     NaN policy per spec) for `S₁ = 0` degenerate windows.
   - Tests: MIMO `freq_map`, multi-output time-series, 3-D `L = 1`, MIMO
     plotting smoke tests, both languages.
4. **#139 — model-order Hankel lag (both).** Drop the lag-0 IFFT sample
   (align with `lti_freq_io`), lift the `r/2` gap-search cap (spec says
   `k = 1..r−1`), unify the DC-bin extrapolation across the two IFFT consumers
   (pick one, document it in SPEC §8.12.12), document or remove the
   time-series fallback. Add a strictly-proper-with-finite-zero test plant
   (the family current tests are blind to) and a MIMO block-Hankel test.
   **Numerics change** → regenerate `reference_model_order.json`.
5. **#140 — multi-trajectory fit and residual (both).** Compute NRMSE per
   trajectory and average (SPEC §15.3); define and document what
   `predicted`/`measured`/`residual` contain for multi-traj input (proposal:
   per-trajectory arrays, with the scalar `fit` the average). Fix `residual`:
   per-trajectory whiteness with an L-adjusted bound (or per-trajectory
   correlograms), cross-correlate against each trajectory's own input. Add
   multi-traj tests in both languages, including the mirror-image degenerate
   pair from the issue. **Numerics change** → regenerate
   `reference_compare.json`.

### Phase 3 — Convention decisions: spec first, then both ports

Each item starts with a SPEC.md edit (reviewed as its own commit inside the
PR), then MATLAB, then Python, then regenerated references. Recommended
decisions are stated; the spec commit is the place to veto them.

1. **#143 + #141 in one PR — degenerate-excitation handling.**
   - *Decision:* adopt §3.3's sentinel — `σ = Inf`, `Ĝ = NaN` (affected row
     only), plus warning — as the single convention for BT, BTFDR, ETFE, and
     the Welch inner path. Resolve the §10.2 vs §10.3 contradiction by adding
     an absolute excitation floor (input-variance check) so the constant-input
     contract is enforceable; specify NaN-row (not whole-matrix) substitution;
     specify the MATLAB/NumPy `max(NaN, ·)` clamp divergence away by defining
     clamp-before-compare semantics.
   - *Implementation:* one shared helper per language (`private/` /
     `_internal/`) implementing rcond guard, NaN substitution, PSD clamp of
     `Φ̂_v` (§2.7), and warnings; all four estimator paths call it. This is
     what closes #141's BTFDR gaps (missing singularity guard, missing PSD
     clamp, silent `M_k` clamping → warn/error per §10.1, §5.4 short-data
     default).
   - *Tests:* fix the two tests that pin the current violation
     (`test_uncertainty.py::test_zero_coherence_clamp`, MATLAB
     `test_sidUncertainty.m` Test 8); add degenerate-data tests (constant u,
     zero u, collinear MIMO inputs, u = y); add the free oracle
     "BTFDR with constant `R = 2π/M` ≡ BT with window `M`" (§5.1).
2. **#142 — Welch one-sided scaling.**
   - *Decision:* no Nyquist doubling (align §6.5 with §7.3 and MathWorks
     convention); state the BT↔Welch two-sided/one-sided scale relation
     explicitly in §6.6 (do not silently renormalize BT); fix §6.5's variance
     text to `2Φ²/ν` (code is right); restrict `ν ≈ 1.8J` to Hann/50% overlap
     and use the general Welch DoF formula otherwise.
   - *Tests:* quantitative cross-check of the Welch path against
     scipy/tfestimate at machine precision including the Nyquist bin;
     spectrogram Nyquist-bin check.
   - **Numerics change** → regenerate any affected freq_map vectors.
3. **#137 — COSMIC posterior scaling.** (Requires Phase 1 #145a merged.)
   - *Spec:* §8.9.2 Step 1 becomes `S(k) = N·S_scaled(k)` with off-diagonal
     couplings `(N·λ_k)²` (equivalently `P(k) = P_scaled(k)/N`); update
     `spec/cosmic/uncertainty_derivation.md` §3.4/§5.1 (`A = VᵀV + FᵀΥF` with
     the estimator's Υ); reconcile §8.12.2 vs §8.12.3 (state explicitly that
     the user-facing λ has effective weight N·λ on unscaled residuals, and
     make the reported cost J the same objective the algorithm descends).
   - *Code:* `ltv_uncertainty_backward_pass` (both languages); DoF hat-trace
     accordingly; `evaluateFullCost`/Python equivalent in the IO path for the
     cost-weight reconciliation.
   - *Verification:* regenerate `reference_cosmic_internals.json` (now
     correctly generated, thanks to Phase 1); tightened Monte-Carlo
     calibration test at mid-λ (the regime that exposed the bug: reported
     variance 1.86× at λ=1e2 in the reproduction); Python gains its first
     `uncertainty=True` tests.

### Phase 4 — Larger reworks (depend on Phase 3)

1. **#138 — Output-COSMIC trust region rewrite.** After #137 (same files,
   and the reported-cost fix changes what "monotone J" means).
   - Rewrite to the spec's two-level structure: outer μ schedule, inner
     alternating loop run to convergence with its own budget (§8.12.4 /
     `spec/cosmic/output.md` §4.3/§7).
   - Resolve the interpolation target in the spec: `μ·I` per §8.12.4 vs the
     implemented `μ·A₀` (Ho-Kalman init) — pick one, document the rationale.
   - Reject path: revert to best iterate and **terminate** per step 5.
   - Also: non-vacuous convergence test, warning on hitting `max_iter`
     without convergence, fix MATLAB `MaxIter=0` crash, convergence
     denominator behavior for `J < 1`.
   - *Fallback:* if the rewrite stalls, ship a small PR that disables
     `trust_region` (error or warn+ignore) — the issue shows 'off' beats the
     current implementation by two orders of magnitude — and land the rewrite
     separately.
2. **#144 — stabilization + IO-frozen contract.**
   - `_stabilize`: reflect eigenvalues in real Schur form, or guard
     `cond(V)` and fall back with a warning (defective-A repro: 3×3 Jordan
     block → entries of 3e14). Guard `_ho_kalman`'s `√σ` (error on requesting
     order above numerical rank).
   - *Spec decision:* define the IO-frozen contract (§8.12.11): either
     `ltv_disc_frozen` honors `result.h` — returning `H(e^{jω}I−A)⁻¹B` with
     uncertainty propagated through the left-multiplication — or it rejects
     IO results with a clear error. Recommendation: honor `H`; that is what
     the documented §8.12.11 comparison workflow needs.
   - *Tests:* repeated/defective-eigenvalue stabilization; `lti_freq_io`
     accuracy at `H ≠ I` (currently only shape-checked); `ltv_disc_frozen`
     with an IO result.

### Phase 5 — Close the loop (#145d–e)

1. **#145d — missing reference vectors:** multi-trajectory estimation, ETFE
   time-series periodogram, Welch `freq_map`, var-len COSMIC/StateEst,
   auto-λ (L-curve), `sidLTVdiscTune`, vector-R BTFDR, COSMIC
   `AStd/BStd/Σ̂/ν`. Several double as regression pins for Phases 2–4, so any
   that are ready earlier should land with their phase's PR.
2. **#145e — metadata and docs:** generator version/date/seed in every JSON;
   `testdata/README.md` corrections (`function_name` key, Python validator
   already in CI, per-file inventory); CI trigger for
   `testdata/generate_reference.m` path in `tests.yml`.

---

## 3. Working method

- **One branch + PR per numbered item** (Phase 3 items include their spec
  edit as a distinct commit). Suggested execution order:
  `145a-c → 136 → 134 → 135 → 139 → 140 → 143+141 → 142 → 137 → 138 → 144 → 145d-e`.
- **Full local CI before every PR:** MATLAB tests + examples, Python tests
  (incl. cross-validation) + notebooks, `validate_reference.m`, MISS_HIT,
  ruff, both header checks. All verified runnable locally.
- **Reference regeneration:** where a fix legitimately changes numerics
  (#139, #140, #142, #137, #145a), regenerate the affected
  `testdata/*.json` locally via `generate_reference.m` and include them in
  the PR, flagged in the PR description. Local MATLAB is R2024a vs CI's
  R2025a; the post-merge CI regeneration commit is the authoritative one and
  may produce tiny diffs.
- **Issue hygiene:** each PR closes its issue via `Closes #NNN`; issues that
  land in a combined PR (#141/#143) are both referenced. Deviations from
  this plan get a comment on the issue explaining why.

## 4. Dependency graph

```
#145a ──► #137 ──► #138
#143 ═══ #141 (one PR)
(everything else independent; #136/#134/#135/#139/#140 can proceed in
parallel with Phase 1 as they touch disjoint files)
```

## 5. Out of scope

Open issues below #134 (seed-adoption epic #111–#117, coverage issues
#119–#124, #131) are untouched by this plan except where a Phase 2–5 test
gap happens to overlap a coverage issue; overlaps will be noted on those
issues as they close incidentally.
