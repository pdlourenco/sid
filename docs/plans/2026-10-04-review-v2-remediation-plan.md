# Remediation plan for the 2026-10-04 repo-wide review (v2)

**Date:** 2026-10-04
**Source:** [`docs/analyses/2026-10-04-repo-wide-review-v2.md`](../analyses/2026-10-04-repo-wide-review-v2.md) (anchored at `1785dd3`). Finding ids below (`§3.2`, `F1`, `U1`, `C3`, …) refer to that document.
**Status:** proposed — every item here is a recommendation for the maintainer to accept, reorder or drop (`CLAUDE.md` §4: recommend, don't decide). Items marked **[decision]** change a contract or lock in a trade-off and need an explicit go-ahead (and, where non-obvious, an ADR) before implementation.

The ordering principle is the same as the July plan's, with one addition learned from this cycle: **repair the gates before adding the tests that rely on them** (Phase A), because the July remediation added 10 cross-vectors and ~120 tests into an infrastructure that, as §4.1 shows, could not have rejected several of the errors it was meant to catch.

---

## 0. Issue inventory (proposed)

One issue per row; the PR column is the suggested grouping (one branch + PR per group, disjoint files per the July working method).

| Id | Finding | Severity | Ports | Spec edit | PR group |
|---|---|---|---|---|---|
| A1 | Required checks: tests + cross-validation + docs (§4.2) | gate | CI | ADR-0005 amendment | A-ci |
| A2 | Regeneration `[skip ci]` → validated before landing; provenance ancestry/`unknown` rejected; reproducibility policy (§4.1 F1, F5) | gate | CI, testdata | ADR-0002 amendment | A-ci |
| A3 | Octave validator NaN-blind; `Inf`/`NaN` on-disk ambiguity (F3) | gate | testdata | ADR-0002 amendment **[decision]** | A-vectors |
| A4 | Python consumer structural field loop; complex tolerance key parity; SV truncation (F2, F4, F8a) | gate | python | — | A-vectors |
| A5 | Generator-side atol-floor gate; default-tolerance audit (rule 1, F8b) | gate | testdata | — | A-vectors |
| A6 | Re-derive the three 1 % vectors: converged IO case, var-len at 1e-6, frozen-of-IO from the converged case (F6, F7, §3.1) | gate | testdata | — | A-vectors |
| A7 | Spec-derived oracle job over `testdata/*.json` (§4.1 "both ports wrong") | gate | CI, testdata | — | A-oracle |
| A8 | Missing vector / missing JSON = failure, not skip (§3.6) | gate | python, testdata, CI | — | A-ci |
| A9 | Minimum-dependency leg; pin Octave channel; `filterwarnings = error`; `local-ci` runs the docs API gate (§4.2) | gate | CI | — | A-ci |
| B1 | `test_sidModelOrder.m` Test 5 bound restored (§3.2) | test | matlab | — | B-relax |
| B2 | MC calibration: campaign in a scheduled job; gate band re-centred; MATLAB Test 13 made a verifier (§3.3) | test | both | — | B-relax |
| B3 | Tautological / smoke tests rewritten or deleted (§3.4 table) | test | both | — | B-relax |
| B4 | Deterministic checks tightened to 1e-10-class (§3.5) | test | both | — | B-relax |
| B5 | Silent skips → explicit skips with a tally; `warning('off')` restored via `onCleanup`; runner workspace isolation (§3.6, §4.3) | test | matlab, python | — | B-relax |
| B6 | `Verified by:` over-claims corrected (§3.7) | spec | spec | yes | B-relax |
| C1 | Free oracles: MIMO static-gain plant, L-copies ensemble identity, OLS at N = 2, FFT = direct, `q ≥ 2` everywhere (§4.5) | coverage | both | — | C-oracles |
| C2 | Every warning identifier and error code asserted once per port (§4.3, §4.4) | coverage | both | — | C-ids |
| C3 | Never-called options and input modes (§4.3, §4.4, §4.5) | coverage | both | — | C-options |
| C4 | Cross-port symmetry: port the MATLAB-only and Python-only suites (§4.3) | coverage | both | — | C-symmetry |
| C5 | `util_msd*` vs `EXAMPLES.md` §1/§2 test points; positivity checks (§5.3) | coverage | both | EXAMPLES digits | C-helpers |
| D1 | `freq_map` MISO crash; 3-D `L = 1` in `freq_etfe`/Welch map/`detrend` (F1, F2, U5) | **critical** | python | — | D-crash |
| D2 | SS `residual`/`compare` sizing, `u = None`, `MaxLag` validation (U1, U4) | **critical** | both | §14.5 note | D-utils |
| D3 | Single-trajectory MIMO ETFE (F3) **[decision]** | major | both | §4.1 | D-spec |
| D4 | ETFE degenerate `Φ̂_v` convention (F4); §2.7 substitution rule (SP2) **[decision]** | major | both | §2.7, §4 | D-spec |
| D5 | `freq_domain_sim` DC / out-of-grid rule (U2) **[decision]** | major | both | §14.2, §15.2 | D-spec |
| D6 | §14 multi-trajectory contract; freq-domain 3-D divergence (U3) **[decision]** | major | both | §14 | D-spec |
| D7 | Python warning identifiers (`SidWarning` + code) (F5) **[decision]** | major | python | python guide | D-warn |
| D8 | Python parity: `show_confidence`, `map_plot` on spectrograms, `model_order(plot=)` (F6, F7, §2.2 #26) | major | python | — | D-parity |
| D9 | `sidSpectrogram` cell/list input (F8); `S₁ = 0` guards in MATLAB (F9) | major | both | — | D-shapes |
| D10 | COSMIC-half defects (§5.2) | — | — | — | D-cosmic |
| D11 | Minor conformance items F10–F20, U6–U8, `p_cov` naming, `ltv_state_est` shape (§5.0, §5.1, §5.3) | minor | both | §9 field table | D-minor |
| E1 | Spec hygiene: de-narrate, version + changelog, `output.md` §7 reconciled, frozen-TF I/O table, §8.2 options, §8.4.3 options, `Θ_k`, `N_eff`, `Method` values, §1 over-promise, §8.5/§8.12.9 `Cost` (§6) | spec | spec | yes | E-spec |
| E2 | `EXAMPLES.md`: stale branch note, row 7, noise wording, digits (§5.3 E2–E4) | spec | spec | yes | E-spec |
| E3 | "not yet in SPEC.md" headers + header-checker section-existence rule (§2.2 #33, P5) | hygiene | both, CI | — | E-headers |
| F1 | Docs site factual fixes (§7 items 1–8) | docs | docsite | — | F-docs |
| F2 | Dev-tracking leakage: curated changelog/contributing pages; API-page changelog blocks; spec narration (§7, §6) | docs | docsite, spec | — | F-docs |
| F3 | Rendering: notebook MathJax, RST directives, `sidResultTypes`, hook off-site rewrite, SUMMARY pages, TOC anchors (§7) | docs | docsite, scripts | — | F-docs |
| F4 | Missing content: estimator guide, diagnostics FAQ, citation + license, version, PyPI/`sid-toolbox` note (§7) | docs | docsite | — | F-content |
| G1 | Packaging: setuptools floor, `py.typed`, sdist contents, PyPI release, `CITATION.cff` (P1–P2) | release | python | — | G-release |
| G2 | #204 release suppression; #199–#203 | release | CI, docs | — | G-release |
| H | Publication track (§8) | — | — | — | H-papers |

---

## 1. Phase A — make the gates able to bite

Exit criterion: a deliberate 2×-tolerance mutation of *any* stored field in *any* vector fails both consumers; a regeneration commit cannot land on `main` without both validators running on it; a red test job blocks a normal merge.

1. **A1 — required checks (ADR-0005 amendment).** Add `Octave tests & examples`, `MATLAB tests & examples`, `Tests (Python 3.12)` (one matrix leg), `Validate Python port`, `Validate Octave against MATLAB reference`, and `Build site (strict)` (#199) to `.github/rulesets/main.json`. All six now start on every PR (`tests.yml` is unfiltered; `cross-validate.yml` and `python-tests.yml` are path-filtered — make them always-start with an in-job relevance gate like `tests.yml`, otherwise a required check stays pending on an unrelated PR, which is the exact failure ADR-0005 avoided). Keep the admin bypass. Fix the stale "the required check" comment in `tests.yml`.
2. **A2 — regeneration lands validated (ADR-0002 amendment).** Replace the `[skip ci]` commit with either (a) a bot PR carrying the regenerated vectors, which runs every validator before merge, or (b) a `workflow_run`-triggered validation of the regen commit that reverts it on failure. Make the provenance gate reject `git_sha: "unknown"` and verify the SHA is an ancestor of HEAD. **Reproducibility policy [decision]:** §4.1 F1 shows the pinned R2025a runner is not byte-reproducible (ULP churn in 20 files, 4.5e-5 in `ltv_io`). Options: (i) stop auto-committing — regenerate on PRs, `git diff --exit-code` against the committed vectors *within tolerance* (a numeric diff, not bytes), and commit only on a semantic change; (ii) keep auto-commit but only when a numeric diff exceeds the stored tolerance; (iii) accept churn and drop the README rule. Recommended: (i) — it makes "payloads change only by regeneration" true and kills the provenance-only re-stamp commits. Closes #204 as a side effect if the bot stops committing on `main`.
3. **A3 — NaN-capable comparison [decision].** `validate_reference.m:340` must treat `NaN`/`Inf` explicitly: equal iff both NaN, or both Inf with the same sign; otherwise fail. On disk, `jsonencode` writes both as `null`; options: (i) encode sentinels as strings (`"NaN"`, `"Inf"`, `"-Inf"`) in `writeJSON` and decode them in both consumers, (ii) store a parallel `<field>_mask`. Recommended: (i) — one convention, readable, and it unblocks degenerate-input vectors (A6/F9). Amend ADR-0002 and `testdata/README.md`.
4. **A4 — Python structural gate.** In `test_cross_validation.py`, replace the per-file hand lists with one loop that asserts every key in `ref["output"]` was compared (a `COMPARED` set per test, checked at teardown), mirroring `validate_reference.m:282–303`; route *every* comparison through `_tol()`, resolve complex keys exactly as Octave does (`<name>_real/_imag` → `<name>_rel/_atol`, and make the generator emit that form for `Phi_xx`/`Phi_xz`/`DFT`), fail on an unknown tolerance key, and compare singular values at full length. Add a `--testdata` pytest option so mutation checks can be first-class (F10) and add three permanent mutation tests (rtol, atol, missing field) against a fixture copy.
5. **A5 — atol floors enforced.** In `generate_reference.m`, a `writeJSON` pre-check that refuses to write any output field whose `min|value| < 1e-12·max|value|` (or that contains an exact zero) without a `<field>_atol`; the same check as a test in both consumers over the committed files. Audit the six default-rtol fields and give them explicit entries.
6. **A6 — the three 1 % vectors.** Regenerate `reference_ltv_io` from a case that *converges* (raise `MaxIter`, or pick a better-conditioned plant; assert `Iterations < MaxIter` in the generator and store `Iterations` and the `Cost` history), then set its tolerance from the measured cross-engine error with a floor (expect 1e-6-class if converged; if the EM path is intrinsically engine-sensitive, store the per-iteration cost so the gate pins the *trajectory*, and record the measured drift as the rationale in the generator comment). `reference_ltv_cosmic_varlen` is a closed-form solve: 1e-6 rel + 1e-10 atol like `reference_ltv_cosmic`. `reference_frozen_of_io` follows its parent. Add the IO uncertainty fields (`P`, `AStd`, `BStd`) to the IO vector.
7. **A7 — a third job that needs no port.** `testdata/check_payloads.py` (CI job `Reference payload oracles`): for each vector with an analytic plant (`model_order`, `ltv_frozen`, `lti_freq_io*`, `test_msd`, the new C1 vectors) assert the stored output against the closed form within the spec's own bound; for all vectors assert invariants — coherence ∈ [0, 1], PSD ≥ 0 where the spec guarantees it, `Φ̂_xz = conj(Φ̂_zx)`, `P(k)` symmetric PD, IO `Cost` non-increasing, `util_msd.Ad == expm(A·Ts)`. This is the only mechanism that catches a vector both ports match wrongly (ADR-0001's stated gap).
8. **A8 — no silent skips.** `test_cross_validation.py::_load` raises; `python-tests.yml` drops the `if ls` guard; `validate_reference.m` errors on an empty directory; add a test that the set of `reference_*.json` equals the set the consumers know (closes rule 3 mechanically).
9. **A9 — matrix and hygiene.** Add a `min-deps` leg (`numpy==1.22.*`, `scipy==1.8.*` on Python 3.10, or raise the floors in `pyproject.toml` to what is actually supported — **[decision]**); pin the Octave snap channel or vendor a version check that fails below 11; `filterwarnings = ["error"]` in `pyproject.toml` with explicit `pytest.warns` where a warning is the point; `scripts/local-ci` runs `scripts/check_python_api_pages.py`; fork-PR MATLAB job tolerates a missing App token.

## 2. Phase B — the relaxed, vacuous and over-claiming tests

Exit criterion: every test in §3.4's table either asserts the property its name/docstring claims or is deleted; every tolerance on a deterministic computation is ≤ 1e-8 unless a comment states the measured error and why; no `Verified by:` cites a test that does not assert the rule.

1. **B1** — `test_sidModelOrder.m` Test 5: restore `n_thresh <= 6` (or assert `== 2` on the analytic plant that Test 15 already builds) and delete the `#139` marker.
2. **B2** — calibration: run `TestUncertaintyCampaign` in a weekly scheduled workflow (`SID_MC_CAMPAIGN=1`), with the λ grid extended to 1e0 (the point that fails the current gate band on the high side); re-centre the gate band on the measured value (1.07 at λ = 1e2 — update the docstring) with an *asymmetric* band that rejects understatement harder than overstatement (e.g. [0.9, 1.3]); make MATLAB Test 13 estimate `Σ` (drop `'NoiseCov', Sigma_true`), tighten `[0.3, 3.0]` to a band the pre-#137 1.64 fails, and add the same mid-λ point.
3. **B3** — rewrite: `test_convergence` / `test_monotone_cost` on a fixture that converges (assert `iterations < max_iter` and no warning, then a separate test that the cap *does* warn); Welch `Inf` sentinel on an input that reaches the branch (constant `u` segment); `sid:noConsistentLambda` read back; MATLAB Test 6 asserts the OLS closed form; `test_sidLTIfreqIO.m` Test 15 split into "raises `orderExceedsRank`" and "defective matrix bounded"; `eig_err`/`obs_err` bounds replaced by the actual recovery error at a stated tolerance; whiteness/independence pass flags asserted on white vs AR(1) residuals in both ports; "helps" tests assert `err_L < err_1` on a seeded case where it is true; delete the `isfinite`-only and `in-range`-only assertions or replace them with the formula.
4. **B4** — `test_dft`/`test_sidDFT` 5 % → 1e-12 and fix the comment; detrend absolute tolerances → 1e-10; state-est/frozen/compare noiseless cases → 1e-8; `test_ltv_disc.py:131` ramp test asserts the *slope* (tracking) not a 0.3 band; ETFE `y = 2u` exact.
5. **B5** — MATLAB: explicit `SKIPPED` tally in `runAllTests.m`; `onCleanup` restore around every `warning('off', …)`; `clear`/`rng`/`close all` between files, or run each file in a function scope; drop `test_template.m` from discovery. Python: replace `pytest.raises(Exception)` and exception tuples with the code.
6. **B6** — correct the `Verified by:` lines listed in §3.7.

## 3. Phase C — coverage, symmetric across ports

Exit criterion: §4.5's matrix has no ✗ in a unit column for a rule the spec states; every warning identifier / error code in each port is asserted at least once; every option of every public function is called at least once per port.

1. **C1 — free oracles** (both ports, one PR each): MIMO 2×2 static-gain plant → `Ĝ == gain` to 1e-10 for BT, BTFDR, ETFE (`L ≥ 2`), `freq_map`; `L` identical copies of one trajectory reproduce the single-trajectory estimate exactly (every multi-trajectory path); N = 2 OLS closed form for COSMIC; `q = 2` and `q = 3` in every LTV test fixture (parametrise the shared simulator) and in at least one COSMIC cross-vector; a strictly-proper plant with a finite zero for `lti_freq_io` at `H ≠ I`.
2. **C2 — identifiers**: one test per `warning('sid:…')` site and one per `SidError` code (49 codes; a parametrised table test per port is enough), asserting the identifier, not the prefix.
3. **C3 — options and modes**: the lists in §4.3/§4.4 (plot options, `covariance_mode` full/isotropic with a full ≠ diagonal assertion, user `noise_cov` / `R` / `Q`, per-step λ, `lambda_grid`, `max_lag`, `initial_state`, 3-D `L = 1` for every function, list input for BT/ETFE/BTFDR/`ltv_disc`/`ltv_disc_io`/`lti_freq_io` with the trim warning asserted, MISO/SIMO `freq_map`, multi-output time series everywhere).
4. **C4 — symmetry**: port `test_sidLTVdiscUncertainty.m`, `test_sidLTVdiscVarLen.m`, the FD-Jacobian check and the ±5 % 1/√L ratio to Python; port the five `TestFreqMapVariableLength` cases, the Welch-vs-scipy oracle (as a rect-window periodogram oracle), zero-input / `u = y` / complex-`u`, defective stabilization, MIMO block Hankel, degenerate window, MIMO noise-spectrum plot to MATLAB.
5. **C5 — helpers**: assert every `EXAMPLES.md` §1 / §2.2.3 / §2.3.4 test point in both ports once E2 extends the printed digits (or lowers the demanded tolerance to what the digits support); enforce §2.1.4 positivity; fix `util_msd_ltv.m` `(1×N)` mass.

## 4. Phase D — conformance defects, spec first where a decision is needed

Order within the phase: crashes (no decision) → the four spec decisions (one PR each, spec commit first) → parity → minor.

1. **D1 (Python, no decision)** — `freq_map` MISO storage (`ny == 1 and nu == 1` for the ravel branch, or allocate per `(ny, nu)`); 3-D `L = 1` in `dft.py` (squeeze the trailing axis once, at validation) and the Welch inner path; `detrend` 3-D branch on `x.ndim`. Tests from C3 pin them.
2. **D2 (both, no decision)** — SS `residual`/`compare`: size from the data, require `N_data ≤ DataLength` with `sid:dimMismatch`/`dim_mismatch` and a clear message, typed error for `u = []` on SS models; `1 ≤ MaxLag ≤ N − 1`.
3. **D3 [decision]** — single-trajectory MIMO ETFE: (i) error `sid:etfeNeedsTrajectories` unless `L ≥ n_u`, (ii) fall back to a per-input SIMO ratio with a documented bias, (iii) leave NaN but document it. Recommended: (i) — the estimator is undefined, and a silent all-NaN result with a "some frequencies" warning is the worst option.
4. **D4 [decision]** — §2.7 at degenerate frequencies: one convention for `Φ̂_v` and `γ̂²` when the guard or the excitation check fires, applied by every estimator including ETFE. Recommended: NaN for both (what §2.7 already says "is produced"), because `Φ̂_y` with coherence 0 is a *different* estimate presented in the noise-spectrum field. Then fix `sidFreqETFE.m:253` and `degenerate.py` / `sidRegularizeResponse.m`, and store a degenerate-input cross-vector (unblocked by A3).
5. **D5 [decision]** — `freq_domain_sim`: specify DC/out-of-grid handling. Recommended: hold `Ĝ(ω₁)` down to DC (zero-order extrapolation, consistent with `lti_freq_io`'s DC convention being an extrapolation too), force the Nyquist bin real, and state in §14.3 that whiteness on non-detrended data tests the model *and* the detrending.
6. **D6 [decision]** — §14 multi-trajectory: write the pooled-bound contract both ports already implement; for frequency-domain models with 3-D data, either reject in both ports or simulate per trajectory — recommended: reject with `sid:dimMismatch` until a use case needs it.
7. **D7 [decision]** — Python `SidWarning(UserWarning)` with a `code` attribute mapped 1:1 to the MATLAB identifiers (`window_reduced` ↔ `sid:windowReduced`, …), documented in `python/CONTRIBUTING.md`, asserted by C2.
8. **D8 (Python)** — `show_confidence` in `bode_plot`/`spectrum_plot` (and document it in the MATLAB headers); `map_plot(plot_type='spectrum')` accepts `SpectrogramResult`; `model_order(plot=)`.
9. **D9** — `sidSpectrogram` cell input (trim-to-shortest + `sid:trimmedTrajectories`, as §2.3) and a Python list branch that does not silently mis-dispatch; `S₁ ≤ 0`/NaN guards in `sidFreqMap.m`/`sidSpectrogram.m` and the `window_length = 1` case in both.
10. **D10 — COSMIC half (§5.2).** In order: **C3** (no decision — Python catches the singular pivot and warns-and-returns non-finite like MATLAB, or both ports raise a typed `sid:singularLbd` error: **[decision]**, recommended warn-and-return to match §8.3.4, plus a test built from the collinear-state fixture); **C2 [decision]** — restore the relative test, with a documented absolute floor only if the maintainer wants one (then it is `max(|J|, ε_abs)` with `ε_abs` in the spec, not 1); **C4** — divide by `N` (numerics change → regenerate `reference_ltv_tune`, and add a unit oracle with a constructed error so the vector is no longer the only verifier); **C1 [decision, ADR-0004 revisit]** — options: (i) include the cliff ratio when `σ_{L+1}` is a numerical zero (exact-data case) while keeping the exclusion when the tail is resolvable, (ii) a data-aware floor from the artifact level (`|Re(2Ĝ₁ − Ĝ₂)|` scale), (iii) replace the DC extrapolation by a one-sided fit that removes the artifact and then include the cliff; recommended (iii) + (i), with the 3rd-order plant at `nf = 512` and a MIMO `n = 3` plant as the new oracles and the cross-vector regenerated; **C5** — the final μ = 0 pass warns when it caps (it *is* the base alternation); **C6** — one time-step convention for frequency-mode tuning (0-based, per §6.7), pinned by a frequency-mode cross-vector; **C7–C16** — typed validation errors (`bad_lambda`, `bad_noise_cov` with symmetry/PSD, `R` SPD, `lambda_grid` sorted/positive/≥ 3, `MaxLag`), copy `noise_cov`, clamp `P` diagonals before `sqrt`, specify or remove the second DoF fallback, fast-path `cost`/`iterations` consistency, `ltv_state_est` input immutability and 2-D `A`, `blk_tri_solve` check ordering, `Horizon` clamp warning, delete the orphan expression, specify the `lti_freq_io` trajectory discard/trim and the frequency-mode tuning extras (or remove them), unify the condition-number metric (**[decision]**: `rcond`-style 1-norm in both, or 2-norm in both), align `tooFewSV`/`too_short`, strip the `#137` docstring citation.
11. **D11** — F10–F20, U6–U8; rename `p_cov` → `p` and `complex_stft` → `complex` (or amend §9's "purely syntactic" claim with an explicit exceptions table — **[decision]**, a public API change); `ltv_state_est` returns the shape it was given; add `ResidualResult` extras to §14.6 or drop them; `ltv_state_est`/`detrend` list input errors typed.

## 5. Phase E — the specification as a contract, not a changelog

1. **E1** — one PR, spec-only, no behaviour change: remove the "Correction (…)", "previously", "earlier drafts", "supersedes" narration and the issue numbers from normative text (move the history to a `spec/CHANGELOG.md` or to the ADRs); bump `SPEC.md` to a version that reflects 0.2.0 and add the versioning rule `EXAMPLES.md` §6 already has; reconcile `spec/cosmic/output.md` §7/§8/§10 and Appendix B with §4 (delete or clearly mark Appendix B as an alternative never used); add the `sidLTVdiscFrozen` inputs/outputs table; add `Uncertainty`/`NoiseCov`/`CovarianceMode` to §8.2 and `Method`/`SegmentLength`/`ConsistencyThreshold`/`CoherenceThreshold` to §8.4.3; fix `Θ_k`; define `N_eff`; fix the `Method` enumeration and either amend it to the Python values or change the Python values (**[decision]**, API change); scope §1's "all functions" sentence; give the IO `Cost` field its own name or state the overload; put `L` into §3.4/§3.5; state the ETFE `WindowSize`; resolve SP1 (ETFE and the `Inf` sentinel), SP3 (BT short-data default), SP5 (non-negativity claim), SP8 (Welch `σ_G`, MIMO NaN, DC bin); disambiguate §14.2.
2. **E2** — `EXAMPLES.md`: delete the branch note and the "Python is the reference implementation" sentence (ADR-0001), fix §3.6 row 7, write noise levels as std, state the chirp formula as implemented, extend the tabulated digits (or lower the demanded tolerances), bump the version.
3. **E3** — replace the twelve "not yet in SPEC.md" placeholders with the section numbers; make the header checkers verify that every cited `§x.y` exists in `SPEC.md` and reject "not yet in SPEC" (both already pass with 0 dangling citations).

## 6. Phase F — documentation site

1. **F1** — the eight factual fixes in §7 (compatibility floor first).
2. **F2** — stop including `CHANGELOG.md` and `CONTRIBUTING.md` verbatim: a curated `about/changelog.md` without issue numbers and a user-facing `about/contributing.md` (how to file issues, run tests, add a port) that links to the full engineering guide on GitHub; strip the generated "Changelog — 2026-04-08" blocks from API pages (a mkdocstrings filter or a docstring section the gen script drops); the spec narration goes with E1.
3. **F3** — MathJax: either let `mkdocs-jupyter` render math (set `include_source`/`ignore_h1_titles` and add `arithmatex` to the notebook pipeline, or configure `tex.inlineMath`/`displayMath` and drop the `ignoreHtmlClass` restriction for `.jupyter-wrapper`) and verify in a browser; convert `.. math::`/`.. [1]` in the two docstrings to plain numpy-style text; special-case `sidResultTypes.m` in `build_matlab_api.py` (render as a reference page, not a function) or move its content into the spec include; fix `rewrite_external_links.py` to leave links that resolve inside `docsite/`; `not_in_nav`/`exclude_docs` the four `SUMMARY.md`; remove backticks from notebook headings; pin `requirements-docs.txt`; add ruff to it for signature formatting; make `build_matlab_api.py` fail on an empty header and derive `PYTHON_EQUIV` from the roadmap catalogue instead of a hand list.
4. **F4** — new pages: "Choosing an estimator" (BT vs BTFDR vs ETFE vs `freq_map`; COSMIC vs Output-COSMIC; when to tune λ and how), "Diagnostics and warnings" (the contract identifiers, what NaN/`Inf` mean, what to do), "Cite" (`CITATION.cff` + the two arXiv references), a license page, the version in the header, the PyPI/`sid-toolbox` note, a note that `util_msd.py` ships beside the notebooks, examples for `lti_freq_io`/`ltv_state_est` (also fixes the "every estimator ships with an example" claim); then #201 (`mike`) and #200 (lychee).

## 7. Phase G — release mechanics

1. **G1** — `setuptools >= 77`; `py.typed`; sdist includes `conftest.py`, `examples/util_msd.py` and the vectors the shipped tests need (or stop shipping tests); first PyPI release of `sid-toolbox` (trusted publishing from `release.yml` on the `-python` tag); `CITATION.cff`; sync the Python-version claims (3.10–3.14).
2. **G2** — #204 (choose option 1: real releases from `workflow_dispatch`, and A2 removes the bot commit from `main` anyway), #203, #202.

## 8. Phase H — publication track (after A–C are green)

Order by readiness (analysis §8): (1) JOSS/SoftwareX paper for the toolbox — needs G1 and a green, required CI; (2) the Bayesian-uncertainty note — needs the MC campaign running (B2), a sandwich/bootstrap comparison, and the Comet Interceptor demonstration; (3) Output-COSMIC — needs E1's `output.md` reconciliation, A6's converged vector, a trust-region characterisation study, and experiments vs EM-LTV / windowed subspace; (4) λ-tuning by consistency — needs a correctness oracle for §8.11.2 (C3) and benchmarks vs L-curve / validation; (5) online COSMIC — needs an implementation; (6) the methodology experience report — needs the maintainer's decision on discussing the agent workflow, and can cite this cycle's data.

---

## 9. Suggested execution order and dependencies

```
A1 A8 A9 ──────────────────────────────┐
A3 ──► A4 A5 ──► A6 ──► A7            │   (gates first; A3 unblocks degenerate vectors)
B1 B4 B5 B6 (independent)              │
B2 (needs A9's scheduled-job pattern)  ├──► C1 C2 C3 C4 C5 (coverage, both ports, one PR per row)
D1 D2 (no decision; can start day 1)   │
D3 D4 D5 D6 D7 (spec commit, then both ports, one PR each; D4 then adds the degenerate vector)
D8 D9 D10 D11 ─────────────────────────┘
E1 E2 E3 (spec-only PRs; E1 before D3–D6 land so they edit a de-narrated spec)
F1 (day 1) → F2 F3 → F4
G1 G2 → H
```

Working method as in July: one branch + PR per numbered item, spec commit first inside the PR where a decision is involved, full `scripts/local-ci` before every PR, cross-vectors regenerated only where numerics legitimately change and flagged in the PR, every PR closes its issue, every deviation from this plan noted on the issue.
