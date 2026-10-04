# sid v0.2.0 — Repo-wide review v2 (software, testing facilities, documentation site)

**Date:** 2026-10-04
**Reviewed at:** `1785dd3bc` (`main` HEAD at review time — the post-#110 reference-data re-stamp; `v0.2.0-matlab` / `v0.2.0-python` tags point at `63b606f`). File:line anchors refer to this commit.
**Scope:** `spec/`, `python/` (sources, tests, examples, packaging), `matlab/` (sources, tests, examples), `testdata/` (generator, both validators, all 32 vectors), `.github/` (workflows, ruleset, header checkers), `docs/` (ADRs, plans, analyses), `docsite/` + `mkdocs.yml` + `scripts/` (the public site), the READMEs and release notes.
**Method:** eight independent review passes run in parallel by reviewer agents, each seeded with [`docs/REVIEW_CONTEXT.md`](../REVIEW_CONTEXT.md) and the relevant spec sections — (1) git-history archaeology of every test / gate relaxation over all 306 commits, (2) the Python test suite, (3) the MATLAB/Octave test suite, (4) the cross-vector infrastructure, with *mutation testing* of both consumers, (5) the documentation site, built with `mkdocs build --strict` including notebook execution, (6) spec-to-implementation conformance of the frequency-domain half (§1–§7, §9–§11), (7) of the COSMIC half (§8 + `spec/cosmic/*`), (8) of the utilities, example suite and packaging (§13–§15, `EXAMPLES.md`) — plus first-hand numerical probes by the lead reviewer where a claim mattered (calibration ratios, multi-input COSMIC, regeneration drift). Every finding below was re-verified against the tree at the anchor commit before being recorded.
**Toolchain:** Python 3.11.15, NumPy 2.4.6, SciPy 1.17.1, pytest 9.1.1, mkdocs 1.6.1. The Python suite passes at HEAD (415 passed, 4 skipped). **No MATLAB engine was available, and no Octave while the review passes ran**: every MATLAB finding in §2–§7 is a static reading of the `.m` source and is labelled as such. After the passes completed, GNU Octave 8.4 (apt) was installed and `runAllTests.m`, `validate_reference.m` and `runAllExamples.m` were executed under it; those results are recorded in the companion plan's §10.1 (374/380 test cases and 32/32 vectors pass; the six failing test files and 11/12 examples are blocked by an Octave-8.4 figure-rendering defect on the image, not by sid). The site builds clean with `--strict` (all 11 notebooks execute).
**Baseline:** the previous review [`2026-07-12-repo-wide-review.md`](2026-07-12-repo-wide-review.md) and its remediation plan [`2026-07-22-review-remediation-plan.md`](../plans/2026-07-22-review-remediation-plan.md). This review deliberately re-checks what that cycle claimed to close (§2) before looking for new defects.

Severity scale (same as the July review): **critical** = wrong results or crash on documented usage; **major** = spec violation, silent numerical error in realistic use, or a gate that cannot do its job; **minor** = contract / edge-case / documentation defect; nit = cosmetic.

Companion plan: [`docs/plans/2026-10-04-review-v2-remediation-plan.md`](../plans/2026-10-04-review-v2-remediation-plan.md).

---

## 1. Executive summary

**The mathematics is sound; the verification apparatus is not yet trustworthy enough to prove it.** Every numerical core re-derived in this review — the Blackman-Tukey machinery, the §3 variance formulas, Welch and spectrogram scaling, the COSMIC closed form, the corrected posterior covariance and degrees of freedom, the RTS smoother including ragged trajectories, the frozen-transfer-function variance, the real-Schur stabilization, the two-level trust region — matches an independent dense or closed-form oracle to machine precision in Python and is implemented identically in MATLAB. The July remediation genuinely closed every critical finding and nine of ten majors.

What this review adds is mostly about the **right side of the V**:

1. **The cross-language gate cannot do the job the project assigns to it.** The canonical generator is not reproducible (a regeneration commit with no code change rewrote 21 of 32 payloads, one by 4.5e-5), and CI commits such churn under `[skip ci]` so no validator runs on it; the Python consumer silently ignores stored fields (a mutated frequency grid and a zeroed DFT both pass); the Octave validator cannot fail on NaN; the three vectors guarding Output-COSMIC, variable-length COSMIC and frozen-of-IO sit at 1 % with no recorded rationale, and one of them pins a *non-converged* iterate; and no mechanism at all can catch a vector both ports match wrongly. Tests and cross-validation are not required status checks — only lint is. (§3.1, §4.1, §4.2)
2. **Tests were not loosened to make PRs green — but a number of them cannot fail for the reason they exist.** One bound was widened to absorb a since-fixed bug and never restored; both Monte-Carlo calibration bands would have passed the very bug §8.9 was corrected for, and the Python gate's documented ratio (0.89) is not what the code produces (1.07); the Output-COSMIC "convergence" test runs on a fixture that never converges; the Welch `Inf`-sentinel test never reaches the sentinel; whiteness/independence pass flags are asserted nowhere in either port; eleven deterministic linear-algebra checks run at 1–50 % tolerance where 1e-12 holds; the only MIMO numerical check in the MATLAB suite is toolbox-gated and silently counted as passed when skipped. (§3)
3. **Functionality not covered is systematic, not scattered:** no test or vector in either port uses more than one input (`q = 1`) for any LTV function; MIMO numerics for every frequency estimator are shape-only; 35 of 49 Python error codes and 16 of 21 MATLAB warning identifiers are never asserted; the two suites are asymmetric in ~40 named cases; the example-suite helpers are never checked against the spec's own tabulated test points. (§4.3–§4.5)
4. **New defects on documented inputs** — 3 critical (Python `freq_map` crashes on MISO data; 3-D single-trajectory input still crashes ETFE and the Welch map; state-space `residual`/`compare` crash or silently truncate on validation data of a different length), 17 major (among them the same 3-D crash in `detrend`), some forty minor — plus three places where both ports agree on a number the spec does not sanction: the model-order gap rule of ADR-0004 cannot detect an exactly rank-`n` Hankel and returns the wrong order on an exact 3rd-order plant; the Output-COSMIC convergence test is absolute below `J = 1`; the validation loss divides by `N + 1`. Roughly thirty of the July minor findings are unchanged. (§2.2, §5)
5. **The specification is drifting from contract toward changelog**: it narrates its own corrections, cites issue numbers in normative text, carries a stale version header, contradicts itself in `output.md` (two initialisations, two trust-region anchors), has no I/O table for one catalogued function, and its `Verified by:` annotations over-claim in at least nine places — all of which is published verbatim on the new documentation site. (§6)
6. **The documentation site builds clean and is structurally good**, but its MATLAB install page states a platform floor two major versions below the real one, four API-index sentences are false, the Python uncertainty page uses MATLAB field names, notebook math is almost certainly unrendered, and the changelog/contributing includes leak the entire engineering process onto user pages. (§7)
7. **Six things are publishable**, in readiness order: the toolbox (JOSS), closed-form Bayesian uncertainty for COSMIC, Output-COSMIC, λ-tuning by parametric/non-parametric consistency, online COSMIC, and the spec-as-contract/agent-driven methodology as an experience report. None has appeared in the literature. (§8)

The plan (§9, companion document) repairs the gates first, then the tests that lean on them, then coverage, then the conformance defects spec-first, then spec and site hygiene, then release mechanics and the publication track.

---

## 2. What the July remediation actually closed

The July review filed 12 issues (#134–#145) from 5 critical, 10 major and 43 minor findings. All 12 issues are closed and the tracker holds only six open docs/CI follow-ups (#199–#204). Re-checking each finding against the tree, rather than the tracker, gives a less tidy picture.

### 2.1 Critical and major findings

| July finding | Status at `1785dd3` | Evidence |
|---|---|---|
| C1 var-len smoother terminal block | **Fixed, both ports, oracle-pinned** | `ltv_state_est.py:213–230`, `sidLTVStateEst.m:150–165`; dense-LSQ oracle `test_ltv_state_est.py:248` (atol 1e-10), `test_sidLTVStateEst.m:367–426` |
| C2 multi-output time-series `freq_map` crash | **Fixed, both** | storage branch now `ny == 1` only; `test_freq_map.py` multi-output TS case; no MATLAB twin (§4.3) |
| C3 3-D `L = 1` crashes | **Fixed, Python** (`cov.py:78`, `spectrogram.py:344`); covered for `freq_bt`/`freq_map`/`spectrogram`, **not** for `freq_btfdr`/`freq_etfe`/`ltv_disc` | §4.4 |
| C4 complex `u` cast | **Fixed** (`validate_data.py`), `test_error_complex_u`; MATLAB has no complex-`u` test | §4.3 |
| C5 MIMO plotting | **Fixed**, smoke-tested (`test_plotting.py:203,287`); `sidMapPlot` accepts spectrogram results | — |
| S1 posterior scaling | **Fixed, spec-first**, dense-Hessian oracles in both ports; **the calibration tests that guard it are weaker than they look** (§3.4) | — |
| S2 trust region | **Fixed** (two-level schedule); the only tests are behaviour pins at a tuned budget (§3.5) | `test_sidLTVdiscIO.m:387–409`, `test_ltv_disc_io.py:351–410` |
| S3 model-order lag-0 | **Fixed** (+ ADR-0004 gap convention); **a test bound widened to absorb S3 was never restored** (§3.3) | `test_sidModelOrder.m:98–104` |
| S4 multi-traj fit | **Fixed, both**, mirror-image test present | — |
| S5 BTFDR guards | **Fixed** via shared `degenerate.py` / `sidRegularizeResponse.m` | — |
| S6 Welch scaling | **Fixed** (Nyquist un-doubled; 2× BT↔Welch documented) | — |
| S7 σ conventions | **Fixed, spec-first** (`Inf` sentinel, absolute excitation check); the two tests pinning the violation were corrected, not deleted | `test_uncertainty.py:99` |
| S8 `_stabilize` blow-up | **Fixed** (real-Schur), defective-matrix test in Python; MATLAB counterpart is an "errors **or** finite" test (§3.5) | `test_sidLTIfreqIO.m:468–479` |
| S9 frozen-of-IO | **Fixed, spec-first**, cross-vector added — **but that vector pins a non-converged iterate** (§4.1 F7) | — |
| S10 cross-validation blind spots | **Partially**: the `P` mis-generation, the four orphans and the divergent tolerances were fixed; the mechanism that was supposed to make the vectors trustworthy has five new holes (§4.1) | — |

### 2.2 Minor findings still open

The following July minor findings are **unchanged** at the anchor commit (numbering as in the July report):

- **#6 / #7** — ETFE `WindowSize` reported as `N`; `Method` strings. Python returns `'freq_bt'`, `'ltv_disc'`, … while SPEC §9 enumerates `'sidFreqBT'`, … (`_results.py:59–60`, `SPEC.md:1860`). Neither the spec nor the port was changed; the §9 PascalCase note covers field *names*, not values.
- **#26** — `model_order` has no `'Plot'` option in Python (`model_order.py` signature) though SPEC §8.12.12 lists it.
- **#32** — `'ShowConfidence'` (SPEC §11.3) is absent from `bode_plot` / `spectrum_plot` in Python, and present but **undocumented** in the MATLAB headers (`sidBodePlot.m:18–23` lists five options; `defs.ShowConfidence` at `:56`).
- **#33** — twelve headers/docstrings still say "not yet in SPEC.md" for behaviour that *is* specified: `detrend` (§13), `residual` (§14), `compare` (§15), `bode_plot` / `spectrum_plot` (§11), `validate_data` (§10.1), in both ports (`residual.py:89`, `compare.py:84`, `detrend.py:121`, `bode_plot.py:88`, `spectrum_plot.py:76`, `_internal/validate_data.py:88`; `sidCompare.m:44`, `sidDetrend.m:40`, `sidResidual.m:41`, `sidBodePlot.m:36`, `sidSpectrumPlot.m:33`, `private/sidValidateData.m:47`). Both header checkers accept the placeholder — they test only that the word `SPEC.md`/`SPECIFICATION` appears (`check_python_headers.py:231–238`, `check_headers.py:216–218`) — so the traceability gate is vacuous for these six functions per port. mkdocstrings publishes the Python placeholders on the API pages.
- **#35** — `conftest.load_reference` and the `rng` fixture are dead (`python/tests/conftest.py:34–55`); `test_cross_validation.py:32–41` redefines its own loader. `python/CONTRIBUTING.md:393–406` documents the dead fixtures as the shared ones.
- **#37** — `Θ_k = D(k)ᵀ X'(k)ᵀ` at `SPEC.md:968` is still dimensionally impossible (`(p+q)×L · p×L`); both ports compute `D(k)ᵀ X'(k)`.
- **#38** — the §8.2 input table (`SPEC.md:914–920`) still omits `Uncertainty`, `NoiseCov`, `CovarianceMode` which §8.5/§8.9 presuppose and both ports accept.
- **#39 (second half)** — `uncertainty_derivation.md:598–602` still says covariance mode `Full` is the "Default"; SPEC §8.9.3 and both ports default to `diagonal`.
- **#40** — `EXAMPLES.md:12–16` still says the spec was "authored against the Python example suite on branch `claude/smd-util-helpers`" and that "the Python port is the v1.0.0 reference implementation of this spec" — a stale work-branch reference and a sentence that contradicts ADR-0001 (no port is a reference implementation).
- **#41** — `testdata/README.md:16` documents the schema key as `function` (files use `function_name`); `:57–59` still says "currently Octave; Python and Julia will be added in future phases" although the Python consumer has existed since April.

Closed minor findings (verified): #1 ETFE NaN propagation, #3 BT time-series clamp, #10 Welch tiny-segment, #13 MIMO whole-slice NaN, #15 NaN clamp semantics (spec'd in §2.7), #16 singular pivot, #19 `_frequency_tune` N, #25 DC extrapolation unified, #36 dead `lpv_extension_theory.md` reference, #42 `tests.yml` path trigger. Findings #2, #4, #5, #8, #9, #11, #12, #14, #17, #18, #20–#24, #27–#31, #34 are re-assessed in §5 (conformance) and §4 (coverage) below.

---

## 3. Tests that were wrongfully relaxed, weakened, or never tightened

This is the question the review was asked to focus on. The short answer: **nothing in the remediation window (PRs #134–#198) was loosened to make a PR green** — every July tolerance edit added an absolute floor with a stated reason, replaced a weaker assertion with a stronger one, or corrected a test that pinned a spec violation. The relaxations that matter are older, or are *bounds that were never tightened once the defect they absorbed was fixed*, or are tests whose stated purpose their assertion cannot deliver.

### 3.1 History: the three 1 % cross-vectors have no recorded rationale

| Commit | Date | What | Verdict |
|---|---|---|---|
| `8a1c436` | 2026-04-08 | `reference_ltv_cosmic` `A_rel`/`B_rel` **1e-10 → 1e-6** ("Relax … for cross-engine validation", no measurement) | Questionable: the same day `2b71650` measured the actual cross-engine error at ~4e-10 *absolute*; an atol floor (what #145c eventually did) was the right fix, not a 10⁴× rtol. Harmless today only because 1e-6 is still tight for O(1) entries. |
| `da5a625` | 2026-04-08 | `reference_ltv_io` `A_rel`/`B_rel`/`Cost_rel` **born at 1e-2** | **No justification recorded.** 10⁴× looser than every other LTV vector. The only rationale on record is `testdata/README.md:34` ("the LTV-IO solve drifts ~1e-4"), written three months later. |
| `d870645` (#185) | 2026-07-25 | `reference_ltv_cosmic_varlen` A/B/Cost **born at 1e-2**, `atol 0` on A/B | **No justification recorded.** Equal-length COSMIC on the *same closed-form solver* is pinned at 1e-6; the ragged path inherited the IO figure by copy. This is the exact code area of July C1. |
| `528a242` (#190) | 2026-07-26 | `reference_frozen_of_io` `Response_rel` **born at 1e-2** | Coherent with its IO parent, but the parent is unjustified (see F7, §4.1: the pinned IO iterate is non-converged). |

Mutation test (§4.1): scaling **all** of `reference_ltv_io.A` by 1.005 passes the Python consumer. These three vectors cannot adjudicate a 0.5 % regression in Output-COSMIC, variable-length COSMIC, or the frozen-of-IO contract — the three paths with the most remediation churn.

The April loosening `af95b60`/`9a9da93` (`reference_internals` `C_rel` 1e-10 → 1e-8, `C_atol` 1e-10) and every July atol floor (`29c30e9`, `cc54447`, `954f952`, `2ae1457`) are justified and measured. The #145c switch from hardcoded to JSON tolerances made nothing looser on the Python side and tightened six fields; the only Octave-side loosening was the `atol = 1e-10` floor on `reference_ltv_cosmic` that fixes the July S10(e′) flake.

### 3.2 A bound widened to absorb a bug, never restored after the fix

`matlab/tests/test_sidModelOrder.m:98–104` (Test 5, threshold method on a true 2nd-order plant):

```matlab
% ... The threshold method's overcounting on noisy tails is part of the
% deferred noise-floor work (#139).
assert(n_thresh >= 1 && n_thresh <= 12, ...
```

The bound was widened from `<= 6` to `<= 12` during #157 to absorb July S3; #139 is **closed** (lag-0 fix landed; #160/ADR-0004 settled the floor). The marker now points at a closed issue — exactly what ADR-0003 rule 2 forbids — and a threshold method returning `n = 12` for a 2nd-order plant passes. This is the one unambiguous "relaxed and forgotten" test in the tree.

### 3.3 Monte-Carlo calibration: both bands are too wide to verify §8.9

- `test_sidLTVdiscUncertainty.m:302–362` (Test 13, `nMC = 200`, λ = 1): asserts every ratio in **[0.3, 3.0]** and passes the *true* `Σ` through `'NoiseCov'`, so it tests neither the noise-covariance estimate nor the DoF correction. The July S1 bug over-stated σ by ≤ 1.9× at mid-λ — inside this band. It is cited in SPEC §8 `Verified by:` as the MATLAB verifier of §8.9; it would not have caught the defect the section was corrected for.
- `python/tests/test_ltv_uncertainty_calibration.py:131–143`: band **[0.6, 1.2]** on the median ratio at λ = 1e2. The docstring claims "new code: median ratio ≈ 0.89"; **measured at HEAD the ratio is 1.070** (same seeds, same call), so the docstring is stale and the stated "would be decoration if it passed both" argument was never re-checked. The lower bound admits a 40 % *understatement* of σ (the anti-conservative direction) without failing. A finer probe across λ (first-hand, 120 trials, seeds as in the test):

  | λ | reported/empirical | with true Σ | Σ̂/Σ |
  |---|---|---|---|
  | 1e0 | 1.33 | 1.34 | 1.00 |
  | 1e1 | 1.18 | 1.16 | 1.03 |
  | 1e2 | 1.07 | 1.05 | 1.05 |
  | 1e3–1e5 | 1.05 | 1.02 | 1.05 |

  The ratio is ≥ 1 everywhere — consistent with `uncertainty_derivation.md` §4.2 (Bayesian ≽ sandwich) — and the 1e0 point is *outside* the gate band on the high side. The band is centred on neither theory nor measurement.
- The 4-λ campaign (`:146–156`) is `skipif(not SID_MC_CAMPAIGN)`; nothing in CI sets it, so the 500-trial sweep adopted under ADR-0014 has never run anywhere. The only MC evidence in the gate is the single λ = 1e2 point.

### 3.4 Tests whose assertion cannot deliver their stated purpose (both ports)

Class 5 (tautological) and class 4 (smoke) tests that *claim* to verify a contract rule. The most consequential, with the rule they are cited for:

| Test | Claims | Actually asserts | Verifier of |
|---|---|---|---|
| `test_ltv_disc_io.py:287` `test_convergence` | convergence | `iterations <= max_iter`; the fixture (seed 200, `H=[1 0]`) **hits the cap** 50/50 with the `notConverged` warning uncaught and 44 % mean-A error | §8.12.3 |
| `test_ltv_disc_io.py:139` `test_monotone_cost` | monotone J | runs on the same non-converging fixture | §8.12.5 |
| `test_freq_map.py:674–683` | `σ = Inf` sentinel on Welch degenerate input | "no NaN"; measured 0 Inf, 0 NaN, min coherence 1.7e-5 — the sentinel branch is never reached | §6.5, `Verified by:` at `SPEC.md:759` |
| `test_sidLTVdiscTune.m:184–195` (Test 10) | fallback warning `sid:noConsistentLambda` | `lastwarn('')` cleared and never read; asserts a fraction ∈ [0,1] | §8.11.2 |
| `test_sidLTVdiscUncertainty.m:119–137` (Test 6) | "N = 2 OLS closed-form check" (credited in July §5) | `eig(P) > 0` only | §8.9 |
| `test_sidLTIfreqIO.m:468–479` (Test 15) | order-above-rank / defective stabilization | passes if it errors **or** returns finite | §8.13.1 |
| `test_sidLTVdiscIO.m:458, :492` | eigenvalue / observability recovery | `eig_err < 1.0`, `obs_err < 1.0` (100 % relative error passes) | §8.12 |
| `test_reviewFixes.m:57–97` (BUG-2) | FD Jacobian | `isfinite` of the FD derivative | §8.11.1 |
| `test_sidResidual.m:6–40` (Tests 1–2) | whiteness pass/fail | field presence; comment admits `WhitenessPass` is not asserted. `whiteness_pass` / `independence_pass` are asserted **nowhere in either port** | §14.3–14.4 |
| `test_ltv_disc_tune.py:302–303, :336, :406–408` | tuning selects a λ | `grid[0] <= best <= grid[-1]`; fraction ∈ [0,1] | §8.4.3, §8.11 |
| `test_ltv_disc.py:255–262` | cost triple | `cost[0] == cost[1] + cost[2]`, which `ltv_evaluate_cost.py:107` computes | §8.3.3 |
| `test_ltv_disc.py:176`, `test_freq_etfe.py:254`, `test_sidLTIfreqIO.m:170`, `test_sidLTVdiscIO.m:308,341`, `test_sidFreqETFE.m:172` | multi-trajectory / R-weighting *helps* | `err_L <= 1.5–3 × err_1` — the opposite passes | §1, §8.12.10 |
| `test_sidModelOrder.m:146–156` (Test 8) | `'Plot'` option | `try/catch` rethrows only if the message lacks "figure"/"display"/"gnuplot" | §8.12.12 |
| `test_sidModelOrder.m:240–257` (Test 13) | cap on noise-only data | `n < nSV` | §8.12.12 |
| `test_util_msd.py:198–199` | — | `ratio = … if … else 0; assert isfinite(ratio)` | `EXAMPLES.md` §2 |

`test_sidLTVdiscIO.m` Test 11 (`:387–409`) and `test_ltv_disc_io.py:351–410` are the only verifiers of the §8.12.4 trust region. Their own comment says `cost_tr <= cost_off` "is NOT a general invariant" and that the `0.5×` margin and `MaxIter = 40` were chosen because "at this budget TR beats off by ~2.5 decades": they are revert-checks against the pre-#138 loop, not contract tests, and the µ-schedule clauses the spec marks normative (best-iterate-per-stage, guarded final µ = 0 pass) are `manual`.

### 3.5 Deterministic linear algebra at statistical tolerances

Measured error in parentheses where the lead reviewer or a pass re-ran the case:

- `test_dft.py:47` / `test_sidDFT.m:39` — FFT fast path vs direct DFT on the default grid at **5 %** (actual 2.2e-14); the MATLAB comment "tolerance needs to be moderate due to FFT binning/interpolation" is wrong — the fast path lands exactly on the grid by construction (§2.5.1), and `test_sidWindowedDFT.m:22` holds the same machinery to 1e-10. Open since the July review named it.
- `test_detrend.py:39,67,105,163,196` and `test_sidDetrend.m:17–19,44,74–75,157` — absolute tolerances **0.1–5** on noiseless polynomial fits that LS reproduces to 1e-13.
- `test_ltv_state_est.py:115,121,160` (1 %, 5 %, 1e-4 on noiseless data; actual 1e-13) and `test_sidLTVStateEst.m:40,89` (1e-6, 1e-4) while `:422` / `test_ltv_state_est.py:283` prove 1e-10 holds.
- `test_ltv_disc_frozen.py:128` / `test_sidLTVdiscFrozen.m:76` — frozen TF vs analytic on noiseless data at **0.02** (actual 6e-8).
- `test_ltv_disc.py:131` — `atol = 0.3` on an A-entry ramp 0.5 → 0.9: a constant 0.7 passes; measured max error 0.247. The test cannot distinguish LTV tracking from LTI averaging — i.e. it cannot tell whether λ is doing anything.
- `test_compare.py:67,199,218` — `fit > 90` for a noiseless perfect model (99.99998); `test_sidCompare.m:71,92,180` `Fit > 50` / `< 50`.
- `test_freq_etfe.py:177` — 0.01 for `y = 2u`, where the ETFE is exactly 2 at every bin.
- `test_compareMultiTraj.m:111` — ETFE vs `etfe(merge(...))` at **50 %** max relative error; the comment concedes the estimators differ, which makes the oracle a non-oracle.
- `test_sidLTIfreqIO.m:285` `mp_err < 0.5`; `test_sidLTVdiscUncertainty.m:88` `P` vs `(DᵀD)⁻¹` at 0.1 where ~1e-5 holds at λ = 1e-8.

### 3.6 Silent skips and swallowed failures

- `test_cross_validation.py:35–41` `_load()` → `pytest.skip` when a JSON is missing, and `python-tests.yml:45–51` wraps the whole leg in `if ls testdata/*.json`; `validate_reference.m:36–42` returns "SKIP" on an empty directory. A deleted or renamed vector becomes 34 skipped tests, not a failure. The April bootstrap rationale no longer applies with 32 committed vectors.
- Five MATLAB toolbox-comparison files (`test_compare{Etfe,MultiTraj,Spa,Spafdr,Welch}.m`) early-`return` when the toolbox is absent; `runAllTests.m` counts them as **passed files with zero cases** and prints no "skipped" tally — on the Octave leg they are always silent.
- `test_sidFreqBTFDR.m:183–185`, `test_sidFreqMap.m:326–328,342–344` turn `warning('off','all')` with no restore on error; a throw leaves warnings off for every later file, and the `lastwarn`-based ETFE Tests 16/17 then fail spuriously (their own comment at `:185`).
- `test_sidModelOrder.m:146–156` catch that converts any display-related failure into a pass.
- `test_ltv_disc_io.py:166` `if hasattr(cost, "__len__") and len(cost) >= 2:` guards the only assertions of its test; `test_compare.py:176–178` asserts `A or B`.
- The 4-case MC campaign (§3.3) is skipped on every run everywhere.

### 3.7 `Verified by:` annotations that over-claim

- `SPEC.md:478` — ETFE `ResponseStd`/`NoiseSpectrumStd = NaN` "held by the response cross-vector's NaN std fields": `reference_siso_etfe.json` stores **no** std field at all.
- `SPEC.md:1827` — cites MATLAB Test 15 (exact-Hessian oracle: real) *and* implicitly the MC Test 13 (§3.3: not a verifier).
- `SPEC.md:1833` — the `none` tag for a trust-region cross-vector is rationalised as "threshold-dependent benefit, brittle"; a vector pins *agreement*, not benefit (see §4.1 F9).
- `SPEC.md:295, :1920` — `none` for degenerate-input cross-vectors "are not in a stored vector": the real blocker is that the Octave validator cannot compare NaN (§4.1 F3). Presenting a validator limitation as a design choice hides the one regime where the ports are known to diverge.
- `SPEC.md:1826` — `sid:singularLbd` "not reachable from normal public-API inputs": it is reachable from the always-tried λ = 1e-3 grid point on rank-deficient data (July minor #16), which is why it was made a diagnostic.

---

## 4. Testing facilities: gaps

### 4.1 Cross-vector infrastructure (ADR-0002) — rule-by-rule

All claims below were established by scripts (inventory of all 113 stored output fields; mutation of copies of the vectors with the Python consumer redirected to them) — not by reading prose.

| ADR-0002 rule | Verdict | Evidence |
|---|---|---|
| 1. Every field has an absolute floor | **Not enforced** | 78 of 113 stored output fields have no `<field>_atol` (default 0); 20 of those contain exact zeros or ≤1e-15 entries (`freqmap_bt.Response_imag` 14 zeros, `siso_btfdr.Response_imag` min 2.4e-16, `internals.DFT_imag` min 2.4e-16). No gate checks that a floor exists; the rule is a review convention. The 1e-16 flips between CI runs (F1) are exactly the S10(e′) class. |
| 2. Stored tolerances authoritative, one default | **Partially** | Defaults agree (rtol 1e-6 / atol 0, `validate_reference.m:320,327`, `test_cross_validation.py:70–71`). But `test_cross_validation.py:367,390` bypass `_tol()` (rtol-only on the complex modulus), the two consumers resolve tolerance keys differently for every complex field other than `Response` (F4), `model_order.n` is `==` in Python and rtol 1e-6 in Octave, and `_tol` silently falls back to the default on a misspelled key. |
| 3. No orphan artifacts | **Enforced today, not mechanically** | All 32 JSONs are loaded by Python and dispatched by `validate_reference.m`; no test asserts "every `reference_*.json` is consumed", so the next orphan is caught only by review. |
| 4. Payloads change only by regeneration; provenance | **Partially, and self-defeating** | Provenance present in all 32, SHA exists. But both gates accept `git_sha: "unknown"` (`generate_reference.m:994–995` fallback; `validate_reference.m:57–58` checks `isfield` only; `test_cross_validation.py:1521` truthiness), never check ancestry, and **F1** below shows the rule's own red flag firing on `main` with nothing able to act on it. |
| 5. Structural gates over outcome tests | **Octave yes, Python no** | `validate_reference.m:282–303` iterates every stored field and fails on a missing one. Python has no such loop (F2). Nothing anywhere would catch a vector both ports match wrongly. |

**F1 (major). The canonical generator is not reproducible, and CI auto-commits the churn unvalidated.** `testdata/README.md:31–37` forbids committing engine-ULP churn and rule 4 calls byte-for-byte regeneration "a checkable invariant". Numerically diffing every regeneration commit (first-hand): `5fd6cbe` (2026-07-29) changed payloads in **21 of 32 files** with *only three documentation commits* in range (`5df968d..5fd6cbe^`: #198's `.md` edits, nothing under `matlab/` or `testdata/`). Deltas: `ltv_io.A` **4.5e-5 rel / 3.9e-4 abs**, `ltv_io.Cost` 3.9e-5, `frozen_of_io` 2e-10, everything else 1e-16…1e-11 (`model_order.SingularValues` tail flips at 2e-15, i.e. 48 % relative on O(ε) entries). The same pinned R2025a runner produces different bytes run-to-run; `tests.yml:173–191` commits them under `[skip ci]`, so neither `cross-validate.yml` nor `python-tests.yml` ran on the new payloads at landing. Of the five regeneration commits since provenance landed with #172, three are provenance-only re-stamps (the SHA is HEAD at generation time, so every push touching `matlab/`/`testdata/` manufactures a commit) — the mechanism produces the churn the policy forbids.

**F2 (major). The Python consumer silently ignores stored fields.** Fields stored but never compared by `test_cross_validation.py`: `btfdr_vecres` {Frequency, Coherence}; `multitraj_bt` {Frequency, Coherence}; `freqmap_welch` {Time, Frequency}; `internals` {DFT_real, DFT_imag}; `freqmap_bt`, `ltv_frozen`, `mimo_bt`, `siso_btfdr`, `siso_etfe`, `timeseries_bt` {Frequency}. Mutation: `mimo_bt.Frequency[5] = 999` and `internals.DFT_real := 0` → **63 passed, no failure**. `sid_dft` therefore has zero Python cross-coverage, and the frequency grid under `SampleTime` is pinned only for three files.

**F3 (major, latent). The Octave validator cannot fail on NaN.** `validate_reference.m:340` is `any(absDiff > thresh)`; `NaN > x` is false, so an actual NaN at any element — or an expected NaN, since `jsondecode` maps `null` → NaN — **passes**. Python (`assert_allclose`) fails number-vs-NaN and passes NaN-vs-NaN. `jsonencode` writes both NaN and Inf as `null`, so the §3.3 `Inf` sentinel is indistinguishable from NaN on disk. No committed payload holds a `null` today, so the hole is latent — but it is precisely the degenerate-input regime (July minor #15: MATLAB `max(NaN,0) = 0` vs NumPy) where the ports are known to diverge, and it is the real reason the three `none` cross-vector tags exist.

**F4 (minor). Complex-field tolerance resolution diverges.** Octave special-cases only `Response_real/_imag → Response_rel/_atol` (`validate_reference.m:289–294`); any other `*_real/_imag` field looks up `<name>_real_rel` / `<name>_imag_rel`. The generator stores `Phi_xx_imag_atol = 1e-14` (`generate_reference.m:164`), so Octave checks `|Δimag| ≤ 1e-14` while Python (`test_cross_validation.py:367`) checks `|Δcomplex| ≤ 1e-10·|Φ_xx|` and never reads the stored atol.

**F5 (minor). Provenance gate weaker than stated.** Beyond the `"unknown"` acceptance: on PRs the MATLAB job regenerates (`tests.yml:167–171`) but neither commits nor `git diff --exit-code`s, so a hand-edited tolerance block that disagrees with the generator is not detected on the PR; on `main` the regeneration silently overwrites it.

**F6 (major). The 1 % gates are unfloored and cannot adjudicate.** `ltv_io` A/B/Cost rtol 1e-2 atol 0 (`generate_reference.m:571`); `ltv_cosmic_varlen` A/B rtol 1e-2 atol 0 (`:808`); `frozen_of_io` Response 1e-2 (`:929`). Mutation: `ltv_io.A ·= 1.005` passes. See §3.1.

**F7 (major). `reference_frozen_of_io` pins a non-converged iterate.** Running the consumer emits `ltv_disc_io: alternating loop did not converge in 50 iterations`. The vector pins the state after `MaxIter` truncation — which is *why* it needs 1 % — and any change to iteration arithmetic or stopping logic moves it. The IO uncertainty fields (`P`, `AStd`) are stored in no vector.

**F8 (minor).** (a) `test_cross_validation.py:847` truncates actual singular values to `len(expected)` — extra SVs pass in Python, fail in Octave. (b) Six fields rely on the undocumented default rtol (five `Frequency`, `model_order.n`). (c) `SPEC.md:478` claims a NaN-std cross-vector that does not exist.

**F9 (major, by consequence).** The three `none` tags (`SPEC.md:295, :1833, :1920`): the degenerate-input rationale does not hold (nothing prevents storing a constant-`u` input; the blocker is F3); the trust-region rationale partially holds (benefit is threshold-dependent, but a vector pins agreement at a fixed schedule and the 1e-2 precedent already exists for the same loop; a vector pinning `Cost` history and `Iterations` is feasible).

**F10 (minor).** The loader is a module-level `TESTDATA` constant; monkeypatching it from a scratch module works but there is no fixture / `--testdata` option, so "does the gate bite" tests are not first-class.

**What would catch "both ports wrong":** nothing today. It needs spec-derived oracle assertions on the *payloads* themselves, run as a third job over `testdata/*.json` without any port: analytic-plant vectors where true `G(e^{jω})` is known (asserted within the §3 variance bound), invariants (coherence ∈ [0,1], PSD ≥ 0, `Φ_xz = conj(Φ_zx)`, `util_msd.Ad == expm(A·Ts)`, `P(k)` symmetric PD, `Cost` monotone for IO). A Julia port would add independence, not an oracle.

### 4.2 CI gates

- **Required checks are lint only.** `.github/rulesets/main.json:27–34` requires `MATLAB/Octave lint` and `Python lint`; `Python Tests`, `Octave/MATLAB tests & examples`, `Cross-Language Validation` and `Docs` are advisory. `tests.yml:3–10` describes the test job as "the required check" — it is not. ADR-0005 documents this as a first step (path-filtered checks wedge PRs); `tests.yml` has since been made path-unfiltered precisely so it can be required, and the ruleset was not updated. #199 (docs) is the only open follow-up. A red test job does not block merge; `RepositoryRole 5` bypasses always.
- **`[skip ci]` on reference regeneration** (`tests.yml:184`) means regenerated vectors are never validated by any consumer at landing (F1), and tag pushes on a bot commit silently produce no release (#204, open). The MATLAB job never runs `validate_reference.m` itself.
- `cross-validate.yml` path filters: a `spec/`-only change triggers nothing, by design; `python-tests.yml` also runs the Python consumer on 3.10–3.14, so the matrix is covered — by the other workflow.
- **Minimum dependencies are never tested** (`python-tests.yml:40` installs latest NumPy/SciPy; the floor `numpy >= 1.22, scipy >= 1.8` is a claim). July §5 named this; still open.
- Octave is `snap install octave` on an unpinned channel (`tests.yml:66–70`, `cross-validate.yml:72–74`); only the version is echoed. The 11+ floor is echoed, not enforced: a runner shipping apt Octave 8.4 would run below the floor unnoticed. (The July S10(e′) 8.4 failure itself no longer recurs now that `reference_ltv_cosmic` has an atol floor — verified under apt Octave 8.4, plan §10.1.)
- Fork PRs: the MATLAB job mints an App token first (`tests.yml:104–109`); without the secrets the step fails before any test runs.
- Both header checkers (`.github/scripts/check_headers.py`, `check_python_headers.py`) validate section presence and order only; neither validates that a cited `SPEC.md §` exists, and both accept "not yet in SPEC.md" (§2.2 #33). `scripts/local-ci` does not run `scripts/check_python_api_pages.py` (the docs API gate), so a missing API stub surfaces only in the docs workflow.
- `pyproject.toml [tool.pytest.ini_options]` has no `filterwarnings = error`: nine warnings escape the suite, including the `notConverged` warning from the vacuous convergence test (§3.4) and the tuning fallback warning in a test *named* `test_fallback_strict_threshold` that never asserts it.

### 4.3 MATLAB/Octave suite (static reading; 41 files, ~464 numbered cases)

**Runner** (`runAllTests.m:62–87`): a file that throws before its first assertion is counted as one failed case and the runner exits non-zero — correct. But `run()` executes scripts in the runner's own workspace, so any test can read or overwrite `runner__nPassed` / `runner__*`; variables, RNG state (two `test_sidDetrend.m` tests have no `rng`) and warning state leak across files; `test_template.m` and `example_template.m` are discovered and counted; there is no per-case timeout.

**Warning identifiers** (21 `warning('sid:…')` sites): asserted by id in a test — `constantInput`, `singularPhiU` (ETFE path only), `deadInputChannel` (BT only), `stabilized`, `singularLbd` (`sidLTVcosmicSolve` only). **Never asserted:** `windowReduced` (BT: suppressed at `test_reviewFixes.m:319`, the clamp checked as `<= 10` not `== 10`; BTFDR: none — Python has `test_window_reduced_warns`), `shortData` (§10.1 N < 10), `trimmedTrajectories` (§2.3), `notConverged` (Test 24 must fire it; never read), `noConsistentLambda`, `illConditioned`, `detrendOrderReduced`, `allZeroSV`, `preconditionUnsupported` (always suppressed), the BTFDR-path `deadInputChannel`/`singularPhiU`, `singularLbd` in `sidLTVblkTriSolve.m:83`.

**Error identifiers never asserted:** BT/ETFE/BTFDR `badTs`, `badFreqs`; BT `badWindowSize` (prefix only); every `sidSpectrogram` id (`complexData`, `nonFinite`, `invalidWindowLength`, `invalidSampleTime`, `invalidNFFT`, `tooFewSegments`, `windowSizeMismatch`, `invalidWindow`); every `sidFreqMap` id (11); `sidDetrend` (4); `sidResidual`/`sidCompare` `badModel`; `sidModelOrder` `badHorizon`, `badThreshold`, `tooShort`; `sidLTIfreqIO` `tooShort`/`badInput`/`dimMismatch` (prefix only), `orderExceedsRank` (either-or); `sidLTVdiscIO` missing-`Lambda`, `R` mismatch, cell-path mismatch; **all 16** `sidLTVStateEst` sites; `sidLTVdiscFrozen` `badTimeSteps`; `sidLTVdiscTune` `badMethod`; `sidLTVdisc` wrong-length λ vector; `sidValidateData` `trajMismatch` and the cell path; `sidMapPlot` `noCoherence`.

**Options never exercised:** `sidResidual 'MaxLag'`; `sidLTVdiscTune 'CoherenceThreshold'`, `'Algorithm'`, frequency-mode `'Precondition'`; `sidLTVdisc 'LambdaGrid'` without `'auto'`, explicit `'NoiseCov','estimate'`; `sidLTVdiscIO` per-step `'Lambda'` vector; `sidFreqMap 'Window'` numeric vector, `'Frequencies'` with `'welch'`; `sidBodePlot/sidSpectrumPlot 'ShowConfidence'`, `'Color'`, `'LineWidth'`; `sidSpectrogramPlot`/`sidMapPlot 'Axes'`.

**Input-shape modes:** 3-D with `L = 1` — BT only; **variable-length cell input to `sidFreqMap` has no test** (the §6.2 per-segment filtering at `sidFreqMap.m:95,306–320` is exercised only by Python's five cases); `sidSpectrogram` has no `iscell` path at all (§7.2 requires one — see §5); zero input (§10.3 `realmin` branch), `u = y`, complex `u` — untested (Python has all three). **Every LTV test in both ports uses `q = 1`**: the `B = permute(C(p+1:end,:,:),[2 1 3])` extraction (`sidLTVdisc.m:171,454`, `sidLTVdiscIO.m:515`) and every downstream consumer (`AStd`/`BStd` indexing, frozen TF, `sidCompare`, `sidLTVStateEst`) is never checked numerically with `q ≥ 2`, in either port, and no cross-vector has `q ≥ 2` (§4.5). A first-hand Python probe at `q = 2, 3` found the stack correct to 1e-14; the MATLAB side is unverified.

**Numerical rules with no verifier:** §2.6/§3.6 MIMO `Ĝ = Φ_yu Φ_u⁻¹` and `σ_Gij` on a known plant (the only MIMO numeric check is the toolbox-gated `test_compareSpafdr.m` Test 6 at 15 %; `test_crossMethod.m:153–154` builds a static 2×2 gain plant whose Ĝ must equal `[1 0.3; 0.5 1]` and asserts only `size`); §10.2 clamps; §8.9.3 noise-covariance/DoF formulas (structure only); §8.4.2 corner; multi-trajectory ensemble averaging (free oracle unused: L identical copies of one trajectory must reproduce the single-trajectory estimate exactly); §8.11 consistency score; `sidModelOrder` MIMO (Test 4 accepts n ∈ [2,4]); `sidResidual`/`sidCompare` on an IO (`H ≠ I`) model; MIMO Bode/spectrum plots.

**Cross-port asymmetry** (the ports are supposed to be verified independently against the same contract, so one-sided coverage is itself a gap): MATLAB-only — `test_sidLTVdiscUncertainty.m` (covariance modes, user `NoiseCov`, `badNoiseCov`, var-len + uncertainty, DoF, isotropic), `test_sidLTVdiscVarLen.m` (no Python list-input unit test for `ltv_disc`; cross-vector at 1e-2 only), the FD-Jacobian check, the ±5 % 1/√L ratio, `test_validation.m`, `test_crossMethod.m`, the five toolbox oracles, IO Tests 9/17/18/20/21/25, `sidLTVStateEst` 6/8/9/10, `sidLTIfreqIO` 9–12, `sidModelOrder` 11–14. Python-only — `TestFreqMapVariableLength` (5), `test_mimo_freq_map`, Welch-vs-scipy bit-exact, `test_welch_is_twice_bt_off_nyquist`, zero-input / `u = y` / complex-`u`, BTFDR partial degeneracy and `windowReduced`, `test_stabilize_defective_is_bounded`, `test_order_exceeds_rank_errors`, `test_accuracy_partial_observation`, `test_mimo_block_hankel`, `test_time_series_compare`, `test_degenerate_window_raises`, `test_mimo_noise_spectrum`, IO `test_uncertainty_fields`.

### 4.4 Python suite (28 files; 415 passed, 4 skipped; 86 % line coverage)

Per-module coverage, lowest first, with the uncovered ranges that matter:

```
lti_freq_io.py        71%  233-292  entire list / var-len path; 558-562 H-basis warning; 609-621
residual.py           72%  384-430  plot=True; 297-305
freq_etfe.py          74%  301-320  multi-output time series; 416-427, 479-488 MIMO multi-traj
ltv_disc_io.py        75%  472-525  entire var-len path; 570-577, 586-596 R / covariance_mode validation
compare.py            76%  318-338  plot=True
spectrogram_plot.py   76%  105-114  channel=
model_order.py        78%  258-271  all-SV-zero warning
map_plot.py           80%  error codes at 110/117/133/149/171
ltv_disc.py           83%  543-555  user noise_cov; 389-456 dim errors
estimate_noise_cov.py 85%  105, 123-125, 134
validate_data.py      86%  105-170  list-path errors; 186/286 N<10 warning; 215 trim warning
freq_bt.py            93%  199-209  bad_ts / window reduction / bad_freqs
```

**Error codes:** 49 distinct `SidError` codes are raised; **14 are asserted** by a test. Never asserted include `dim_mismatch` (39 raise sites), `bad_ts`, `bad_freqs`, `bad_window_size`, `bad_noise_cov`, `bad_cov_mode`, `bad_algorithm`, `bad_method`, `bad_model`, `bad_horizon`, `bad_threshold`, `bad_time_steps`, `bad_order`, `bad_segment_length`, `traj_mismatch`, `no_active_traj`, `too_few_segments`, `too_few_sub_segments`, `sub_segment_too_long`, `invalid_*` (window/window_size/window_length/sample_time/sub_overlap/sub_segment_length/channel/result/plot_type), `degenerate_window`, `window_size_mismatch`, `no_response`, `no_coherence`. Several tests use `pytest.raises(Exception)` or a five-exception tuple without checking the code (`test_plotting.py:117–119,270–272`, `test_residual.py:171–181`, `test_lti_freq_io.py:151–164`).

**Warnings:** 8 of 21 `warnings.warn` sites asserted. Unasserted: `freq_bt.py:203` window reduction (uncovered), `validate_data.py:186,286` N < 10 (uncovered), `:215` trimming (only "not trimmed" is asserted), `ltv_blk_tri_solve.py:116`, `lti_freq_io.py:558`, `ltv_disc_io.py:371` (`notConverged`), `ltv_disc_tune.py:413`, `model_order.py:259,271`, `detrend.py:197`.

**Options never called:** `ltv_disc(noise_cov=<matrix>)`, `covariance_mode='isotropic'` (`'full'` appears once as a fixture with no full ≠ diagonal assertion); `ltv_disc_io(R=, tolerance=, covariance_mode=, uncertainty=False)`; `ltv_disc_tune(algorithm=, coherence_threshold=)`; `compare(initial_state=, plot=True)`; `residual(plot=True)`; `bode_plot(color=, line_width=)` and any MIMO result (`bode_plot.py:124–135`); `spectrum_plot(frequency_unit='Hz', ax=)`; `map_plot(frequency_unit=, clim=, ax=)`, SISO `plot_type='phase'`; `spectrogram_plot(channel=, clim=, ax=)`.

**Input modes:** list/var-len input to `freq_bt`, `freq_btfdr`, `freq_etfe` — zero unit calls (the §2.3 trim-to-shortest rule and `trimmedTrajectories` warning are untested); `ltv_disc` var-len — unit-untested, cross-validated at 1e-2 only; `ltv_disc_io` and `lti_freq_io` var-len — whole path uncovered; 3-D `L = 1` for `freq_btfdr`/`freq_etfe`/`ltv_disc` — none; ETFE MIMO multi-trajectory — none; `sid_cov` ensemble 1/(L·N) path (§2.3) — no 3-D case in `test_cov.py`.

**Numerics with no test:** MIMO `Ĝ = Φ_yu Φ_u⁻¹` on a known plant for any estimator (`test_freq_bt.py:114`, `test_freq_btfdr.py:169`, `test_freq_etfe.py:203`, `test_freq_map.py:620` are shape-only); §2.7 MIMO PSD eigenvalue clamp and NaN preservation (`degenerate.py:229` uncovered); §10.2 γ² > 1 clamp (no > 1 case constructed); §8.9 DoF formula and `'estimate'` vs supplied Σ (`a_std/b_std/noise_cov/DoF` are MATLAB-pinned only); §8.12 `R` weighting; the N = 2 OLS closed form; `ltv_disc_tune` frequency method (all four tests are class 3/5 — no correctness oracle for §8.11.2).

**Infrastructure:** dead `conftest` fixtures (§2.2 #35); the LTI simulation loop is copied 34 times across 7 test files; the two trust-region tests re-run the identical 11–14 s fit; no pytest markers (campaign gating is an env var; CI selects the cross-validation leg by `-k` name match); `test_plotting.py:17–19` imports matplotlib unconditionally (collection error, not skip, without it); nbmake runs all 11 notebooks on every push for five Python versions (execution only — no output assertions).

`test_util_msd.py` **never asserts the `EXAMPLES.md` §1 normative discretization test points** (Plant A/B ≥ 1e-9, Plant C ≥ 1e-6, Plant E impulse ≥ 1e-6); its SDOF case uses `c = 0.5` (Plant A is `c = 2.0`) against `scipy.linalg.expm`, and the §2.2.4 "bit-identical" LTI collapse is checked at `atol = 1e-14` instead of `array_equal`. The only spec-parameter check is `reference_test_msd.json` — MATLAB agreement, not spec conformance (principle 5). First-hand: `util_msd` reproduces Plant A to 5e-11 and Plant E to every printed digit; Plant C as printed in the spec carries only four significant digits on its smallest entry (`0.00001056`), so a literal test needs per-entry tolerances or the spec's digits extended.

### 4.5 Functionality not fully covered — consolidated matrix

Legend: ✅ numeric test against an oracle · ◐ shape/smoke/loose only · ✗ none · 1e-2 = cross-vector at rtol 1e-2.

| Feature (spec §) | unit(M) | unit(Py) | cross-vector |
|---|---|---|---|
| MIMO `Ĝ`, `Φ̂_v` numerics, BT (§2.6–2.7) | ◐ (toolbox-gated, 15 %) | ◐ | ✅ `mimo_bt` (2×2) |
| MIMO ETFE / BTFDR / `freq_map` | ◐ | ◐ | ✗ |
| MIMO variance diag. approx. (§3.6) | ✗ | ✗ | ✗ |
| Multi-trajectory ETFE / BTFDR / `freq_map` / spectrogram numerics | ◐ | ◐ | ✗ (only `multitraj_bt`) |
| Var-len list input, BT/ETFE/BTFDR trim + warning (§2.3) | ◐ (trim, no warning) | ✗ | ✗ |
| Var-len `freq_map` per-segment filtering (§6.2) | ✗ | ✅ | ✗ |
| Cell input to `sidSpectrogram` (§7.2) | ✗ (unsupported) | — | ✗ |
| Custom `Frequencies`, `SampleTime ≠ 1` in freq estimators | ◐ | ◐ | ✗ |
| Welch `SubSegmentLength`/`SubOverlap`/`Window`/`NFFT` | ◐ | ✅ (scipy) | ✗ |
| Degenerate inputs (§10.2–10.3) | ✅ | ✅ | ✗ (`none`, F3) |
| §2.7 NaN-preserving clamp, MIMO PSD projection | ✗ | ✗ | ✗ |
| COSMIC `q ≥ 2` (any LTV function) | ✗ | ✗ (probe only) | ✗ |
| `Lambda = 'auto'` L-curve corner (§8.4.2) | ◐ interior + 15 % | ◐ | ✗ |
| Per-step λ vector | ◐ | ◐ | ✗ |
| Uncertainty: `CovarianceMode` full/isotropic, user `NoiseCov`, DoF (§8.9.3) | ◐ structure | ✗ | ✅ `cosmic_uncertainty` (diag only) |
| Var-len COSMIC (§8.8) | ✅ dense oracle | ✗ | 1e-2 |
| `sidLTVdiscTune` validation loss formula (§8.4.3) | ◐ | ◐ | ✅ `ltv_tune` |
| `sidLTVdiscTune` frequency method, consistency score (§8.11.2) | ◐ | ◐ (class 3/5) | ✗ |
| Output-COSMIC converged case, `R` weighting, per-step λ, var-len | ◐ | ✗ | 1e-2 (non-converged, F7) |
| Output-COSMIC uncertainty outputs | ◐ | ◐ | ✗ |
| Trust region (§8.12.4) normative clauses | ◐ pin | ◐ pin | ✗ (`none`) |
| `sidLTVStateEst` user `R`, `Q` | ◐ | ◐ | ✗ |
| `sidLTIfreqIO` `Horizon`, `MaxStabilize`, stabilization path | ◐ | ✅ defective | ✗ |
| `sidModelOrder` `Threshold`, `Horizon`, MIMO | ◐ (n ∈ [1,12], [2,4]) | ✅ MIMO | ✗ |
| `sidDetrend` `Order` 0/2, `SegmentLength`, multi-traj | ◐ loose | ◐ loose | ✗ |
| `sidResidual` pass/fail flags, `MaxLag`, time-series, state-space, multi-traj | ✗ flags | ✗ flags | ✗ |
| `sidCompare` freq-domain model, `InitialState`, multi-traj | ◐ | ◐ | ✗ |
| Plotting options (`ShowConfidence`, `Color`, `LineWidth`, `Axes`, `clim`, `channel`), MIMO Bode | ✗ | ✗ | — |
| `util_msd*` vs `EXAMPLES.md` §1 test points | ✗ | ✗ | ✅ (agreement only) |

---

## 5. Bugs and spec-conformance defects found in this review

Each item names the port(s), the spec rule, and — for Python — a reproduction run at the anchor commit. MATLAB items are static readings.

### 5.1 Frequency-domain half (§1–§7, §9–§11)

The core math is correct: biased covariance, Hann / `C_W`, the K-strided FFT buffer, the direct DFT with `Rneg`, SISO/MIMO `Ĝ`, every §3 variance formula including `L·N` and the MIMO diagonal approximation, ETFE conventions, the BTFDR ≡ BT oracle, segmentation and time vectors, Welch scaling including the Nyquist bin (bit-exact vs `scipy.signal.welch`, rel. err. 7e-16), spectrogram PSD, and the ensemble `Complex` field were all re-verified in Python to machine precision against independent computations, and the MATLAB code implements the same formulas. The §10.3 excitation check, the §2.6/§2.7 regularisation and the §3.3 `Inf` sentinel are genuinely shared through `degenerate.py` / `sidRegularizeResponse.m`. The defects are at the shape/edge contract and in cross-port asymmetry.

| # | Sev. | Defect | Where | Rule |
|---|---|---|---|---|
| F1 | **critical** | **`freq_map` crashes on MISO data (`n_y = 1`, `n_u ≥ 2`), both algorithms** — the MIMO branch allocates `(nf, K, 1, 1)` but `_store_segment` branches on `ny == 1` and ravels. Repro: `freq_map(randn(1000), randn(1000, 2), segment_length=500)` → `ValueError: could not broadcast (128,) into (128,1,1)`. MATLAB works (trailing singletons drop). No MISO `freq_map` test in either port. | Py `freq_map.py:573–589, 603–605` | §6.1, §6.8 |
| F2 | **critical** | **3-D single-trajectory `(N, ch, 1)` input still crashes `freq_etfe` and Welch `freq_map`** (residual of July C3; the fix covered bt/map-bt/spectrogram only). Repro: `freq_etfe(y[:,None,None], u[:,None,None])` → `TypeError: only 0-dimensional arrays…`; `freq_map(…, algorithm='welch')` → `ValueError: non-broadcastable output operand (128,) … (128,128)`. | Py `dft.py:77–91`, `freq_map.py:108–118` | §1 |
| F3 | **major** | **Single-trajectory MIMO/MISO ETFE (`n_u ≥ 2`) is 100 % NaN by construction, both ports.** `Φ̂_u(ω) = U(ω)U(ω)ᴴ` is rank 1 for `L = 1`, so the §2.6 `cond > 1/ε` guard fires at every bin and the "near-singular at *some* frequencies" warning misleads. Repro: `freq_etfe(randn(N,2), randn(N,2))` → NaN fraction 1.00; `L = 3` → 0.00. §4.1 defines only the SISO ratio and the multi-trajectory H1 estimator; neither port errors, documents, or requires `L ≥ n_u`. The existing "MIMO" ETFE tests are SIMO (`n_u = 1`), shape-only. | `freq_etfe.py:407–413`, `sidFreqETFE.m:264–271` | §4.1 |
| F4 | **major** | **ETFE noise spectrum at degenerate frequencies: three conventions.** MATLAB `max(PhiV, 0)` turns NaN into 0 — the exact hazard §2.7 names, and §2.6 says the clamp is "shared verbatim by … `sidFreqETFE`"; Python preserves NaN; BT/BTFDR/Welch substitute `Φ̂_y` with coherence 0 (spec-silent — §2.7 says NaN "is produced"). Constant-`u` ETFE: `Φ̂_v` all-NaN (Py) vs all-0 (M). | `sidFreqETFE.m:253`, `freq_etfe.py:396`, `degenerate.py:318,354`, `sidRegularizeResponse.m:86,113` | §2.7 |
| F5 | major | **Python warnings carry no stable identifier.** MATLAB uses `sid:windowReduced`, `sid:trimmedTrajectories`, `sid:shortData`, `sid:constantInput`, `sid:singularPhiU`, `sid:deadInputChannel`; Python emits bare `UserWarning` text, no `SidWarning` class, no code attribute, nothing in `python/CONTRIBUTING.md`; tests match on free text. Error codes, by contrast, are stable (`SidError.code` ↔ `sid:*`). Principle 8 is unmet on the Python side; the spec names identifiers Python cannot be checked against. | `freq_bt.py:203`, `validate_data.py:186,215,286`, `degenerate.py:328,362`, `freq_btfdr.py:253`, `freq_etfe.py:67–77`, … | §10, principle 8 |
| F6 | major | `bode_plot` / `spectrum_plot` have no `show_confidence` (July #32, Python side; MATLAB has it undocumented). | `bode_plot.py:24–32`, `spectrum_plot.py:23–31` | §11.3 |
| F7 | major | `map_plot` rejects spectrogram results (`SidError(invalid_result)`); MATLAB implements §7.5. | `map_plot.py:109` | §7.5 |
| F8 | major | **`sidSpectrogram` cell/list input (mandated by §2.3/§7.2):** MATLAB rejects a cell with a misleading `sid:complexData` (July #14, unchanged); Python **silently** treats a list of 300 trajectories of 600 samples as `N = 300` samples × 600 channels with `num_trajectories = 1`. | `sidSpectrogram.m:92–104`, `spectrogram.py:243–257` | §7.2 |
| F9 | major | **Welch / spectrogram divide by `S₁ = 0`** for a symmetric Hann with `L_sub` or `L = 2` in MATLAB (silent NaN/Inf); Python fixed (July #10, Python only). Both ports: `spectrogram(window_length=1, window='hann')` → `0/0` window → silent all-NaN (the `S₁ ≤ 0` guard sees NaN, not 0). | `sidFreqMap.m:430`, `sidSpectrogram.m:155`; `spectrogram.py:65`, `sidSpectrogram.m:240` | §6.5, §7.3 |
| F10 | minor | ETFE boxcar smoothing propagates one NaN bin to `S − 1` neighbours (July #1, both). Repro: one NaN bin, `S = 3` → `[1, 1, nan, nan, nan]`. | `freq_etfe.py:50–54`, `sidFreqETFE.m:388–392` | §4.2 |
| F11 | minor | BT time-series `Φ̂_y` is not clamped (July #3, both): `freq_bt(sin(0.7t))` → min −0.18 for BT, BTFDR and map-bt; and SPEC §2.3/§6.6 still *claim* the biased estimator guarantees non-negativity (true only for Bartlett). | `freq_bt.py`, `SPEC.md:115, :665` | §2.3, §2.7 |
| F12 | minor | `y` 2-D + `u` 3-D (`L > 1`) passes validation and dies in `cov` (July #4, both): raw `ValueError: matmul…` (Py, all five estimators); `sidCov.m:69` dimension error (M). | `cov.py:90`, `sidCov.m:69` | §10.1 |
| F13 | minor | list-`y` + `u = None/[]` (time-series multi-trajectory) → `SidError(size_mismatch)` in Python for `freq_bt`/`freq_btfdr`/`freq_etfe`/uniform `freq_map` (July #5); MATLAB handles `[]`. | `validate_data.py:232` | §1 |
| F14 | minor | ETFE `WindowSize` reported as `N` (spec equates ETFE to BT with `M = N − 1`; July #6, both, spec silent). | `freq_etfe.py:513`, `sidFreqETFE.m:375` | §4.1, §9 |
| F15 | minor | §3.4/§3.5 variance formulas have no `L`; both ports use `N_eff = L·N` (July #8, half-fixed: §3.3 has it). | `uncertainty.py:84,90`, `sidUncertainty.m:68,76` | §3.4–3.5 |
| F16 | minor | Welch MIMO `σ_G` is all-NaN in both ports and the spec is silent (July #11; the PSD projection half is fixed). | `freq_map.py:199–200`, `sidFreqMap.m:553–554` | §6.5 |
| F17 | minor | BTFDR with `N ∈ {2, 3}` → `M_k = 1` after the `⌊N/2⌋` clamp, contradicting "M < 2 is invalid"; `hann_win(1)` runs (both). | `freq_btfdr.py:258`, `sidFreqBTFDR.m:158` | §5.2, §10.1 |
| F18 | minor | Spec-silent formulas/thresholds both ports agree on: Welch SISO `σ_G = |Ĝ|√((1−γ²)/(γ²ν))` (§6.5 gives no `Ĝ` variance; Bendat–Piersol's complex-total form is √2 larger); `ν = max(2, 1.8J)` floor; `SegmentLength ≥ 4`; MIMO condition metric is 2-norm `cond` in Python vs 1-norm `rcond` estimate in MATLAB (can disagree on borderline bins). | `freq_map.py:163,186,402`, `sidFreqMap.m:125,524,543`, `degenerate.py:351`, `sidRegularizeResponse.m:111` | §6.5, §6.2, §2.6 |
| F19 | minor | Default `M = min(⌊N/10⌋, 30) = 1` for `10 ≤ N ≤ 19` raises `bad_window_size` in both ports — literal §2.1/§10.1, but §5.4 claims BT has a short-data floor of 2 that it does not have. | `freq_bt.py:191,200`, `sidFreqBT.m:122,132` | §2.1, §5.4 |
| F20 | nit | Welch map title: MATLAB `sprintf('%d', [])` truncates "M="; Python prints `M=None` (July #34, both). `freq_map` re-emits the degenerate-input warning once per segment. | `sidMapPlot.m:173`, `map_plot.py:208` | — |

Fixed and verified from the July list in this area: C2, C4, C5 (MIMO), S5, S6, S7, minor #2, #10 (Python), #13, #15 (BT/BTFDR/Welch — not ETFE, F4).

### 5.2 COSMIC / LTV half (§8, `spec/cosmic/*`)

The core math is correct. A dense brute-force solve of the full block-tridiagonal normal equations (N = 6, p = 2, q = 1, L = 3, per-step λ) matches `ltv_disc` on `C` (2e-16), the cost triple (3e-17), `P(k) = [(N·A_scaled)⁻¹]_kk` (4e-16), the hat-trace `ν` (2e-15), `Σ̂` in all three modes (2e-16) and the `AStd`/`BStd` index map (4e-16; the transposed map differs by 0.52, so the test is discriminating); variable-length lists `[6, 4, 3]` to 1e-14; `ltv_state_est` vs a dense `J_state` minimiser with non-diagonal `Q`, `R ≠ I`, uniform and ragged, to 5e-13; `ltv_disc_frozen` std vs a central-difference Jacobian `J(Σ⊗P)Jᴴ` to 5e-11 including `H ≠ I`; the real-Schur stabilization preserves arguments and bounds a Jordan block; the trust region runs exactly the normative budget `50·(⌈log₂10⁶⌉ + 2) = 1100` iterations and lowers J 9× on the hard case; μ = 0 is monotone. The defects are in the model-order rule, validation, diagnostics, the convergence test, and cross-port edge behaviour.

| # | Sev. | Defect | Where | Rule |
|---|---|---|---|---|
| C1 | **major** | **The gap rule of §8.12.12 / ADR-0004 cannot return `n` for an exactly rank-`n` Hankel, and its answer depends on the DC-extrapolation artifact.** Step 5b sets `K = min(L − 1, ⌊m/2⌋)`; on exact data `L = n`, so `k = n` is excluded and the search returns ≤ n − 1 unless an artifact singular value lifts `L`. First-hand on exact samples: `(z+0.5)/(z²−1.2728z+0.81)` → `n = 2` only because σ₃ = 1.3e-5 (artifact) exists; the 3rd-order plant `(z+.3)(z−.2)/((z−.5)(z−.7)(z−.9))` → **`n = 2` at `nf = 512`** (artifact 1.6e-2 vs weakest mode 0.27), `n = 3` at `nf = 4096`; the pass's 2×2 MIMO `n = 3` example returned 4 (not reproduced on a different MIMO plant, which returned 3 — the outcome is plant- and grid-dependent, which is the point). The spec's rationale ("the lag-1 reconstruction already suppresses artifacts to near-ε", `SPEC.md:1662`) is false; `test_sidModelOrder.m` Test 15 and `test_finite_zero_second_order` pass *because* of the artifact, and `test_mimo_block_hankel` sidesteps the gap method. The cross-vector pins the two ports to the same wrong rule (principle 5). ADR-0004's premise needs revisiting — a **[decision]**. | `model_order.py:284–293`, `sidModelOrder.m:242–251` | §8.12.12 |
| C2 | **major** | **Convergence test is absolute for J < 1** (July #22, unfixed, both): `|ΔJ|/max(|J_prev|, 1)` instead of §8.12.3's `|ΔJ|/|J|`. First-hand: data scaled so `J ≈ 4e-6` → "converged" after 3 iterations with relative change **0.19** (tolerance 1e-6), no warning. The in-code comment defends it as needed "on ordinary data" — a contract change encoded silently. | `ltv_disc_io.py:357`, `sidLTVdiscIO.m:290` | §8.12.3 |
| C3 | **major** | **Exactly singular pivot crashes Python from the default public API** (July #16): `ltv_disc(X, U)` with collinear states and zero input → `LinAlgError: Singular matrix` at the always-tried λ = 1e-3 grid point; MATLAB warns and returns Inf/NaN. Neither matches §8.3.4 ("still returns a result"), and `SPEC.md:1826`'s "not reachable from normal public-API inputs" is false. | `ltv_cosmic_solve.py:104,118`, `sidLTVcosmicSolve.m:70–78` | §8.3.4 |
| C4 | **major** | **Validation loss divides by `N + 1`** (July #18, unfixed, both): `np.mean` over the `N + 1` rows including the zero `k = 0` row, so an error of 3 at step `N` alone gives `√(9/(N+1))` where §8.4.3 says `√(9/N)`. `reference_ltv_tune` pins `AllLosses` across ports — the shared drift is locked in by the gate (principle 5). | `ltv_disc_tune.py:479`, `sidLTVdiscTune.m:346` | §8.4.3 |
| C5 | major | **The guarded final μ = 0 pass can cap out silently.** §8.12.4 step 5 "guarantees the returned iterate is a fixed point of the μ = 0 alternation" and the base alternation must warn on the cap; with `trust_region = 1` all 22 stages capped at 50, the final pass capped, no warning, and one more μ = 0 sweep from the returned iterate lowers J a further 1.4e-3 relative — not a fixed point. | `ltv_disc_io.py:393–395`, `sidLTVdiscIO.m:203–211` | §8.12.4 |
| C6 | major | **Frequency-mode tuning evaluates the frozen TF one step apart across ports:** MATLAB uses `round(Time/Ts)` as a **1-based** `TimeSteps`, Python as **0-based**; §6.7's `t_k/Ts` is a 0-based sample position. No cross-vector covers frequency mode, so the asymmetry is invisible to the gate. | `sidLTVdiscTune.m:236`, `ltv_disc_tune.py:358` | §8.11, §6.7 |
| C7 | minor | `lambda_` string other than `'auto'`: Python raw `ValueError`; MATLAB a 1-char string is `isscalar` and is **silently coerced to its ASCII code** (`'x'` → λ = 120). `lambda_grid`: negative values accepted, < 3 points returns `grid[0]` silently, empty → raw `ValueError`, an unsorted grid yields a different corner (1.3e-2 vs 1.8e14 verified). (July #17, both.) | `ltv_disc.py:258–260,617–620`, `sidLTVdisc.m:144–150,358,439–443` | §8.4.2 |
| C8 | minor | `_frequency_tune` derives `N` from the first trajectory only (July #19, both); validation trajectories longer than training → raw `IndexError`; wrong `p` → raw `ValueError`. The #189 1-D/2-D acceptance *is* fixed and tested in both ports. | `ltv_disc_tune.py:320`, `sidLTVdiscTune.m:189` | §8.4.3 |
| C9 | minor | User `noise_cov` aliased (`shares_memory` True) and never checked symmetric/PSD (July #20, both): `diag(1, −1)` accepted → NaN `a_std` with only a NumPy `RuntimeWarning`. `R` never SPD-checked in `ltv_disc_io` (−1 accepted; singular → raw `LinAlgError`). | `ltv_disc.py:517`, `sidLTVdisc.m:377–389`; `ltv_disc_io.py:586–591`, `sidLTVdiscIO.m:485–487` | §8.9.3, §8.12.8 |
| C10 | minor | Spec-silent second DoF fallback `max(total_obs, 1)`; `extract_std` takes `sqrt` of tiny negative `P` diagonals → NaN stds (verified at λ = 1e15 on rank-deficient data). (July #21, both.) | `estimate_noise_cov.py:122–125`, `sidEstimateNoiseCov.m:97–102`; `extract_std.py:86–89`, `sidExtractStd.m:59–63` | §8.9.3–8.9.4 |
| C11 | minor | Fast path (`rank(H) = n`) returns `iterations = 0` with `cost` of length 1, against §8.12.9's `(n_iter × 1)` (July #22, both). | `ltv_disc_io.py:276–277`, `sidLTVdiscIO.m:145–147` | §8.12.9 |
| C12 | minor | `ltv_state_est` rebinds the caller's list elements (July #23); 2-D `A` → raw `IndexError` in Python while MATLAB runs it as `N = 1`; list `A` → `AttributeError`. | `ltv_state_est.py:269,310–311` | §8.12.13 |
| C13 | minor | `blk_tri_solve` conditioning check lags use (block 0 used before checked, last block never checked; full-SVD `cond` per block). (July #24, both.) | `ltv_blk_tri_solve.py:111–114`, `sidLTVblkTriSolve.m:78–81` | §8.12.14 |
| C14 | minor | `Horizon` silently clamped in `lti_freq_io` and `model_order` (July #25); orphan expression `nf - 1` at `lti_freq_io.py:392`; the time-series fallback of `model_order` (a spectrum's Hankel rank as plant order) remains spec-silent (July #26). `'Plot'` still missing in Python (§2.2 #26). | `lti_freq_io.py:163–164,392`, `model_order.py:147–151,221–222` | §8.12.12, §8.13.2 |
| C15 | minor | Spec-silent behaviour both ports agree on: `lti_freq_io` list input *discards* trajectories shorter than `⌈2N_max/3⌉` and trims the rest (§8.13.1 says "average across trajectories"); frequency-mode tuning adds `segment_length = min(N/4, 256)`, a `coherence_threshold = 0.3` mask, `≥` instead of §8.11.2's `>`, and an argmax fallback with warning — none in §8.11. | `lti_freq_io.py:273–292`, `sidLTIfreqIO.m:182–202`; `ltv_disc_tune.py`, `sidLTVdiscTune.m` | §8.13.1, §8.11 |
| C16 | minor | Cross-port edge asymmetries: `trust_region='bogus'` raw `ValueError` (Py) vs `sid:badInput` (M); 1-D `Y` to `ltv_disc_io` raw `IndexError` (Py) vs fine (M); every rcond threshold is a 1-norm LU estimate in MATLAB vs `1/np.linalg.cond` (2-norm SVD) in Python — can disagree near `eps` / `1e3·eps`; Ho-Kalman too-few-SV id `sid:tooFewSV` vs Python `too_short`. `q = 0`, `p = 1`, `N = 2` work in Python; `N = 1` → `too_short`. | various | — |
| C17 | nit | `ltv_uncertainty_backward_pass.py:70` docstring cites "issue #137" (rendered on the API reference; principle 9). | — | — |

Fixed and verified from the July list in this area: C1, S1, S2 (with the weakness that `test_trust_region_mu_advances_and_terminates` asserts `iterations > max_iter`, which the mandatory final μ = 0 pass satisfies even on an immediate reject — it does not prove μ advanced), S3 (but see C1 above), S8 (Python-pinned only — MATLAB has no defective/integrator test; `test_sidLTIfreqIO.m` Test 14 checks the warning on a diagonalizable matrix, and `SPEC.md:1836`'s "integrator revert-check" claim is false for MATLAB), S9, minor #25 (DC convention unified and in the spec), `sid:orderExceedsRank`, `sid:notConverged`, `TrustRegion ∈ [0, 1]`, `MaxIter = 0` rejected, #189.

### 5.3 Utilities (§13–§15), example suite (`EXAMPLES.md`), packaging

Verified correct in Python: `detrend` equals per-segment `np.polyfit` exactly with `x == x_d + trend` to 0.0 and the partial last segment handled; `r_ee` equals `np.correlate/N` to 0.0; an AR(1) residual fails whiteness on 50/50 seeds; the state-space residual uses *measured* `x(k)`; `compare.fit` equals a hand NRMSE per channel, free-run from `x(0)`, `initial_state` honoured; the §15.5 multi-trajectory shapes, NaN-skip rule and mirror-pair case; default `M_test = min(25, N//5)`.

| # | Sev. | Defect | Where | Rule |
|---|---|---|---|---|
| U1 | **critical** | **State-space `residual` / `compare` size everything from `model.data_length` (July #27, unfixed, both ports).** Validation data shorter than the model → `IndexError`; longer → `residual` **silently** evaluates the first `N_model` steps (verified: 80-sample data, 50-row residual), `compare` → `ValueError`; `u = None` (allowed by §14.5) → `TypeError`. §15.7 explicitly shows `sidCompare(ltv, X_val, U_val)`. | `residual.py:282,300–302`, `compare.py:203,233–239`, `sidResidual.m:221,246–248`, `sidCompare.m:138,171,178` | §14.5, §15.7 |
| U2 | **major** | **`freq_domain_sim` zeroes DC and every FFT bin below the model grid** (July #30, unfixed, both): `y ≡ u`, `Ĝ ≡ 1`, `mean(u) = 3` → `whiteness_pass = False`, `r_ee(1..3) ≈ 0.999`, residual mean 2.99. With the default BT grid (`ω₁ = π/128`) and `N = 2000`, DC plus 7 bins are zeroed, so any non-detrended record fails §14.3 regardless of model quality and §15 `Predicted` is high-passed; §14.2/§15.2 say nothing about interpolation, DC or out-of-grid bins. The even-`N` Nyquist bin keeps a complex value that `real(ifft)` discards. | `freq_domain_sim.py:88,105`, `sidFreqDomainSim.m:56,66` | §14.2, §15.2 |
| U3 | major | **§14 multi-trajectory behaviour is spec-silent and the ports diverge for frequency-domain models.** State-space: per-trajectory `Residual (N×n_y×L)`, own-input cross-correlation, bound `2.58/√(L·N)`, `DataLength = L·N` — implemented identically in both (July #28 fixed *in code*) but §14.5/§14.6 still say `(N×n_y)`, `2.58/√N`, `DataLength = N`. Frequency-domain model with 3-D `y, u`: Python raises `AxisError`; MATLAB (static) silently predicts from trajectory 1 and broadcasts the subtraction over `L`. | `residual.py:357` → `freq_domain_sim.py:93`; `sidFreqDomainSim.m:66`, `sidResidual.m:288` | §14, principle 3 |
| U4 | major | `MaxLag` unvalidated (July #29, both): `max_lag = 0` → vacuous `whiteness_pass = True` (verified); `−1` → `ValueError: negative dimensions`; `> N` → raw `matmul` error. | `residual.py:132,155`, `sidResidual.m:93,117` | §14.5 |
| U5 | major | **`detrend` crashes on `(N, n_ch, 1)`** (same class as July C3, missed by #135): `ValueError: could not broadcast (1000,) into (1000,1)`; the docstring advertises 3-D input. MATLAB unaffected. | `detrend.py:151–154,180–183,208–211` | §1, §13.4 |
| U6 | minor | Variable-length list/cell input rejected by all three utilities with raw/misleading errors (Python raw `ValueError`; MATLAB `isreal(cell)` → `sid:complexData`). §1's "all sid functions accept cell arrays" is the inconsistent party. | `detrend.py:133`, `sidDetrend.m:62` | §1, §13.4 |
| U7 | minor | `ResidualResult` drift vs §14.6: `CrossCorr` is `(2M+1 × n_y·n_u)` not `(2M+1 × 1)`; `AutoCorrAll`, `WhitenessPassAll`, `IndependencePassAll` are unspecified extras (both ports). | `sidResidual.m:188–199`, `residual.py:229–240` | §14.6 |
| U8 | minor | `compare`/`residual` silently FFT-pad/truncate `u` to `len(y)` (fit 56 % returned for 1000 vs 900 samples); `order = 2.0` rejected in Python, accepted in MATLAB; both ports fit polynomials on raw `t = 0..N−1` (`cond(Vandermonde) = 1.8e25` at `d = 5`, `N = 1e5` — NumPy's column scaling hides it, MATLAB `polyfit` will warn); the family-wise whiteness pass rate for truly white residuals at `M = 25` is 0.785 (≈ 0.99²⁵), while examples and docstrings say "99 % confidence level". | — | §13–14 |

**Example suite vs `EXAMPLES.md` (conformance by value, both ports):** all 11 examples exist, are auto-discovered, seeded with the same seeds, import helpers by sibling path, and the Python notebooks execute (`pytest --nbmake`: 11 passed in 107 s). The MATLAB examples were not executed by the conformance pass (static reading; `runAllExamples.m` does fail the job on an erroring example); under apt Octave 8.4 only 1 of 12 runs, blocked by the figure-rendering defect recorded in the plan's §10.1. Deviations:

- **E1 (must-fix, Python)** — `example_freq_map.ipynb` Duffing cell uses `k_cubic = 5e4`; §1.5/§3.7 mandate `1e5` (MATLAB `exampleFreqMap.m:202` is `1e5`). A Plant-parameter MUST violated in the port §6.3 calls the reference.
- **E2 (spec defect, both)** — §3.6 row 7 mandates printing "MIMO response_std contains NaN: True", but SPEC §3.6 specifies the diagonal-approximation MIMO variance and both ports implement it: Python prints `… NaN: False` under a markdown cell claiming NaN; `exampleMIMO.m:171–174` says "the Python port returns NaN" (stale).
- **E3 (spec wording)** — "ε ~ N(0, 2·10⁻⁴)" reads as a variance; both ports use std 2e-4 (also 5e-4, 1e-4). §3.8's chirp "instantaneous frequency f₀ + (f₁−f₀)·t/(2·T_end)" is the phase-rate coefficient; both ports implement `f₀ + (f₁−f₀)t/T_end`, which is what the prose intends.
- **E4** — `python/examples/README.md` claims `example_output_cosmic` demonstrates `compare`/`residual`; it does not. `example_spectrogram` has no RNG seed in either port (§4.4 MUST, technically unmet; the chirp is deterministic). Neither `lti_freq_io`/`sidLTIfreqIO` nor `ltv_state_est`/`sidLTVStateEst` appears in any example.
- **Helpers (§2):** signatures, defaults, `expm` ZOH, RK4 tableau and shapes match in both ports; Python reproduces Plant A/B to 2e-11, the LTV test points to 4e-9, the impulse point to 2e-11, the LTI collapse bit-identically. **No unit test in either port checks any §1 / §2.2.3 / §2.3.4 tabulated test point**; `reference_test_msd.json` is an off-catalog `n = 3` plant covering `util_msd` only. The tabulated tolerances are unreachable from the tables as printed (Plant C to 8 significant digits yet 1e-6 demanded: measured 4.4e-4; Plant D `Ad` to 8 decimals with 1e-9 demanded: 8.9e-7; §2.2.3's 10-digit table cannot support 1e-9: 4.2e-9). §2.1.4 positivity constraints are unenforced in both ports (negative `m`, `k`, `c`, `Ts` accepted silently). `util_msd_ltv.m:153` misclassifies a `(1×N)` time-varying mass for `n = 1` as `N` masses (Python accepts it).

**Packaging and hygiene:**

- **P1 (July #31, unfixed)** — `pyproject.toml` `license = "MIT"` (PEP 639) with `requires = ["setuptools >= 68"]`: `pip wheel --no-build-isolation` fails with setuptools 68.x (`project.license must be valid exactly by one definition`); CI passes only via isolated builds pulling latest. Needs `setuptools >= 77`.
- **P2** — the sdist ships `tests/test_*.py` but not `tests/conftest.py`, `examples/util_msd.py` or `testdata/`, so shipped tests cannot run; the wheel has no `py.typed` (annotations unusable downstream); `sid-toolbox` vs `import sid` is stated nowhere user-facing; no `CITATION.cff`; the package is not on PyPI although the site says "install".
- **P3** — `CONTRIBUTING.md:313` and `python/RELEASE_NOTES.md:285` say Python 3.10–3.13; CI tests 3.10–3.14. No minimum-dependency leg.
- **P4** — `sidInstall.m` adds `matlab/sid` only; `sidInstall; exampleSISO` from elsewhere fails with "util_msd undefined" (README does not say to `cd`). `matlab/RELEASE_NOTES.md:11–13` / `python/RELEASE_NOTES.md:11–13` cite `#137`, `#138`, `#144` (user-facing; the CHANGELOG already carries the trail).
- **P5** — both header checkers: the spec rule is a *warning* that tests for a heading / substring only (§2.2 #33); a cross-check finds 52 distinct `SPEC.md §x.y` citations with 0 dangling, so a section-existence rule plus a ban on "not yet in SPEC" would pass today and close the hole. `private/sidFreqDomainSim.m` and `sidResultTypes.m` lack `SPECIFICATION:`; `freq_domain_sim.py:60` says "not a standalone section" though §15's `Verified by:` names it.

### 5.0 Found by the lead reviewer's own probes

- **`P` → `p_cov` breaks the documented field mapping** (minor). SPEC §9 and `docs/roadmap.md:29–41` state the MATLAB→Python field mapping is "purely syntactic" PascalCase→snake_case. `LTVResult.p_cov` / `LTVIOResult.p_cov` (`_results.py:186,298`) and `SpectrogramResult.complex_stft` (`:82`) are not that mapping (`P` → `p`, `Complex` → `complex`); `ResidualResult.auto_corr_all`, `whiteness_pass_all`, `independence_pass_all` (`:372–393`) and `LTVIOResult.output_dim` (`:281`) are fields the spec does not define (§14.6, §8.12.9). The docs site's uncertainty concept page uses the MATLAB names for the Python API as a consequence (§7).
- **`ltv_state_est` returns a trailing singleton for single-trajectory input** (minor). `ltv_state_est(Y2d, U2d, …)` returns `(N+1, n, 1)`, while every other function in the port returns the 2-D shape it was given (and MATLAB drops the trailing singleton). Every Python test that calls it carries a six-line squeeze idiom to cope (`test_ltv_state_est.py`, six occurrences).
- **No I/O contract for `sidLTVdiscFrozen`** (spec gap, minor). `SPEC.md` has no inputs/outputs table for the catalogued function: `'Frequencies'`, `'TimeSteps'`/`time_steps`, `'SampleTime'` and the output struct (`Frequency`, `FrequencyHz`, `TimeSteps`, `Response (nf×p×q×nk)`, `ResponseStd`, `SampleTime`, `Method`) exist only in code (`ltv_disc_frozen.py`, `sidLTVdiscFrozen.m`, `_results.py:325–354`). §8.11 gives the formula and §8.12.11 one usage line.
- **Multi-input COSMIC is correct in Python** (verified, not a bug): `q = 2, 3` on a noiseless 3-state plant recovers `A`, `B` to 1e-14, uncertainty shapes `(3, q, N)` / `(3+q, 3+q, N)` finite, frozen TF to 1e-14 vs analytic, `compare` fit 100 %, `ltv_state_est` to 1e-15, IO fast path identical to `ltv_disc`, partial-observation IO reaches 3e-3 / 2e-2 relative frozen-TF error. Recorded here because no test or vector covers it (§4.5).
- `sidLTIfreqIO.m:29` documents `'Horizon'` default `min(floor(nf/3), 50)` "where nf is the number of frequency bins"; the code uses `N_imp` (`:91–95`), the impulse-response length. Python documents `N_imp` correctly.
- `python/README.md:45` ships a placeholder line in the Quick Start (`y = np.convolve(u, [1.0], mode="full")[:N]  # placeholder; use scipy.signal.lfilter`) immediately overwritten by the real `lfilter` call.

---

## 6. Specification defects

Items in §2.2 (#37, #38, #39, #40) plus:

- **`spec/cosmic/output.md` contradicts itself on initialisation and the trust-region anchor.** §4.1 and §4.3 (and SPEC §8.12.3–8.12.4) specify the LTI initialisation via `sidLTIfreqIO` and interpolation toward `A₀`; the "Algorithm Summary" §7 (`output.md:205–216`) still prescribes the composite `A = I` solve of Appendix B as step 1 and `Ã = (1−μ)A + μ·I` in step 3, §8 costs "Initialisation: composite blocks", and §10 (`:242`) argues for "the `A = I` initialisation proposed here". The normative text and its own summary describe two different algorithms; Appendix B is dead weight the implementation never uses.
- **The spec is a changelog in places.** `SPEC.md` carries "Correction (was a scaling bug)" (`:1229`), "Correction (issue #138)" (`:1536`), "(Earlier drafts interpolated toward `μ·I` …)" (`:1520`), "previously implied" (`:642`), "Previously the reported `J` …" (`:1478`), "(This whole-slice behaviour supersedes the earlier 'affected row' wording)" (`:241`), issue numbers `#113`, `#120`, `#121`, `#137`, `#138`, `#144`, `#145d`, `#158`, `#189` inside normative sections. The spec is included verbatim on the public site (`docsite/spec/*` are `include-markdown` of `spec/*.md`), so principle 9 is violated on the page users read the contract from; and a contract that narrates its own history is harder to port (the next port has to decide which sentence is binding). `spec/cosmic/output.md:113` and `uncertainty_derivation.md:437` do the same.
- **Version/date header is stale.** `SPEC.md:3–4` says "Version 1.0.0, Date 2026-04-04"; the document has since changed normative behaviour in §2.6, §2.7, §3.3, §5.2, §6.5, §6.6, §8.9.2, §8.12.2, §8.12.4, §8.12.11, §8.12.12, §8.13.1, §10.3, §15.5 (the 0.2.0 release) with no version bump and no spec changelog. `EXAMPLES.md` has versioning rules (§6) that `SPEC.md` lacks.
- `SPEC.md:1860` — `Method` is `'sidFreqBT'`, …, "`'sidFreqMap'`, or `'welch'`": `'welch'` is an `Algorithm` value, not a `Method`; and the Python values differ anyway (§2.2 #7).
- `SPEC.md:389` uses `N_eff` without defining it (§3.6).
- `SPEC.md:59` — "All `sid` functions accept multiple independent trajectories": `sidDetrend`, `sidResidual` (time-series), `sidModelOrder` do not take trajectories; the sentence over-promises.
- §8.5 defines `Cost` as `(1 × 3)` and §8.12.9 as `(n_iter × 1)` for a result that "extends the standard `sidLTVdisc` output struct" — the same field name with two incompatible semantics; `sidResidual`/`sidCompare` must special-case it.
- §8.4.3 lists only `'LambdaGrid'` and `'Algorithm'`; both ports also accept `method`/`'Method'` (`validation` | `frequency`), `segment_length`, `consistency_threshold` (0.90), `coherence_threshold` (0.3), `precondition` — the §8.11.2 thresholds are hard-coded in the spec text but exposed as spec-silent options.
- §14.2 state-space residual: `x̂(k+1) = A(k) x̂(k) + …`, `e(k) = x(k+1) − x̂(k+1)` is ambiguous between a one-step-ahead residual (`x̂(k) := x(k)`) and a simulation error (propagated `x̂`); the two differ materially for an unstable or lightly damped plant.
- §12 and §9's "`Verified by:` manual" on bibliographic references and metadata fields are noise that dilutes the tag.
- §8.12.2 says the reported IO cost equals "`N` times §8.3.3's `f(C)`"; it is `2N·f(C)` because §8.3.3 carries a ½ (verified exactly on the fast path). §8.4.2 writes `F_j`, `R_j` without the ½ that `automatic_tuning.md` §2.1 and both ports use (July #39, first half — corner invariant, contract self-inconsistent). `uncertainty_derivation.md` §1 uses `X'(k) ∈ ℝ^{p×L}` while SPEC §8.3.2 uses `ℝ^{L×p}`; its §7 Step 2 `ν` formula lacks the `N·` factor §8.9.3 and the code use; §8.9.2 writes `H⁻¹` for the Hessian inverse while `H` is the observation matrix elsewhere in §8.
- §8.12.12's rationale for the machine-ε floor is false and step 5b makes exact rank-`n` undetectable (§5.2 C1) — ADR-0004's premise needs revisiting.
- §8.12.9 `Cost (n_iter × 1)` vs the fast path's `Iterations = 0`, `Cost` length 1 (§5.2 C11).
- Frequency-domain spec items found by the §5.1 pass: §2.6/§3.3 list `sidFreqETFE` among the estimators sharing the `Inf` sentinel while §4.5 sets ETFE std to NaN (SP1); §2.7 says NaN "is produced" at degenerate frequencies but every BT-family estimator substitutes `Φ̂_y` with coherence 0 (SP2); §5.4 claims a BT short-data floor that does not exist (SP3); §2.3/§6.6 claim non-negativity the Hann lag window does not guarantee (SP5); §3.4/§3.5 omit `L` (SP6); §4.1 is undefined for single-trajectory MIMO ETFE, §4.2 silent on NaN bins, §4 silent on ETFE `WindowSize`/`Coherence` (SP7); §6.5/§6.6 have no Welch `σ_G` formula and no statement that Welch MIMO `σ_G` is NaN, §6.10's compatibility claim omits the DC bin (SP8); §7.2/§2.3 mandate cell input to `sidSpectrogram` that neither port supports (SP9); §10.3 is silent on `Φ̂_v`/`γ̂²` when the excitation check fires (SP12); `SPEC.md:409` cites "collinear-MIMO cases in … `test_sidFreqBT.m`" that do not exist (SP11).
- `Verified by:` over-claims in §8: `SPEC.md:1824` "`test_sidLTVdiscVarLen.m` (dense-LSQ oracle)" — the file has no oracle (Tests 1–7 are shape/recovery); "`unit(Py)` `test_ltv_disc.py` var-len cases" — none exist; `:1836` "integrator revert-check" — no MATLAB defective-stabilization test; `:1826` "not reachable from normal public-API inputs" — reachable (§5.2 C3).

---

## 7. Documentation site (`docsite/`, `mkdocs.yml`, `scripts/`)

**Build:** `mkdocs build --strict` with notebook execution exits 0 in 127 s with zero warnings; all 11 notebooks execute; 89 pages, 22 MB. `gh-pages` was last deployed 2026-07-30 from the #110 merge. Internal links: 0 broken. The site is in good structural shape; the defects are content.

**Factual errors on user-facing pages (fix first):**

1. `docsite/getting-started/install-matlab.md:3–4, 36–39` — "Tested on **MATLAB R2016b+** and **GNU Octave 8+**" and a compatibility table with those floors. Every other source (`README.md:33`, `matlab/README.md:5–6`, `matlab/RELEASE_NOTES.md:35–36`, `tests.yml`) says **R2024a+ / Octave 11+**, and CI explicitly rejects apt's Octave 8 as below the floor. The highest-impact error on the site.
2. `docsite/concepts/uncertainty.md:19–23` — "`Uncertainty=True` (Python)" and fields `AStd`, `BStd`, `P`: the Python kwarg is `uncertainty=True` and the fields are `a_std`, `b_std`, `p_cov` (§5.0).
3. `docsite/api/python/index.md:21` — "`lti_freq_io` — LTI frequency response from input/output data" is backwards; it returns an `(A0, B0)` Ho-Kalman realization from I/O data.
4. `docsite/api/index.md:53` — "Every function returns a frozen dataclass (Python) or struct (MATLAB)": false for `model_order` (tuple), `detrend` (tuple), `lti_freq_io` (tuple), `ltv_state_est` (ndarray/list), `ltv_disc_tune` (tuple), the plot functions (dict/handles).
5. `docsite/examples/index.md:3` — "Every public estimator ships with a runnable example": `lti_freq_io`/`sidLTIfreqIO` and `ltv_state_est`/`sidLTVStateEst` appear in no notebook or MATLAB example.
6. `docsite/index.md:46` — "confidence bands for all estimation functions": ETFE returns NaN std by contract (§4.5).
7. `docsite/getting-started/quick-start.md:3–5` — "Both produce numerically identical results": the two snippets use different RNGs.
8. Nowhere on the site: the package is **not on PyPI** and the distribution name `sid-toolbox` differs from the import name `sid`; the landing card gives `pip install -e "./python[plot]"` with no clone step.

**Dev-tracking leakage (principle 9):** `about/changelog.md` includes `CHANGELOG.md` verbatim (20 `#NNN` references and `ADR-0002` render; the page intro precedes the included H1); `about/contributing.md` includes `CONTRIBUTING.md` verbatim ("Pre-push self-review (agent convention)", "reviewer subagent", two `CLAUDE.md` links, `#145`, `#172`, `#113`, `#160`, three ADRs, `.github/rulesets/main.json`) — internal engineering process on the user site, contradicting the `mkdocs.yml:146–150` comment that publishing from `docs/`-type content is a deliberate decision; every generated Python API page ends with a "Changelog — 2026-04-08: First version (Python port)" block; the spec includes carry the issue numbers listed in §6.

**Rendering defects:** `.. math::` and `.. [1]` RST directives in `ltv_disc.py` / `lti_freq_io.py` docstrings render literally on the API pages (mkdocstrings is in numpy mode without RST processing); "Specification: SPEC.md S2" (ASCII `S` for `§`, unlinked); signatures unformatted (no ruff/black in `requirements-docs.txt`); `api/python/results.md` hand headings duplicate mkdocstrings' own `show_root_heading` so every dataclass appears twice in the TOC; the generated `api/matlab/sidResultTypes` page treats a documentation-only file as a function and collapses its tables into run-on prose, and it is listed in the function index. **Notebook math is almost certainly unrendered:** `docsite/javascripts/mathjax.js:9–10` restricts MathJax 3 to `arithmatex` spans, nbconvert emits MathJax-2 config that v3 ignores, and the built `example_output_cosmic` HTML contains raw `$$x(k+1) = …$$` with zero `arithmatex` wrappers (reasoned from the built HTML, not browser-verified; 5 of 11 notebooks use `$` math). Twelve right-hand-TOC anchors on notebook pages are broken because mkdocs-jupyter slugs headings containing backticks differently from the ids it emits.

**Generator / gate findings:** `docsite/hooks/rewrite_external_links.py:67–78` rewrites *every* resolvable relative link to GitHub, including links to files that are themselves published on the site (its docstring at `:10–12` says the opposite) — `spec/cosmic/online-recursion` links `uncertainty_derivation.md` off-site instead of `../uncertainty-derivation/`, and `SPEC.md` links from the examples-spec and changelog pages go to GitHub. `scripts/build_matlab_api.py`: a malformed or empty header yields a title-only page silently; `DENYLIST = {"sidInstall"}` is dead (the file lives in `matlab/`, not `matlab/sid/`); `PYTHON_EQUIV` (`:310–330`) and `NOTEBOOK_TITLES` (`link_notebook_examples.py:19–31`) are hand-kept manifests — a new MATLAB function silently loses its Python cross-link, and nothing checks MATLAB↔Python 1:1 pairing or that `api/index.md` lists new MATLAB functions (the Python side *is* gated by `check_python_api_pages.py`, which runs inside the build and would catch a missing stub). Four `SUMMARY` pages are built, indexed and searchable. `requirements-docs.txt` is unpinned. Notebooks are committed without outputs, so figures exist only after build-time execution (every PR re-executes all 11; GitHub's notebook viewer shows no plots). The spec pages are includes (no drift risk); the edit button on them points at the include stub, not `spec/SPEC.md`.

**Missing content a toolbox user would expect:** an estimator-selection guide (the concept pages are 12–40-line link lists: no "BT vs BTFDR vs ETFE vs `freq_map`; when COSMIC vs Output-COSMIC"); a FAQ/troubleshooting page for the contract diagnostics (`NaN`/`Inf` sentinels, `sid:windowReduced`, `sid:notConverged`, `sid:stabilized`, `orderExceedsRank` — undocumented on-site outside docstrings); citation (no `CITATION.cff`; the two arXiv references appear nowhere on the site); a license page; a version number anywhere on the site; a changelog without dev references; versioned docs (#201 open); per-function MATLAB result-struct docs (the function index promises "struct layouts documented inline on each function page" and never links the struct reference); a note that `util_msd.py` is a sibling helper the notebooks need. Present and fine: search, dark mode, responsive theme, 404 page, sitemap, Binder badges (URLs correct, `main` branch, `postBuild` executable).

---

## 8. Potentially publishable parts

A literature check (arXiv/Google Scholar, 2026-10) finds the base COSMIC algorithm (Carvalho, Soares, Lourenço, Ventura, arXiv:2112.04355) and the control-application paper (Łaszkiewicz et al., arXiv:2509.13531) but **no publication** of the extensions this repository already specifies and implements. Ranked by readiness:

1. **The toolbox itself — JOSS / SoftwareX paper.** sid meets the JOSS bar (research software, open license, tests, documentation site, examples, ≥ 6 months of development, a statement of need: a dependency-free, cross-validated MATLAB/Octave + Python system-identification toolbox covering non-parametric *and* LTV state-space identification with uncertainty). Missing before submission: `CITATION.cff`, a `paper.md` + `paper.bib`, a PyPI release (the docs say "install" but the package is not published), a corrected compatibility floor on the site, and the test-gate repairs in §4 (a JOSS reviewer runs the tests and reads CI). Effort: small once Phase A of the plan lands.
2. **Closed-form Bayesian uncertainty for regularized LTV identification** (`spec/cosmic/uncertainty_derivation.md`, SPEC §8.9). Contribution: matrix-normal posterior `Σ ⊗ P(k)` for COSMIC, O(N(p+q)³) extraction of the diagonal blocks of the inverse block-tridiagonal Hessian via left/right Schur complements (= RTS smoother covariance in parameter space), exact hat-trace degrees of freedom for the noise-covariance estimate, the scaling subtlety (`N·λ` effective prior) that the July review found and fixed, and propagation to frozen transfer-function bands via the rank-1 Jacobian (§8.11.1). Format: a 6-page conference paper (CDC / ECC / IFAC SYSID) or a short Automatica/IEEE CSL note. Missing: the MC calibration campaign (exists, never run — §3.3), a comparison against the sandwich/frequentist covariance and against bootstrap, a real-data demonstration (the Comet Interceptor example from the 2022 paper is the obvious vehicle), and a decision on how to present the Bayesian-vs-frequentist caveats of `uncertainty_derivation.md` §4/§10.
3. **Output-COSMIC** (`spec/cosmic/output.md`, SPEC §8.12–8.13). Contribution: LTV identification from partial observations by alternating an RTS state smoother with the COSMIC closed form (a MAP-EM with a smoothness prior), Ho-Kalman initialisation from a Blackman-Tukey estimate transformed to the `H`-basis, the two-level trust-region homotopy toward the LTI initialisation, and the Grippo–Sciandrone stationarity argument. Format: full paper (Automatica / IEEE TAC technical note / IJACSP). Missing: identifiability treatment beyond the similarity-ambiguity remark (§8.12.7), systematic experiments against EM-LTV and windowed subspace methods, characterisation of when the trust region helps (today "threshold-dependent", §3.4, and the stored vector is non-converged, §4.1 F7), and the convergence rate / local-minimum study the repository's own tests hint at (`test_convergence` fixture never converges). The spec's own §7 summary must first be reconciled with §4 (§6).
4. **λ selection by parametric/non-parametric consistency** (SPEC §8.11, `ltv_disc_tune(method='frequency')`; `automatic_tuning.md` Method 3 spectral pre-scan). Contribution: tune the regularisation of a parametric LTV model by χ²-testing its frozen transfer-function bands against non-parametric `sidFreqMap` bands; the spectral pre-scan for a *time-varying* λ_k schedule is explicitly flagged "conjecture, unvalidated" (`automatic_tuning.md` §4.10, §6). Format: SYSID/CDC paper once validated. Missing: everything empirical — there is no correctness oracle in the tests (§4.4), no benchmark against L-curve and validation tuning, and no implementation of Method 3 at all.
5. **Online/recursive COSMIC as a Kalman smoother in parameter space** (`spec/cosmic/online_recursion.md`). Theory is written (information-form filter, windowed smoother, λ selection from innovation consistency); deferred to v2 with nothing implemented. Publishable with an implementation and a streaming use case.
6. **The development methodology** — spec-as-contract polyglot numerical software, cross-language reference vectors as contract artifacts with drift hardening, ADRs, and parallel agent-driven porting under a review/remediation cycle. This repository is a measured case study: two full reviews, 12 remediation issues, five regeneration-drift incidents (§4.1 F1), a documented bug class the vectors cannot catch (joint drift), and governance docs (`CLAUDE.md`, `REVIEW_CONTEXT.md`, `DISCIPLINE_ADOPTION.md`). Format: experience report (ICSE-SEIP / FSE industry / IEEE Software / JOSS "software engineering" track / Journal of Open Research Software). Missing: an honest write-up of what did *not* work (§3–§4 here), and permission to discuss the agent workflow.

Not publishable as novel: BTFDR (equivalent to MathWorks `spafdr`), the Welch/BT scaling reconciliation, the model-order gap convention (ADR-0004; a paragraph in paper 1 at most).

---

## 9. Plan

The sequenced remediation plan is in [`docs/plans/2026-10-04-review-v2-remediation-plan.md`](../plans/2026-10-04-review-v2-remediation-plan.md). In one paragraph: **(A)** make the gates able to bite before trusting them again — required checks, `[skip ci]` regeneration validated, NaN-capable Octave comparison, Python structural field loop, atol floors enforced by a generator-side gate, the three 1 % vectors re-derived (converged IO case; var-len at 1e-6), spec-derived oracle job over the payloads, missing-vector = failure; **(B)** tighten or delete the relaxed/vacuous tests enumerated in §3 and add the free oracles (MIMO static-gain plant, L-copies ensemble identity, OLS at N = 2, FFT = direct at 1e-12, `q ≥ 2` everywhere) in both ports; **(C)** close the coverage matrix of §4.5 symmetrically, starting with warnings/error identifiers and the never-called options; **(D)** fix the conformance defects of §5 spec-first; **(E)** spec hygiene (§6) — de-narrate, version, reconcile `output.md`, add the frozen-TF I/O table; **(F)** docs site content fixes (§7); **(G)** the publication track (§8) once A–C are green.
