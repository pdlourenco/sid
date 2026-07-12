# sid v0.1 — Full Codebase Review

**Date:** 2026-07-12
**Scope:** `spec/` (SPEC.md, EXAMPLES.md, cosmic/*.md), `python/` (all sources + tests), `matlab/` (all sources + tests), `testdata/`, `docs/`, `.github/`.
**Method:** every source file was read in full and checked against the mathematical specification. On top of static review, the highest-risk math was verified **independently and numerically**:

- COSMIC solution, cost triple, posterior blocks `P(k)`, `AStd`/`BStd` extraction, and DoF/noise-covariance formulas cross-checked against a dense brute-force solve of the full block-tridiagonal system (per-step λ, N=6, p=2, q=1, L=4) — all match to machine precision.
- `ltv_state_est` cross-checked against a dense LSQ solve of `J_state` with general `R`, `Q` — matches to 1e-14 for uniform horizons.
- BT FFT fast path vs direct DFT vs a hand-written correlogram at single frequencies — all match to ~1e-15, including the `K>1` buffer regime (M=150).
- `spectrogram` vs `scipy.signal.spectrogram(mode='psd', scaling='density')` — bit-exact (max rel err 7e-16).
- Welch `freq_map` transfer function and coherence vs `scipy.signal` — machine precision; the sid Welch estimate matches the **true** plant phase (a naive scipy `csd` comparison has the conjugation backwards).
- `C_W = 3M/4` closed form verified numerically for M ∈ {2, 5, 30, 101}; the algebra in SPEC §3.1 (Σcos = −1, Σcos² = M/2) checks out.
- Monte-Carlo calibration + exact sandwich-covariance computation for the COSMIC posterior (finding S1 below).
- Both test suites pass at HEAD: Python 329/329; MATLAB/Octave 39/39 files, 402/402 cases (executed under GNU Octave 8.4.0). The two critical findings shared with MATLAB (C1, C2) were reproduced by executing the MATLAB implementation under Octave as well, not only the Python port.
- `testdata/validate_reference.m` executed under the same Octave: 21 of 22 reference files pass; `reference_ltv_cosmic.json` fails on a near-zero element — a validator tolerance-design flaw, see S10(e′).

Severity scale: **critical** = wrong results or crash on documented usage; **major** = spec violation / silent numerical error in realistic use; **minor** = contract or edge-case defect; nit = cosmetic.

Cross-language notes: the Python port mirrors the MATLAB reference closely, so most defects exist **in both implementations**; file references below give Python first, MATLAB second where applicable.

---

## 1. Executive summary

The mathematical core of this library is in very good shape. The genuinely hard parts — biased covariance estimation, the K-strided FFT lag-buffer with wrapped negative lags, cross-covariance conjugation at negative lags, the COSMIC forward/backward recursions and λ-indexing, the left/right Schur uncertainty passes, the RTS-smoother normal equations (uniform horizons), Ho-Kalman realization, and the frozen-TF rank-1 uncertainty propagation — are all **correct**, independently verified, and in several cases protected by exemplary regression tests (the FFT large-M suite, the finite-difference Jacobian check, the MATLAB Monte-Carlo uncertainty test).

The defects cluster in three places:

1. **The edges of the input-shape contract** — four reproducible crashes on documented inputs, all silent in the test suites (multi-output time-series maps; 3-D single-trajectory arrays; MIMO plotting; degenerate Welch segments).
2. **Silently wrong numbers in less-traveled numerical paths** — the variable-length RTS smoother (wrong terminal block), `model_order`'s lag-0 off-by-one, `compare`'s average-then-fit, the COSMIC posterior scaling (a spec-level derivation error), the non-functional Output-COSMIC trust region, and the eigen-reconstruction blowup in `_stabilize`.
3. **Spec-internal inconsistencies faithfully inherited by both ports** — degenerate-coherence σ conventions, Welch one-sided scaling (Nyquist / factor 2), the §10.3 constant-input contract that the §10.2 relative threshold can never satisfy, and the λ-vs-N·λ objective mismatch in Output-COSMIC.

The test suites are strong on happy-path math and historical regressions, and systematically weak on: MIMO numerics, multi-trajectory semantics, variable-length inputs, warnings, degenerate excitation, and the entire COSMIC uncertainty path in Python (zero direct tests; the one stored cross-language reference for `P` is generated incorrectly).

---

## 2. Critical findings

### C1. Variable-length `ltv_state_est` uses the wrong terminal diagonal block — silently wrong states  *(both languages)*
`python/sid/ltv_state_est.py:213` (blocks built at 164–174); `matlab/sid/sidLTVStateEst.m:94–104, 148–154`.
The shared diagonal blocks are built once for the maximal horizon N; a trajectory with `Nl < N` gets `S_blk[Nl]` as its terminal block, which is the *interior* block `HᵀR⁻¹H + Q⁻¹ + A(Nl)ᵀQ⁻¹A(Nl)`. Per `spec/cosmic/output.md` App. A.1 the terminal block must be `HᵀR⁻¹H + Q⁻¹` (no dynamics constraint leaves the last state). The RHS *is* built with the correct terminal form, so the solved system corresponds to a phantom constraint `x(Nl+1)=0, u(Nl)=0`, biasing `x̂(Nl)` toward `null(A(Nl))`.
**Verified:** N=8, Nl=5 — short-trajectory states off by **2.4 max-abs** vs the dense LSQ minimizer; the full-length trajectory matches at 1.3e-14. Poisons every variable-length `ltv_disc_io` run (state step no longer minimizes J; the §8.12.5 monotone-decrease guarantee is void). No test uses list inputs to `ltv_state_est`.

### C2. `freq_map` crashes on multi-output time-series data  *(both languages, both algorithms)*
`python/sid/freq_map.py:588`; `matlab/sid/sidFreqMap.m:379–381`.
The storage branch `if ny == 1 or is_time_series:` ravels an `(nf, ny, ny)` per-segment noise spectrum into a slice of incompatible shape. **Verified:** `freq_map(np.random.randn(1000, 2), None, segment_length=256)` → `ValueError: could not broadcast (512,) into (128,2,2)`. SPEC §6.1 explicitly supports this mode. Fix: condition on `ny == 1` only.

### C3. 3-D single-trajectory input `(N, ch, 1)` crashes `freq_bt` / `freq_map` / `spectrogram`  *(Python)*
Root cause `python/sid/_internal/cov.py:77–81` (the `L == 1` branch does 2-D algebra on a 3-D array); separate instance at `python/sid/spectrogram.py:334–338` (branch on `n_traj > 1` instead of `x.ndim == 3`).
**Verified:** `freq_bt(randn(500,1,1), randn(500,1,1))` and `spectrogram(randn(1000,1,1), window_length=128)` both raise `ValueError`. Docstrings document `(N, p, L)` with no exclusion of L=1; programmatic slicing `arr[..., :1]` hits this immediately.

### C4. Complex `u` is silently cast to real instead of raising  *(Python)*
`python/sid/_internal/validate_data.py:234` (dead checks at 250–254, 267–271).
`np.asarray(u, dtype=np.float64)` runs before the complex check, discarding the imaginary part with only a NumPy `ComplexWarning`; the subsequent `iscomplexobj` check is unreachable. SPEC §10.1 mandates an error; the docstrings promise `SidError('complex_data')`. **Verified.** (The `y` path is protected; the list-input path is correct.) Silent data corruption for any user passing analytic signals.

### C5. MIMO plotting is broken end-to-end  *(Python; MATLAB partially)*
- `python/sid/map_plot.py:138–179`: 4-D MIMO `(nf, K, ny, nu)` response reduced only one axis → `pcolormesh` rejects it. **Verified crash** on a 2×2 map, all plot types.
- `python/sid/spectrum_plot.py:104–122`: 3-D `(nf, ny, ny)` noise spectrum reduced to 2-D → tuple-unpack `ValueError`. **Verified.**
- `matlab/sid/sidMapPlot.m:64–67`: additionally rejects `sidSpectrogram` results outright, contradicting SPEC §7.5, and its `'spectrum'` branch reads a field spectrogram results don't have.

---

## 3. Major findings

### S1. Spec-level: COSMIC posterior covariance is inconsistent with the estimator's own scaling — reported stds too wide at practical λ  *(spec + both languages)*
SPEC §8.9.2 Step 1 (`S(k) = N(S_scaled(k) − reg(k)) + reg(k)`), implemented verbatim in `python/sid/_internal/ltv_uncertainty_backward_pass.py:94–103` and `matlab/sid/private/sidLTVuncertaintyBackwardPass.m`.
The `1/√N` data scaling (§8.3.2) makes the returned MAP mean the minimizer of `‖unscaled residuals‖² + N·λ‖ΔC‖²` — the effective prior weight is **N·λ** (that is the stated purpose of the convention). The reconstructed Hessian instead uses weight **λ**, so the reported `P(k)` is the posterior of a *different* estimator than the returned mean.
**Verified exactly:** with the estimator's true (linear-in-noise) sampling covariance computed in closed form, the self-consistent precision `VᵀV + Nλ FᵀF` tracks the sandwich covariance almost perfectly across λ, while the spec formula overstates variance by **1.86×** at λ=1e2 (example problem), converging only in the λ→0 and λ→∞ limits. Monte-Carlo calibration (200 trials) confirms reported `AStd` exceeds empirical spread. The mid-λ regime — where the L-curve typically lands — is the worst case. Downstream, the DoF hat-trace also overcounts effective parameters (inflating `Σ̂`, conservative).
**Fix belongs in the spec first** (Step 1 should scale the reg terms by N too, equivalently `P = P_scaled/N`, and the off-diagonal Schur couplings become `(Nλ)²`), then both ports and `uncertainty_derivation.md` §3.4/§5.1.
Related manifestation (MATLAB `sidLTVdiscIO.m:441–487`, Python equivalent): Output-COSMIC's reported/convergence-tested cost `J` uses weight λ while the embedded COSMIC step minimizes with effective weight N·λ — so the documented monotone decrease of `Cost` is not actually guaranteed, and the user's λ is N× stronger than §8.12.2 implies. The spec (§8.12.2 vs §8.12.3) is self-contradictory here.

### S2. Output-COSMIC trust region is non-functional as specified and empirically harmful  *(Python verified; same structure in MATLAB)*
`python/sid/ltv_disc_io.py:300–376`; `matlab/sid/sidLTVdiscIO.m:196–213`.
(a) Spec §8.12.4 prescribes an outer μ-loop, each stage running the inner alternating loop **to convergence**; the code fuses them into a single `max_iter` budget and only halves μ when the inner test happens to fire — at μ=1 the iteration is not coordinate descent on any single objective, so it may never fire. **Verified:** standard partial-observation example, `trust_region=1.0, max_iter=100`: μ stays 1 for all 100 iterations, cost bottoms at 19.0 then **rises**; `trust_region='off'` reaches J=0.235 — two orders of magnitude better.
(b) The interpolation target is the LTI initialization `A₀`, not the identity the spec specifies (`Ã = (1−μ)A + μI`).
(c) The reject path sets μ=0 and keeps iterating instead of "revert and terminate" (both languages).
Recommendation: rewrite to the spec's two-level loop with a per-μ inner budget, or remove the parameter; no test exercises `trust_region` at all.

### S3. `model_order` builds the Hankel from lag 0 — order overestimated by one for generic strictly proper plants  *(both languages)*
`python/sid/model_order.py:189–223`; `matlab/sid/sidModelOrder.m:157–159`.
SPEC §8.12.12 defines the Hankel from `g(1)…g(2r−1)` (`g(k) = HA^{k−1}B`); the code keeps the IFFT's lag-0 sample, so the Hankel has the McMillan degree of `z⁻¹G(z)` = n+1 unless G has a zero at the origin. `lti_freq_io` correctly discards lag 0 — the two consumers disagree.
**Verified:** exact samples of `G(z) = (z+0.5)/(z²−1.2728z+0.81)` (true n=2) → estimated **n=3**, spurious structural singular value 0.41. Every existing test plant has numerator `z^k` or is biproper — precisely the family that masks the bug. The documented workflow (§8.12.12: model_order → H → `ltv_disc_io`) hands an inflated state dimension to the identifier for typical physical plants.

### S4. `compare` multi-trajectory fit averages signals, not fits  *(both languages)*
`python/sid/compare.py:200–228`; `matlab/sid/sidCompare.m:128–160`.
SPEC §15.3: "fit is computed per trajectory and averaged." Code averages `y` and `ŷ` across trajectories, then computes one NRMSE on the ensemble means — systematically optimistic (independent errors cancel), and degenerate: **verified** that two mirror-image trajectories (`[X, −X]`, `[U, −U]`) give `measured ≡ 0` and `fit = NaN` where the spec answer is ≈100% per trajectory.

### S5. BTFDR MIMO path skips the singularity guard and PSD clamp entirely  *(both languages)*
`python/sid/freq_btfdr.py:390–399`; `matlab/sid/sidFreqBTFDR.m:262–269`.
§5.1 says BTFDR is "identical to sidFreqBT except" the window mapping; §2.6/§2.7 require the rcond→NaN+warning guard and the negative-eigenvalue zeroing of `Φ̂_v`. **Verified:** exactly collinear inputs → raw `LinAlgError: Singular matrix` (Python); near-singular inputs return garbage silently; the returned MIMO noise spectrum had **eigenvalue −0.14 / diagonal −0.18** on ordinary 2×2 data. `freq_bt` implements both protections; BTFDR has neither. Also (minor): BTFDR silently clamps `M_k` to `⌊N/2⌋` without the §10.1-mandated warning, silently raises `M_k<2` to 2, and its short-data default deviates from §5.4/`freq_bt` behavior for N∈[10,19].

### S6. Welch one-sided scaling: Nyquist bin doubled; BT↔Welch differ by an exact factor 2  *(both languages; partially a spec bug)*
`python/sid/freq_map.py:84–85, 135–138`; MATLAB equivalent.
(a) The one-sided doubling is applied to **all** bins including Nyquist. **Verified:** sid Welch noise spectrum = exactly 2× the scipy/`cpsd` reference at ω=π, machine-precision equal elsewhere; the package's own `spectrogram` correctly excludes Nyquist (bit-exact vs scipy). SPEC §6.5 says "excluding DC" while §7.3 says "excluding DC and Nyquist" — the code faithfully implements each section, making the two estimators mutually inconsistent.
(b) The Welch factor 2 makes `noise_spectrum` from `algorithm='welch'` exactly **2×** that from `'bt'` on the same data (BT's correlogram is a two-sided convention) — a +3 dB step when switching algorithms, contradicting §6.6's "same convention" claim. Spec must pick a convention; code follows.
(c) `noise_spectrum_std` uses `Var = 2Φ²/ν` (correct χ²) where SPEC §6.5 says `Φ²/ν` — a √2 spec error, code right.
(d) `ν ≈ 1.8J` is applied for *any* overlap and window, though Harris's value is specific to 50 % Hann.

### S7. Degenerate-coherence/excitation σ handling: three conventions, none matching the spec, one cemented by tests  *(both languages)*
SPEC §3.3/§10.2: `γ² < ε` or `Φ̂_u ≈ 0` ⇒ `σ_G = Inf` + NaN response + warning. Reality: BT floors coherence at 1e-10 → huge **finite** σ (`python/sid/_internal/uncertainty.py:98–101`; `matlab/.../sidUncertainty.m:87–90`), and NaN (not Inf) where G is NaN; BTFDR returns Inf on one branch and the floor on another; the Welch inner path writes NaN and never warns; ETFE's SISO branch NaNs without warning. `test_uncertainty.py::test_zero_coherence_clamp` and MATLAB `test_sidUncertainty` Test 8 **assert the non-spec behavior**. Additionally the §10.2 threshold is *relative* (`< 1e-10·max|Φ_u|`), which means the §10.3 "constant input ⇒ all-NaN + warning" contract can never fire: **verified** `freq_bt(randn(500), ones(500))` returns |G| ranging 0.04–109 with zero NaNs and zero warnings (the biased covariance of a constant is `(N−τ)/N ≠ 0`). §10.2 and §10.3 are mutually inconsistent as written; an absolute floor or input-variance check is needed.

### S8. `lti_freq_io` eigen-reconstruction explodes for defective A  *(both languages)*
`python/sid/lti_freq_io.py:576–592`; `matlab/sid/sidLTIfreqIO.m:362–387`.
Stabilization reflects eigenvalues then reconstructs `A = V·diag(λ)·V⁻¹` with no `cond(V)` guard. **Verified:** a 3×3 Jordan block with λ=1.5 produces entries of magnitude **3e14**. Realistic trigger: integrator chains (repeated eigenvalue 1, defective). Needs a Schur-form reflection or a conditioning fallback.

### S9. `ltv_disc_frozen` of an Output-COSMIC result silently returns the state response, not `H(zI−A)⁻¹B`  *(both languages; spec gap)*
`python/sid/ltv_disc_frozen.py:182–184`. Implements §8.11.1 literally (no H), **verified** — but SPEC §8.12.11 explicitly pipes an `sidLTVdiscIO` result into it for Bode comparison against output-based BT bands, which is dimensionally wrong for H ≠ I (n×q state response vs p_y×q). The uncertainty propagation itself (`Var(G_ab) = (vᴴP(k)v)(rₐΣrₐᴴ)`) is exactly right, conjugation included. Spec must define the IO-frozen contract; code should honor `result.h` or reject IO results.

### S10. Testdata / cross-validation gaps that hide all of the above  *(infrastructure)*
- `testdata/generate_reference.m:489–490` divides `S` by N **again** before calling the uncertainty backward pass (production passes `S` directly), so the stored reference `P` in `reference_cosmic_internals.json` does not match what either implementation produces — and the Python cross-check reads only the cost triple from that file, so nothing notices.
- Four reference files (`reference_detrend.json`, `reference_freq_domain_sim.json`, `reference_ltv_io.json`, `reference_uncertainty.json`) are consumed by **no** Python or MATLAB test — the Octave `validate_reference.m` glob is their only consumer. The Python ports of those modules are never cross-checked against MATLAB.
- Python hardcodes tolerances (sometimes 100× looser, sometimes 10⁴× tighter) instead of reading the JSON `tolerance` blocks; Octave defaults absent tolerances to 1e-6 while Python uses 1e-10/1e-12 for the same fields — the "cross-language equivalence contract" is enforced differently per language.
- (e′) **Observed in practice:** the validator uses `atol = 0` for every field, so near-zero elements are held to a pure-relative check that any BLAS/rounding difference between the generating engine and the validating engine can trip. Executing `validate_reference.m` under GNU Octave 8.4.0 against the committed (MATLAB R2025a-generated) vectors: 21/22 pass, but `reference_ltv_cosmic.json` **fails** on `A` element 22 — |diff| = 2.17e-10 absolute on an expected value of 1.77e-4 (machine noise for this computation), vs threshold `1e-6 × 1.77e-4 = 1.77e-10`. The gate is therefore environment-flaky (it passes on CI's Octave build), while a real absolute error up to `1e-6·|element|` on O(1) elements passes. Fix: give every comparison an absolute floor (e.g. `atol = rtol·‖expected‖∞` per field, or explicit `*_atol` keys as `reference_test_msd.json` already uses for `Bd`), honored identically by both consumers.
- No reference vectors exist for: multi-trajectory anything, ETFE time-series periodogram, Welch `freq_map`, var-len COSMIC/StateEst, auto-λ, `ltv_disc_tune`, vector-R BTFDR, or COSMIC `AStd/BStd/Σ̂/ν`.
- No JSON carries generator version/date/seed metadata, so staleness is undetectable from file contents.

---

## 4. Minor findings (condensed; each verified by the reviewing agent, spot-checked where marked ✓)

**Frequency-domain estimators**
1. ETFE boxcar smoothing propagates one NaN bin to S−1 neighbors (`freq_etfe.py:50–54`; MATLAB same). Use `nanmean`/masking.
2. ETFE never warns on near-singular Φ_u in SISO/single-input branches (`freq_etfe.py:313–339, 400–405`).
3. BT time-series `noise_spectrum` is not clamped ≥0; the Hann lag kernel has negative sidelobes — **✓ verified min −0.25** on a sinusoid (`freq_bt.py:236`). Also falsifies SPEC §2.3's claim that the biased estimator "guarantees" non-negativity (true only for Bartlett).
4. `y` 2-D + `u` 3-D passes validation, dies deep in `cov` (`validate_data.py:274–280`).
5. `validate_data`: list-`y` + `u=[]` mis-detected then rejected instead of time-series mode (`validate_data.py:102, 232`).
6. ETFE `window_size` reported as N; spec equates ETFE to BT with M=N−1 (both languages; pinned by tests).
7. `method` strings (`'freq_bt'` …) don't match §9's enumerated values (`'sidFreqBT'` …); the PascalCase note covers field *names*, not values.
8. Multi-trajectory `Var{Φ̂_v}` uses `N_eff = L·N` — statistically sensible but unauthorized by §3.4/§3.5 (doc gap).
9. `sid_dft` FFT path would silently evaluate at rounded bins if ever called off-grid (currently unreachable); dead `squeeze_output` code.
10. Welch: NFFT default computed before `Lsub` validation → raw `math domain error` for tiny segments (`freq_map.py:444`); symmetric Hann/Hamming with `Lsub=2` (or spectrogram `window_length≤2`) gives `S₁=0` → silent all-NaN.
11. Welch MIMO: no uncertainty (all-NaN std) and no PSD projection — undocumented second-class path (`freq_map.py:176–190`).
12. Welch omits the DC bin `tfestimate` returns (deliberate, but §6.10's compat claim should say so).
13. MATLAB `sidFreqBT` MIMO sets the whole `G(ω_k)` to NaN; spec says the affected row (`sidFreqBT.m:219`).
14. MATLAB `sidSpectrogram` rejects cell-array input required by §2.3/§7.2 with a misleading `sid:complexData` error (`sidSpectrogram.m:92–108`).
15. NaN-clamp semantics differ across languages: MATLAB `max(NaN,0)=0`, NumPy `maximum(NaN,0)=NaN` — a cross-port divergence on degenerate data that no reference vector pins.

**COSMIC / LTV stack**
16. Exactly singular Λ pivot crashes with raw `LinAlgError` where SPEC §8.3.4 promises warn-and-return (`ltv_cosmic_solve.py:118–122`; MATLAB warns+returns Inf/NaN — also a cross-port divergence). Reachable from the always-tried λ=1e-3 auto-grid point on rank-deficient data. ✓
17. Invalid λ string → raw `ValueError` not `SidError('bad_lambda')` (`ltv_disc.py:521`); `lambda_grid` never validated (negative candidates accepted); degenerate L-curve (grid < 3 points, flat curve) silently returns `grid[0]`; curvature finite differences assume uniform log-spacing and an unsorted grid is not sorted.
18. `ltv_disc_tune` validation loss divides by N+1 including the zero k=0 row — `√(N/(N+1))` bias vs §8.4.3 (argmin unaffected; both languages).
19. `_frequency_tune` derives N from the first trajectory only — var-len input with a short first trajectory mis-scores every λ (`ltv_disc_tune.py:312–351`); no train/val shape validation (raw `ValueError`/`IndexError` on mismatch).
20. User-supplied `noise_cov` neither copied (aliasing, ✓ `shares_memory` true) nor validated for symmetry/PSD.
21. Silent double DoF fallback (`max(total_obs, 1)`) beyond the spec's single fallback; `extract_std` takes `sqrt` without clamping tiny negative P diagonals (NaN stds near singularity). Both languages.
22. `ltv_disc_io`: convergence denominator `max(|J_prev|, 1.0)` is an absolute test for J<1 (spec: relative); non-convergence at `max_iter` is silent (no warning, no flag); fast path returns `iterations=0` but `cost` length 1; `trust_region` value not validated to [0,1]; R/Q never checked SPD; MATLAB `MaxIter=0` crashes on undefined variables.
23. `ltv_state_est` mutates the caller's input lists (`ltv_state_est.py:296–297`); 2-D A raises bare `IndexError`.
24. `blk_tri_solve` pivot-conditioning warning checks lag the actual use (block 0 used before checked; last block never checked); `np.linalg.cond` (full SVD) per block is needlessly expensive.
25. `lti_freq_io`: DC-bin extrapolation differs between the two consumers (`2G(ω₁)−G(ω₂)` vs `G(ω₁)`) and injects a uniform bias (~4e-4 measured) into every Markov parameter — spec gap; `_ho_kalman` divides by `√σ` without a positivity guard; user `horizon` silently clamped; orphan dead expression at `lti_freq_io.py:391`.
26. `model_order` gap search truncated to the first half of the singular values (`min(·, r/2)`) — orders > r/2 undetectable regardless of data (both languages); undocumented time-series fallback treats a spectrum's Hankel rank as plant order (no spec basis); spec'd `'Plot'` option missing in Python.

**Analysis / plotting / packaging**
27. State-space `residual`/`compare` size everything from `model.data_length` — validation data of different length dies with bare `IndexError` (§15.7 explicitly validates on held-out data); ss-`residual` with `u=None` raises `TypeError` though §14.5 allows it.
28. Multi-trajectory `residual` averages residuals across trajectories (shrinking variance L× against an unadjusted `2.58/√N` bound — anticonservative) and cross-correlates against trajectory 1's input only (both languages).
29. `residual` doesn't validate `max_lag` (0 → vacuous pass; negative → bare `ValueError`).
30. `freq_domain_sim` zeroes DC and all bins below the model grid → predictions are exactly zero-mean and high-passed; whiteness tests fail on any non-detrended data regardless of model quality (documented internally, not in `residual`'s docstring); even-N Nyquist bin left complex and truncated by `np.real` (harmless, implicit).
31. `pyproject.toml`: `license = "MIT"` (PEP 639) requires setuptools ≥ 77 but `requires = ["setuptools >= 68"]` — ✓ verified build failure with setuptools 68; CI passes only via latest isolated builds.
32. Missing `show_confidence` option (§11.3) — conflated with `confidence=0`. Band math itself is spec-conformant.
33. Six docstrings claim "not yet in SPEC.md" for features that are specified (§10.1, §11, §13–§15) — defeats the header-checker's traceability goal.
34. MATLAB `sidMapPlot` Welch title `sprintf` runs out of arguments (`M=` truncation); Python renders `M=None`.
35. `conftest.load_reference` is dead code (test file defines its own `_load`).

**Spec / docs / examples / CI**
36. SPEC §8.14 references `spec/lpv_extension_theory.md`, which does not exist. ✓
37. SPEC.md:899 `Θ_k = D(k)ᵀX'(k)ᵀ` is dimensionally impossible (should be `D(k)ᵀX'(k)`; code is correct).
38. SPEC §8.2's input table omits the `Uncertainty` / `NoiseCov` / `CovarianceMode` options that §8.5/§8.9 presuppose.
39. ½-factor mismatch between SPEC §8.4.2 and `automatic_tuning.md` §2.1 (code follows the latter; corner unaffected); `uncertainty_derivation.md` says covariance-mode default is `full`, SPEC and code say `diagonal`.
40. EXAMPLES.md:13 references a stale work branch (`claude/smd-util-helpers`). All 11 mandated examples exist in both suites with the required helpers (conformance of §3 parameter values not exhaustively re-verified — see §6).
41. `testdata/README.md`: schema key documented as `function` but files use `function_name`; "16 significant digits" is actually 17; "Python … will be added" is stale (it exists); no per-file inventory.
42. MATLAB `tests.yml` triggers only on `matlab/**.m` paths — edits to `testdata/generate_reference.m` alone don't regenerate references.
43. `docs/TODO.md` is fully struck-through (fine); roadmap current; binder config fine.

---

## 5. Test-coverage assessment

**Strengths.** The FFT fast-path regression suites in both languages (large-M envelope, FFT-vs-direct at 1e-10) are exemplary; hand-computed covariance/window values are pinned; MATLAB has a Monte-Carlo COSMIC-uncertainty calibration test, an N=2 OLS closed-form check, and a finite-difference frozen-TF Jacobian check; Python's variable-length `freq_map` per-segment filtering test is unusually strong; reference-vector design is RNG-safe (inputs stored, outputs recomputed).

**Systematic gaps** (each maps to at least one finding above):
- **Zero Python tests** for COSMIC `uncertainty=True`, `covariance_mode`, `noise_cov`, or list/var-len inputs — the entire §8.9 path rides on one mis-generated reference JSON that isn't even read. Port MATLAB's `test_sidLTVdiscUncertainty` (14 tests) and `test_sidLTVdiscVarLen` (7 tests).
- No MIMO numerical-correctness tests anywhere (shapes only): no known-2×2-plant check of `Ĝ = Φ_yu Φ_u⁻¹`, no PSD-clamp test, no MIMO plotting smoke test, no MIMO `freq_map` test at all.
- No degenerate-excitation tests: nothing exercises the singularity guards, NaN substitution, σ=Inf, or the §10.3 constant-input row; two tests (one per language) actively pin the coherence-floor spec violation.
- No warning assertions at all (window reduction, trimming, N<10, near-singular).
- No BTFDR-vs-BT identity test at constant `R = 2π/M` (§5.1's free oracle — passes at 1e-14 when run manually).
- No `trust_region` test; `test_convergence` is vacuous (`iterations <= max_iter` is true by construction); no smoother-vs-dense-reference test (loose 0.05–0.1 heuristics would not catch C1 even with list inputs).
- `model_order` plants all have zeros at the origin or are biproper — the one family blind to S3; no MIMO Hankel test; `lti_freq_io` accuracy tested only at H=I.
- Tolerance sloppiness: `test_dft.py` accepts 5 % where 1e-13 is expected (Python and MATLAB `test_sidDFT` Test 3 alike).
- CI never tests the declared minimum dependencies (numpy 1.22 / scipy 1.8).

---

## 6. Independently verified as correct (highlights)

Biased covariances (incl. ensemble 1/(L·N) and negative-lag transposition); Hann window and `C_W = 3M/4`; FFT fast-path buffer sizing/wrapping/bin extraction in both languages (hand-traced + machine-precision cross-checks, K>1 regime included); direct-DFT path incl. cross-covariance conjugation; SISO/MIMO `Ĝ` algebra; coherence and all §3 variance formulas; ETFE conventions (1-based DFT phase, H1 averaging, periodogram 1/N); BTFDR ≡ BT at constant R; outer segmentation and time vectors (`freq_map`/`spectrogram` aligned, matching scipy exactly); spectrogram PSD scaling (bit-exact vs scipy); Welch transfer/coherence (machine precision vs scipy, correct phase); COSMIC data matrices, λ-indexing across all four recursions, cost triple, forward/backward solve (dense-verified), left/right Schur extraction (dense-verified *given the spec's matrix* — see S1), DoF/noise-cov formulas, AStd/BStd transposition; RTS normal equations for uniform horizons (hand-derived, dense-verified with general R, Q); generic `blk_tri_solve` on non-uniform blocks; Ho-Kalman factors and H-basis transform; frozen-TF value and rank-1 uncertainty propagation; `detrend` (incl. partial last segment); residual whiteness/independence lag conventions and bounds; NRMSE fit formula; `freq_domain_sim` conjugate symmetry for even/odd N; §9 output-struct field mapping; the full 22-file reference-vector generate/validate loop.

**Residual limitations of this review:** `spec/cosmic/online_recursion.md` was only skimmed (feature deferred to v2); EXAMPLES.md §3's per-example parameter values were checked for presence/structure, not re-derived line-by-line; the MATLAB implementation was executed under GNU Octave 8.4.0 (suite + validator + finding reproductions), not under MATLAB proper.

---

## 7. Recommended priorities

1. **Correctness of shipped numbers:** C1 (var-len smoother terminal block), S3 (model_order lag-0), S4 (compare fit), S5 (BTFDR MIMO guards), S6a (Welch Nyquist).
2. **Crashes on documented inputs:** C2, C3, C4, C5, F-minor 10 (all small, local patches).
3. **Spec decisions, then code:** S1 (posterior scaling — the most consequential; also fixes the IO cost inconsistency), S6b/c (Welch convention), S7 (degenerate-σ convention — pick Inf per spec or amend spec; delete the two tests pinning the violation), S9 (IO-frozen contract), §10.2-vs-§10.3 threshold.
4. **Feature quarantine:** S2 — disable or rewrite `trust_region` before anyone relies on it.
5. **Test/infra debt:** S10 (fix `generate_reference.m` P-generation, wire the four orphaned JSONs, unify tolerance handling), port the MATLAB uncertainty/var-len suites to Python, add MIMO/degenerate/warning tests, add the BTFDR≡BT and Welch-vs-scipy oracles.
