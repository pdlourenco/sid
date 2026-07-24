# Model-order gap search — noise floor, cliff, and cap

**Date:** 2026-07-24
**Analyzed at:** `7fb3760` (`main` HEAD at analysis time). File:line anchors below refer to this commit.
**Scope:** `python/sid/model_order.py`, `matlab/sid/sidModelOrder.m`, `matlab/tests/test_sidModelOrder.m`, `spec/SPEC.md` §8.12.12.
**Method:** static reading of both ports plus independent numerical reproduction under the `.venv` Python (NumPy/SciPy). The proposed unified rule was re-implemented standalone and run against the *actual* `model_order` SVDs, so every `n` below is computed, not argued. Feeds: SPEC §8.12.12 (contract), [ADR-0004](../decisions/ADR-0004-model-order-gap-convention.md) (decision), issue #160 (register).

Motivated by issue #160, which split this cluster out of #139 (PR #157, the lag-0 Hankel fix). The gap method (`'Threshold'` unset) is the default; the threshold method is unaffected throughout.

---

## 1. Summary

The gap search carries two undocumented guards — a noise floor and a `floor(m/2)`
cap — that SPEC §8.12.12 does not mention (the spec says `k = 1..r-1`, which is
un-runnable on estimated data). Three coupled problems were reported. On the
post-#157 tree, measurement shows:

1. **Cross-language divergence (real, the live bug).** MATLAB includes the
   resolvable→floor cliff ratio and Python excludes it, so the default gap method
   returns different `n` on identical singular values for low-order plants.
2. **Artifact singular values (largely resolved upstream).** The lag-1
   reconstruction in #139/PR #157 shrank finite-grid/DC artifacts to near-`ε`; the
   machine-eps floor now recovers the true order on exact plants without a
   data-aware floor.
3. **The `floor(m/2)` cap (still load-bearing).** It is the only guard that bounds
   `n` on no-cliff, noise-dominated data.

The consequence — recorded in ADR-0004 — is the opposite of #160's first sketch:
**keep the machine-eps floor, exclude the cliff (align MATLAB to Python), and
retain the cap.** No data-aware floor; no lifted cap.

## 2. The three problems, measured

Gap logic at [model_order.py:275-293](../../python/sid/model_order.py) and
[sidModelOrder.m:233-250](../../matlab/sid/sidModelOrder.m); both carry the
`issue #160` NOTE.

### 2.1 Cross-language divergence (Problem 2)

The floor identifies the last resolvable singular value `σ_L`. MATLAB
(`lastSig = find(..., 1, 'last')`, 1-based) sets the ratio bound to `L` and so
includes `σ_L/σ_{L+1}`; Python (`last_sig`, 0-based) sets it to `L-1` and excludes
it. That cliff divides a resolvable value by a numerical zero, so it is an
enormous, meaningless ratio that dominates any argmax.

Exact plant `G(z) = (z+0.5)/(z²−1.2728z+0.81)` (strictly proper, one finite zero,
true `n=2`), grid `ω = linspace(π/200, π, 200)`, horizon 50 — singular values
`[6.05, 4.58, 2.76e-3, ~1e-15, …]` (`σ₃/σ₁ = 4.6e-4`):

| implementation | n | why |
|---|---|---|
| current Python | **2** | excludes cliff → argmax is `σ₂/σ₃` |
| current MATLAB | **3** | includes `σ₃/σ₄` (≈`3e12`) → argmax there |

Same `Σ`, different `n`. This is the divergence #160 names.

### 2.2 Artifact singular values (Problem 3) — mitigated by PR #157

On the **pre-#157** reconstruction (DC ≈ `Re(Ĝ(ω₁))`, lag-0 retained) the same
plant produced artifacts at `σ₃/σ₁ ≈ 6.5e-2` and `σ₄/σ₁ ≈ 1.7e-5` — far above the
machine-eps floor (`√m·ε ≈ 1.6e-15`), so the gap search latched onto them and
over-estimated the order. That is what motivated #160's data-aware floor.

On the **post-#157** tree (DC ≈ `Re(2Ĝ(ω₁) − Ĝ(ω₂))`, lag-0 discarded) the
artifact collapses to `σ₃/σ₁ = 4.6e-4` and the next value is already `O(ε)`. The
machine-eps floor now recovers the true order directly:

| plant (exact) | machine-eps floor, cliff excluded, capped | `ρ = 1e-3` floor |
|---|---|---|
| finite-zero (true 2) | **2** ✓ | **1** ✗ (cuts real `σ₂`) |
| third-order (true 3) | **3** ✓ | **2** ✗ |

A `1e-3` relative floor is now strictly harmful: it removes a small-but-real
singular value and, combined with cliff exclusion, drops the order by one.

### 2.3 The cap (Problem 1) — orthogonal to the floor

White-noise BT estimate (`rng 1300`, `N=2000`, `sidFreqBT` window 60), `m=42`:

- capped (`floor(m/2)=21`): `n=19` — passes Test 13's `n < nSV`
  ([test_sidModelOrder.m:248](../../matlab/tests/test_sidModelOrder.m)).
- uncapped: `n=40`.

All 42 values sit above the machine-eps floor *and* above a `1e-3` floor — a
noise-dominated spectrum decays slowly with **no** resolvable/floor boundary, so
no floor engages. Only the cap prevents `n ≈ m`. Floor and cap therefore guard
different regimes and neither replaces the other.

## 3. The proposed rule, verified

Language-neutral: `L =` #{`σ_k > σ_1·√m·ε`}; `K = min(L−1, floor(m/2))`, `K ≥ 1`;
`n = argmax_{k=1..K} σ_k/σ_{k+1}`. Because `L` and `K` are magnitude counts, the
ports agree by construction.

Noisy BT estimates confirm the floor choice is inert where it should be — for a
noisy 2nd-order plant across 6 seeds, a machine-eps, `1e-8`, and `1e-3` floor all
return `n=2`; a noisy 4th-order Butterworth returns `n=3` under all three (a BT
identifiability limit at the tested SNR/window, floor-independent). The floor only
ever *differs* on near-exact data, and there the larger floors regress (§2.2).

## 4. Cross-language reference vector

To pin the agreement, a strictly-proper finite-zero plant (the exemplar #160
recommends over the current biproper `generate_reference.m` test-case 7):

```
G(z) = (z + 0.5) / (z² − 1.2728 z + 0.81)
ω = linspace(π/64, π, 64),  horizon = 21
σ = [6.047747334326, 4.577697211291, 2.755498562572e-3, 9.805e-16, …]
expected n = 2   (both languages, after the MATLAB cliff fix)
```

Under the current code Python returns 2 and MATLAB returns 3; after aligning
MATLAB to the count-based rule both return 2. The exact stored vector/grid is
finalized when `testdata` references are regenerated (ties into the #145d/e
test-case-7 swap).

## 5. What lands

Per ADR-0004 and CONTRIBUTING (spec first, ADR + spec in one PR):

- **SPEC §8.12.12** — replace the `k = 1..r-1` step 5 with the count-based rule
  (floor, cliff exclusion, cap) and its rationale. *(done in this branch)*
- **MATLAB** — one-line fix: bound the ratio search by `L-1` (`lastSig-1`), not
  `lastSig`; remove the `issue #160` NOTE (ADR-0003).
- **Python** — already excludes the cliff (`last_sig == L-1`); reframe the comment,
  remove the NOTE.
- **Reference vector** — add the finite-zero plant above to both suites; assert
  identical `n=2`.
- **issue #160** — close as resolved, noting the "lift the cap" direction is
  superseded (ADR-0004 alternatives).
