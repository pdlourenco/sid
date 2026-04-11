# sid — Examples Specification

**Version:** 1.0.1
**Date:** 2026-04-12
**Reference:** Companion to [`SPEC.md`](SPEC.md). Where `SPEC.md` defines
the binding behavior of the algorithmic functions, this document defines
the binding structure of the example suite shipped with every language
port.

---

> **Implementation status.** This specification is authored against the
> Python example suite on branch `claude/smd-util-helpers` (Phase B). The
> MATLAB port, Julia port, and any future language port MUST conform to
> this document. The Python port is the v1.0.0 reference implementation
> of this spec.

---

## 0. Preamble

### 0.1 Scope

This document is the **source of truth** for the 11 examples that ship
with the sid toolbox in `{lang}/examples/`. It specifies:

- The physical plants the examples use and their parameters (§1).
- The mathematical definition of the `util_msd*` helper functions that
  build the plants (§2).
- For each example: the goal, the plant, the excitation, the required
  pedagogical sections, the `sid.*` function calls that MUST appear, and
  the plots and prints that MUST be produced (§3).
- The conventions every port must follow (§4).
- The definition of "equal" across languages (§5).
- Versioning and change control rules (§6).

### 0.2 Requirement levels

This document uses the RFC 2119 / BCP 14 key words **MUST**, **MUST NOT**,
**SHOULD**, **SHOULD NOT**, and **MAY**. MUST statements are binding —
any implementation that violates one is non-conformant. SHOULD statements
are recommendations: violations must be justified and documented. MAY is
purely informational.

Throughout this document, MUST statements apply to **every language port**.

### 0.3 Binding vs advisory aspects

| Aspect of an example | Requirement level |
|---|---|
| Plant parameters (m, k, c, Ts, N) | MUST match §1 |
| Plant topology (n-mass chain, SDOF, Duffing) | MUST |
| Ordered list of pedagogical sections | MUST |
| Section topic (what each section teaches) | MUST |
| Exact section-header prose | SHOULD |
| Exact markdown / block-comment body prose | SHOULD |
| `sid.*` function calls invoked in each section | MUST |
| Options passed to those calls (window size, grid, λ, etc.) | MUST |
| Plot kinds (Bode, spectrum, map, spectrogram, line chart) | MUST |
| Plot title *topic* and axis label *quantities* | MUST |
| Plot title/label exact wording | SHOULD |
| Plot cosmetics (colors, line widths, figure sizes, fonts) | MAY |
| Printed line *semantics* (what each print conveys) | MUST |
| Exact printed string format | SHOULD |
| Helper function mathematical definition | MUST |
| Helper function names per language | SHOULD (see §2) |
| RNG seed choice | MAY |
| Numerical values produced at runtime | MAY |

In short: a port MUST match the pedagogical skeleton and the `sid.*`
call graph; it MAY use language-idiomatic phrasing, styling, and RNG
seeds.

### 0.4 Why a unified spring-mass-damper plant family?

Every example in this specification is built on a common family of
spring-mass-damper (SMD) plants, ranging from a 1-DoF oscillator to a
3-mass chain, with variants that are linear, linearly time-varying,
and nonlinear (Duffing hardening). Using one physical plant family
across the tutorial suite keeps the narrative concrete and cumulative:
a reader working through the examples sees the same physical systems
again and again, picking up new analysis tools on familiar territory.

Arbitrary `lfilter` coefficients or handcrafted state-space matrices
are explicitly rejected — they have no physical interpretation, their
pole locations are chosen numerically rather than from first
principles, and they make the pedagogical content artificial.

### 0.5 Document structure

| Section | Contents |
|---|---|
| §1 | Plant catalog (five plants A–E) |
| §2 | Helper API (`util_msd*` mathematical contract) |
| §3 | Example catalog (eleven examples) |
| §4 | Conventions (auto-discovery, self-containment, RNG, naming) |
| §5 | Cross-language "equivalence" (what "equal" means) |
| §6 | Versioning and change control |
| §7 | References |

---

## 1. Plant Catalog

Every example in this specification uses one or more of the five plants
defined below. Each plant has a fixed set of physical parameters; the
parameters MUST NOT be changed by a language port without a corresponding
version bump of this specification (see §6).

The plants are referred to by their catalog letter (Plant A, Plant B, …)
throughout §3. When an example refers to a plant it MUST use the
parameters defined here.

### 1.1 Plant A — 1-DoF spring-mass-damper (SDOF)

**Physical story.** A single mass connected to a fixed wall through a
linear spring and a viscous damper — the canonical textbook oscillator.
A force is applied at the mass and the measured output is the mass
position.

**Continuous-time dynamics.**

```
m x''(t) + c x'(t) + k x(t) = u(t)
```

**State vector.** `[x; v]` with `x` position and `v = x'` velocity.

**Parameters.**

| Symbol | Value | Units | Derived |
|---|---|---|---|
| `m`  | `1.0`   | kg       | |
| `k`  | `100.0` | N/m      | `ω_n = sqrt(k/m) = 10 rad/s ≈ 1.59 Hz` |
| `c`  | `2.0`   | N·s/m    | `ζ = c / (2 sqrt(k m)) = 0.1`, `Q = 5` |
| `F`  | `[[1]]` | (n×q)    | force applied directly at the single mass |
| `Ts` | `0.01`  | s        | `fs = 100 Hz`, Nyquist `= π/Ts ≈ 314 rad/s` |

**Discretization testpoint.** Under `util_msd(m, k, c, F, Ts)` (see §2.1):

```
Ad ≈ [[ 0.9950372995,  0.0098841706],
      [-0.9884170600,  0.9752689583]]

Bd ≈ [[ 4.9627005463e-05],
      [ 9.8841705996e-03]]
```

A conforming implementation MUST reproduce these to at least `1e-9`
relative tolerance for this exact `(m, k, c, F, Ts)` tuple.

**Used by.** §3.1 `siso`, §3.2 `etfe`, §3.5 `method_comparison`, §3.9
`ltv_disc` (main LTV section, plus the Duffing linearization reuses
Plant A's `m`, `k_lin`, `c`).

---

### 1.2 Plant B — 2-mass chain

**Physical story.** Two equal masses in series, anchored to a wall
through the first spring. Naturally exhibits a 2×2 force-to-position
MIMO response with cross-coupling through the second spring.

**Topology.** `wall──k₁──m₁──k₂──m₂`

**Continuous-time dynamics.**

```
M x''(t) + C x'(t) + K x(t) = F u(t)
```

with `M = diag(m)` and tridiagonal `K`, `C` built from `k_spring`,
`c_damp` per §2.1.

**State vector.** `[x₁, x₂, v₁, v₂]`.

**Parameters.**

| Symbol | Value | Units | Notes |
|---|---|---|---|
| `m` | `[1, 1]`    | kg    | equal masses |
| `k` | `[100, 80]` | N/m   | `k[0]` wall→m₁, `k[1]` m₁→m₂ |
| `c` | `[2, 2]`    | N·s/m | |
| `Ts`| `0.01`      | s     | |

**Normal modes.** `ω_modes ≈ [5.97, 14.98] rad/s` (obtained as square
roots of the sorted eigenvalues of `M⁻¹K`).

**Input distribution matrix `F`.** Depends on example:

- For `mimo` (§3.6), `coherence` (§3.4) sections that use a second
  input for an unmodeled disturbance: `F = I₂` (two independent force
  channels).
- For `freq_map` (§3.7), `multi_trajectory` (§3.10 sections 1–2, 4),
  `output_cosmic` (§3.11): `F = [[1], [0]]` (force at mass 1 only).

**Discretization testpoint.** Under `util_msd(m=[1,1], k=[100,80],
c=[2,2], F=[[1],[0]], Ts=0.01)`:

```
Ad ≈ [[ 9.9116061487e-01,  3.9125769853e-03,  9.7740440835e-03,  1.1093285197e-04],
      [ 3.9060177716e-03,  9.9605789532e-01,  1.1093285197e-04,  9.8882565423e-03],
      [-1.7504533069e+00,  7.7304889853e-01,  9.5228630424e-01,  2.3238799448e-02],
      [ 7.7109261003e-01, -7.8218589523e-01,  2.3238799448e-02,  9.7650324794e-01]]

Bd ≈ [[4.9268081419e-05],
      [3.6086910233e-07],
      [9.7740440835e-03],
      [1.1093285197e-04]]
```

A conforming implementation MUST reproduce these to at least `1e-9`
relative tolerance.

**Used by.** §3.4 `coherence`, §3.6 `mimo`, §3.7 `freq_map`, §3.10
`multi_trajectory` (sections 1, 2, 4), §3.11 `output_cosmic`.

---

### 1.3 Plant C — 3-mass chain

**Physical story.** Three masses in series, anchored to a wall through
the first spring. Three well-separated bending-mode-like resonances in
the bottom decade of the spectrum; the canonical plant for
demonstrating frequency-dependent resolution.

**Topology.** `wall──k₁──m₁──k₂──m₂──k₃──m₃`

**Parameters.**

| Symbol | Value | Units | Notes |
|---|---|---|---|
| `m` | `[1, 1, 1]`     | kg    | equal masses |
| `k` | `[300, 200, 100]` | N/m | stiffness *decreases* toward the free end |
| `c` | `[8, 8, 8]`     | N·s/m | |
| `F` | `[[1], [0], [0]]` | (3×1) | force at mass 1 only |
| `Ts`| `0.01`          | s     | |

**Normal modes.** `ω_modes ≈ [6.45, 15.15, 25.08] rad/s`.

**Physical interpretation.** A slender cantilever-like structure with
diminishing stiffness toward the tip — the three modes are close enough
to challenge a short-window BT estimator but far enough apart to be
resolved by a long window or by frequency-dependent resolution.

**Discretization testpoint.** Under `util_msd`:

```
Ad ≈ [[ 0.97662704,  0.00906461,  0.00012791,  0.00917488,  0.00038846,  0.00001056],
      [ 0.00882163,  0.98617894,  0.00460777,  0.00038846,  0.00921581,  0.00038430],
      [ 0.00024562,  0.00448490,  0.99526154,  0.00001056,  0.00038430,  0.00960491],
      [-4.50974587,  1.71949205,  0.03779075,  0.83293672,  0.07633266,  0.00315117],
      [ 1.64892960, -2.64861857,  0.88315012,  0.07633266,  0.84490819,  0.07525978],
      [ 0.07158313,  0.84731120, -0.92206093,  0.00315117,  0.07525978,  0.92149667]]

Bd ≈ [[4.7268e-05], [1.3055e-06], [2.6458e-08],
      [9.1749e-03], [3.8846e-04], [1.0555e-05]]
```

Trailing digits shown to the spec's working precision. A conforming
implementation MUST reproduce these to at least `1e-6` relative
tolerance.

**Used by.** §3.3 `freq_dep_res`.

---

### 1.4 Plant D — SDOF for spectrogram (high-frequency SDOF)

**Physical story.** A 1-DoF oscillator scaled to a structural-acoustic
frequency range so it can be driven by an audible chirp and analyzed
with a short-time FFT. Models a stiff mechanical resonator
(`ω_n ≈ 32 Hz`) with moderate damping — the sort of plant you would
hit with an impact hammer in a lab modal test.

**Parameters.**

| Symbol | Value | Units | Derived |
|---|---|---|---|
| `m` | `1.0`     | kg     | |
| `k` | `4·10⁴`   | N/m    | `ω_n = 200 rad/s ≈ 31.83 Hz` |
| `c` | `20.0`    | N·s/m  | `ζ = 0.05`, `Q = 10` |
| `F` | `[[1]]`   | (1×1)  | force at the mass |
| `Ts`| `0.001`   | s      | `fs = 1000 Hz`, Nyquist `= 500 Hz` |

**Discretization testpoint.** Under `util_msd`:

```
Ad ≈ [[ 0.98019872,  0.00098348],
      [-39.33916499,  0.96052913]]

Bd ≈ [[ 4.9503206884e-07],
      [ 9.8347912463e-04]]
```

Relative tolerance `1e-9`.

**Used by.** §3.8 `spectrogram`, §3.10 `multi_trajectory` (section 3 —
chirp-in-noise).

---

### 1.5 Plant E — Duffing hardening SDOF (nonlinear)

**Physical story.** A 1-DoF oscillator with a **cubic stiffness** term
superposed on the linear spring. Physically models a hardening rubber
mount, a geometrically-nonlinear thin plate, or a pre-stressed
structure whose restoring force grows faster than linearly with
displacement.

**Continuous-time dynamics.**

```
m x''(t) + c x'(t) + k_lin x(t) + k_cub x(t)³ = u(t)
```

**State vector.** `[x; v]`.

**Parameters.**

| Symbol  | Value   | Units   | Notes |
|---|---|---|---|
| `m`     | `1.0`   | kg      | same mass as Plant A |
| `k_lin` | `100.0` | N/m     | linear component (Plant A's spring) |
| `k_cub` | `1·10⁵` | N/m³    | hardening; dominates at `|x| ≳ 0.03 m` |
| `c`     | `2.0`   | N·s/m   | linear damping (Plant A's damper) |
| `F`     | `[[1]]` | (1×1)   | force at the mass |
| `Ts`    | `0.01`  | s       | same as Plant A |

**Effective natural frequency** (small-signal linearization around
displacement `x`):

```
ω_eff(x) = sqrt( (k_lin + 3 k_cub x²) / m )
```

At `x = 0`: `ω_eff = 10 rad/s` (Plant A's baseline).
At `x = 0.05 m`: `ω_eff ≈ 21.8 rad/s` (a ~2× frequency shift).

**Integration.** This plant cannot be discretized with exact ZOH
(`expm`) because it is nonlinear. It MUST be simulated with the RK4
integrator specified in §2.3, operating on the continuous-time RHS.

**Integration testpoint.** With `m=1`, `k_lin=100`, `k_cub=1e5`,
`c=2`, `F=[[1]]`, `Ts=0.01`, `substeps=1`, initial state `x0 = 0`,
and impulse input `u = [1, 0, 0, 0, 0]` (one step):

```
x_traj[1] ≈ [4.9626666641e-05, 9.8841633080e-03]
x_traj[2] ≈ [1.4707706719e-04, 9.5906651394e-03]
x_traj[3] ≈ [2.4114284377e-04, 9.2080973779e-03]
x_traj[4] ≈ [3.3096036948e-04, 8.7419986134e-03]
x_traj[5] ≈ [4.1572503242e-04, 8.1986213117e-03]
```

Relative tolerance `1e-6`. A conforming implementation MUST also
satisfy the **reduction-to-linear** test: setting `k_cub = 0` and
calling the RK4 integrator with `substeps ≥ 4` MUST reproduce the
exact ZOH result from §2.1 on the same input to within `1e-8`
absolute (fifth-order RK4 truncation error at `Ts = 0.01`).

**Used by.** §3.7 `freq_map` (Duffing section), §3.9 `ltv_disc`
(Duffing linearization section).

---

## 2. Helper API (`util_msd*`)

Every language port MUST provide three helper functions that build the
plants defined in §1 from their physical parameters. The mathematical
definitions below are binding. Function *names* follow the local
language's naming convention (`util_msd` in Python, `util_msd.m` in
MATLAB, etc.); what matters is that exactly three functions with the
semantics defined here are available to the example suite.

The three helpers are:

- **§2.1** The LTI helper — builds `(Ad, Bd)` for a time-invariant
  n-mass chain via exact zero-order-hold discretization.
- **§2.2** The LTV helper — builds per-step `(Ad(k), Bd(k))` stacks for
  chains whose parameters vary in time.
- **§2.3** The nonlinear simulator — integrates an n-mass chain with
  Duffing-style cubic stiffness via fixed-step RK4.

### 2.1 LTI helper

#### 2.1.1 Construction of `K` and `C`

Given physical parameters `m`, `k_spring`, `c_damp` (all length-`n`
vectors with positive entries), construct `n × n` tridiagonal matrices:

```
K[i, i]   = k_spring[i] + k_spring[i+1]     for i = 0, …, n-2
K[n-1, n-1] = k_spring[n-1]
K[i, i+1] = K[i+1, i] = -k_spring[i+1]      for i = 0, …, n-2
```

with `C` built by the same pattern from `c_damp`. This is the standard
chain topology `wall──k₁──m₁──k₂──m₂──…──kₙ──mₙ`, where `k_spring[0]`
is the wall-to-mass-1 spring and `k_spring[i]` (for `i ≥ 1`) connects
mass `i` to mass `i+1`.

#### 2.1.2 Continuous-time state-space model

With `M = diag(m)`, `M⁻¹ = diag(1/m)`, and the force distribution
matrix `F` of shape `n × q`:

```
Ac = [[  0_{n×n},      I_{n×n} ],
      [ -M⁻¹ K,       -M⁻¹ C   ]]    (shape 2n × 2n)

Bc = [[  0_{n×q}              ],
      [  M⁻¹ F                ]]    (shape 2n × q)
```

The state vector is `[x₁, x₂, …, xₙ, v₁, v₂, …, vₙ]`: positions first,
velocities second.

#### 2.1.3 Exact zero-order-hold discretization

For sample time `Ts > 0`:

```
Ad = expm(Ac · Ts)
Bd = Ac⁻¹ · (Ad - I_{2n}) · Bc
```

The implementation MUST use a library-quality matrix exponential
(`expm` in MATLAB/Octave/SciPy). Taylor-series or Pade truncations
below 6th order are non-conformant because the discretization testpoints
in §1 are tabulated to `1e-9` relative tolerance.

#### 2.1.4 Interface (binding)

Inputs:

| Name | Shape | Type | Constraint |
|---|---|---|---|
| `m` | `(n,)` | real | entries > 0 |
| `k_spring` | `(n,)` | real | entries > 0 |
| `c_damp` | `(n,)` | real | entries ≥ 0 |
| `F` | `(n, q)` | real | — |
| `Ts` | scalar | real | `Ts > 0` |

Outputs:

| Name | Shape | Type |
|---|---|---|
| `Ad` | `(2n, 2n)` | real |
| `Bd` | `(2n, q)` | real |

The function MUST reject mismatched shapes (e.g., `len(k_spring) ≠
len(m)`) with an error. A one-column `F` passed as a `(n,)` vector
MAY be auto-reshaped.

#### 2.1.5 Numerical testpoints

These are the same testpoints tabulated in §1. A conforming LTI helper
MUST reproduce them:

| Plant | Tuple | Tolerance |
|---|---|---|
| Plant A | `util_msd([1], [100], [2], [[1]], 0.01)` | rel 1e-9 |
| Plant B | `util_msd([1,1], [100,80], [2,2], [[1],[0]], 0.01)` | rel 1e-9 |
| Plant C | `util_msd([1,1,1], [300,200,100], [8,8,8], [[1],[0],[0]], 0.01)` | rel 1e-6 |
| Plant D | `util_msd([1], [4e4], [20], [[1]], 0.001)` | rel 1e-9 |

Ports SHOULD also include `n = 2` and `n = 5` regression tests (not
tabulated here) whose expected outputs are computed from the
continuous closed form against the modal decomposition.

---

### 2.2 LTV helper

#### 2.2.1 Contract

The LTV helper builds a stack of per-step `(Ad(k), Bd(k))` matrices for
an n-mass chain whose parameters vary over `N` sample steps. It MUST
be semantically equivalent to calling the LTI helper (§2.1) once per
time index:

```
for k in 0..N-1:
    Ad_stack[:, :, k], Bd_stack[:, :, k] = util_msd(
        m(k), k_spring(k), c_damp(k), F(k), Ts
    )
```

A conforming implementation MAY take a fast path when all of `m`,
`k_spring`, `c_damp`, `F` are time-invariant (call `util_msd` once and
broadcast the result along a newly introduced time axis).

#### 2.2.2 Interface (binding)

Inputs:

| Name | Allowed shapes | Notes |
|---|---|---|
| `m` | `(n,)` or `(n, N)` | constant or per-step |
| `k_spring` | `(n,)` or `(n, N)` | constant or per-step |
| `c_damp` | `(n,)` or `(n, N)` | constant or per-step |
| `F` | `(n, q)` or `(n, q, N)` | constant or per-step |
| `Ts` | scalar | `Ts > 0` |
| `N` (optional) | integer | required only when every other input is constant |

**`N` inference.** The function MUST infer `N` from the first time-
varying input it sees. If every input is constant, the function MUST
require the caller to pass `N` explicitly; in that case the result is
the LTI pair replicated `N` times along the new axis.

**Dimension consistency.** The function MUST raise an error if the
time axes of two different time-varying inputs disagree (e.g.,
`k_spring` has 80 columns but `c_damp` has 100).

Outputs:

| Name | Shape |
|---|---|
| `Ad_stack` | `(2n, 2n, N)` |
| `Bd_stack` | `(2n, q, N)` |

#### 2.2.3 Testpoint

With `m = [1]`, `c_damp = [2]`, `F = [[1]]`, `Ts = 0.01`, and
`k_spring` a single-row `(1, 10)` array built as
`linspace(200, 50, 10)` (1-DoF with `k` ramping 200 → 50 over 10
samples):

```
Ad_stack[:, :, 0] ≈ [[ 0.9900828576,  0.0098676943],
                     [-1.9735388683,  0.9703474690]]

Ad_stack[:, :, -1] ≈ [[ 0.9975176169,  0.0098924149],
                      [-0.4946207456,  0.9777327870]]
```

Relative tolerance `1e-9`.

#### 2.2.4 LTI collapse check

Calling the LTV helper with 1-D `m`, `k_spring`, `c_damp`, 2-D `F`,
and `N = 50` MUST produce an `Ad_stack` whose every slice
`Ad_stack[:, :, k]` is bit-identical (not merely numerically equal)
to `util_msd(m, k_spring, c_damp, F, Ts)`. A conforming implementation
SHOULD take the fast path in this case.

---

### 2.3 Nonlinear simulator (Duffing-style cubic stiffness)

#### 2.3.1 Continuous-time RHS

Given a state vector `state = [x; v]` of length `2n` with
`x = state[0:n]` (positions) and `v = state[n:2n]` (velocities), and
an instantaneous input vector `u_k` of length `q`, the RHS is:

```
net_force[0]     = -( k_lin[0] · x[0] + k_cub[0] · x[0]³ + c_damp[0] · v[0] )

for i in 1..n-1:
    Δx = x[i] - x[i-1]
    Δv = v[i] - v[i-1]
    f_i = k_lin[i] · Δx + k_cub[i] · Δx³ + c_damp[i] · Δv
    net_force[i-1] += f_i
    net_force[i]   -= f_i

net_force += F · u_k                         (external force injection)

acc = net_force / m                          (element-wise)

d state / dt = [v; acc]
```

This implements the chain topology of §1 with a cubic stiffness
contribution on each spring. Setting `k_cub = zeros(n)` reduces the
RHS exactly to the linear dynamics of §2.1.2.

#### 2.3.2 Integration — fixed-step RK4 with zero-order-hold input

Let `h = Ts / substeps`, where `substeps` is an integer ≥ 1.

For each sample step `k = 0, 1, …, N-1`, hold the input `u[k]`
constant over the interval `[k·Ts, (k+1)·Ts)` and take `substeps`
classic RK4 steps of size `h`:

```
for _ in 1..substeps:
    k₁ = rhs(state,              u[k])
    k₂ = rhs(state + 0.5·h·k₁,   u[k])
    k₃ = rhs(state + 0.5·h·k₂,   u[k])
    k₄ = rhs(state + h·k₃,       u[k])
    state = state + (h/6) · (k₁ + 2 k₂ + 2 k₃ + k₄)
```

The `rhs` function is the one defined in §2.3.1. The same `u[k]` is
used for all four stage evaluations (zero-order hold).

The implementation MUST use this exact tableau. Variable-step integrators
(`ode45`, `RK45` with adaptive step) are non-conformant because their
outputs would differ from the tabulated testpoint in §1.5.

#### 2.3.3 Interface (binding)

Inputs:

| Name | Shape | Notes |
|---|---|---|
| `m`       | `(n,)`    | entries > 0 |
| `k_lin`   | `(n,)`    | linear spring constants |
| `k_cub`   | `(n,)`    | cubic spring constants; zeros → linear plant |
| `c_damp`  | `(n,)`    | entries ≥ 0 |
| `F`       | `(n, q)`  | force distribution |
| `Ts`      | scalar    | `Ts > 0` |
| `u`       | `(N, q)`  | input time series |
| `x0`      | `(2n,)`   | initial state (MAY default to zeros) |
| `substeps` | integer  | RK4 sub-steps per `Ts` (MAY default to 1) |

Output:

| Name | Shape |
|---|---|
| `x_traj` | `(N + 1, 2n)` |

`x_traj[0]` MUST equal `x0`; `x_traj[k+1]` is the state after stepping
through `u[k]` for one sample period.

#### 2.3.4 Numerical testpoint — impulse response

With `m = [1]`, `k_lin = [100]`, `k_cub = [1e5]`, `c_damp = [2]`,
`F = [[1]]`, `Ts = 0.01`, `u = [[1], [0], [0], [0], [0]]`, `x0 = 0`,
`substeps = 1`:

```
x_traj[1] ≈ [4.9626666641e-05, 9.8841633080e-03]
x_traj[2] ≈ [1.4707706719e-04, 9.5906651394e-03]
x_traj[3] ≈ [2.4114284377e-04, 9.2080973779e-03]
x_traj[4] ≈ [3.3096036948e-04, 8.7419986134e-03]
x_traj[5] ≈ [4.1572503242e-04, 8.1986213117e-03]
```

Relative tolerance `1e-6`.

#### 2.3.5 Reduction-to-linear check

Setting `k_cub = zeros(n)` and calling the nonlinear simulator with
`substeps ≥ 4` MUST reproduce the ZOH result of §2.1 on the same
linear `(m, k_lin, c_damp, F, Ts, u)` tuple to within `1e-8` absolute
tolerance. This catches tableau errors and chain-assembly mismatches
between the LTI and nonlinear helpers.

---

## 3. Example Catalog

This section specifies the eleven examples that MUST ship with every
language port. Each example has a **canonical identifier** (e.g.,
`siso`, `freq_dep_res`) that is language-neutral. The file-name mapping
for each identifier is in §4.5.

Per §0.3, the structure below is binding: the plant, excitation,
section order, `sid.*` calls, and plot/print outputs. The exact prose
of markdown cells or block comments is advisory.

Every example MUST:

- Start with a title cell or comment block carrying the canonical
  identifier and a one-paragraph description of the goal.
- Import the local `sid` module and the local `util_msd*` helpers
  (see §2). Additional imports (numpy/scipy, plotting library,
  language-idiomatic tools) are allowed.
- Generate its own data via the helper API — no external data files.
- Produce the required plots and prints.
- Be reproducible within a single language invocation: fix an RNG seed
  at the top of the data-generation section.

### 3.1 `siso` — SISO frequency response with Blackman-Tukey

**Goal.** Estimate the frequency response of a physical SDOF oscillator
from a noisy input/output record using `freq_bt`, and walk through the
standard validation workflow (window-size comparison, detrending,
residual analysis, fit comparison, time-series mode).

**Plant.** Plant A (§1.1) with `N = 2048` samples.

**Excitation.** Unit-variance Gaussian white force `u[k] ~ N(0, 1)`.
Measured output `y[k] = x₁[k] + ε[k]` where `ε[k] ~ N(0, 2·10⁻⁴)` is
additive sensor noise.

**Required sections** (binding order, SHOULD use titles near the
wording in the first column):

| # | Section title (SHOULD) | Topic (MUST) |
|---|---|---|
| 1 | Generate test data | Build `(Ad, Bd)` via §2.1, simulate the state recursion, measure position with added sensor noise. |
| 2 | Estimate frequency response using Blackman-Tukey | Call `freq_bt` with a dense custom frequency grid and a window size large enough to resolve the narrow resonance. |
| 3 | Plot Bode diagram | Display magnitude and phase of the estimate with confidence bands. |
| 4 | Plot noise spectrum | Display the BT noise spectrum estimate. |
| 5 | Compare different window sizes | Demonstrate how the Hann window length trades bias for variance; show that a short window smears the narrow resonance. |
| 6 | Preprocessing: detrend data before estimation | Corrupt the signals with a linear drift and DC offset; show that `detrend` removes the low-frequency bias. |
| 7 | Model validation: residual analysis | Call `residual` and print whiteness/independence pass/fail. Ports MUST include a narrative note that the non-parametric BT estimator has finite-window bias on narrow resonances, so the whiteness test is *expected* to FAIL on this plant — the test is reporting unmodelled residual structure, not a bug in the example. Readers should treat the FAIL as the pedagogical point of the section. |
| 8 | Model validation: compare predicted vs measured | Call `compare` and print NRMSE fit. |
| 9 | Time-series mode (no input) | Re-simulate the plant and hand only the output to `freq_bt` (with `u = None`); show the resonance still appears in the output spectrum. |

**Required `sid.*` invocations.**

| # | Function | Options (MUST) |
|---|---|---|
| 2 | `freq_bt` | `window_size = 200`, custom frequency grid `linspace(0.005, π, 512)`, `sample_time = Ts` |
| 3 | `bode_plot` | on the result from #2 |
| 4 | `spectrum_plot` | on the result from #2 |
| 5 | `freq_bt` | three calls with `window_size ∈ {50, 100, 300}`, same frequency grid as #2 |
| 6 | `detrend` | applied to both the drifted output and the offset input; then `freq_bt` on the detrended pair with the same options as #2 |
| 7 | `residual` | on the #2 result with `(y, u)` |
| 8 | `compare` | on the #2 result with `(y, u)` |
| 9 | `freq_bt` | called with the output only (time-series mode), `window_size = 200`, same frequency grid |

**Required outputs.**

| Kind | Section | Binding content |
|---|---|---|
| Plot (Bode 2-panel) | 3 | Title conveys "Bode diagram / freq_bt estimate"; x-axis frequency, y-axes magnitude (dB) and phase (rad or deg); confidence band shown. |
| Plot (spectrum 1-panel) | 4 | Title conveys "noise spectrum"; x-axis frequency, y-axis power (dB). |
| Plot (window comparison, magnitude) | 5 | Title conveys "effect of window size on resonance resolution"; three curves labeled with `M` values; a dotted/dashed vertical reference line at the true `ω_n = 10 rad/s`. |
| Print | 6 | Two lines: "without detrend: max |G| at low freq = …" and "with detrend: max |G| at low freq = …" (semantics; exact formatting is advisory). |
| Print | 7 | Two lines: "whiteness test: PASS/FAIL", "independence test: PASS/FAIL". |
| Print | 8 | One line: "NRMSE fit: …%". |
| Plot (spectrum) | 9 | Title conveys "SDOF output spectrum (time-series mode)". |

---

### 3.2 `etfe` — Empirical transfer function estimate

**Goal.** Estimate the frequency response of Plant A via `freq_etfe`
and explore the smoothing vs resolution trade-off. Contrast with the
BT estimator on the same plant (§3.1).

**Plant.** Plant A (§1.1) with `N = 2048` samples. Same `(m, k, c)` as
§3.1. This is intentional: the pair of examples forms a BT-vs-ETFE
comparison on one physical system.

**Excitation.** Same as §3.1.

**Required sections.**

| # | Section title (SHOULD) | Topic (MUST) |
|---|---|---|
| 1 | Generate test data | Simulate Plant A; compute and stash the exact discrete transfer function `G(e^{jω}) = C (e^{jω} I − Ad)⁻¹ Bd` for later overlays. |
| 2 | Basic ETFE (no smoothing) | Call `freq_etfe` with default options; plot the result. Note that `response_std` is `NaN` for ETFE. |
| 3 | Effect of smoothing | Call `freq_etfe` with three smoothing levels; overlay against the true TF reference. |
| 4 | Known FIR system: pure delay | On a synthetic `y[k] = u[k-1]` pair, call `freq_etfe` and verify |G| ≈ 1 and phase ≈ −ω. |
| 5 | Time-series mode: periodogram | With `u = None` on a fresh simulation of Plant A, show the output periodogram. |
| 6 | Custom frequency grid and Hz display | Call `freq_etfe` with a log-spaced grid and `bode_plot` with `frequency_unit='Hz'`. |

**Required `sid.*` invocations.**

| # | Function | Options |
|---|---|---|
| 2 | `freq_etfe` | default (no smoothing), `sample_time = Ts` |
| 2 | `bode_plot` | on the #2 result |
| 3 | `freq_etfe` | three calls with `smoothing ∈ {1, 11, 21}`, `sample_time = Ts` |
| 4 | `freq_etfe` | on the one-sample-delay pair (sanity check; no plant simulation) |
| 5 | `freq_etfe` | with only output data (`u = None`) |
| 5 | `spectrum_plot` | on the #5 result |
| 6 | `freq_etfe` | `smoothing = 11`, custom log-spaced frequency grid covering `[0.005, π]` rad/sample |
| 6 | `bode_plot` | with `frequency_unit = 'Hz'` |

**Required outputs.**

| Kind | Section | Binding content |
|---|---|---|
| Plot (Bode) | 2 | Title conveys "ETFE raw (no smoothing)". |
| Plot (magnitude overlay) | 3 | Three ETFE curves at different smoothing levels plus the true-TF dashed reference; x-axis rad/sample, y-axis dB. |
| Plot (2-panel magnitude + phase) | 4 | Top: `|G| ≈ 1` for the pure delay; bottom: phase vs `−ω` reference. |
| Plot (spectrum) | 5 | Title conveys "SDOF output periodogram". |
| Plot (Bode) | 6 | Title conveys "ETFE with log frequency grid (Hz)"; x-axis labeled Hz. |

---

### 3.3 `freq_dep_res` — Frequency-dependent resolution

**Goal.** Demonstrate `freq_btfdr` on a multi-mode plant where a short
fixed BT window smears the mode structure but a per-frequency window
resolves it cleanly. Show both the scalar-R and vector-R cases.

**Plant.** Plant C (§1.3) with `N = 6000` samples.

**Excitation.** Unit-variance white force at mass 1. Measurement noise
`ε ~ N(0, 5·10⁻⁴)` on `x₁`.

**Required sections.**

| # | Section title (SHOULD) | Topic (MUST) |
|---|---|---|
| 1 | Generate test data | Simulate Plant C; precompute the true TF at the default BT grid for use as a reference. |
| 2 | Fixed-window `freq_bt`: the resolution-variance trade-off | Overlay `freq_bt` at `M = 15` and `M = 80` against the true TF; show the short window smears modes and the long window is noisier. |
| 3 | Scalar resolution with `freq_btfdr` | Call `freq_btfdr` with a scalar `resolution = 0.2` and plot the Bode. |
| 4 | Per-frequency resolution vector | Build an `R_vec` ramping from fine (0.1) at low frequencies to coarse (1.5) at high frequencies; call `freq_btfdr(..., resolution=R_vec)`; plot `window_size` vs frequency (top) and the magnitude estimate vs true TF (bottom). |
| 5 | Compare BT vs BTFDR side by side | Overlay `freq_bt(M=30)` and `freq_btfdr(R=0.3)` against the true TF. |

**Required `sid.*` invocations.**

| # | Function | Options |
|---|---|---|
| 2 | `freq_bt` | two calls: `window_size=15`, `window_size=80`; `sample_time = Ts` |
| 3 | `freq_btfdr` | `resolution = 0.2`, `sample_time = Ts` |
| 3 | `bode_plot` | on the #3 result |
| 4 | `freq_btfdr` | scalar call first (to get `nf`), then `resolution = linspace(0.1, 1.5, nf)` |
| 5 | `freq_bt` | `window_size = 30` |
| 5 | `freq_btfdr` | `resolution = 0.3` |

**Required outputs.**

| Kind | Section | Binding content |
|---|---|---|
| Plot (magnitude overlay) | 2 | Two BT curves + true-TF dashed reference; title conveys "fixed window: resolution vs variance trade-off". |
| Plot (Bode) | 3 | Title conveys "freq_btfdr with scalar resolution R = 0.2". |
| Plot (2-panel: window size + magnitude) | 4 | Top panel: `window_size` vs frequency (line). Bottom panel: BTFDR magnitude overlaid on true TF. |
| Plot (magnitude overlay) | 5 | BT, BTFDR, true TF on one axes; title conveys "fixed vs frequency-dependent resolution". |

---

### 3.4 `coherence` — Coherence analysis

**Goal.** Use squared coherence `γ²(ω)` to diagnose frequency-local
estimate quality. Show how a colored disturbance entering a plant via
an unmeasured channel depresses coherence in its own spectral band.

**Plant.** Plant B (§1.2) with two input channels:

- Channel 0: commanded force at mass 1, `F_col_0 = [1, 0]`.
- Channel 1: unmeasured disturbance at mass 2, `F_col_1 = [0, 1]`.

`F = [[1, 0], [0, 1]] = I₂`, `N = 4000` samples.

**Excitation.**

- Commanded force `u[k] = 10 · ξ[k]` with `ξ[k] ~ N(0, 1)`.
- Disturbance `d[k] = 0.5 · AR1(e)[k]` where `e[k] ~ N(0, 1)` is white
  noise and `AR1` denotes `y[k] = 0.9 y[k-1] + e[k]` (AR(1) pole at
  0.9, a DC-heavy colored process).

Measured output `y[k] = x₁[k]` (no sensor noise — the coherence drop
must come from the disturbance, not from additive noise).

**Required sections.**

| # | Section title (SHOULD) | Topic (MUST) |
|---|---|---|
| 1 | Generate test data | Build `(Ad, Bd)` with `F = I₂` via §2.1 so the returned `Bd` has both channels. Simulate `x[k+1] = Ad x[k] + Bd[:,0] u[k] + Bd[:,1] d[k]`. |
| 2 | Estimate with `freq_bt` | Call `freq_bt` with the single "known" input channel `u` and the position output `y`; the disturbance is unmodeled. |
| 3 | Plot Bode magnitude and coherence together | 2-panel plot: magnitude on top, squared coherence on bottom with `γ² = 0.5` and `γ² = 0.9` reference lines. |
| 4 | Confidence bands reflect coherence | Call `bode_plot` with `confidence = 2` and `confidence = 3` in side-by-side subplots. |
| 5 | High-disturbance vs low-disturbance comparison | Re-run `freq_bt` with `d` scaled by 0.1× and by 2.0×; overlay the two coherence curves. |
| 6 | Note: ETFE does not provide coherence | Call `freq_etfe` and print a one-line confirmation that `result.coherence is None`. |

**Required `sid.*` invocations.**

| # | Function | Options |
|---|---|---|
| 2 | `freq_bt` | `window_size = 200`, custom grid `linspace(0.01, π, 512)`, `sample_time = Ts` |
| 3 | (plot only) | — |
| 4 | `bode_plot` | two calls: `confidence = 2`, `confidence = 3` |
| 5 | `freq_bt` | two calls, one per disturbance scale |
| 6 | `freq_etfe` | one call; just to inspect `coherence is None` |

**Required outputs.**

| Kind | Section | Binding content |
|---|---|---|
| Plot (2-panel) | 3 | Top: Bode magnitude; bottom: coherence `γ²` with 0.5 and 0.9 horizontal references. |
| Plot (2×2 subplots) | 4 | Four panels forming `confidence=2` Bode and `confidence=3` Bode side by side. |
| Plot (coherence overlay) | 5 | Two curves labeled "low disturbance" and "high disturbance". |
| Print | 6 | One line: "ETFE coherence is None: True". |

---

### 3.5 `method_comparison` — Comparing frequency-domain estimators

**Goal.** Apply all three frequency-domain estimators (`freq_bt`,
`freq_btfdr`, `freq_etfe`) to one plant, compare magnitudes, noise
spectra, and NRMSE fits; end with a summary table of trade-offs.

**Plant.** Plant A (§1.1) with `N = 2048` samples. Same plant as §3.1
and §3.2 — enabling direct BT-vs-ETFE-vs-BTFDR comparison.

**Excitation.** Same as §3.1.

**Required sections.**

| # | Section title (SHOULD) | Topic (MUST) |
|---|---|---|
| 1 | Generate test data | Simulate Plant A. |
| 2 | Estimate with all three methods | Call each of `freq_bt`, `freq_btfdr`, `freq_etfe` on the same custom frequency grid with pedagogically-reasonable options. |
| 3 | Compare Bode magnitude plots | Four curves (`freq_etfe` raw, `freq_etfe` smoothed, `freq_bt`, `freq_btfdr`) plus the true TF overlaid. |
| 4 | Compare noise spectra | Overlay the noise-spectrum estimates from the three methods. |
| 5 | Custom logarithmic frequency grid | Re-run the three methods on a log-spaced grid. |
| 6 | Time-series comparison: periodogram vs smoothed spectrum | Fresh simulation; hand `u = None` to both `freq_bt` and `freq_etfe`. |
| 7 | Model output comparison using `sid.compare` | Call `compare` on each of the non-time-series results and print the NRMSE fits. |
| 8 | Summary of method trade-offs | A markdown table summarizing window-size, uncertainty availability, coherence availability, and best-for notes. |

**Required `sid.*` invocations.**

| # | Function | Options |
|---|---|---|
| 2 | `freq_bt` | `window_size = 200`, custom grid `linspace(0.005, π, 512)`, `sample_time = Ts` |
| 2 | `freq_etfe` | two calls: unsmoothed and `smoothing = 15`, same grid |
| 2 | `freq_btfdr` | `resolution = 0.3`, same grid |
| 3 | (plot only) | — |
| 4 | (plot only) | — |
| 5 | `freq_bt`, `freq_etfe`, `freq_btfdr` | on a log-spaced grid `logspace(log10(0.005), log10(π), 200)` |
| 6 | `freq_bt`, `freq_etfe` | each with `u = None` |
| 7 | `compare` | one call per non-time-series result |

**Required outputs.**

| Kind | Section | Binding content |
|---|---|---|
| Plot (magnitude overlay) | 3 | Four estimator curves + true TF; legend identifies each. |
| Plot (noise-spectrum overlay) | 4 | Three estimator curves. Ports MUST clamp the noise-spectrum array with a small positive floor (e.g. `max(abs(Φ_v), eps)`) before applying `log10` because per §2.7 of `SPEC.md` the noise spectrum is non-negative and hits exact zero at frequencies where the clamp activates. |
| Plot (magnitude overlay) | 5 | Log-grid version of §3 with three estimators + true. |
| Plot (power-spectrum overlay) | 6 | ETFE periodogram (grey) and BT smoothed spectrum (blue). |
| Print | 7 | Four NRMSE lines, one per method (including the second ETFE smoothing variant). |
| Table | 8 | Method trade-off table (pure markdown). |

---

### 3.6 `mimo` — Multi-input multi-output frequency response

**Goal.** Demonstrate MIMO identification with `freq_bt` on a plant
that has natural 2×2 force-to-position structure. Inspect the
`response.shape == (nf, ny, nu)` dimension and the noise spectral
matrix's off-diagonal entries.

**Plant.** Plant B (§1.2). Two sub-sections use two different
input-distribution matrices:

- Sub-section A (2-output, 1-input): `F = [[1], [0]]`, `N = 4000`.
- Sub-section B (2-output, 2-input): `F = I₂`, `N = 4000`.

**Excitation.**

- Sub-section A: white unit-variance force `u` on channel 0.
- Sub-section B: independent white unit-variance forces on both input
  channels.

Measurement noise `ε ~ N(0, 2·10⁻⁴)` on both position outputs.

**Required sections.**

| # | Section title (SHOULD) | Topic (MUST) |
|---|---|---|
| 1 | 2-output, 1-input system | Build the plant with `F = [[1], [0]]`; simulate; call `freq_bt` with a custom grid. Precompute the exact discrete TFs `G₁(z)` (input → x₁) and `G₂(z)` (input → x₂) for reference overlays. |
| 2 | Inspect MIMO result dimensions | Print `result.response.shape`, `result.noise_spectrum.shape`, and the fact that `result.coherence is None` for MIMO. |
| 3 | Plot individual output channels | Two-panel Bode magnitude: `G₁` and `G₂` with BT estimate and dashed true reference on each panel. |
| 4 | Noise spectral matrix | Plot the diagonals `Φ_v[0,0](ω)` and `Φ_v[1,1](ω)` from `result.noise_spectrum`. |
| 5 | 2-output, 2-input system | Rebuild the plant with `F = I₂`; simulate with `(N, 2)` inputs; call `freq_bt`. Print the 2×2 response shape. |
| 6 | Plot the full 2×2 transfer matrix | 2×2 subplot grid of the four magnitudes `G_{ij}(ω)`. |
| 7 | MIMO uncertainty | Print a one-line confirmation that `result.response_std` is all `NaN` for MIMO. |

**Required `sid.*` invocations.**

| # | Function | Options |
|---|---|---|
| 1 | `freq_bt` | `window_size = 200`, custom grid `linspace(0.005, π, 512)`, `sample_time = Ts` |
| 5 | `freq_bt` | on `(N, 2)` input and `(N, 2)` output arrays, same options |

**Required outputs.**

| Kind | Section | Binding content |
|---|---|---|
| Print | 2 | Three lines: response shape, noise-spectrum shape, coherence is None. |
| Plot (2-panel) | 3 | Two Bode magnitudes with estimated vs true reference. |
| Plot (1-panel) | 4 | Noise-spectral-matrix diagonals. |
| Print | 5 | 2×2 MIMO response shape. |
| Plot (2×2 grid) | 6 | Four magnitude plots for `G₁₁, G₁₂, G₂₁, G₂₂`. |
| Print | 7 | One line: "MIMO response_std contains NaN: True". |

---

### 3.7 `freq_map` — Time-varying frequency response maps

**Goal.** Demonstrate `freq_map` on four scenarios: an LTI baseline
(should look stationary), two LTV variants (ramp and step change
built via §2.2), and a nonlinear hardening Duffing oscillator whose
apparent resonance drifts purely from nonlinearity (via §2.3).

**Plants.** Plant B (§1.2) for sections 1–3, 5–7. Plant E (§1.5) for
section 4 (Duffing). Section lengths: `N = 4000` for Plant B LTI,
`N = 4000` for the LTV variants, `N = 4000` for Duffing.

**Excitation.**

- Plant B sections: white unit-variance force at mass 1.
- Duffing section: `u[k] = amp[k] · ξ[k]` where `amp` ramps linearly
  from `0.5` to `10.0` over the record and `ξ[k] ~ N(0, 1)` — a
  ramped-amplitude white force.

Measurement noise `5·10⁻⁴` on the Plant B sections; none on Duffing
(the pedagogical point is purely in the nonlinearity).

**Required sections.**

| # | Section title (SHOULD) | Topic (MUST) |
|---|---|---|
| 1 | LTI baseline: constant 2-mass chain | Build Plant B via §2.1, simulate, call `freq_map`. The magnitude map should look constant along time. |
| 2 | Continuous LTV: ramping first-spring stiffness | Build `k_spring` as a `(2, N)` array with `k[0, :]` a linspace `200 → 20` and `k[1, :] = 80`. Use §2.2 to get the per-step stack; simulate the LTV recursion. Call `freq_map`. |
| 3 | Discrete LTV: step change in stiffness | Same plant but `k[0, :N/2] = 200` and `k[0, N/2:] = 40`; otherwise identical to section 2. |
| 4 | Coherence map | Re-use the section-2 result; plot the coherence map. |
| 5 | BT vs Welch algorithm | Call `freq_map` on the section-2 data with `algorithm='bt'` and `algorithm='welch'` side by side. |
| 6 | Segment length and overlap tuning | Call `freq_map` on the section-2 data with two different `segment_length` choices. |
| 7 | Time-series mode: evolving output spectrum | Fresh simulation of the section-2 LTV plant; hand only the output to `freq_map`; plot the spectrum map. |
| 8 | Duffing hardening oscillator | Use §2.3 with `k_cub = 1e5` and a ramped-amplitude input; call `freq_map` on the response; show the apparent resonance drifting upward over time. Print the small-amplitude and ramped-amplitude effective natural frequencies for context. |

**Required `sid.*` invocations.**

| # | Function | Options |
|---|---|---|
| 1 | `freq_map` | `segment_length = 512`, `overlap = 384`, `sample_time = Ts` |
| 2–3 | `freq_map` | `segment_length = 256`, `overlap = 192`, `sample_time = Ts` |
| 4 | `freq_map` | (reuse result) — just plot coherence |
| 5 | `freq_map` | two calls: `algorithm = 'bt'`, `algorithm = 'welch'`, `segment_length = 256` |
| 6 | `freq_map` | two calls: `segment_length = 128` and `segment_length = 512` |
| 7 | `freq_map` | with `u = None` |
| 8 | `freq_map` | `segment_length = 256`, `overlap = 192`, `sample_time = Ts` |

**Required outputs.**

| Kind | Section | Binding content |
|---|---|---|
| Plot (magnitude map) | 1 | Title conveys "LTI baseline: stationary map". |
| Plot (magnitude map) | 2 | Title conveys "continuous LTV: ramping `k₁`, first mode drifts 7→3 rad/s". |
| Plot (magnitude map) | 3 | Title conveys "discrete LTV: step change at `t = T/2`". |
| Plot (coherence map) | 4 | Title conveys "coherence map". |
| Plot (two magnitude maps) | 5 | Side-by-side BT vs Welch. |
| Plot (two magnitude maps) | 6 | Side-by-side short vs long segments. |
| Plot (spectrum map) | 7 | Title conveys "time-series output spectrum drifts as `k₁` softens". |
| Print | 8 | Two lines reporting early vs late position RMS and their linearized effective natural frequencies. |
| Plot (magnitude map) | 8 | Title conveys "Duffing hardening: apparent resonance rises with amplitude". |

---

### 3.8 `spectrogram` — Short-time FFT spectrogram

**Goal.** Demonstrate `spectrogram` on a physical SDOF plant driven
by a chirp force. The output response spectrogram shows the chirp
track *modulated* by the plant's own resonance — the ridge lights up
when the chirp sweeps through `ω_n`.

**Plant.** Plant D (§1.4) with `N = 5000`, `Ts = 1/Fs`, `Fs = 1000 Hz`.

**Excitation.** Linear chirp force `u[k] = cos(2π · φ(t_k))` with
instantaneous frequency `f(t) = f₀ + (f₁ − f₀) · t / (2 T_end)` and
`f₀ = 20 Hz`, `f₁ = 60 Hz`. This sweeps through `Plant D`'s resonance
at ≈ 31.83 Hz near `t ≈ 1.5 s`.

**Required sections.**

| # | Section title (SHOULD) | Topic (MUST) |
|---|---|---|
| 1 | Chirp force driving the SDOF | Build Plant D via §2.1; simulate the state recursion under the chirp input; compute the spectrogram of the position response. |
| 2 | Window length trade-off | Two `spectrogram` calls with short and long window lengths; side-by-side plot. |
| 3 | Window types | Three `spectrogram` calls with `window ∈ {'hann', 'hamming', 'rect'}` at the same length. |
| 4 | Multi-channel signal | Stack the chirp-response (channel 0) with a constant-tone response (channel 1) driven by `u_tone[k] = cos(2π · 50 · t_k)`. Compute a multi-channel spectrogram and plot both channels. |
| 5 | Log frequency scale and NFFT zero-padding | Single `spectrogram` call with `window_length = 128` and `nfft = 1024`; plot with `frequency_scale = 'log'`. |
| 6 | Accessing raw STFT data | Print the dimensions of `result.time`, `result.frequency`, `result.power`, and `result.complex_stft`. |

**Required `sid.*` invocations.**

| # | Function | Options |
|---|---|---|
| 1 | `spectrogram` | `window_length = 256`, `sample_time = Ts` |
| 2 | `spectrogram` | two calls: `window_length = 64`, `window_length = 512` |
| 3 | `spectrogram` | three calls: `window = 'hann'`, `'hamming'`, `'rect'`; `window_length = 256` |
| 4 | `spectrogram` | on a two-channel signal `(N, 2)` stacked from the chirp response and the 50 Hz tone response |
| 5 | `spectrogram` | `window_length = 128`, `nfft = 1024`, `sample_time = Ts`; plotted with `frequency_scale = 'log'` |

**Required outputs.**

| Kind | Section | Binding content |
|---|---|---|
| Plot (spectrogram) | 1 | Title conveys "SDOF response to chirp force: resonance lights up near ~32 Hz"; x-axis time (s), y-axis frequency (Hz). |
| Plot (two spectrograms) | 2 | Side-by-side short and long window. |
| Plot (three spectrograms) | 3 | Hann / Hamming / Rectangular. |
| Plot (two spectrograms) | 4 | Two channels: chirp response and 50 Hz tone response (the tone response should show a horizontal line at 50 Hz). |
| Plot (spectrogram) | 5 | Log-frequency axis; title conveys "log frequency scale with zero-padding (NFFT = 1024)". |
| Print | 6 | Four lines reporting `len(time)`, `len(frequency)`, `power.shape`, `complex_stft.shape`. |

---

### 3.9 `ltv_disc` — LTV state-space identification (COSMIC)

**Goal.** Walk through the full COSMIC workflow — LTI recovery,
LTV recovery with ramping stiffness, multi-trajectory benefit,
validation lambda tuning, cost decomposition, uncertainty
quantification, frozen transfer function, compare/residual
validation — all on a physical 1-DoF SMD plant. End with a
nonlinear Duffing section that recovers the amplitude-dependent
local linearization.

**Plant.** Plant A (§1.1) for all linear sections. Plant E (§1.5) for
the Duffing section (same `m`, `k_lin`, `c`, `F` as Plant A; adds
`k_cub`).

**Dimensions.** `p = 2` (state = `[x, v]`), `q = 1` (single force
input).

**Removed in v1.0.1.** Earlier drafts included a "Preconditioning
for numerical stability" section and a "Frequency-based lambda
tuning" section. The preconditioning code path is intentionally
disabled in v1.0 of the library (the call site emits a runtime
warning), so demonstrating it produces a noisy example. The
frequency-based tuning method compares the COSMIC frozen TF against
a non-parametric `freq_map` estimate; on physical SMD plants the
velocity state has a band-pass spectrum that `freq_map` cannot
reliably estimate at the low frequencies where the modes live, so
no lambda in any reasonable grid clears the 90% consistency
threshold and the call emits "No lambda achieved 90% consistency".
Both sections were removed to keep this example a clean demo. The
library still supports both features; they are simply not
exercised here.

**Required sections.**

| # | Section title (SHOULD) | Topic (MUST) |
|---|---|---|
| 1 | LTI system recovery | Build `(Ad, Bd)` via §2.1 for Plant A; simulate `L = 10` trajectories of `N = 50` steps with process noise `σ = 0.01`; call `ltv_disc` with `lambda_ = 1e5`; print the true `Ad`, the mean recovered `A(k)` across time, and the Frobenius-norm recovery error. |
| 2 | LTV system: time-varying stiffness | Build an `(n=1, N=250)` `k_spring` array ramping `200 → 50`; use §2.2 for the per-step stack; simulate `L = 30` trajectories with process noise `σ = 0.01`; call `ltv_disc` with `lambda_ = 'auto'`. Plot the recovered `A[1, 0](k)` against the true curve (for a 1-DoF plant with small `Ts`, `Ad[1, 0] ≈ −k(t)·Ts`). The longer record and larger trajectory count are chosen so that the downstream residual whiteness test (section 7) passes cleanly. |
| 3 | Multi-trajectory benefit | Compare COSMIC recovery with `L = 5` vs `L = 30` trajectories on the section-2 LTV plant. |
| 4 | Validation-based lambda tuning | Split `L = 30` trajectories into 20 train + 10 validation; call `ltv_disc_tune` with `method='validation'` and a log-spaced lambda grid; plot the validation-loss curve. |
| 5 | Cost decomposition | Print the three-element `cost = [total, data_fidelity, regularization]`. |
| 6 | Uncertainty quantification | Call `ltv_disc` with `uncertainty = True` and plot `A[1, 0](k) ± 2σ` over the true curve. |
| 7 | Frozen transfer function with `ltv_disc_frozen` | Call `ltv_disc_frozen` on the section-6 result at three time steps `[0, N/2, N-1]`; overlay the three Bode magnitudes. |
| 8 | Model validation with `compare` and `residual` | Call both on the section-2 result; print per-channel NRMSE fits and the whiteness verdict. With `N = 250, L = 30` the whiteness test is expected to PASS. |
| 9 | Weakly-nonlinear Duffing oscillator | Build Plant E; simulate `L = 12` trajectories over `N = 400` steps with a *ramped-amplitude* white input (amplitude profile `linspace(0.5, 8.0, N)`); call `ltv_disc` with a **small manual lambda** (`lambda_ = 0.1`, not `'auto'`, because auto-tuning over-regularises this dataset); plot the recovered `A[1, 0](k)` and overlay the linearized `Ad[1, 0]` reference from §2.1. |

**Required `sid.*` invocations.**

| # | Function | Options |
|---|---|---|
| 1 | `ltv_disc` | `lambda_ = 1e5` |
| 2 | `ltv_disc` | `lambda_ = 'auto'` |
| 3 | `ltv_disc` | two calls on subsets with different `L` |
| 4 | `ltv_disc_tune` | positional args `(X_train, U_train, X_val, U_val)`, `method = 'validation'`, `lambda_grid = logspace(-3, 6, 30)` |
| 6 | `ltv_disc` | `lambda_ = 'auto'`, `uncertainty = True` |
| 7 | `ltv_disc_frozen` | `time_steps = [0, N/2, N-1]` on the uncertainty-enabled result |
| 8 | `compare`, `residual` | each on the section-2 result |
| 9 | `ltv_disc` | `lambda_ = 0.1` (manual) on the Duffing data |

**Required outputs.**

| Kind | Section | Binding content |
|---|---|---|
| Print | 1 | True `Ad`, mean recovered `A`, Frobenius recovery error. |
| Plot (`A[1,0](k)`) | 2 | True curve + COSMIC recovery overlay. |
| Plot (`A[1,0](k)`) | 3 | `L_few`, `L_many`, and true curves. |
| Plot (validation loss curve) | 4 | Loss vs λ on log scale with the minimum marked. |
| Print | 5 | Four-line cost breakdown with a consistency check (total − data − reg ≈ 0). |
| Plot (`A[1,0](k)` with ±2σ band) | 6 | — |
| Plot (Bode magnitudes at 3 time steps) | 7 | — |
| Print | 8 | Per-channel fits and whiteness PASS (expected PASS with the v1.0.1 parameter choice). |
| Plot (2-panel: amplitude profile + recovered `A[1,0](k)`) | 9 | Top panel: excitation amplitude ramp. Bottom panel: recovered `A[1,0](k)` with the small-amplitude linear reference as a horizontal dashed line. |

---

### 3.10 `multi_trajectory` — Multi-trajectory ensemble averaging

**Goal.** Show that ensemble averaging across `L` independent
trajectories reduces variance by `1/L` without sacrificing frequency
resolution. Four sub-sections exercise this benefit across all four
spectral estimators: `freq_bt`, `freq_map`, `spectrogram`, `ltv_disc`.

**Plants.**

- Sections 1, 2, 4: Plant B (§1.2) with `F = [[1], [0]]`.
- Section 3: Plant D (§1.4) — chirp driving the high-frequency SDOF.

**Required sections.**

| # | Section title (SHOULD) | Topic (MUST) |
|---|---|---|
| 1 | LTI ensemble averaging — tighter confidence bands | Simulate `L = 10` trajectories of Plant B (`N = 2000`, measurement noise `5·10⁻⁴`); compare `freq_bt` confidence bands between a single trajectory and the full ensemble. Print the `response_std` ratio and confirm it is approximately `1/sqrt(L)`. |
| 2 | LTV time-varying map — sharper transition detection | Simulate `L = 5` trajectories of Plant B with a step change `k₁: 200 → 50` at `t = N/2`, `N = 4000`; compare `freq_map` between a single trajectory and the ensemble. |
| 3 | Spectrogram averaging — chirp in noise | Simulate `L = 8` trajectories of Plant D driven by the same chirp but with independent `1·10⁻⁴` noise; compare single-trajectory and ensemble spectrograms. |
| 4 | COSMIC + `freq_map` consistency | Reuse the ltv_disc LTV plant (Plant A 1-DoF with ramping stiffness, `N = 80`, `L = 10`); call `ltv_disc(lambda_='auto', uncertainty=True)` and `freq_map` on the same dataset; print that both identifiers use the same `L` trajectories. |

**Required `sid.*` invocations.**

| # | Function | Options |
|---|---|---|
| 1 | `freq_bt` | two calls: single-trajectory and ensemble, `window_size = 80`, `sample_time = Ts` |
| 2 | `freq_map` | two calls: single-trajectory and ensemble, `segment_length = 256`, `sample_time = Ts` |
| 3 | `spectrogram` | two calls: single-trajectory and ensemble, `window_length = 128`, `sample_time = Ts` |
| 4 | `ltv_disc` | `lambda_ = 'auto'`, `uncertainty = True` |
| 4 | `freq_map` | `segment_length = min(N, 30)`, `sample_time = Ts` |

**Required outputs.**

| Kind | Section | Binding content |
|---|---|---|
| Print | 1 | Three lines: single-trajectory max std, ensemble max std, ratio and the expected `1/sqrt(L)`. |
| Plot (two Bode) | 1 | Single-trajectory Bode and ensemble Bode side by side (or stacked). |
| Plot (two magnitude maps) | 2 | Single-trajectory and ensemble `freq_map` magnitude maps side by side. |
| Plot (two spectrograms) | 3 | Single-trajectory and ensemble spectrograms side by side. |
| Print | 4 | Two lines: COSMIC `A(0)` and COSMIC `A(N-1)`; one line confirming `freq_map.num_trajectories` equals `L`. |

---

### 3.11 `output_cosmic` — LTV identification from partial observations

**Goal.** Demonstrate `ltv_disc_io` (Output-COSMIC) on a natural
partial-observation scenario: a 2-mass mechanical plant with position
sensors on both masses but *no* velocity sensors. The hidden-state
dimension is `n = 4`, the measured dimension is `py = 2`, and COSMIC
must infer the hidden velocities along with the dynamics.

**Plant.** Plant B (§1.2) with `F = [[1], [0]]` (force at mass 1
only). Observation matrix `H = [[1, 0, 0, 0], [0, 1, 0, 0]]` (measure
both positions, velocities hidden).

**Trajectories.** `N = 80`, `L = 10`, `u[k] = 5 · ξ[k]` with `ξ ~
N(0, 1)` (scaled up so displacements reach a few centimetres —
displacement magnitudes that let the identification converge). Process
noise `σ_proc = 1·10⁻³`, measurement noise `σ_meas = 1·10⁻⁴`.

**Important framing note (binding).** Output-COSMIC recovers
`(A, B, x)` up to an *unobservable similarity transform* — element-
wise comparison of the recovered `A`/`B` against the simulation's
`Ad`/`Bd` is **not** a meaningful metric. The validation step MUST be
the gauge-invariant **observation reconstruction error**
`‖H · x̂ − y‖ / ‖y‖`.

**Required sections.**

| # | Section title (SHOULD) | Topic (MUST) |
|---|---|---|
| 1 | System setup | Build Plant B; state it is a 4th-order plant with 2 position measurements. Print `Ad` and `H`. |
| 2 | Simulate trajectories | Generate `L = 10` trajectories of length `N = 80` with scaled input and process/measurement noise. |
| 3 | Step 1: estimate frequency response | Call `freq_bt` on the first trajectory (trimming `Y` to match `U`'s length) with `window_size = 20`, `sample_time = Ts`. |
| 4 | Step 2: model-order determination | Call `model_order` on the `freq_bt` result; print the estimate and note that it may overshoot the true order for short lightly-damped records. |
| 5 | Step 3: construct observation matrix | Show `H` and print its shape. |
| 6 | Step 4: identify the LTV model via `ltv_disc_io` | Call `ltv_disc_io(Y, U, H, lambda_=1e5)`; print iterations and final cost. |
| 7 | Convergence history | Plot `cost` vs iteration (semilogy). |
| 8 | State recovery: observed channels vs hidden channels | 2×2 grid: top row shows the two measured positions with true, estimated, and measurement curves; bottom row shows the two hidden velocities with true and estimated curves. |
| 9 | Validation: observation reconstruction error | Compute `H · x̂` for every trajectory and compare against `Y`; print the relative Frobenius error. |
| 10 | Frozen-time inspection of the recovered A and B | Print `A(N/2)` and `B(N/2)` as a sanity check that magnitudes are `O(1)`. |

**Required `sid.*` invocations.**

| # | Function | Options |
|---|---|---|
| 3 | `freq_bt` | `window_size = 20`, `sample_time = Ts` |
| 4 | `model_order` | default |
| 6 | `ltv_disc_io` | `lambda_ = 1e5` |

**Required outputs.**

| Kind | Section | Binding content |
|---|---|---|
| Print | 1 | `Ad` matrix and `H` matrix. |
| Print | 2 | Max absolute output and max absolute hidden-velocity value. |
| Print | 3 | `freq_bt` response shape. |
| Print | 4 | `n_est` estimate with a note that it may overshoot. |
| Print | 5 | `H` matrix. |
| Print | 6 | Iterations and final cost. |
| Plot (convergence) | 7 | Semilog `cost` vs iteration. |
| Plot (2×2 grid) | 8 | Four subplots: `x₁`, `x₂`, `v₁`, `v₂` with true/estimated/measured overlays as appropriate. |
| Print | 9 | One line: "Observation reconstruction error: …". |
| Print | 10 | `A(mid)` and `B(mid)` matrices. |

---

## 4. Conventions

### 4.1 Auto-discovery

Every language port MUST lay out its examples directory so that the
example runner discovers files purely by glob pattern. Hardcoded
manifest lists are non-conformant.

| Language | Directory | Glob pattern |
|---|---|---|
| Python | `python/examples/` | `example_*.ipynb` |
| MATLAB / Octave | `matlab/examples/` | `example*.m` |
| Julia | `julia/examples/` | `example_*.jl` |

The CI test runner (e.g., `pytest --nbmake python/examples/` for
Python, `matlab/examples/runAllExamples.m` for MATLAB) MUST discover
new examples without code changes. This rule mirrors the discovery
convention in `CONTRIBUTING.md` at the repo root.

### 4.2 Self-contained

Every example MUST run top-to-bottom with no external data files. All
simulation data is generated inline via the helper API of §2. Loading
from `testdata/*.json` or any other on-disk fixture is non-conformant.

### 4.3 Outputs cleared

For languages where example files persist runtime outputs (Jupyter
notebooks, R Markdown), every example file MUST be committed with all
outputs cleared. The CI runner validates that the file *executes*
without error; it does not assert output content. MATLAB `.m` scripts,
Julia `.jl` scripts, and similar source-only formats have nothing to
clear.

### 4.4 RNG seeding

Every example MUST fix an RNG seed at the top of its data-generation
section so that a given commit gives deterministic numerical output
within one language.

Seeds are *not* bound across languages. Python's `numpy.random` and
MATLAB's `randn` produce different sequences from the same seed, so
numerical equivalence across languages is not a goal of this spec.
Each port picks seeds that produce visually similar plots; reviewers
compare Python and MATLAB output side by side to catch mismatches,
not by diffing numeric values.

**Recommended convention (SHOULD):** use the same seed value across
languages when possible (e.g., `seed = 42` in both Python's
`default_rng(42)` and MATLAB's `rng(42)`). This makes differences
obvious as RNG differences rather than spec mismatches.

### 4.5 File-name mapping table

Each example's canonical identifier maps to language-specific file
names as follows:

| Identifier | Python (`python/examples/`) | MATLAB (`matlab/examples/`) | Julia (`julia/examples/`) |
|---|---|---|---|
| `siso`              | `example_siso.ipynb`              | `exampleSISO.m`             | `example_siso.jl`             |
| `etfe`              | `example_etfe.ipynb`              | `exampleETFE.m`             | `example_etfe.jl`             |
| `freq_dep_res`      | `example_freq_dep_res.ipynb`      | `exampleFreqDepRes.m`       | `example_freq_dep_res.jl`     |
| `coherence`         | `example_coherence.ipynb`         | `exampleCoherence.m`        | `example_coherence.jl`        |
| `method_comparison` | `example_method_comparison.ipynb` | `exampleMethodComparison.m` | `example_method_comparison.jl`|
| `mimo`              | `example_mimo.ipynb`              | `exampleMIMO.m`             | `example_mimo.jl`             |
| `freq_map`          | `example_freq_map.ipynb`          | `exampleFreqMap.m`          | `example_freq_map.jl`         |
| `spectrogram`       | `example_spectrogram.ipynb`       | `exampleSpectrogram.m`      | `example_spectrogram.jl`      |
| `ltv_disc`          | `example_ltv_disc.ipynb`          | `exampleLTVdisc.m`          | `example_ltv_disc.jl`         |
| `multi_trajectory`  | `example_multi_trajectory.ipynb`  | `exampleMultiTrajectory.m`  | `example_multi_trajectory.jl` |
| `output_cosmic`     | `example_output_cosmic.ipynb`     | `exampleOutputCOSMIC.m`     | `example_output_cosmic.jl`    |

Python pattern: `example_` + snake_case + `.ipynb`. MATLAB pattern:
`example` + camelCase + `.m` (matching the existing MATLAB examples in
the repo). Julia pattern: same as Python but with `.jl` extension.

### 4.6 Helper module placement

The `util_msd*` helpers specified in §2 MUST live **inside the
examples directory** as sibling modules, not in the main package.
This reflects their purpose as example fixtures, not part of the
public `sid` API.

| Language | Location |
|---|---|
| Python | `python/examples/util_msd.py` |
| MATLAB | `matlab/examples/util_msd.m`, `matlab/examples/util_msd_ltv.m`, `matlab/examples/util_msd_nl.m` |
| Julia | `julia/examples/util_msd.jl` (module file) |

Notebooks and scripts MUST import the helpers by sibling-path
reference, not by reaching into the main package's internals.

### 4.7 Example file templates and boilerplate

Every language provides templates for writing new examples. Every
example file in this spec MUST follow its language's template:

- Python: notebook with a title markdown cell, an imports code cell
  including `%matplotlib inline`, and one or more paired
  markdown + code cells per §3 section.
- MATLAB: script with a top `%% exampleName` title block comment, an
  imports block (or a runner-style `runner__nCompleted` counter per
  the MATLAB CONTRIBUTING guide), and `%% Section N — Title` block
  comments separating sections. Matching `matlab/examples/example_template.m`.
- Julia: script or notebook per the local language convention.

### 4.8 Import of the local sid package

Every example MUST exercise functions from the local `sid` package
and MUST NOT call any comparable third-party implementation (e.g.,
`scipy.signal.welch`, `control.freqresp`) in places where a `sid.*`
function would be equivalent. Third-party imports are allowed for
utility tasks (plotting, array construction) but never for the
identification or spectral analysis that is the example's subject.

---

## 5. Cross-Language "Equivalence"

This section defines what it means for two language ports to implement
"the same example". It is the reviewer's contract when approving a new
port.

### 5.1 Structural equivalence (binding)

For every example in §3, two conformant ports MUST share:

1. **Section inventory.** The same ordered list of pedagogical
   sections, identified by topic (not verbatim header prose — see
   §0.3). A reviewer comparing the two implementations side by side
   should be able to say "this section on the left corresponds to
   that section on the right" for every section, without ambiguity.
2. **`sid.*` call graph.** The same sequence of `sid.*` function
   calls in the same sections, with the same binding options (§3.X.5
   tables). A port MAY invoke additional language-idiomatic helpers
   (e.g., a MATLAB port might use `fprintf` where Python uses an
   f-string), but every `sid.*` call listed in the spec MUST be
   present.
3. **Plot kinds.** The same kinds of plots (Bode, spectrum, magnitude
   overlay, time-frequency map, spectrogram, time series, 2×2 grid,
   etc.) in the same sections.
4. **Printed line semantics.** The same pieces of information are
   printed, in the same order. "Whiteness test: PASS" in Python is
   equivalent to "Whiteness test: PASS" in MATLAB even if the exact
   formatting differs.

### 5.2 Numerical divergence (not binding)

RNG draws differ across languages. Numerical outputs (fit
percentages, peak magnitudes, exact curve values, iteration counts,
final cost values, recovered matrix entries) are allowed to differ by
any amount consistent with finite-sample variance.

A language port is NOT required to reproduce the Python port's exact
numerical outputs. Bit-identity is a non-goal.

### 5.3 Visual similarity (recommended)

A reviewer comparing side-by-side plots from two ports of the same
example SHOULD be able to identify which example they are looking at
without reading any code or markdown. The resonance should land at
the same frequency, the confidence bands should have similar widths,
the time-frequency ridges should traverse the same diagonal, etc.

This is a judgment test, not a mechanical one. Mismatches that trace
entirely to RNG differences are acceptable; mismatches that trace to
wrong plant parameters or missing `sid.*` calls are not.

### 5.4 Review checklist (per example)

When reviewing a new language port of the example suite, confirm each
item for each example:

```
For example `<identifier>`:

[ ] Plant parameters match §1 exactly.
[ ] The helper functions from §2 are called with the right shapes.
[ ] The ordered list of pedagogical sections matches §3.X.4.
[ ] Every sid.* call listed in §3.X.5 is present with the specified
    options.
[ ] Every plot listed in §3.X.6 is produced with a title that conveys
    the specified topic and with the specified axis labels.
[ ] Every printed line listed in §3.X.6 is produced with the
    specified semantics.
[ ] The example runs to completion in the language's example runner
    (`pytest --nbmake`, `runAllExamples`, etc.).
[ ] Side-by-side visual comparison against the Python reference shows
    similar plot features (resonance peaks in the same places, bands
    of similar width, time-frequency ridges traversing similar
    paths). Differences consistent with RNG noise are acceptable.
[ ] No external data files are loaded (§4.2).
[ ] Outputs are cleared if the file format persists them (§4.3).
[ ] Auto-discovery picks the file up without manifest changes (§4.1).
```

A port that passes all eleven example checklists is conformant.

### 5.5 What this spec does NOT verify

The following are **not** checked by this specification, even though
they are important in their own right:

- **Cross-language numerical agreement at the algorithm level.** That
  is the job of `SPEC.md` and the reference test vectors in
  `testdata/`. The example suite exercises the algorithms as a user
  would; it is not a verification harness for numerical equivalence.
- **Plot styling (colors, fonts, figure sizes).** Ports are free to
  adopt the idiomatic visual style of their plotting library.
- **Internal variable names inside example code.** A port MAY use
  `N_samples` where the Python reference uses `N`; no renaming test
  is performed.
- **Markdown / comment prose.** As stated in §0.3, narrative prose is
  advisory. Only section *topics* are binding.

---

## 6. Versioning and Change Control

### 6.1 Semantic versioning

This specification uses semantic versioning (`MAJOR.MINOR.PATCH`).
The current version is printed at the top of the document.

| Change type | Version bump | Examples |
|---|---|---|
| Typo fix, clarification, non-binding advisory text | PATCH | Rewording §0.4, fixing a table alignment |
| New example added | MINOR | Adding a 12th example for a new `sid.*` function |
| Plant parameter tweak that preserves the example's pedagogy | MINOR | Changing `c` from 2 to 2.2 because the Bode plot looks nicer |
| New required section inside an existing example | MINOR | Adding a new pedagogical beat to `siso` |
| New binding `sid.*` call inside an existing example | MINOR | Adding a `residual` check where there wasn't one |
| Example removal | MAJOR | Retiring `siso` in favor of something new |
| Plant deletion or renumbering | MAJOR | Removing Plant C |
| Breaking change to a helper API's interface | MAJOR | Changing the argument order of the LTI helper |

### 6.2 Deprecation cycle

When an example is removed, its entry MUST stay in this document with
a "Removed in v1.Y.0" note for at least one full minor cycle before
being deleted. This gives language ports time to delete the
corresponding file during their next catch-up pass.

### 6.3 Authoritative reference implementation

The Python example suite on `main` (or the latest merged PR) is the
reference implementation of this spec at any given commit. If a
conflict arises between this document and the Python port, fix the
spec first (in a dedicated commit) and then update the Python port to
match the spec. Do not silently update the spec to match a Python
drift.

This rule mirrors the "spec is source of truth" principle in the root
`CONTRIBUTING.md`.

### 6.4 Language ports lag the spec

New language ports MAY target any stable minor version of this spec.
A port against v1.0.0 is conformant until v1.0.0 is explicitly
deprecated (major version bump). A port SHOULD update to a newer
minor version when practical, but is not required to chase every
MINOR bump immediately.

### 6.5 Changelog

**v1.0.1 — 2026-04-12.** Cleanup pass motivated by spurious
warnings in the v1.0.0 example suite:

- §3.5 `method_comparison`: clarified that ports MUST clamp the
  noise-spectrum array with a small positive floor before applying
  `log10`, so the noise-spectrum overlay plot does not emit a
  divide-by-zero warning at frequencies where the spectrum has been
  clamped to zero per §2.7 of `SPEC.md`.
- §3.1 `siso`: clarified that the residual whiteness test is
  *expected* to FAIL on the Plant A SDOF. The non-parametric BT
  estimator has finite-window bias that shows up as residual
  auto-correlation; the FAIL is the pedagogical content of the
  section, not a bug.
- §3.9 `ltv_disc`:
  - Removed "Preconditioning for numerical stability" section.
    Preconditioning is intentionally disabled in the v1.0 library
    release, so the call emits a runtime warning with no
    pedagogical payload.
  - Removed "Frequency-based lambda tuning" section.
    `ltv_disc_tune(method='frequency')` relies on a non-parametric
    `freq_map` estimate of each state channel; the velocity
    component of a 1-DoF SMD has a bandpass spectrum that
    `freq_map` cannot reliably estimate, so the method finds 0%
    consistency on any reasonable lambda grid and emits "No lambda
    achieved 90% consistency". The library method is not broken,
    but it is the wrong tool for physical SMD plants.
  - Changed the LTV section's simulation parameters from
    `N = 80, L = 15` to `N = 250, L = 30`. With the longer record
    and larger ensemble, the Section 8 residual whiteness test
    passes cleanly in both Python and MATLAB, so the example can
    print PASS without narration.
  - Renumbered remaining sections from 11 to 9.

**v1.0.0 — 2026-04-11.** Initial release.

---

## 7. References

- [`SPEC.md`](SPEC.md) — Algorithm specification. All `sid.*`
  functions referenced in this document are defined there.
- [`../python/CONTRIBUTING.md`](../python/CONTRIBUTING.md) — Python
  notebook conventions, docstring template, and inline comment style.
- [`../matlab/CONTRIBUTING.md`](../matlab/CONTRIBUTING.md) — MATLAB
  function header standard and template files.
- [`../CONTRIBUTING.md`](../CONTRIBUTING.md) — Root contributing
  guide, including the "spec as source of truth" principle and the
  auto-discovery convention.
- [`../python/examples/README.md`](../python/examples/README.md) —
  User-facing summary of the Python example suite. This file is a
  downstream consumer of the spec; it describes what the spec
  requires in plain-English form for end users.
- [`../matlab/examples/README.md`](../matlab/examples/README.md) —
  The MATLAB counterpart to the Python examples README.
- Ljung, L. *System Identification: Theory for the User*, 2nd ed.,
  Prentice Hall, 1999. Referenced by `SPEC.md` for the underlying
  theory of all spectral and state-space identification methods
  exercised in the examples.

---

**End of specification.**

