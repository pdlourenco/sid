# Cross-Language Implementation Strategy for `spa`

## MATLAB-Primary, Python & Julia Secondary

---

## 1. The Three Viable Architectures

### Option A: C Core Library + Language Wrappers

Write the numerical algorithms once in C, expose them through a shared library (`.so`/`.dll`/`.dylib`), and wrap in each language.

```
┌─────────────────────────────────────────────┐
│             libspa (C99, uses FFTW)         │
│  ┌────────────┬──────────────┬────────────┐ │
│  │ covariance │ windowed_dft │ spa_core   │ │
│  └────────────┴──────────────┴────────────┘ │
└──────────┬──────────┬──────────┬────────────┘
           │          │          │
     ┌─────▼──┐  ┌────▼───┐  ┌──▼────┐
     │  MEX   │  │ cffi / │  │ ccall │
     │wrapper │  │ctypes  │  │wrapper│
     │(MATLAB)│  │(Python)│  │(Julia)│
     └────────┘  └────────┘  └───────┘
```

**How each language calls C:**

| Language | Mechanism | Boilerplate | Build Complexity |
|----------|-----------|-------------|------------------|
| MATLAB   | MEX API (`mex` compiler, `mexFunction` entry point) | Moderate — must marshal `mxArray` ↔ C arrays | Medium |
| Python   | `cffi` or `ctypes` (no compilation needed on Python side) | Low — just declare signatures | Low |
| Julia    | `ccall` (built-in, zero overhead, no glue code) | Minimal — one-liner per function | Very Low |

**Pros:**
- Single source of truth for all numerics — guaranteed bit-identical results across languages.
- C is the universal FFI target; every scientific language can call it.
- High performance (though for this algorithm, performance is rarely the bottleneck).
- Easy to add more language bindings later (R, Fortran, Rust, etc.).

**Cons:**
- You own a C codebase: manual memory management, no exceptions, harder debugging.
- Build system complexity: must compile for each OS/architecture and distribute binaries.
- FFTW dependency adds licensing considerations (GPL, or FFTW's commercial license).
- MEX compilation can be fragile across MATLAB versions.

---

### Option B: Pure Idiomatic Implementations (3 Codebases)

Write clean, native code in each language, sharing only a specification and test vectors.

```
┌──────────────┐  ┌──────────────┐  ┌──────────────┐
│  spa.m       │  │  spa.py      │  │  Spa.jl      │
│  (MATLAB)    │  │  (NumPy)     │  │  (Julia)     │
│  Primary     │  │  Port        │  │  Port        │
└──────┬───────┘  └──────┬───────┘  └──────┬───────┘
       │                 │                 │
       └────────┬────────┴────────┬────────┘
                │                 │
         ┌──────▼──────┐  ┌──────▼──────┐
         │ Shared test │  │   Shared    │
         │  vectors    │  │    spec     │
         │  (.mat/.h5) │  │   (docs)    │
         └─────────────┘  └─────────────┘
```

**Pros:**
- Each language feels completely native and idiomatic.
- No build system / compilation / binary distribution headaches.
- MATLAB users get a pure `.m` file — just add to path and use.
- Easy to maintain for contributors who only know one language.

**Cons:**
- Three codebases to maintain in sync.
- Risk of subtle numerical divergence between implementations.
- Bug fixes must be applied three times.

---

### Option C: Write in One Language, Call from Others

Pick one implementation language and call it from the other two.

**Sub-option C1: Julia as core, called from MATLAB and Python**
- MATLAB calls Julia via `jlcall` (community MEX bridge) or system calls.
- Python calls Julia via `juliacall` / `PyJulia`.
- Julia is fast, high-level, has native FFT, and produces easy-to-read numerical code.
- But: Julia startup time is significant, and the MATLAB↔Julia bridge is fragile.

**Sub-option C2: Python as core, called from MATLAB and Julia**
- MATLAB calls Python via its built-in `py.` interface (since R2014b).
- Julia calls Python via `PyCall.jl`.
- But: adds a Python dependency to MATLAB workflows, which MATLAB users dislike.

**Sub-option C3: MATLAB as core, called from Python and Julia**
- Python calls MATLAB via `matlab.engine` (requires MATLAB license).
- Defeats the purpose of being open source.

**Verdict:** Option C is generally the worst choice here because it chains every user to a specific runtime, and the inter-language bridges add latency, startup costs, and fragility.

---

## 2. Recommendation: Option A (C Core) — With a Pragmatic Twist

For a MATLAB-primary project where cross-language numerical consistency matters, the C core approach is the strongest architecture. But pure C for the entire algorithm is overkill. Here's the refined strategy:

### Use C99 for the Core Compute Kernels

The algorithm has just three compute-intensive operations. Everything else is glue code that belongs in the wrapper layer.

```c
// ── libspa.h ──────────────────────────────────────────────

// 1. Biased cross-covariance estimation
//    Computes R_xy(τ) for τ = 0..max_lag
int spa_covariance(
    const double *x,       // input signal x, length N
    const double *y,       // input signal y, length N  (can equal x for auto-cov)
    int N,                 // number of samples
    int max_lag,           // M
    double *cov_real,      // output: real part of R_xy(0..max_lag)
    double *cov_imag       // output: imag part (zero for real signals)
);

// 2. Windowed spectral estimate at arbitrary frequencies
//    Computes Φ(ω) = Σ R(τ) W(τ) e^{-jωτ} for each frequency
int spa_windowed_spectrum(
    const double *cov_real,    // covariance R(0..M), length M+1
    const double *cov_imag,    // imaginary part
    int M,                     // max lag (window size)
    const double *freqs,       // frequency vector, length n_freq (rad/sample)
    int n_freq,
    int window_type,           // 0=Hann, 1=Hamming
    double *spec_real,         // output: real part of Φ(ω), length n_freq
    double *spec_imag          // output: imag part
);

// 3. Windowed spectral estimate using FFT (fast path for uniform grid)
int spa_windowed_spectrum_fft(
    const double *cov_real,
    const double *cov_imag,
    int M,
    int n_fft,                 // FFT length (>= 2*M+1, power of 2)
    int window_type,
    double *spec_real,         // output: length n_fft
    double *spec_imag
);
```

The key design principles for the C layer:

- **No memory allocation in C.** The caller (MATLAB/Python/Julia) allocates all buffers. This eliminates memory management bugs and lets each language use its native array types.
- **No FFT library dependency in the C core.** The `spa_windowed_spectrum_fft` function takes pre-computed covariance and applies the window + DFT. The actual FFT call happens in the wrapper layer using each language's native FFT (MATLAB's `fft`, numpy's `fft`, Julia's `FFTW`). This sidesteps FFTW licensing entirely.
- **Stateless functions.** No global state, no structs to initialize/destroy. Pure functions operating on arrays.
- **Integer return codes** for error handling (0 = success).

### Actually — Reconsider Whether C Is Needed at All

Let's be honest about the computational profile of this algorithm:

```
Covariance estimation:  O(N × M)  or O(N log N) via FFT
Windowed DFT:           O(M × n_freq) or O(M log M) via FFT  
Transfer function:      O(n_freq) — just division
Uncertainty:            O(n_freq) — just arithmetic
```

With typical values of N = 1000–100000, M = 30, and n_freq = 128, the entire computation finishes in **milliseconds** in any language. The algorithm is not compute-bound. The value is in getting the **mathematics and defaults right**, not in squeezing out cycles.

This tilts the balance toward **Option A-Lite: a specification-driven approach with a C reference implementation but idiomatic primary implementations.**

---

## 3. Recommended Architecture: Spec-Driven with C Reference

```
┌──────────────────────────────────────────────────────────┐
│                    spa-spec (repository)                  │
│                                                          │
│  ┌─────────────┐  ┌──────────────┐  ┌─────────────────┐ │
│  │  SPEC.md    │  │ test_vectors/│  │  libspa_ref/    │ │
│  │  (algorithm │  │  input_*.mat │  │  (C99 reference │ │
│  │  definition)│  │  output_*.mat│  │   implementation│ │
│  │             │  │  (HDF5 also) │  │   for edge cases│ │
│  └─────────────┘  └──────────────┘  │   and validation│ │
│                                     └─────────────────┘ │
│  ┌──────────────────────────────────────────────────────┐│
│  │  Language Packages (each is a standalone package)    ││
│  │                                                      ││
│  │  spa-matlab/          spa-python/       spa-julia/   ││
│  │  ├── +spa/            ├── spa/          ├── src/     ││
│  │  │   ├── spa.m        │   ├── core.py   │   └─Spa.jl││
│  │  │   ├── bode.m       │   ├── plot.py   ├── test/   ││
│  │  │   └── spafdr.m     │   └── ...       └── ...     ││
│  │  ├── test/            ├── tests/                     ││
│  │  │   └── test_spa.m   │   └── test_spa.py            ││
│  │  └── README.md        └── pyproject.toml             ││
│  └──────────────────────────────────────────────────────┘│
└──────────────────────────────────────────────────────────┘
```

### The Spec Document (`SPEC.md`)

This is the most important artifact. It defines:

1. **Exact mathematical formulas** with unambiguous notation.
2. **Default parameter values** (window size, frequency grid, etc.).
3. **Normalization conventions** (this is where cross-language bugs hide).
4. **Edge case behavior** (what happens when N < M, when Φ_u ≈ 0, etc.).
5. **Reference test vectors** with expected outputs to 15 significant digits.

### Test Vectors (`test_vectors/`)

Generated once from the C reference implementation (or from MATLAB itself if you have a license), stored in both `.mat` (v7.3 / HDF5) and `.h5` formats so all three languages can load them.

```
test_vectors/
├── siso_whitenoise.mat     # White noise input, known system
├── siso_colored.mat        # Colored noise, first-order system
├── mimo_2x2.mat            # 2-input 2-output system
├── timeseries_ar1.mat      # No input, AR(1) output
├── short_data.mat          # N=50, edge case
├── custom_freqs.mat        # Non-default frequency grid
└── README.md               # Description of each test case
```

Each file contains: `u`, `y`, `Ts`, `M`, `freqs`, and the expected outputs `G_real`, `G_imag`, `G_std`, `Phi_v`, `Phi_v_std`.

---

## 4. Language-Specific Implementation Details

### MATLAB (Primary Target)

```matlab
function result = spa(y, u, varargin)
%SPA Estimate frequency response via Blackman-Tukey spectral analysis.
%
%   G = spa(y, u)
%   G = spa(y, u, 'WindowSize', M, 'Frequencies', w)
%   G = spa(y, [], ...)    % time series (no input)
%
%   Returns a struct with fields:
%     .Frequency      - frequency vector (rad/sample)
%     .Response       - complex frequency response G(e^{jw})
%     .ResponseStd    - standard deviation of G
%     .NoiseSpectrum  - Phi_v(w)
%     .NoiseSpectrumStd
%     .SampleTime
%     .WindowSize

    p = inputParser;
    addRequired(p, 'y');
    addRequired(p, 'u');
    addParameter(p, 'WindowSize', []);
    addParameter(p, 'Frequencies', []);
    addParameter(p, 'SampleTime', 1.0);
    addParameter(p, 'NumFrequencies', 128);
    parse(p, y, u, varargin{:});

    % ... implementation using MATLAB's built-in fft, xcorr, etc.
end
```

**Key decisions for the MATLAB package:**

- **Pure `.m` files only.** No MEX, no compiled code, no Java. Users just add the folder to their path. This is critical for adoption — MATLAB users expect `addpath` to be sufficient.
- **No toolbox dependencies.** Only core MATLAB functions: `fft`, `ifft`, `conv`, basic array ops. Do not require Signal Processing Toolbox or any other toolbox.
- **Return a struct**, not a custom class (unless you want to implement a full `idfrd`-compatible object, which is a much larger scope). A struct with clear field names is immediately usable.
- **Namespace via `+spa` package folder** so users call `spa.spa(y,u)` or import it. Alternatively, use a flat function name with a distinctive prefix to avoid collisions with the toolbox function.
- **Octave compatibility.** Test against GNU Octave. Avoid MATLAB-only syntax (e.g., string arrays with `"..."` — use `'...'` char vectors). This gives free/open access to users without a MATLAB license.

### Python

```python
import numpy as np
from dataclasses import dataclass

@dataclass
class SpaResult:
    frequencies: np.ndarray
    response: np.ndarray          # complex
    response_std: np.ndarray
    noise_spectrum: np.ndarray
    noise_spectrum_std: np.ndarray
    sample_time: float
    window_size: int

def spa(y, u=None, *, window_size=None, frequencies=None,
        sample_time=1.0, n_frequencies=128):
    """Blackman-Tukey spectral analysis for system identification."""
    ...
```

- Pure Python + NumPy. No compiled extensions needed.
- Publish on PyPI as `pip install openspa` (or similar).
- Include matplotlib-based `bode_plot()` and `spectrum_plot()`.

### Julia

```julia
module OpenSpa

using FFTW, LinearAlgebra

struct SpaResult
    frequencies::Vector{Float64}
    response::Vector{ComplexF64}
    response_std::Vector{Float64}
    noise_spectrum::Vector{Float64}
    noise_spectrum_std::Vector{Float64}
    sample_time::Float64
    window_size::Int
end

function spa(y::AbstractVector, u::Union{AbstractVector,Nothing}=nothing;
             window_size::Union{Int,Nothing}=nothing,
             frequencies::Union{AbstractVector,Nothing}=nothing,
             sample_time::Float64=1.0,
             n_frequencies::Int=128)::SpaResult
    ...
end

end # module
```

- Register as a Julia package.
- Use `FFTW.jl` for FFT (standard in Julia ecosystem).
- Add `Plots.jl` or `Makie.jl` recipe for Bode plots.

---

## 5. The C Reference Implementation

Even though the primary implementations are native, a C reference is valuable for:

1. **Generating authoritative test vectors** that don't depend on any one language's FFT implementation.
2. **Embedding via MEX** if someone needs maximum performance on huge datasets.
3. **Binding to other languages** in the future (R, Rust, etc.).

```
libspa_ref/
├── CMakeLists.txt
├── include/
│   └── spa.h            # Public API (the header shown above)
├── src/
│   ├── covariance.c     # Biased cross-covariance
│   ├── window.c         # Hann/Hamming window generation
│   ├── spectrum.c       # Direct DFT at arbitrary frequencies
│   └── spa.c            # Orchestrator: covariance → window → DFT → ratios
├── test/
│   ├── test_covariance.c
│   └── test_spa.c
├── matlab_mex/
│   └── spa_mex.c        # MEX gateway (optional, for power users)
├── python_cffi/
│   └── spa_cffi_build.py
└── julia_ccall/
    └── libspa.jl
```

**Build with CMake:**
```cmake
cmake_minimum_required(VERSION 3.15)
project(libspa C)

set(CMAKE_C_STANDARD 99)

# No external dependencies! Uses a minimal built-in DFT
# for the reference. Language wrappers use their native FFT.
add_library(spa SHARED
    src/covariance.c
    src/window.c
    src/spectrum.c
    src/spa.c
)

target_include_directories(spa PUBLIC include/)

# Optional: MEX target
if(MATLAB_FOUND)
    matlab_add_mex(NAME spa_mex SRC matlab_mex/spa_mex.c LINK_TO spa)
endif()
```

---

## 6. Cross-Language Consistency Guarantee

The critical question is: **do all three implementations produce identical results?**

### Strategy: Shared Acceptance Tests

```
tests/
├── generate_test_vectors.m    # MATLAB script to create ground truth
├── test_vectors/
│   ├── case_001.h5            # HDF5 (readable by all 3 languages)
│   └── case_002.h5
├── run_matlab_tests.m
├── run_python_tests.py
└── run_julia_tests.jl
```

Each test script:
1. Loads the test vector (input data + expected output).
2. Runs the native `spa` implementation.
3. Asserts results match within tolerance.

**Tolerance levels:**
- `|computed - expected| / |expected| < 1e-10` for the FFT path (default frequencies).
- `|computed - expected| / |expected| < 1e-8` for the direct DFT path (custom frequencies, minor floating-point ordering differences).

### Where Numerical Differences Creep In

| Source | Risk | Mitigation |
|--------|------|------------|
| FFT implementation differences | Low — all use Cooley-Tukey | Specify FFT length and verify against known DFT |
| Covariance computation order | Medium — summation order affects rounding | Specify that summation proceeds `t = 1..N-τ` in ascending order |
| Window function evaluation | Low | Define window values exactly: `W(τ) = 0.5 * (1 + cos(π*τ/M))` |
| Division by small Φ_u | High — regularization choices diverge | Spec defines: if `Φ_u(ω) < ε * max(Φ_u)`, set `G(ω) = NaN` |
| Normalization convention | High — the #1 source of cross-language bugs | Spec defines exact formula with no ambiguity |

---

## 7. Repository Layout

```
openspa/
├── SPEC.md                     # Algorithm specification (normative)
├── LICENSE                     # MIT
├── README.md
│
├── test_vectors/               # Shared ground truth
│   ├── generate.m
│   └── *.h5
│
├── matlab/                     # MATLAB/Octave package
│   ├── +openspa/
│   │   ├── spa.m
│   │   ├── spafdr.m
│   │   ├── bodeplot.m
│   │   └── spectrumplot.m
│   ├── test/
│   │   └── test_spa.m
│   ├── install.m               # Adds package to path
│   └── README.md
│
├── python/                     # Python package
│   ├── src/openspa/
│   │   ├── __init__.py
│   │   ├── core.py
│   │   └── plotting.py
│   ├── tests/
│   ├── pyproject.toml
│   └── README.md
│
├── julia/                      # Julia package
│   ├── src/OpenSpa.jl
│   ├── test/runtests.jl
│   ├── Project.toml
│   └── README.md
│
└── libspa_ref/                 # C reference (optional)
    ├── CMakeLists.txt
    ├── include/spa.h
    ├── src/*.c
    └── bindings/
        ├── spa_mex.c           # MEX gateway
        ├── spa_cffi.py         # Python cffi build
        └── libspa.jl           # Julia ccall wrapper
```

---

## 8. Development Roadmap

| Phase | What | Deliverable |
|-------|------|-------------|
| **1** | Write `SPEC.md` | Complete algorithm specification with formulas, defaults, edge cases, normalization |
| **2** | MATLAB implementation | Pure `.m` files, SISO + time series, tested against MathWorks `spa` |
| **3** | Generate test vectors | HDF5 files from MATLAB, covering all cases |
| **4** | Python port | NumPy implementation passing all test vectors |
| **5** | Julia port | Julia implementation passing all test vectors |
| **6** | MIMO extension | All three languages, new test vectors |
| **7** | Uncertainty estimation | Asymptotic variance formulas in all three |
| **8** | Plotting | Bode + spectrum plots with confidence bands, per language |
| **9** | `spafdr` | Frequency-dependent resolution variant |
| **10** | C reference + MEX | Optional high-performance path |
| **11** | Packaging & docs | PyPI, Julia registry, MATLAB File Exchange, Octave package |

---

## 9. Summary of Recommendation

**Don't write a C library and wrap it.** Instead:

1. **Write a rigorous spec** (`SPEC.md`) that pins down every formula, default, and edge case.
2. **Implement natively in each language** — pure `.m` for MATLAB (your primary audience), pure NumPy for Python, pure Julia for Julia. Each feels native to its users.
3. **Share test vectors** (HDF5 format, readable everywhere) to guarantee cross-language consistency.
4. **Keep a C reference implementation** in the repo for generating authoritative test outputs and for anyone who needs a compiled library, but don't make it the primary codebase.

This gives you the best of both worlds: MATLAB users get a drop-in `.m` file that requires no compilation, no MEX, no external dependencies, and zero friction — while Python and Julia users get equally idiomatic packages that produce bit-compatible results.
