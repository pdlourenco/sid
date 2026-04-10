# Contributing to sid — MATLAB/Octave

This guide covers MATLAB/Octave-specific contribution standards.
For general project guidelines, see the root [CONTRIBUTING.md](../CONTRIBUTING.md).

> **⚠ Read this first.** Before writing or modifying any algorithmic code,
> read the [Specification as Source of Truth](../CONTRIBUTING.md#specification-as-source-of-truth)
> section in the root contributing guide. `spec/SPEC.md` is the binding
> contract for this implementation — the MATLAB code is *not* the ground
> truth. If the current MATLAB behaviour disagrees with the spec, the
> spec wins and the MATLAB code must be fixed.

Please ensure that `matlab/tests/runAllTests.m` passes on both MATLAB and
Octave before submitting — the CI pipeline checks both platforms automatically.

---

## Function Header Standard

Every `.m` function file (public and private) **must** follow the canonical
header template below. This ensures consistency, enables MATLAB `help` to
display useful documentation, and keeps a clear link between code and the
algorithm specification.

### Canonical Template

```matlab
function [out1, out2] = sidFunctionName(in1, in2, varargin)
% SIDFUNCTIONNAME Brief one-line description.
%
%   out = sidFunctionName(in1, in2)
%   out = sidFunctionName(in1, in2, 'Option', value)
%   out = sidFunctionName(in1, in2, posArg)
%
%   Extended description paragraph(s). What the function does, context,
%   and any important notes.
%
%   INPUTS:
%     in1 - Description, (dimension) type. Constraints.
%     in2 - Description. Use [] for alternative.
%
%   NAME-VALUE OPTIONS:
%     'OptionName' - Description. Default: value.
%
%   OUTPUTS:
%     out1 - Description, (dimension) type.
%     out2 - Description.
%
%   EXAMPLES:
%     % Basic usage
%     result = sidFunctionName(x, y);
%
%   ALGORITHM:
%     1. Step description.
%     2. Step description.
%
%   REFERENCES:
%     Author, "Title", Publisher, Year. Sections X.Y.
%
%   SPECIFICATION:
%     SPEC.md §X.Y — Section Title
%
%   See also: sidRelated1, sidRelated2
%
%   Changelog:
%   YYYY-MM-DD: Description by Author Name.
%
%  -----------------------------------------------------------------------
%   Copyright (c) 2026 Pedro Lourenço, All rights reserved.
%   This code is released under the MIT License. See LICENSE file in the
%   project root for full license information.
%
%   This function is part of the Open Source System Identification
%   Toolbox (SID).
%   For full documentation and examples, visit
%   https://github.com/pdlourenco/sid-matlab
%  -----------------------------------------------------------------------
```

### Section Order (fixed)

Sections must appear in this exact order:

1. `% FUNCTIONNAME` — brief one-line description (ALL CAPS function name)
2. **Usage signatures** — all calling forms, including positional variants
3. **Extended description** — one or more paragraphs
4. **`INPUTS:`** — bullet list of required positional arguments
5. **`NAME-VALUE OPTIONS:`** — bullet list of optional name-value pairs
6. **`OUTPUTS:`** — bullet list of return values
7. **`EXAMPLES:`** — runnable code snippets with inline comments
8. **`ALGORITHM:`** — numbered steps describing the computational approach
9. **`REFERENCES:`** — academic citations (author, title, publisher, year)
10. **`SPECIFICATION:`** — cross-reference to `SPEC.md` section
11. **`See also:`** — comma-separated list of related functions
12. **`Changelog:`** — entries in `YYYY-MM-DD: Description by Author.` format
13. **Copyright block** — MIT license notice with dashed separators

### Required vs Optional Sections

| Section | Required? |
|---------|-----------|
| One-line description | Always |
| Usage signatures | Always |
| Extended description | Always (can be brief for simple helpers) |
| INPUTS | Always (unless the function takes no arguments) |
| NAME-VALUE OPTIONS | Only if the function accepts name-value pairs |
| OUTPUTS | Always (unless the function returns nothing) |
| EXAMPLES | Always |
| ALGORITHM | Only for non-trivial algorithms |
| REFERENCES | Only when citing academic papers |
| SPECIFICATION | When a corresponding `SPEC.md` section exists |
| See also | Always |
| Changelog | Always |
| Copyright block | Always |

### Result Struct Definitions

All result struct types are centrally defined in
[`sid/sidResultTypes.m`](sid/sidResultTypes.m). When adding or modifying a
function that produces a result struct:

1. **Update `sidResultTypes.m`** to reflect any new or changed fields.
2. **Keep the function's `OUTPUTS:` section in sync** with the central
   definition — both must list the same fields with matching dimensions.
3. **Reference the central definition** from any function that *consumes* a
   result struct (e.g., plotting functions, `sidCompare`, `sidResidual`).
   Use the format: `(see sidResultTypes §N)` where N is the section number.

This mirrors the Python approach where result types are formally defined as
frozen dataclasses in `python/sid/_results.py`.

### Key Rules

- **Section headings are always PLURAL**: `INPUTS:`, `OUTPUTS:`, `EXAMPLES:`,
  `REFERENCES:`. Exception: `See also:` and `Changelog:` keep their
  traditional casing.
- **Positional calling forms** go in the usage signatures at the top — do not
  create a separate "POSITIONAL SYNTAX" section.
- **SPEC.md cross-reference**: if the function implements or relates to a
  section of `SPEC.md`, add a `SPECIFICATION:` entry pointing to it
  (e.g., `SPEC.md §2 — Blackman-Tukey Spectral Analysis`).
- **Copyright block** uses dashed separators (`% -------...`) and is always
  the last part of the header comment.
- **Changelog dates** use ISO 8601 format (`YYYY-MM-DD`).
- **`See also:`** appears exactly once per file.
- **MATLAB `help` compatibility**: all documentation sections must appear
  *before* the copyright block, inside a single contiguous comment block.
  Do not place documentation after the copyright separator.

---

## Naming Conventions

Functions follow the pattern:

```
sid [Domain] [Method/Variant]
 │     │          │
 │     │          └── BT, BTFDR, ETFE, ARX, N4SID, AR, ...
 │     │
 │     └──────────── Freq, TF, SS, TS, LTV, ...
 │
 └─────────────────── system identification (root)
```

Examples: `sidFreqBT`, `sidLTVdisc`, `sidBodePlot`, `sidModelOrder`.

Internal helper functions live in the `sid/private/` directory and use the same
`sid` prefix with camelCase (e.g., `sidCov`, `sidValidateData`,
`sidLTVcosmicSolve`). MATLAB's `private` directory mechanism makes these
functions available to all functions in `sid/` without requiring them on the
path, and invisible to end users.

---

## Code Style

| Rule | Value |
|------|-------|
| Indentation | 4 spaces (no tabs) |
| Line length | 100 characters max |
| Line endings | LF |
| Charset | UTF-8 |
| Semicolons | Required (suppress console output) |
| Changelog dates | ISO 8601 (`YYYY-MM-DD`) |

See `.editorconfig` and `miss_hit.cfg` for automated enforcement.

### Inline Comments

Code comments within function bodies should make the mathematical intent
clear and link back to the specification. Follow these guidelines:

**Section separators.** Use `% ---- Name ----` to mark major computational
phases. These should correspond to steps in the function's ALGORITHM header
section:

```matlab
% ---- Build data matrices (SPEC.md §8.3.2) ----
[D, Xp] = sidLTVbuildDataMatrices(X, U);
```

**SPEC.md cross-references.** When a code block implements a specific
equation or algorithm step from SPEC.md, cite the section number:

```matlab
% Schur complement forward pass (SPEC.md §8.3.4, Eq. 8.3):
%   Lambda(k) = S(k) - lambda(k-1)^2 * Lambda(k-1)^{-1}
Lbd(:,:,k) = S(:,:,k) - lambda(k-1)^2 * (Lbd(:,:,k-1) \ I);
```

**Mathematical steps.** Annotate non-obvious operations — matrix
inversions, Schur complements, spectral transformations, and
regularization terms. Write the formula in comment notation before
the code that implements it:

```matlab
% G(w) = Phi_yu(w) / Phi_u(w)  — transfer function estimate
G = Phi_yu ./ Phi_u;

% Phi_v(w) = Phi_y(w) - |Phi_yu(w)|^2 / Phi_u(w)  — noise spectrum
Phi_v = Phi_y - abs(Phi_yu).^2 ./ Phi_u;
```

**Variable-to-notation mapping.** When a variable name differs from the
mathematical notation in SPEC.md, state the correspondence on first use:

```matlab
% Lbd corresponds to Lambda_k in SPEC.md §8.3 (forward Schur complement)
Lbd = zeros(d, d, N);
```

**Dimensions.** Annotate array dimensions on the line that creates or
returns them, using trailing comments:

```matlab
Ryy = sidCov(y, y, M);   % (M+1 x ny x ny) biased auto-covariance
```

**What not to comment.** Do not comment self-explanatory operations
(loop counters, standard input parsing, trivial assignments). Focus
comments on *why*, not *what*:

```matlab
% Bad: increment k by 1
k = k + 1;

% Good: skip the first segment (it has incomplete overlap)
k = k + 1;
```

---

## Testing

- All tests: `matlab/tests/runAllTests.m`
- All examples: `matlab/examples/runAllExamples.m`
- Both must pass on **MATLAB R2016b+** and **GNU Octave 8.0+**
- CI runs lint (`miss_hit`) and tests on both platforms automatically

### Auto-discovery

The test and example runners **discover files by naming convention** — there
is no manifest to edit. To add a new test or example, simply create a file
that follows the naming pattern:

| Directory | Pattern | Example |
|-----------|---------|---------|
| `matlab/tests/` | `test_*.m` | `test_sidNewFeature.m` |
| `matlab/examples/` | `example*.m` | `exampleNewFeature.m` |

The runners use `dir('test_*.m')` and `dir('example*.m')` respectively.
Files are sorted alphabetically and executed in that order. **Do not add
entries to a hardcoded list** — it defeats the purpose of auto-discovery
and creates merge conflicts.

This convention applies across all language ports (Python, Julia). See the
root [`CONTRIBUTING.md`](../CONTRIBUTING.md) for the cross-language naming
table.

---

## Example and Test Scripts

Example scripts (`examples/`) and test scripts (`tests/`) use `%%` section
markers and do not require the full function header template. They should
have a brief `%%` title and description at the top.

### Templates

Copy the appropriate template when creating new examples or tests:

| Type | Template | Key variable |
|------|----------|-------------|
| Test | [`tests/test_template.m`](tests/test_template.m) | `runner__nPassed` |
| Example | [`examples/example_template.m`](examples/example_template.m) | `runner__nCompleted` |

The templates include the instrumentation variables that the runners use to
track per-file progress (section counts for examples, assertion counts for
tests). Following the template ensures your file integrates with the CI
summary output automatically.
