# Contributing to sid-matlab

Contributions are welcome via issues and pull requests. Please ensure that
`tests/runAllTests.m` passes on both MATLAB and Octave before submitting —
the CI pipeline checks both platforms automatically.

---

## Function Header Standard

Every `.m` function file in this repository (public and internal) **must**
follow the canonical header template below. This ensures consistency,
enables MATLAB `help` to display useful documentation, and keeps a clear
link between code and the algorithm specification.

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

Internal helper functions live in the `internal/` directory and use the same
`sid` prefix with camelCase (e.g., `sidCov`, `sidValidateData`,
`sidLTVcosmicSolve`).

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

---

## Testing

- All tests: `tests/runAllTests.m`
- All examples: `examples/runAllExamples.m`
- Both must pass on **MATLAB R2016b+** and **GNU Octave 8.0+**
- CI runs lint (`miss_hit`) and tests on both platforms automatically

---

## Example and Test Scripts

Example scripts (`examples/`) and test scripts (`tests/`) use `%%` section
markers and do not require the full function header template. They should
have a brief `%%` title and description at the top.
