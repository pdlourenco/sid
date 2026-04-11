# Contributing to sid

Contributions are welcome via issues and pull requests.

## Project Structure

sid is a multi-language system identification toolbox built around a shared
mathematical specification:

- [`spec/`](spec/) — Algorithm specification and mathematical derivations
  (single source of truth for all implementations)
- [`testdata/`](testdata/) — Cross-language reference test vectors (JSON; see [README](testdata/README.md) for format)
- [`matlab/`](matlab/) — MATLAB/Octave implementation (stable)
- [`python/`](python/) — Python implementation (stable)
- [`julia/`](julia/) — Julia implementation (planned)

## Specification as Source of Truth

**[`spec/SPEC.md`](spec/SPEC.md) is the binding contract for every
implementation.** This applies equally to human contributors and AI coding
agents. Read it before touching any algorithmic code.

### Core rules

1. **The spec defines required behaviour.** Defaults, edge cases, error
   conditions, output struct fields, normalization conventions, and
   numerical diagnostics are all part of the contract — not just the core
   formulas. If a requirement is in `SPEC.md`, every implementation must
   satisfy it.

2. **Implementations conform to the spec, not to each other.** MATLAB is
   *not* a ground truth. "The MATLAB version does it this way" is not a
   valid justification — if MATLAB and the spec disagree, MATLAB is wrong.
   Cross-language numerical equivalence is a *consequence* of every
   implementation independently satisfying the spec, not a goal pursued
   by copying one implementation into another.

3. **Fix the spec first when it is ambiguous or wrong.** If you find that
   the spec is silent, contradictory, or disagrees with a clearly correct
   algorithm, **update `SPEC.md` first** (in a dedicated commit, ideally
   with reviewer approval), *then* update the implementations to match.
   Never silently encode behaviour that the spec does not describe — the
   next language port will have no way to recover the same decision.

4. **Shared helpers can cause shared drift.** A bug in a helper like
   `sidValidateData` / `validate_data` can make *every* downstream caller
   silently violate the spec in the same way, which will not be caught by
   cross-validation tests (because the MATLAB and Python CI both exhibit
   the same drift). When you touch a shared helper, audit each caller
   against the relevant spec section.

5. **Cross-language reference vectors are a check, not a proof.** Tests
   in [`testdata/`](testdata/) verify that MATLAB and Python agree
   numerically on specific input data. They do *not* prove either
   implementation satisfies the spec. New features must come with a
   direct spec-to-implementation read-through, not only reference
   vector comparisons.

### Workflow for algorithmic changes

When adding or modifying an algorithm, a function default, an output
field, or any edge-case behaviour:

1. **Read the relevant `SPEC.md` section end to end** before touching
   code. Note every normative statement (defaults, bounds, NaN handling,
   regularization thresholds, output field names and shapes, warnings).
2. **If the spec does not cover what you need**, update `SPEC.md` in a
   dedicated commit before writing code. Include the rationale in the
   commit message.
3. **Implement in every maintained language** (MATLAB and Python today),
   referencing the same spec section. Use the `SPEC.md §X.Y` comment
   convention to mark each step.
4. **Cite the spec in the PR description.** State which sections the
   change implements or modifies, and call out any spec updates.
5. **Write tests against spec requirements**, not against the current
   output of a reference implementation. A test of the form "assert the
   result equals what MATLAB returned today" does not detect joint drift.

### Checklist for reviewers (and self-review)

- [ ] Is every new default/bound/threshold covered by `SPEC.md`?
- [ ] Do all touched functions cite the relevant `SPEC.md §` in comments?
- [ ] If behaviour was added or changed, was `SPEC.md` updated first?
- [ ] Are the MATLAB and Python behaviours derived independently from the
      spec, rather than ported copy-by-copy from one to the other?
- [ ] Are tests written against the spec's requirements, not against the
      current output of the other language?

## General Guidelines

- The spec rules above apply to every implementation. The per-language
  guides below cover language-specific style, naming, and testing.
- Cross-language test vectors in [`testdata/`](testdata/) ensure numerical
  consistency across implementations. New algorithms should include
  reference vectors.
- The project is MIT-licensed. See [`LICENSE`](LICENSE).

## Code Style

The root [`.editorconfig`](.editorconfig) enforces basic formatting rules
(UTF-8, LF line endings, trailing whitespace). Language-specific linting
and style rules are documented in each language's contributing guide.

## Language-Specific Guidelines

Each language has its own contributing guide with conventions for naming,
documentation, code style, and testing:

- **MATLAB/Octave**: [`matlab/CONTRIBUTING.md`](matlab/CONTRIBUTING.md)
- **Python**: [`python/CONTRIBUTING.md`](python/CONTRIBUTING.md)

## Test and Example Auto-Discovery

Test and example runners in every language **discover files by naming
convention** — there is no hardcoded manifest to maintain. To add a test
or example, create a file matching the pattern for that language:

| Language | Tests | Examples |
|----------|-------|----------|
| MATLAB/Octave | `test_*.m` | `example*.m` |
| Python | `test_*.py` | `example_*.py` |
| Julia | `test_*.jl` | `example_*.jl` |

Runners sort discovered files alphabetically and execute them in order.
**Do not maintain hardcoded file lists** — auto-discovery prevents the
common failure mode where a new test exists but is never executed because
it was not added to a manifest.

### Templates

Each language provides template files for tests and examples. Copy the
template when creating a new file — it includes the runner instrumentation
variables that enable per-file progress tracking in CI output.

| Language | Test template | Example template |
|----------|--------------|-----------------|
| MATLAB/Octave | `matlab/tests/test_template.m` | `matlab/examples/example_template.m` |

Each language's contributing guide documents the discovery mechanism and
templates in detail. When starting a new language port, implement the runner
with auto-discovery from day one and provide starter templates.

## CI

CI workflows run per-language:

- **MATLAB Tests** — MATLAB and GNU Octave test suites
- **MATLAB Lint** — MISS_HIT style/lint checks and function header validation
- **Python Lint** — ruff style/lint checks and docstring validation
- **Python Tests** — pytest on Python 3.10–3.13
- **Cross-Language Validation** — reference test vector consistency

All checks must pass before merging.
