# Contributing to sid

Contributions are welcome via issues and pull requests.

## Project Structure

sid is a multi-language system identification toolbox built around a shared
mathematical specification:

- [`spec/`](spec/) — Algorithm specification and mathematical derivations
  (single source of truth for all implementations)
- [`testdata/`](testdata/) — Cross-language reference test vectors (JSON; see [README](testdata/README.md) for format)
- [`matlab/`](matlab/) — MATLAB/Octave implementation (stable)
- [`python/`](python/) — Python implementation (planned)
- [`julia/`](julia/) — Julia implementation (planned)

## General Guidelines

- All implementations follow the algorithms defined in [`spec/SPEC.md`](spec/SPEC.md).
  When in doubt, the spec is authoritative.
- Cross-language test vectors in [`testdata/`](testdata/) ensure numerical
  consistency across implementations. New algorithms should include reference
  vectors.
- The project is MIT-licensed. See [`LICENSE`](LICENSE).

## Code Style

The root [`.editorconfig`](.editorconfig) enforces basic formatting rules
(UTF-8, LF line endings, trailing whitespace). Language-specific linting
and style rules are documented in each language's contributing guide.

## Language-Specific Guidelines

Each language has its own contributing guide with conventions for naming,
documentation, code style, and testing:

- **MATLAB/Octave**: [`matlab/CONTRIBUTING.md`](matlab/CONTRIBUTING.md)

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

- **Tests** — MATLAB and GNU Octave test suites
- **Lint** — MISS_HIT style and lint checks for MATLAB code
- **Cross-Language Validation** — reference test vector consistency

All checks must pass before merging.
