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

## CI

CI workflows run per-language:

- **Tests** — MATLAB and GNU Octave test suites
- **Lint** — MISS_HIT style and lint checks for MATLAB code
- **Cross-Language Validation** — reference test vector consistency

All checks must pass before merging.
