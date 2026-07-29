"""Fail the docs build if a public Python function has no API page.

Runs as a `mkdocs-gen-files` script. The Python API reference is hand-authored
(one small stub per function, each carrying a curated cross-link to its MATLAB
counterpart) rather than generated, so nothing otherwise stops a newly exported
function from silently missing the site. This is the structural gate that closes
that hole: it derives the expected page set from `sid.__all__` -- the same
auto-discovery-over-manifests principle the test and example suites follow -- and
raises on any mismatch in either direction.

Both directions matter. A missing page means an undocumented public function; an
orphan page means the site advertises something the package no longer exports.

Standalone: `python scripts/check_python_api_pages.py` (exit 1 on mismatch).
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
API_DIR = REPO_ROOT / "docsite" / "api" / "python"
SUMMARY = API_DIR / "SUMMARY.md"

# Exported names that are deliberately not one-page-per-name: the version string,
# and the result dataclasses + exception which share `results.md`.
NON_FUNCTION_EXPORTS = {
    "__version__",
    "CompareResult",
    "FreqMapResult",
    "FreqResult",
    "FrozenResult",
    "LTVIOResult",
    "LTVResult",
    "ResidualResult",
    "SidError",
    "SpectrogramResult",
}

# Pages that exist for reasons other than a single exported function.
NON_FUNCTION_PAGES = {"index", "results", "SUMMARY"}


def public_functions() -> set[str]:
    """Names in `sid.__all__` that should each have their own API page."""
    sys.path.insert(0, str(REPO_ROOT / "python"))
    import sid

    return {name for name in sid.__all__ if name not in NON_FUNCTION_EXPORTS}


def existing_pages() -> set[str]:
    return {p.stem for p in API_DIR.glob("*.md") if p.stem not in NON_FUNCTION_PAGES}


def check() -> list[str]:
    """Return a list of problems; empty means the API reference is complete."""
    expected = public_functions()
    actual = existing_pages()
    problems = []

    for name in sorted(expected - actual):
        problems.append(
            f"sid.{name} is exported in __all__ but has no page at "
            f"docsite/api/python/{name}.md -- add the stub and list it in SUMMARY.md"
        )
    for name in sorted(actual - expected):
        problems.append(
            f"docsite/api/python/{name}.md documents '{name}', which sid.__all__ "
            f"does not export -- remove the page or re-export the function"
        )

    # literate-nav renders only what SUMMARY.md lists, so a stub absent from it
    # is invisible on the site even though the file exists.
    if SUMMARY.exists():
        summary = SUMMARY.read_text(encoding="utf-8")
        for name in sorted(expected & actual):
            if f"({name}.md)" not in summary:
                problems.append(
                    f"docsite/api/python/{name}.md exists but SUMMARY.md does not "
                    f"link it -- the page would not appear in the site nav"
                )
    else:
        problems.append(f"missing nav file: {SUMMARY}")

    return problems


def main() -> None:
    problems = check()
    if problems:
        raise SystemExit(
            "Python API reference is out of sync with sid.__all__:\n  - "
            + "\n  - ".join(problems)
        )


if __name__ == "__main__":
    main()
    print("Python API reference is in sync with sid.__all__")
else:
    # Imported by mkdocs-gen-files during the build: fail the build on mismatch.
    main()
