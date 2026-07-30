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

import inspect
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
API_DIR = REPO_ROOT / "docsite" / "api" / "python"
SUMMARY = API_DIR / "SUMMARY.md"
RESULTS_PAGE = API_DIR / "results.md"
FUNCTION_INDEX = REPO_ROOT / "docsite" / "api" / "index.md"
PYTHON_OVERVIEW = API_DIR / "index.md"

# Pages that exist for reasons other than a single exported function.
NON_FUNCTION_PAGES = {"index", "results", "SUMMARY"}


def _sid():
    sys.path.insert(0, str(REPO_ROOT / "python"))
    import sid

    return sid


def partition_exports() -> tuple[set[str], set[str]]:
    """Split `sid.__all__` into (functions, everything else).

    Determined by introspection, not a hand-kept list: a maintained list is the
    same hardcoded-manifest failure mode this gate exists to remove -- adding a
    result type to `__all__` would either break every docs build or silently
    exempt itself, depending on whether someone remembered to update it.
    """
    sid = _sid()
    functions, others = set(), set()
    for name in sid.__all__:
        if name.startswith("__"):
            continue  # dunder metadata (e.g. __version__), not documentable
        (functions if inspect.isfunction(getattr(sid, name)) else others).add(name)
    return functions, others


def existing_pages() -> set[str]:
    return {p.stem for p in API_DIR.glob("*.md") if p.stem not in NON_FUNCTION_PAGES}


def check() -> list[str]:
    """Return a list of problems; empty means the API reference is complete."""
    expected, others = partition_exports()
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

    # Non-function exports (result dataclasses, SidError) share results.md.
    if RESULTS_PAGE.exists():
        results = RESULTS_PAGE.read_text(encoding="utf-8")
        for name in sorted(others):
            if name not in results:
                problems.append(
                    f"sid.{name} is exported in __all__ but is not documented in "
                    f"docsite/api/python/results.md"
                )
    else:
        problems.append(f"missing results page: {RESULTS_PAGE}")

    # Hand-written index tables are manifests too: a stub can exist and be in the
    # nav yet still be missing from these, and a missing table row is invisible to
    # the strict build.
    #
    # Match each page's link form, not a bare name: the two indexes link with
    # different prefixes, and a bare substring test is satisfied by a longer
    # sibling -- a dropped `spectrogram` row would go unnoticed because
    # `spectrogram_plot` contains it.
    for page, link in ((FUNCTION_INDEX, "(python/{name}.md)"), (PYTHON_OVERVIEW, "({name}.md)")):
        if not page.exists():
            problems.append(f"missing index page: {page}")
            continue
        text = page.read_text(encoding="utf-8")
        for name in sorted(expected):
            if link.format(name=name) not in text:
                problems.append(
                    f"sid.{name} is missing from "
                    f"{page.relative_to(REPO_ROOT).as_posix()} "
                    f"(expected a link {link.format(name=name)})"
                )

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
