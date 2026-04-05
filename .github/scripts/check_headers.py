#!/usr/bin/env python3
"""Validate MATLAB function headers against the CONTRIBUTING.md standard.

Checks library functions in the repo root and internal/ directory.
Skips examples/ and tests/ (scripts, not functions).

Exit codes:
    0 — all checks pass
    1 — one or more hard errors found
"""

import glob
import os
import re
import sys

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Directories to scan (relative to repo root)
SCAN_DIRS = ["matlab/sid", "matlab/sid/internal"]

# Files exempt from SPECIFICATION check
SPEC_EXEMPT = set()

# Files exempt from INPUTS/OUTPUTS (no arguments or no return value)
INPUTS_EXEMPT = set()
OUTPUTS_EXEMPT = set()

# Canonical section order — entries that may appear in the header.
# Each tuple: (regex pattern, label, is_required)
# "required" means the section must appear in every function header
# (subject to per-file exemptions above).
SECTION_DEFS = [
    (r"^%\s+INPUTS(?:\s+\([^)]*\))?:",  "INPUTS",            True),
    (r"^%\s+NAME-VALUE OPTIONS:",        "NAME-VALUE OPTIONS", False),
    (r"^%\s+OUTPUTS:",                   "OUTPUTS",           True),
    (r"^%\s+EXAMPLES:",                  "EXAMPLES",          True),
    (r"^%\s+ALGORITHM:",                 "ALGORITHM",         False),
    (r"^%\s+REFERENCES:",                "REFERENCES",        False),
    (r"^%\s+SPECIFICATION:",             "SPECIFICATION",     False),
    (r"^%\s+See also:",                  "See also",          True),
    (r"^%\s+Changelog:",                 "Changelog",         True),
]

# Singular headings that violate the plural rule
SINGULAR_PATTERNS = [
    (r"^%\s+INPUT:",     "INPUT:",     "INPUTS:"),
    (r"^%\s+OUTPUT:",    "OUTPUT:",    "OUTPUTS:"),
    (r"^%\s+Example:",   "Example:",   "EXAMPLES:"),
    (r"^%\s+REFERENCE:", "REFERENCE:", "REFERENCES:"),
]

COPYRIGHT_MARKER = re.compile(r"^%\s*-{10,}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_library_m_files(repo_root):
    """Return sorted list of .m files in SCAN_DIRS (skip examples/tests)."""
    files = []
    for d in SCAN_DIRS:
        pattern = os.path.join(repo_root, d, "*.m")
        files.extend(glob.glob(pattern))
    return sorted(files)


def extract_header(filepath):
    """Return (lines, header_end_line) for the header comment block.

    The header is the contiguous block of '%' comment lines after the
    function signature line.  Blank lines within the comment block are
    included.
    """
    lines = []
    header_lines = []
    in_header = False
    header_end = 0

    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    in_signature = False
    for i, line in enumerate(lines):
        stripped = line.rstrip("\n")
        if not in_header:
            # Skip leading function signature line(s), including continuations
            if stripped.startswith("function ") or stripped.startswith("function\t"):
                in_signature = True
                if not stripped.rstrip().endswith("..."):
                    in_signature = False
                continue
            if in_signature:
                # Continuation of function signature
                if not stripped.rstrip().endswith("..."):
                    in_signature = False
                continue
            if stripped.startswith("%"):
                in_header = True
                header_lines.append((i + 1, stripped))  # 1-based line number
                header_end = i + 1
            else:
                break
        else:
            if stripped.startswith("%") or stripped.strip() == "":
                header_lines.append((i + 1, stripped))
                if stripped.startswith("%"):
                    header_end = i + 1
            else:
                break

    return header_lines, header_end


def check_file(filepath, repo_root):
    """Check a single file. Return list of (level, line, message)."""
    errors = []   # (line, message) — hard errors
    warnings = []  # (line, message) — soft warnings
    relpath = os.path.relpath(filepath, repo_root)
    basename = os.path.basename(filepath)

    header, _ = extract_header(filepath)
    if not header:
        errors.append((1, "No header comment block found"))
        return errors, warnings

    header_text = [line for _, line in header]

    # --- Rule 1: No singular headings ---
    for lineno, line in header:
        for pattern, singular, plural in SINGULAR_PATTERNS:
            if re.match(pattern, line):
                errors.append(
                    (lineno, f"Singular heading '{singular}' — use '{plural}'")
                )

    # --- Rule 2: Required sections present ---
    # Find copyright separator position
    copyright_line = None
    for lineno, line in header:
        if COPYRIGHT_MARKER.match(line):
            copyright_line = lineno
            break

    section_positions = {}  # label -> first line number
    for lineno, line in header:
        for pattern, label, _ in SECTION_DEFS:
            if re.match(pattern, line):
                if label not in section_positions:
                    section_positions[label] = lineno

    for pattern, label, required in SECTION_DEFS:
        if required and label not in section_positions:
            # Per-file exemptions
            if label == "INPUTS" and basename in INPUTS_EXEMPT:
                continue
            if label == "OUTPUTS" and basename in OUTPUTS_EXEMPT:
                continue
            errors.append((1, f"Missing required section: {label}"))

    # Check Copyright block
    has_copyright = any("Copyright" in line for _, line in header)
    if not has_copyright:
        errors.append((1, "Missing Copyright block"))

    has_changelog = "Changelog" in section_positions
    if not has_changelog:
        errors.append((1, "Missing Changelog section"))

    # --- Rule 3: Exactly one See also ---
    see_also_count = sum(
        1 for _, line in header if re.match(r"^%\s+See also:", line)
    )
    if see_also_count > 1:
        errors.append(
            (1, f"Duplicate 'See also:' — found {see_also_count}, expected 1")
        )

    # --- Rule 4: No documentation after copyright separator ---
    if copyright_line is not None:
        doc_sections_after = []
        for lineno, line in header:
            if lineno <= copyright_line:
                continue
            for pattern, label, _ in SECTION_DEFS:
                if re.match(pattern, line):
                    doc_sections_after.append((lineno, label))
        for lineno, label in doc_sections_after:
            errors.append(
                (lineno, f"'{label}' appears after copyright separator (line {copyright_line})")
            )

    # --- Rule 5: Section order ---
    canonical_order = [label for _, label, _ in SECTION_DEFS]
    found_order = []
    for label in canonical_order:
        if label in section_positions:
            found_order.append((section_positions[label], label))
    found_order.sort(key=lambda x: x[0])

    for i in range(len(found_order) - 1):
        lineno_a, label_a = found_order[i]
        lineno_b, label_b = found_order[i + 1]
        idx_a = canonical_order.index(label_a)
        idx_b = canonical_order.index(label_b)
        if idx_a > idx_b:
            errors.append(
                (lineno_b,
                 f"Section order violation: '{label_b}' (line {lineno_b}) "
                 f"appears before '{label_a}' (line {lineno_a})")
            )

    # --- Rule 6 (warning): Missing SPECIFICATION ---
    if basename not in SPEC_EXEMPT and "SPECIFICATION" not in section_positions:
        warnings.append((1, "Missing SPECIFICATION section (SPEC.md reference)"))

    return errors, warnings


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # Determine repo root (script lives in .github/scripts/)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.normpath(os.path.join(script_dir, "..", ".."))

    files = find_library_m_files(repo_root)
    if not files:
        print("No .m files found — nothing to check.")
        return 0

    total_errors = 0
    total_warnings = 0

    for filepath in files:
        relpath = os.path.relpath(filepath, repo_root)
        file_errors, file_warnings = check_file(filepath, repo_root)

        for lineno, msg in file_errors:
            print(f"::error file={relpath},line={lineno}::{msg}")
            total_errors += 1
        for lineno, msg in file_warnings:
            print(f"::warning file={relpath},line={lineno}::{msg}")
            total_warnings += 1

    print()
    print(f"Header check: {len(files)} file(s) analysed, "
          f"{total_errors} error(s), {total_warnings} warning(s)")

    return 1 if total_errors > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
