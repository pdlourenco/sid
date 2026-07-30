"""Generate MkDocs API reference pages from MATLAB H1 headers.

Runs as a `mkdocs-gen-files` script: walks every public `sid*.m` file under
`matlab/sid/`, parses the leading comment block, and writes one Markdown
page per function into the virtual docs tree at `api/matlab/`.

Also generates `api/matlab/index.md` (alphabetical catalog) and
`api/matlab/SUMMARY.md` (literate-nav listing) so the nav populates
automatically.

The script can be executed directly for standalone debugging:
    python scripts/build_matlab_api.py /tmp/out
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = REPO_ROOT / "matlab" / "sid"

# Section keywords recognised in the H1 block. Order matters only for output
# ordering when we emit; section detection itself is a set lookup.
SECTIONS: dict[str, str] = {
    "INPUTS": "Inputs",
    "NAME-VALUE OPTIONS": "Name-value options",
    "OUTPUTS": "Outputs",
    "EXAMPLES": "Examples",
    "ALGORITHM": "Algorithm",
    "FIT METRIC": "Fit metric",
    "FREQUENCY METHOD OPTIONS": "Frequency method options",
    "VALIDATION METHOD": "Validation method",
    "REFERENCES": "References",
    "SPECIFICATION": "Specification",
    "SEE ALSO": "See also",
    "CHANGELOG": "Changelog",
}

SECTION_HEAD_RE = re.compile(
    r"^\s*(?P<name>[A-Z][A-Z0-9 \-/]+?):\s*(?P<rest>.*)$"
)

# Conservative math whitelist — wrap only well-known forms. Single combined
# regex with longest alternatives first prevents nested re-wrapping.
MATH_RE = re.compile(
    r"G\(e\^\{[^}]+\}\)"
    r"|e\^\{-?j[^}]*\}"
    r"|[A-Z][A-Za-z]*\^\{-?\w+\}"
    r"|Phi_(?:yu|uy|y|u|v)(?:\([^)]+\))?"
    r"|x\(k\+1\)"
    r"|G\(w(?:_k|,\s*k)?\)"
)

# Functions to expose. Filename glob is `sid*.m` minus this denylist.
DENYLIST = {"sidInstall"}


def strip_comment(line: str) -> str:
    """Strip a single leading '%' plus at most one space, preserving indent."""
    if not line.startswith("%"):
        return line
    body = line[1:]
    if body.startswith(" "):
        body = body[1:]
    return body.rstrip()


SEPARATOR_RE = re.compile(r"^\s*-{5,}\s*$")


def extract_header_lines(text: str) -> list[str]:
    """Return the contiguous comment block immediately after `function ...`.

    Stops at the first horizontal separator (``-----`` line) so the copyright
    footer is not parsed as Changelog content.
    """
    lines = text.splitlines()
    if not lines:
        return []
    out: list[str] = []
    in_block = False
    for raw in lines[1:]:  # skip `function ...` declaration
        if raw.lstrip().startswith("%"):
            in_block = True
            stripped = strip_comment(raw.lstrip())
            if SEPARATOR_RE.match(stripped):
                break
            out.append(stripped)
        elif in_block:
            break
    while out and not out[0].strip():
        out.pop(0)
    while out and not out[-1].strip():
        out.pop()
    return out


def parse_header(header: list[str]) -> dict:
    """Split header into title, lede, and named sections."""
    title_line = header[0].strip() if header else ""
    title_match = re.match(r"^([A-Z][A-Z0-9]+)\s+(.*)$", title_line)
    title_brief = title_match.group(2) if title_match else title_line

    sections: dict[str, list[str]] = {"_lede": []}
    current = "_lede"

    for line in header[1:]:
        stripped = line.strip()
        m = SECTION_HEAD_RE.match(line) if stripped else None
        key = m.group("name").upper() if m else None
        if key in SECTIONS:
            current = key
            sections[current] = []
            rest = m.group("rest").strip()
            if rest:
                sections[current].append(rest)
            continue
        if stripped.lower().startswith("see also"):
            current = "SEE ALSO"
            sections[current] = []
            after = stripped[len("see also") :].lstrip(": ").strip()
            if after:
                sections[current].append(after)
            continue
        if stripped.lower().startswith("changelog"):
            current = "CHANGELOG"
            sections[current] = []
            after = stripped[len("changelog") :].lstrip(": ").strip()
            if after:
                sections[current].append(after)
            continue
        sections.setdefault(current, []).append(line)

    return {"brief": title_brief, "sections": sections}


def wrap_math(text: str) -> str:
    """Wrap whitelisted math fragments in `$...$`, outside code fences."""
    out_lines = []
    in_fence = False
    for line in text.splitlines():
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
            out_lines.append(line)
            continue
        if in_fence:
            out_lines.append(line)
            continue
        replaced = MATH_RE.sub(lambda m: f"${m.group(0)}$", line)
        out_lines.append(replaced)
    return "\n".join(out_lines)


def render_lede(lines: list[str], fn_name: str) -> str:
    """Render the pre-section block: call signatures + summary prose."""
    paragraphs: list[list[str]] = [[]]
    for line in lines:
        if not line.strip():
            if paragraphs[-1]:
                paragraphs.append([])
            continue
        paragraphs[-1].append(line)
    paragraphs = [p for p in paragraphs if p]

    out = []
    for para in paragraphs:
        if all(re.match(r"\s*(?:\[[^\]]+\]|\w+)\s*=\s*" + fn_name, p) or p.lstrip().startswith(fn_name) for p in para):
            out.append("```matlab\n" + "\n".join(p.strip() for p in para) + "\n```")
        else:
            out.append(" ".join(p.strip() for p in para))
    return "\n\n".join(out)


def render_io_section(lines: list[str], heading: str) -> str:
    """Render INPUTS / NAME-VALUE OPTIONS / OUTPUTS as a 2-column table."""
    rows: list[tuple[str, list[str]]] = []
    for line in lines:
        if not line.strip():
            continue
        m = re.match(r"^\s{0,6}([^\s\-][^-]*?)\s+-\s+(.*)$", line)
        if m and rows is not None:
            name = m.group(1).strip()
            desc = m.group(2).strip()
            rows.append((name, [desc]))
        elif rows:
            rows[-1][1].append(line.strip())
        else:
            rows.append(("", [line.strip()]))

    if not rows:
        return f"## {heading}\n"

    out = [f"## {heading}\n", "| Name | Description |", "|---|---|"]
    for name, desc_lines in rows:
        desc = "<br>".join(d for d in desc_lines if d)
        desc = desc.replace("|", "\\|")
        out.append(f"| `{name}` | {desc} |" if name else f"| | {desc} |")
    return "\n".join(out) + "\n"


def render_examples(lines: list[str]) -> str:
    """Render EXAMPLES: each blank-line group becomes a matlab fenced block."""
    groups: list[list[str]] = [[]]
    for line in lines:
        if not line.strip():
            if groups[-1]:
                groups.append([])
            continue
        groups[-1].append(line.rstrip())
    groups = [g for g in groups if g]

    out = ["## Examples\n"]
    for group in groups:
        common = min((len(line) - len(line.lstrip(" ")) for line in group if line.strip()), default=0)
        body = "\n".join(line[common:] if len(line) >= common else line for line in group)
        out.append("```matlab\n" + body + "\n```")
    return "\n\n".join(out) + "\n"


LIST_ITEM_RE = re.compile(r"^\s*(?:[-*]|\d+[.)])\s+")


def render_paragraph_section(lines: list[str], heading: str) -> str:
    """Render free-form prose, preserving numbered / bulleted list structure."""
    body_lines = [line.rstrip() for line in lines]
    while body_lines and not body_lines[0].strip():
        body_lines.pop(0)
    while body_lines and not body_lines[-1].strip():
        body_lines.pop()
    if not body_lines:
        return ""

    blocks: list[dict] = []
    for line in body_lines:
        stripped = line.strip()
        if not stripped:
            if blocks and blocks[-1]["kind"] != "break":
                blocks.append({"kind": "break"})
            continue
        if LIST_ITEM_RE.match(line):
            blocks.append({"kind": "item", "text": stripped})
        elif blocks and blocks[-1]["kind"] == "item":
            blocks[-1]["text"] += " " + stripped
        elif blocks and blocks[-1]["kind"] == "para":
            blocks[-1]["text"] += " " + stripped
        else:
            blocks.append({"kind": "para", "text": stripped})

    out: list[str] = [f"## {heading}", ""]
    prev = None
    for b in blocks:
        if b["kind"] == "break":
            if prev and prev != "break":
                out.append("")
            prev = "break"
            continue
        if b["kind"] == "item":
            out.append(b["text"])
        else:
            if prev == "item":
                out.append("")
            out.append(b["text"])
        prev = b["kind"]
    return "\n".join(out).rstrip() + "\n"


def render_see_also(lines: list[str], public_names: set[str]) -> str:
    """Render `See also` as inline links, skipping references to private helpers."""
    text = " ".join(line.strip() for line in lines if line.strip())
    if not text:
        return ""
    text = text.rstrip(".")
    parts = [p.strip() for p in re.split(r"[,;]", text) if p.strip()]
    links = []
    for p in parts:
        m = re.match(r"^(sid[A-Za-z0-9]+)", p)
        if m and m.group(1) in public_names:
            fn = m.group(1)
            links.append(f"[`{fn}`]({fn}.md)")
        elif m:
            links.append(f"`{m.group(1)}`")
        else:
            links.append(p)
    return "## See also\n\n" + ", ".join(links) + "\n"


DATE_ENTRY_RE = re.compile(r"^\d{4}-\d{2}-\d{2}\b")


def render_changelog(lines: list[str]) -> str:
    """Render Changelog as a bullet list. Each `YYYY-MM-DD:` line is a new item."""
    items: list[str] = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if DATE_ENTRY_RE.match(stripped) or stripped.startswith(("-", "*")):
            items.append(stripped.lstrip("-* ").strip())
        elif items:
            items[-1] += " " + stripped
        else:
            items.append(stripped)
    if not items:
        return ""
    return "## Changelog\n\n" + "\n".join(f"- {item}" for item in items) + "\n"


PYTHON_EQUIV: dict[str, str] = {
    "sidFreqBT": "freq_bt",
    "sidFreqBTFDR": "freq_btfdr",
    "sidFreqETFE": "freq_etfe",
    "sidFreqMap": "freq_map",
    "sidSpectrogram": "spectrogram",
    "sidLTVdisc": "ltv_disc",
    "sidLTVdiscIO": "ltv_disc_io",
    "sidLTVdiscTune": "ltv_disc_tune",
    "sidLTVdiscFrozen": "ltv_disc_frozen",
    "sidLTIfreqIO": "lti_freq_io",
    "sidLTVStateEst": "ltv_state_est",
    "sidModelOrder": "model_order",
    "sidDetrend": "detrend",
    "sidResidual": "residual",
    "sidCompare": "compare",
    "sidBodePlot": "bode_plot",
    "sidSpectrumPlot": "spectrum_plot",
    "sidMapPlot": "map_plot",
    "sidSpectrogramPlot": "spectrogram_plot",
}


def render_page(fn_name: str, parsed: dict, public_names: set[str]) -> str:
    brief = parsed["brief"]
    sections = parsed["sections"]

    parts: list[str] = [f"# {fn_name}\n"]
    py_equiv = PYTHON_EQUIV.get(fn_name)
    if py_equiv:
        parts.append(f"> Python equivalent: [`sid.{py_equiv}`](../python/{py_equiv}.md)\n")
    if brief:
        parts.append(brief + "\n")

    if sections.get("_lede"):
        parts.append(render_lede(sections["_lede"], fn_name))

    render_order = [
        ("INPUTS", "Inputs", render_io_section),
        ("NAME-VALUE OPTIONS", "Name-value options", render_io_section),
        ("FREQUENCY METHOD OPTIONS", "Frequency method options", render_io_section),
        ("VALIDATION METHOD", "Validation method", render_paragraph_section),
        ("OUTPUTS", "Outputs", render_io_section),
        ("EXAMPLES", None, render_examples),
        ("ALGORITHM", "Algorithm", render_paragraph_section),
        ("FIT METRIC", "Fit metric", render_paragraph_section),
        ("REFERENCES", "References", render_paragraph_section),
        ("SPECIFICATION", "Specification", render_paragraph_section),
        ("SEE ALSO", None, render_see_also),
        ("CHANGELOG", None, render_changelog),
    ]
    for key, heading, renderer in render_order:
        if key not in sections:
            continue
        if key == "SEE ALSO":
            rendered = renderer(sections[key], public_names)
        elif heading is None:
            rendered = renderer(sections[key])
        else:
            rendered = renderer(sections[key], heading)
        if rendered.strip():
            parts.append(rendered)

    body = "\n\n".join(p.rstrip() for p in parts if p.strip())
    return wrap_math(body) + "\n"


def discover_files() -> list[Path]:
    return sorted(
        p
        for p in MATLAB_DIR.glob("sid*.m")
        if p.stem not in DENYLIST
    )


def build_catalog(files: list[Path]) -> list[tuple[str, str]]:
    catalog: list[tuple[str, str]] = []
    for path in files:
        text = path.read_text(encoding="utf-8", errors="replace")
        header = extract_header_lines(text)
        parsed = parse_header(header)
        catalog.append((path.stem, parsed["brief"]))
    return catalog


def render_index(catalog: list[tuple[str, str]]) -> str:
    out = [
        "# MATLAB/Octave API",
        "",
        "Reference pages are auto-generated from each function's H1 header.",
        "See the [function index](../index.md) for the Python ↔ MATLAB pairing.",
        "",
        "| Function | Description |",
        "|---|---|",
    ]
    for name, brief in catalog:
        out.append(f"| [`{name}`]({name}.md) | {brief} |")
    return "\n".join(out) + "\n"


def render_summary(catalog: list[tuple[str, str]]) -> str:
    lines = ["* [Overview](index.md)"]
    for name, _ in catalog:
        lines.append(f"* [{name}]({name}.md)")
    return "\n".join(lines) + "\n"


def main_genfiles() -> None:
    import mkdocs_gen_files

    files = discover_files()
    public_names = {p.stem for p in files}
    catalog: list[tuple[str, str]] = []
    for path in files:
        text = path.read_text(encoding="utf-8", errors="replace")
        header = extract_header_lines(text)
        parsed = parse_header(header)
        page = render_page(path.stem, parsed, public_names)
        with mkdocs_gen_files.open(f"api/matlab/{path.stem}.md", "w") as fp:
            fp.write(page)
        catalog.append((path.stem, parsed["brief"]))
    with mkdocs_gen_files.open("api/matlab/index.md", "w") as fp:
        fp.write(render_index(catalog))
    with mkdocs_gen_files.open("api/matlab/SUMMARY.md", "w") as fp:
        fp.write(render_summary(catalog))


def main_standalone(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    files = discover_files()
    public_names = {p.stem for p in files}
    catalog: list[tuple[str, str]] = []
    for path in files:
        text = path.read_text(encoding="utf-8", errors="replace")
        header = extract_header_lines(text)
        parsed = parse_header(header)
        page = render_page(path.stem, parsed, public_names)
        (out_dir / f"{path.stem}.md").write_text(page, encoding="utf-8")
        catalog.append((path.stem, parsed["brief"]))
    (out_dir / "index.md").write_text(render_index(catalog), encoding="utf-8")
    (out_dir / "SUMMARY.md").write_text(render_summary(catalog), encoding="utf-8")
    print(f"Wrote {len(files)} function pages + index + SUMMARY to {out_dir}")


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO_ROOT / "_matlab_api_out"
    main_standalone(out)
else:
    main_genfiles()
