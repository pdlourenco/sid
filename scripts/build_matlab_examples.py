"""Generate MkDocs example pages from MATLAB example scripts.

Each `matlab/examples/example*.m` becomes a page that shows the leading
comment block as prose and the full script source in a `matlab` code fence.
Utility plant files (`util_msd*.m`) and the runner are skipped.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = REPO_ROOT / "matlab" / "examples"

# Fallback only. `repo_source_prefix` prefers mkdocs.yml's `repo_url` so the repo
# location is declared once, in the site config.
DEFAULT_REPO_PREFIX = "https://github.com/pdlourenco/sid/blob/main"


def repo_source_prefix() -> str:
    """`<repo_url>/blob/main`, from the mkdocs config when reachable."""
    try:
        import mkdocs_gen_files

        repo_url = (mkdocs_gen_files.config or {}).get("repo_url")
    except Exception:  # noqa: BLE001 - best-effort lookup; the fallback is correct
        # Deliberately broad. Outside a mkdocs build `mkdocs_gen_files.config`
        # reaches for the active build context and loads mkdocs.yml, which raises
        # ConfigurationError when the cwd is not the project root -- not just
        # ImportError/AttributeError. This script must stay runnable standalone
        # from any directory, so any failure here falls back to the constant.
        repo_url = None
    if not repo_url:
        return DEFAULT_REPO_PREFIX
    return str(repo_url).rstrip("/") + "/blob/main"

SKIP_PREFIXES = ("util_", "runAllExamples", "example_template")


def extract_intro(text: str) -> tuple[str, str]:
    """Split the leading comment block (intro) from the rest of the source."""
    lines = text.splitlines()
    intro_lines: list[str] = []
    i = 0
    while i < len(lines):
        stripped = lines[i].lstrip()
        if not stripped:
            if intro_lines:
                break
            i += 1
            continue
        if stripped.startswith("%"):
            body = stripped[1:].lstrip("%").lstrip()
            intro_lines.append(body.rstrip())
            i += 1
        else:
            break
    title = ""
    if intro_lines and " - " in intro_lines[0]:
        title = intro_lines[0].split(" - ", 1)[1].strip()
        intro_lines = intro_lines[1:]
    intro = "\n".join(intro_lines).strip()
    return title, intro


def discover_examples() -> list[Path]:
    paths = sorted(EXAMPLES_DIR.glob("*.m"))
    return [p for p in paths if not p.stem.startswith(SKIP_PREFIXES)]


def render_page(path: Path, repo_prefix: str | None = None) -> str:
    repo_prefix = repo_prefix or DEFAULT_REPO_PREFIX
    text = path.read_text(encoding="utf-8", errors="replace")
    title, intro = extract_intro(text)
    heading = title or path.stem
    parts = [f"# {heading}", ""]
    if intro:
        parts.extend([intro, ""])
    src_url = f"{repo_prefix}/matlab/examples/{path.name}"
    parts.append(f"Source: [`matlab/examples/{path.name}`]({src_url})")
    parts.append("")
    parts.append("```matlab")
    parts.append(text.rstrip())
    parts.append("```")
    return "\n".join(parts) + "\n"


def render_index(catalog: list[tuple[str, str, str]]) -> str:
    out = [
        "# MATLAB examples",
        "",
        "Self-contained scripts that exercise the public API on physical",
        "spring-mass-damper plants. Run any example from a MATLAB or Octave",
        "session after adding `matlab/sid/` to the path.",
        "",
        "| Example | Description |",
        "|---|---|",
    ]
    for stem, heading, intro_first_line in catalog:
        desc = heading if heading != stem else (intro_first_line or "—")
        out.append(f"| [`{stem}`]({stem}.md) | {desc} |")
    return "\n".join(out) + "\n"


def render_summary(catalog: list[tuple[str, str, str]]) -> str:
    lines = ["* [Overview](index.md)"]
    for stem, _, _ in catalog:
        lines.append(f"* [{stem}]({stem}.md)")
    return "\n".join(lines) + "\n"


def first_sentence(intro: str) -> str:
    if not intro:
        return ""
    one_line = re.sub(r"\s+", " ", intro).strip()
    match = re.match(r"^(.+?[.!?])(?:\s|$)", one_line)
    return match.group(1) if match else one_line[:120]


def _build(write):
    paths = discover_examples()
    repo_prefix = repo_source_prefix()
    catalog: list[tuple[str, str, str]] = []
    for path in paths:
        text = path.read_text(encoding="utf-8", errors="replace")
        title, intro = extract_intro(text)
        write(f"examples/matlab/{path.stem}.md", render_page(path, repo_prefix))
        catalog.append((path.stem, title, first_sentence(intro)))
    write("examples/matlab/index.md", render_index(catalog))
    write("examples/matlab/SUMMARY.md", render_summary(catalog))
    return len(paths)


def main_genfiles() -> None:
    import mkdocs_gen_files

    def write(rel_path: str, content: str) -> None:
        with mkdocs_gen_files.open(rel_path, "w") as fp:
            fp.write(content)

    _build(write)


def main_standalone(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    def write(rel_path: str, content: str) -> None:
        full = out_dir / Path(rel_path).name
        full.write_text(content, encoding="utf-8")

    n = _build(write)
    print(f"Wrote {n} example pages + index + SUMMARY to {out_dir}")


if __name__ == "__main__":
    import sys

    out = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO_ROOT / "_matlab_examples_out"
    main_standalone(out)
else:
    main_genfiles()
