"""Copy example notebooks + util_msd.py into the virtual docs tree.

`mkdocs-jupyter` renders notebooks in place; for `import util_msd` to resolve
during build-time execution, the helper must sit next to the notebooks. We
write both as virtual files via `mkdocs-gen-files` — no on-disk copies, no
commits.

Also emits `examples/python/index.md` (catalog) and `examples/python/SUMMARY.md`
(literate-nav listing).
"""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
NB_DIR = REPO_ROOT / "python" / "examples"

NOTEBOOK_TITLES: dict[str, str] = {
    "example_siso": "SISO frequency response (1-DoF SMD)",
    "example_etfe": "Empirical transfer function estimate",
    "example_freq_dep_res": "Frequency-dependent resolution (3-mass chain)",
    "example_coherence": "Coherence with colored disturbance",
    "example_method_comparison": "BT vs BTFDR vs ETFE",
    "example_mimo": "MIMO frequency response (2-mass chain)",
    "example_freq_map": "Time-varying frequency map",
    "example_spectrogram": "Spectrogram of a chirp-driven SDOF",
    "example_ltv_disc": "LTV state-space identification (COSMIC)",
    "example_multi_trajectory": "Multi-trajectory ensemble averaging",
    "example_output_cosmic": "Output-only COSMIC (partial observations)",
}


def discover_notebooks() -> list[Path]:
    return sorted(NB_DIR.glob("example_*.ipynb"))


def render_index(notebooks: list[Path]) -> str:
    out = [
        "# Python examples",
        "",
        "Jupyter notebooks built on the spring-mass-damper plants defined in",
        "[`spec/EXAMPLES.md`](../../spec/examples-spec.md). Every notebook is",
        "executed at site build time; the outputs you see below are produced",
        "by the current `sid` codebase.",
        "",
        "Each page also carries a Binder badge — click to launch a runnable",
        "copy in your browser.",
        "",
        "| Notebook | Description |",
        "|---|---|",
    ]
    for path in notebooks:
        title = NOTEBOOK_TITLES.get(path.stem, path.stem.replace("_", " "))
        out.append(f"| [`{path.stem}`]({path.name}) | {title} |")
    return "\n".join(out) + "\n"


def render_summary(notebooks: list[Path]) -> str:
    lines = ["* [Overview](index.md)"]
    for path in notebooks:
        title = NOTEBOOK_TITLES.get(path.stem, path.stem)
        lines.append(f"* [{title}]({path.name})")
    return "\n".join(lines) + "\n"


def main_genfiles() -> None:
    import mkdocs_gen_files

    notebooks = discover_notebooks()

    for path in notebooks:
        with mkdocs_gen_files.open(f"examples/python/{path.name}", "wb") as fp:
            fp.write(path.read_bytes())

    util_path = NB_DIR / "util_msd.py"
    if util_path.exists():
        with mkdocs_gen_files.open("examples/python/util_msd.py", "wb") as fp:
            fp.write(util_path.read_bytes())

    with mkdocs_gen_files.open("examples/python/index.md", "w") as fp:
        fp.write(render_index(notebooks))
    with mkdocs_gen_files.open("examples/python/SUMMARY.md", "w") as fp:
        fp.write(render_summary(notebooks))


if __name__ != "__main__":
    main_genfiles()
