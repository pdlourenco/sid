"""Prepend a Binder launch badge to each rendered example notebook page.

`mkdocs-jupyter` converts notebooks straight to HTML and bypasses the
markdown-processing pipeline, so we hook into `on_page_content` (which
operates on the rendered HTML) rather than `on_page_markdown`.

The badge URL targets `python/examples/<name>.ipynb` on `main` — i.e.
the notebook's canonical source location, not the docs-side copy.
"""

from __future__ import annotations

BINDER_BASE = "https://mybinder.org/v2/gh/pdlourenco/sid/main"


def on_page_content(html, page, config, files):
    src = page.file.src_path.replace("\\", "/")
    if not (src.startswith("examples/python/example_") and src.endswith(".ipynb")):
        return html
    nb_name = src.rsplit("/", 1)[-1]
    badge = (
        f'<p><a href="{BINDER_BASE}?labpath=python%2Fexamples%2F{nb_name}" '
        f'target="_blank" rel="noopener">'
        f'<img alt="Launch on Binder" src="https://mybinder.org/badge_logo.svg">'
        f'</a></p>\n'
    )
    return badge + html
