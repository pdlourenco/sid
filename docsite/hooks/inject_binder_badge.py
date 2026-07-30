"""Prepend a Binder launch badge to each rendered example notebook page.

`mkdocs-jupyter` converts notebooks straight to HTML and bypasses the
markdown-processing pipeline, so we hook into `on_page_content` (which
operates on the rendered HTML) rather than `on_page_markdown`.

The badge URL targets `python/examples/<name>.ipynb` on `main` — i.e.
the notebook's canonical source location, not the docs-side copy.
"""

from __future__ import annotations

# Fallback only. `_binder_base` prefers mkdocs.yml's `repo_url` so the repo
# location is declared once, in the site config.
DEFAULT_BINDER_BASE = "https://mybinder.org/v2/gh/pdlourenco/sid/main"


def _binder_base(config) -> str:
    """`mybinder.org/v2/gh/<owner>/<repo>/main`, derived from `repo_url`."""
    repo_url = (config or {}).get("repo_url") or ""
    prefix = "https://github.com/"
    if not repo_url.startswith(prefix):
        return DEFAULT_BINDER_BASE
    slug = repo_url[len(prefix):].strip("/").removesuffix(".git")
    if slug.count("/") != 1:
        return DEFAULT_BINDER_BASE
    return f"https://mybinder.org/v2/gh/{slug}/main"


def on_page_content(html, page, config, files):
    src = page.file.src_path.replace("\\", "/")
    if not (src.startswith("examples/python/example_") and src.endswith(".ipynb")):
        return html
    nb_name = src.rsplit("/", 1)[-1]
    badge = (
        f'<p><a href="{_binder_base(config)}?labpath=python%2Fexamples%2F{nb_name}" '
        f'target="_blank" rel="noopener">'
        f'<img alt="Launch on Binder" src="https://mybinder.org/badge_logo.svg">'
        f'</a></p>\n'
    )
    return badge + html
