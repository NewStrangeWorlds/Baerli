# Sphinx configuration for the Bärli documentation.
#
# Build (from the repository root):
#     sphinx-build -b html docs_src docs_test
# or `make html` inside docs_src.

import datetime

project = "Bärli"
author = "Daniel Kitzmann"
copyright = f"2024-{datetime.date.today().year}, {author}"
release = "dev"
version = "dev"

extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx.ext.intersphinx",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The pages are written in reStructuredText only.
source_suffix = {".rst": "restructuredtext"}
master_doc = "index"
language = "en"

# cross references to the Python standard library and the scientific stack
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
}
intersphinx_timeout = 5

todo_include_todos = False

# default role so that `name` renders as inline code
default_role = "code"
highlight_language = "python"
pygments_style = "friendly"
pygments_dark_style = "monokai"

numfig = True
math_numfig = True
math_eqref_format = "Eq. {number}"
numfig_format = {"figure": "Fig. %s", "table": "Table %s", "code-block": "Listing %s"}

# MathJax: a few macros used throughout the theory chapters
mathjax3_config = {
    "tex": {
        "macros": {
            "dd": r"\mathrm{d}",
            "vec": [r"\boldsymbol{#1}", 1],
            "Fnet": r"F_{\mathrm{net}}",
            "Fint": r"F_{\mathrm{int}}",
            "Teff": r"T_{\mathrm{eff}}",
            "Tint": r"T_{\mathrm{int}}",
            "micron": r"\mu\mathrm{m}",
        },
        "tags": "ams",
    }
}

# -- HTML output ---------------------------------------------------------------------------------

html_theme = "furo"
html_title = "Bärli documentation"
html_short_title = "Bärli"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_show_sourcelink = False
html_copy_source = False


html_theme_options = {
    "light_logo": "baerli-badge-light.png",
    "dark_logo": "baerli-badge-dark.png",
    "light_css_variables": {
        "color-brand-primary": "#8b4513",
        "color-brand-content": "#a0522d",
        "font-stack--monospace": "'JetBrains Mono', 'Fira Mono', Menlo, Consolas, monospace",
    },
    "dark_css_variables": {
        "color-brand-primary": "#e6a86b",
        "color-brand-content": "#f0b884",
    },
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
    "top_of_page_buttons": [],
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/NewStrangeWorlds",
            "html": "",
            "class": "fa-brands fa-github",
        },
    ],
}
