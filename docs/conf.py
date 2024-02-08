from pathlib import Path

HERE = Path(__file__).parents[1].resolve()

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "tidytcells"
copyright = "2022, Yuta Nagano"
author = "Yuta Nagano"
release = (HERE / "VERSION.txt").read_text(encoding="utf-8")

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
autosummary_imported_members = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_title = "tidytcells"
html_theme_options = {
    "repository_url": "https://github.com/yutanagano/tidytcells",
    "path_to_docs": "docs",
    "use_repository_button": True,
    "use_issues_button": True
}
# html_static_path = ['_static']
