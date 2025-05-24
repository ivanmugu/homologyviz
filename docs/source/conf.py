from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "HomologyViz"
copyright = "2025, Iván Muñoz Gutiérrez"
author = "Iván Muñoz Gutiérrez"
release = "0.1.12"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # for Numpy/Google-style docstrings
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",  # Optional: for autosummary tables
    "sphinx_autodoc_typehints",  # If installed, for cleaner type hint formatting
]

autosummary_generate = True

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = "alabaster"
html_theme = "furo"
html_static_path = ["_static"]
html_title = "HomologyViz Documentation"
html_logo = "_static/logo.png"  # Make sure this file exists
# html_favicon = "_static/favicon.ico"

toc_object_entries_show_parents = "hide"

language = "en"
