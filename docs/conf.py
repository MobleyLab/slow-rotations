# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'slow-rotations'
copyright = '2025, Meghan Osato, Travis Dabbous'
author = 'Meghan Osato, Travis Dabbous'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

master_doc = "index"

extensions = [
    "myst_parser",            # Markdown support
    "sphinx.ext.autodoc",     # Optional, for docstrings
    "sphinx_autodoc_typehints",
    "sphinx.ext.napoleon"
]

# Tell Sphinx which file types to read
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

templates_path = ['_templates']
exclude_patterns = []

import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

autodoc_mock_imports = [
    "MDAnalysisTests",
    "openff_toolkit",
    "openff",
    "rdkit",
    "MDAnalysis",
    "parmed",
]



