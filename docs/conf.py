# conf.py
import os
import sys
sys.path.insert(0, os.path.abspath('.'))

# -- Project information -----------------------------------------------------
project = 'vOMIX-snakemake'
copyright = '2026, Ho Lab, HKU'
author = 'Ho Lab, HKU'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# Enable MyST for Markdown processing and standard Sphinx extensions
extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode'
]

# Paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# Use the official ReadTheDocs responsive theme
html_theme = 'sphinx_rtd_theme'
html_static_path = []

# Configure MyST Parser to handle common Markdown flavors (like tables)
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
]
