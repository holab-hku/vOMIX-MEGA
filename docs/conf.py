# docs/conf.py
import os
import sys

# 1. Project Information
project = 'vOMIX-snakemake'
copyright = '2026, Ho Lab, HKU'
author = 'Ho Lab, HKU'
release = '1.0.0'

# 2. Extensions Setup
# These modules allow Sphinx to read Markdown and format code blocks properly
extensions = [
    'myst_parser',          # Enables Markdown (.md) support
    'sphinx.ext.autodoc',   # Core documentation generator
    'sphinx.ext.viewcode'   # Adds links to source code
]

# 3. Theme Customization
# Sets the classic, mobile-friendly ReadTheDocs sidebar layout
html_theme = 'sphinx_rtd_theme'

# 4. Markdown Settings
# Configures MyST to support advanced GitHub Markdown features like tables and colon blocks
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
]

# 5. Build Exclusions
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
