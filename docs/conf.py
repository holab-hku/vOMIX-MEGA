# docs/conf.py
import os
import sys

# 1. Project Information
project = 'vOMIX-snakemake'
copyright = '2026, Ho Lab, HKU'
author = 'Erfan Shekarriz, HKU'
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
html_theme = "furo"
pygments_style = "tango"
pygments_dark_style = "tango"
html_theme_options = {
    "sidebar_hide_name": False,  # Optional: keeps your project title visible
    "navigation_with_keys": True,
}

# Add this so MyST automatically generates tracking anchors for markdown headers (up to H3)
myst_heading_anchors = 3


# 4. Markdown Settings
# Configures MyST to support advanced GitHub Markdown features like tables and colon blocks
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
]

# 5. Build Exclusions
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
