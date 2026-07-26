# docs/conf.py
import os
import sys

# 1. Project Information
project = "vomix-snakemake"
copyright = "2026, Ho Lab, HKU"
author = "Erfan Shekarriz, HKU"
version = "local-dev"
release = "local-dev"

# Check if building on Read the Docs
if os.environ.get("READTHEDOCS") == "True":
    rtd_version = os.environ.get("READTHEDOCS_VERSION")
    if rtd_version:
        version = rtd_version
        release = rtd_version


# 2. Extensions Setup
# These modules allow Sphinx to read Markdown and format code blocks properly
extensions = [
    "myst_parser",  # Enables Markdown (.md) support
    "sphinx.ext.autodoc",  # Core documentation generator
    "sphinx.ext.viewcode",  # Adds links to source code
    "sphinx_design",  # Allows code tabs and other pretty things
]

# 3. Theme Customization
# Sets the classic, mobile-friendly ReadTheDocs sidebar layout
html_theme = "furo"
pygments_style = "tango"
pygments_dark_style = "tango"
html_logo = "./_static/logo_vomix-snakemake.svg"
html_theme_options = {
    "sidebar_hide_name": True,  # Optional: keeps your project title visible
    "navigation_with_keys": True,
    "collapse_navigation": True,
    "navigation_depth": 3,
}

# Add this so MyST automatically generates tracking anchors for markdown headers (up to H3)
myst_heading_anchors = 2


# 4. Markdown Settings
# Configures MyST to support advanced GitHub Markdown features like tables and colon blocks
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
]

# 5. Build Exclusions
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
