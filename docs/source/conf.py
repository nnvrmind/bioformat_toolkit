import os
import sys

sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
)

project = "BioFormat Toolkit"
author = "Команда BioFormat Toolkit"
language = "ru"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

autodoc_typehints = "description"
autodoc_member_order = "bysource"

html_theme = "alabaster"

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
