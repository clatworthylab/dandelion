# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from sphinx.deprecation import RemovedInSphinx40Warning

with open("README.md", "r") as readme_file:
    readme = readme_file.read()


sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath('.'))

import dandelion

# -- Project information -----------------------------------------------------

project = 'dandelion'
copyright = '2020, zktuong'
author = 'zktuong'

# The full version, including alpha/beta/rc tags
# release = __version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    "sphinx.ext.intersphinx",
    'sphinx.ext.autosummary',
    "sphinx_autodoc_typehints",
    'sphinx_rtd_theme',
    'nbsphinx',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'recommonmark'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', "**.ipynb_checkpoints"]

nitpicky = True  # Warn about broken links
needs_sphinx = "2.0"  # Nicer param docs
nitpick_ignore = [('py:class', 'type')]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_logo = "notebooks/img/dandelion_logo.png"
html_favicon = "notebooks/img/dandelion_logo.png"

html_theme_options = {'logo_only': True}

master_doc = 'index'

napoleon_use_param = False
autodoc_member_order = 'bysource'
autosummary_generate = True