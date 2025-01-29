# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'celltag_tools'
copyright = '2025, Kunal Jindal'
author = 'Kunal Jindal'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
             'nbsphinx',
             'myst_nb']

templates_path = ['_templates']
exclude_patterns = []
autodoc_mock_imports = ["igraph"]

nitpicky = True


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
autodoc_member_order = "bysource"
master_doc = 'index'


import os
import sys
sys.path.insert(0, os.path.abspath('../'))
html_theme = 'sphinx_rtd_theme'

#disable tutorial NB execution
nbsphinx_execute = 'never'


