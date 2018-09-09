# -*- coding: utf-8 -*-

import sphinx_bootstrap_theme


# -- Project information -----------------------------------------------------

project = 'MPSoC Symmetry Reduction Library'
copyright = '2018, Timo Nicolai'
author = 'Timo Nicolai'

version = '1.0'
release = '1.0.0 alpha'


# -- General configuration ---------------------------------------------------

master_doc = 'index'
source_suffix = '.rst'
exclude_patterns = []

language = None
pygments_style = 'sphinx'

templates_path = ['_templates']

extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages',
]

mathjax_path = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js'


# -- Options for HTML output -------------------------------------------------

html_theme = 'bootstrap'
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()
html_static_path = ['_static']
htmlhelp_basename = 'MPSoCSymmetryReductionLibrarydoc'
