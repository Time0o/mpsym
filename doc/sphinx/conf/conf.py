# -*- coding: utf-8 -*-


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

extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages',
]

mathjax_path = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js'


# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinxdoc'
htmlhelp_basename = 'MPSoCSymmetryReductionLibrarydoc'
