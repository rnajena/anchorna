# Configuration file for the Sphinx documentation builder.
# Build doc locally with with
# sphinx-build -anE -c . src _build

import sys
import os

os.environ['SPHINX_BUILD'] = '1'


def version():
    from pathlib import Path
    import re
    init_path = Path(__file__).parent.parent / 'anchorna/__init__.py'
    with open(init_path) as f:
        content = f.read()
    match = re.search("""__version__\s*=\s*['"]([^\s]+)['"]""", content)
    if match:
        return match.group(1)


project = 'AnchoRNA'
copyright = '2024, Tom Eulenfeld'
author = 'Tom Eulenfeld'
version = version()

### General configuration
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
default_role = 'py:obj'
templates_path = ['_templates']
exclude_patterns = ['_build']

# Add any Sphinx extension module names here
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx.ext.intersphinx',
              'sphinx.ext.viewcode',
              ]

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'private-members': True,
    'show-inheritance': True,
}

autodoc_mock_imports = ['sugar', 'tqdm']

html_theme = 'furo'
html_title = f'AnchoRNA <br>v{version} <br>documentation'
html_theme_options = {
    'footer_icons' : [],
    'top_of_page_buttons': [],
}

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'sugar': ('https://rnajena-sugar.readthedocs.io/en/latest', None)
    }
