# Configuration file for the Sphinx documentation builder.
# Build doc locally with
# sphinx-build -anE -c . src _build

def parse_version():
    from pathlib import Path
    import re
    init_path = Path(__file__).parent.parent / 'anchorna/__init__.py'
    with open(init_path) as f:
        content = f.read()
    regex = r"""__version__\s*=\s*['"]([^\s]+)['"]"""
    match = re.search(regex, content)
    if match:
        return match.group(1)


def parse_imports():
    from pathlib import Path
    import re
    root_path = Path(__file__).parent.parent / 'anchorna/'
    regex = r'^\s*import\s+(\w+)(?:.\w+)*|^\s*from\s+(\w+)(?:.\w+)*\s+import'
    modules = set()
    for fname in root_path.glob('**/*.py'):
        with open(fname) as f:
            content = f.read()
        for match in re.finditer(regex, content, flags=re.MULTILINE):
            modules |= set(match.groups())
    modules -= {None, 'anchorna'}
    return sorted(modules)


project = 'AnchoRNA'
copyright = '2024, Tom Eulenfeld'
author = 'Tom Eulenfeld'
version = parse_version()

default_role = 'py:obj'
templates_path = ['_templates']
exclude_patterns = ['_build']

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
autodoc_mock_imports = parse_imports()
print(f'set {autodoc_mock_imports=}')

html_title = f'AnchoRNA <br>v{version} <br>documentation'
html_theme = 'furo'
html_theme_options = {
    'footer_icons' : [],
    'top_of_page_buttons': [],
}

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'sugar': ('https://rnajena-sugar.readthedocs.io/en/latest', None)
    }
