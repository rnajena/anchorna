# (C) 2024, Tom Eulenfeld, MIT license

from setuptools import setup


def clean_readme():
    with open('README.md') as f:
        return '\n'.join(line for line in f.read().splitlines()
                         if not line.startswith('[![') and
                         not line.startswith('<img'))


setup(long_description=clean_readme(),
      long_description_content_type='text/markdown')
