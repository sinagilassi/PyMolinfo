from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

APP_NAME = 'PyMolinfo'
AUTHOR = 'Sina Gilassi'
VERSION = '1.8.0'
LICENSE = 'MIT'
DESCRIPTION = 'PyMolinfo provides comprehensive molecular information and analysis.'
LONG_DESCRIPTION = 'PyMolinfo is a Python package designed for advanced molecular analysis by converting molecular structures into graph representations'

# Setting up
setup(
    name=APP_NAME,
    version=VERSION,
    author=AUTHOR,
    author_email="<sina.gilassi@gmail.com>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(exclude=['tests', '*.tests', '*.tests.*']),
    license=LICENSE,
    install_requires=['pandas', 'pillow', 'requests',
                      'urllib3', 'matplotlib', 'PubChemQuery', 'numpy', 'plotly', 'networkx'],
    keywords=['python', 'chemistry', 'chemistry-visualization',
              'PyMolinfo', 'molecular-graph'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ],
    python_requires='>=3.10',
)
