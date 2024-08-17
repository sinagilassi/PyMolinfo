from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

APP_NAME = 'Molinfo'
<<<<<<< HEAD
VERSION = '1.0.0'
=======
VERSION = '1.2.0'
>>>>>>> c4f5d6f634eb8bb37285ee86c97872be334d0984
DESCRIPTION = 'Molinfo provides comprehensive molecular information and analysis.'
LONG_DESCRIPTION = 'Molinfo is a Python package designed for advanced molecular analysis by converting molecular structures into graph representations'

# Setting up
setup(
    name=APP_NAME,
    version=VERSION,
    author="Sina Gilassi",
    author_email="<sina.gilassi@gmail.com>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(exclude=['tests', '*.tests', '*.tests.*']),
    license='MIT',
    install_requires=['pandas', 'pillow', 'requests',
<<<<<<< HEAD
                      'urllib3', 'matplotlib', 'PubChemQuery', 'numpy', 'plotly'],
    keywords=['python', 'chemistry', 'chemistry-visualization', 'Molinfo'],
=======
                      'urllib3', 'matplotlib', 'PubChemQuery', 'numpy', 'plotly', 'networkx[default]'],
    keywords=['python', 'chemistry', 'chemistry-visualization',
              'Molinfo', 'molecular-graph'],
>>>>>>> c4f5d6f634eb8bb37285ee86c97872be334d0984
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ],
    python_requires='>=3.6',
)
