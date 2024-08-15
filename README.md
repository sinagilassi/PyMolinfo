# Molinfo

![Downloads](https://img.shields.io/pypi/dm/Molinfo) ![PyPI](https://img.shields.io/pypi/v/Molinfo) ![Python Version](https://img.shields.io/pypi/pyversions/Molinfo.svg) ![License](https://img.shields.io/pypi/l/Molinfo)

**MolInfo** is a Python package designed for advanced molecular analysis by converting molecular structures into graph representations. This package enables researchers and chemists to load various molecular file formats, transform them into graphs, and extract valuable information through graph-based methods.

**Features**

* `File Format Support`: Load molecular data from multiple file formats, including SDF and JSON (soon).
* `Graph Conversion`: Transform molecular structures into graph representations for detailed analysis.
* `Functional Group Identification`: Detect and analyze functional groups within the molecular graph.
* `Distance Measurement`: Compute distances between atoms and bonds in the molecular graph.
* `Bond Angle Calculation`: Measure angles between bonds using graph-based methods.

**Getting Started:**

To use Molinfo, simply install the package and import it into your Python script. Refer to the example code snippets above for a quick start.


## Google Colab

You can use the following code to run `Molinfo` in Google Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1rQXg92p_jxviVfKJFf_-1qQwmOgrMLUD?usp=sharing)

## Installation

Install molinfo with pip

```python
  pip install molinfo
```

## Documentation

Import package as:

```python
import molinfo as mi
# check version
print(mi.__version__)
```

## Examples

* Create a graph

```python
# sdf file
sdf_file_name_1 = 'test\Structure2D_COMPOUND_CID_261.sdf'
sdf_file = os.path.join(os.getcwd(), sdf_file_name_1)
# create graph
res = mi.create_graph(sdf_file)
print(type(res))
print(res)
```

* Display a graph:

```python
# visualize compound by sdf file
mi.g3d(sdf_file)
```

* Check the availability of functional groups:

```python
# check functional groups
res = mi.check_functional_group(sdf_file, res_format='dataframe')
print(res)
```

## FAQ

For any question, contact me on [LinkedIn](https://www.linkedin.com/in/sina-gilassi/) 


## Authors

- [@sinagilassi](https://www.github.com/sinagilassi)
