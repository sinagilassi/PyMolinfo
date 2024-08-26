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

## Binder

Test this package by launching our example notebooks on Binder:

| Description | Launch Binder | 
| --- | --- |
| Load a sdf file | [![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/sinagilassi/Molinfo/cb4f3c4b58501786da3dc5a2413a67720f01d579?urlpath=lab%2Ftree%2Fnotebook%2Fdoc-1.ipynb) |
| Visualize a compound | [![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/sinagilassi/Molinfo/cb4f3c4b58501786da3dc5a2413a67720f01d579?urlpath=lab%2Ftree%2Fnotebook%2Fdoc-2.ipynb) |
| Check and count functional groups | [![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/sinagilassi/Molinfo/cb4f3c4b58501786da3dc5a2413a67720f01d579?urlpath=lab%2Ftree%2Fnotebook%2Fdoc-3.ipynb)|
| Create custom functional groups | [![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/sinagilassi/Molinfo/cb4f3c4b58501786da3dc5a2413a67720f01d579?urlpath=lab%2Ftree%2Fnotebook%2Fdoc-4.ipynb)|


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
res, comp1 = mi.check_functional_group(sdf_file, res_format='dataframe')
print(res)
```

* Calculate angle/distance between atoms

```python
# distance matrix
res_distance = comp1.distance_matrix(dataframe=True)
print(res_distance)

# distance between two atoms
distance = comp1.distance_atoms(['O1', 'C2'])
print(distance)

# angle between atoms
angle = comp1.angle_atoms(['O1', 'C2', 'H3'])
print(angle)

# dihedral angle
dihedral = comp1.d_angle_atoms(['H6', 'O1', 'C2', 'H3'])
print(dihedral)
```

* Create custom functional groups:

[`atom1-element`][`atom1-number`][`bond-type`][`atom2-element`][`atom2-number`]

|  Bond Types | Format  | 
|:----------|:----------|
| single bond CC   | C1-C2   | 
| double bond CC   | C1=C2   | 
| triple bond CC   | C1#C2   | 

**How to create a custom functional group?**

|  Name |  Symbol | Format |
|:-----------|:------------:|-------------:|
|  cyanide-1     |     CCN   | ["N1#C2"]      |
| custom_fg      | NCH       | ["N1-C2", "C2-H3"]       |
| NC=O | NC=O | ["N1-C2", "C2=O3"] |

And coded as:

```python
# C1-C2#N3
custom_functional_group = [
    {'cyanide': ["C1-C2", "C2#N3"]},
]

# define different custom functional groups as:
# N#C
# NCH
# NCO
custom_functional_group = [
    {'N#C': ["N1#C2"]},
    {'custom_fg': ["N1-C2", "C2-H3"]},
    {'NC=O': ["N1-C2", "C2=O3"]},
]

# create custom graph
custom_g = mi.create_custom_functional_groups(custom_functional_group)

# visualize custom graph
# custom_g.d("cyanide")

# find custom functional groups in a compound
res = mi.check_functional_group(
    sdf_file, functional_groups=[custom_g])
print(res)
```

## FAQ

For any question, contact me on [LinkedIn](https://www.linkedin.com/in/sina-gilassi/) 


## Authors

- [@sinagilassi](https://www.github.com/sinagilassi)
