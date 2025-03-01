# 🌟 PyMolinfo

![PyMolinfo](https://drive.google.com/uc?export=view&id=1VZz79YQbWXMosfUwGBeTrTdHIiOXt_Ps)

![Downloads](https://img.shields.io/pypi/dm/PyMolinfo) ![PyPI](https://img.shields.io/pypi/v/PyMolinfo) ![Python Version](https://img.shields.io/pypi/pyversions/PyMolinfo.svg) ![License](https://img.shields.io/pypi/l/PyMolinfo)

**PyMolInfo** (previously molinfo) is a Python package designed for advanced molecular analysis by converting molecular structures into graph representations. This package enables researchers and chemists to load various molecular file formats, transform them into graphs, and extract valuable information through graph-based methods.

## ✨ Features

* `File Format Support`: Load molecular data from multiple file formats, including SDF and JSON (soon).
* `Graph Conversion`: Transform molecular structures into graph representations for detailed analysis.
* `Functional Group Identification`: Detect and analyze functional groups within the molecular graph.
* `Distance Measurement`: Compute distances between atoms and bonds in the molecular graph.
* `Bond Angle Calculation`: Measure angles between bonds using graph-based methods.

## 🚀 Getting Started

To use PyMolinfo, simply install the package and import it into your Python script. Refer to the example code snippets above for a quick start.

## 📚 Binder

Test this package by launching our example notebooks on Binder:

- **Load a sdf file**: [![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/sinagilassi/Molinfo/cb4f3c4b58501786da3dc5a2413a67720f01d579?urlpath=lab%2Ftree%2Fnotebook%2Fdoc-1.ipynb)
- **Visualize a compound**: [![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/sinagilassi/Molinfo/cb4f3c4b58501786da3dc5a2413a67720f01d579?urlpath=lab%2Ftree%2Fnotebook%2Fdoc-2.ipynb)
- **Check and count functional groups**: [![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/sinagilassi/Molinfo/cb4f3c4b58501786da3dc5a2413a67720f01d579?urlpath=lab%2Ftree%2Fnotebook%2Fdoc-3.ipynb)
- **Create custom functional groups**: [![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/sinagilassi/Molinfo/cb4f3c4b58501786da3dc5a2413a67720f01d579?urlpath=lab%2Ftree%2Fnotebook%2Fdoc-4.ipynb)

## 🌐 Google Colab

You can use the following code to run `PyMolinfo` in Google Colab:

- **Version 1.6.0**: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1-mkVbXbznEJGeKWdQKtJT8xkWb2Bcvw_?usp=sharing)
- **Version < 1.6.0**: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1rQXg92p_jxviVfKJFf_-1qQwmOgrMLUD?usp=sharing)

## 🛠️ Installation

Install molinfo with pip

```python
  pip install PyMolinfo
```

## 📖 Documentation

Import package as:

```python
import pyMolinfo as mi
# check version
print(mi.__version__)
```

## 💡 Examples

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

## Creating Custom Functional Groups

To create custom functional groups, you need to define the bonds between atoms using the following format:

`[atom1-element][atom1-number][bond-type][atom2-element][atom2-number]`

Here are the formats for different bond types:

- **Single Bond**: Represented as `C1-C2` where `C1` and `C2` are the atoms connected by a single bond.
- **Double Bond**: Represented as `C1=C2` where `C1` and `C2` are the atoms connected by a double bond.
- **Triple Bond**: Represented as `C1#C2` where `C1` and `C2` are the atoms connected by a triple bond.

### Examples

1. **Cyanide Group**: A cyanide group can be represented as `N1#C2`.

```python
custom_functional_group = [
    {'cyanide': ["N1#C2"]},
]
```

2. **Custom Functional Group**: A custom functional group with a single and a double bond can be represented as `N1-C2` and `C2=O3`.

```python
custom_functional_group = [
    {'custom_fg': ["N1-C2", "C2=O3"]},
]
```

3. **Multiple Functional Groups**: You can define multiple functional groups in a list.

```python
custom_functional_group = [
    {'N#C': ["N1#C2"]},
    {'custom_fg': ["N1-C2", "C2-H3"]},
    {'NC=O': ["N1-C2", "C2=O3"]},
]
```

Once you have defined your custom functional groups, you can create and visualize them as follows:

```python
# create custom graph
custom_g = mi.create_custom_functional_groups(custom_functional_group)

# visualize custom graph
# custom_g.d("cyanide")

# find custom functional groups in a compound
res = mi.check_functional_group(
    sdf_file, functional_groups=[custom_g])
print(res)
```

## ❓ FAQ

For any question, contact me on [LinkedIn](https://www.linkedin.com/in/sina-gilassi/)

## 👥 Authors

[@sinagilassi](https://www.github.com/sinagilassi)
