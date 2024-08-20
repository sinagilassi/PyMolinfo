# import packages/modules
from pprint import pprint as pp
import pprint
import os
import molinfo as mi

# check version
# print(mi.__version__)

# sdf file
# Methanol
# Conformer3D_COMPOUND_CID_887
# Benzene
# Conformer3D_COMPOUND_CID_241
# Glycol Dimethacrylate
# Conformer3D_COMPOUND_CID_7355
# 2-(2-Methoxyethoxy)ethanol
# Conformer3D_COMPOUND_CID_8134
# Bromobenzyl cyanide
# Conformer3D_COMPOUND_CID_22044
# Amineptine
# Conformer3D_COMPOUND_CID_34870
# Amidephrine
# Conformer3D_COMPOUND_CID_15010
# N,N-Dimethylformamide
# Conformer3D_COMPOUND_CID_6228
# Acetanilide
# Conformer3D_COMPOUND_CID_904
# Ascorbic Acid
# Conformer3D_COMPOUND_CID_54670067
# Malathion
# Conformer3D_COMPOUND_CID_4004
# Benzbromarone
# Structure2D_COMPOUND_CID_2333
# Butyraldehyde
# Structure2D_COMPOUND_CID_261
sdf_file_name_1 = 'test\Conformer3D_COMPOUND_CID_241.sdf'
sdf_file = os.path.join(os.getcwd(), sdf_file_name_1)

# sdf string
# sdf_string = """
# 241
#   -OEChem-08032421353D

#  12 12  0     0  0  0  0  0  0999 V2000
#    -1.2131   -0.6884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -1.2028    0.7064    0.0001 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -0.0103   -1.3948    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#     0.0104    1.3948   -0.0001 C   0  0  0  0  0  0  0  0  0  0  0  0
#     1.2028   -0.7063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#     1.2131    0.6884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -2.1577   -1.2244    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#    -2.1393    1.2564    0.0001 H   0  0  0  0  0  0  0  0  0  0  0  0
#    -0.0184   -2.4809   -0.0001 H   0  0  0  0  0  0  0  0  0  0  0  0
#     0.0184    2.4808    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#     2.1394   -1.2563    0.0001 H   0  0  0  0  0  0  0  0  0  0  0  0
#     2.1577    1.2245    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#   1  2  2  0  0  0  0
#   1  3  1  0  0  0  0
#   1  7  1  0  0  0  0
#   2  4  1  0  0  0  0
#   2  8  1  0  0  0  0
#   3  5  2  0  0  0  0
#   3  9  1  0  0  0  0
#   4  6  2  0  0  0  0
#   4 10  1  0  0  0  0
#   5  6  1  0  0  0  0
#   5 11  1  0  0  0  0
#   6 12  1  0  0  0  0
# M  END
# > <PUBCHEM_COMPOUND_CID>
# 241

# > <PUBCHEM_CONFORMER_RMSD>
# 0.4

# > <PUBCHEM_CONFORMER_DIVERSEORDER>
# 1

# > <PUBCHEM_MMFF94_PARTIAL_CHARGES>
# 12
# 1 -0.15
# 10 0.15
# 11 0.15
# 12 0.15
# 2 -0.15
# 3 -0.15
# 4 -0.15
# 5 -0.15
# 6 -0.15
# 7 0.15
# 8 0.15
# 9 0.15

# > <PUBCHEM_EFFECTIVE_ROTOR_COUNT>
# 0

# > <PUBCHEM_PHARMACOPHORE_FEATURES>
# 1
# 6 1 2 3 4 5 6 rings

# > <PUBCHEM_HEAVY_ATOM_COUNT>
# 6

# > <PUBCHEM_ATOM_DEF_STEREO_COUNT>
# 0

# > <PUBCHEM_ATOM_UDEF_STEREO_COUNT>
# 0

# > <PUBCHEM_BOND_DEF_STEREO_COUNT>
# 0

# > <PUBCHEM_BOND_UDEF_STEREO_COUNT>
# 0

# > <PUBCHEM_ISOTOPIC_ATOM_COUNT>
# 0

# > <PUBCHEM_COMPONENT_COUNT>
# 1

# > <PUBCHEM_CACTVS_TAUTO_COUNT>
# 1

# > <PUBCHEM_CONFORMER_ID>
# 000000F100000001

# > <PUBCHEM_MMFF94_ENERGY>
# 13.148

# > <PUBCHEM_FEATURE_SELFOVERLAP>
# 5.074

# > <PUBCHEM_SHAPE_FINGERPRINT>
# 16714656 1 18123201108121732318
# 20096714 4 18339082562313515547
# 21015797 1 9222644709913001990
# 21040471 1 18194402190896164549

# > <PUBCHEM_SHAPE_MULTIPOLES>
# 123.48
# 1.59
# 1.59
# 0.62
# 0
# 0
# 0
# 0
# 0
# 0
# 0
# 0
# 0
# 0

# > <PUBCHEM_SHAPE_SELFOVERLAP>
# 251.998

# > <PUBCHEM_SHAPE_VOLUME>
# 67.5

# > <PUBCHEM_COORDINATE_TYPE>
# 2
# 5
# 10

# $$$$
# """

# compound
# comp1 = mi.compound(sdf_file)
# sdf string
# comp1 = mi.compound(sdf_string)
# compound by cid
# comp1 = mi.compound_by_cid(241)
# compound by inchi
# comp1 = mi.compound_by_inchi('InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H')
# log
# print(comp1)
# pp(comp1.atom_bond_block)
# print("-"*100)
# pp(comp1.atom_bond_block_1d)
# print("-"*100)
# pp(comp1.atom_xyz)

# # create graph
# res = mi.create_graph(sdf_file)
# print(type(res))
# print(res)

# visualize compound by sdf file
# mi.g3d(sdf_file, display_bond_length=True)

# visualize compound by inchi
# mi.g3d_by_inchi(
#     'InChI=1S/C14H22O6/c1-11(2)13(15)19-9-7-17-5-6-18-8-10-20-14(16)12(3)4/h1,3,5-10H2,2,4H3')

# network
# res = mi.check_functional_group(sdf_file, functional_groups=['ether'])
# print(res)

# res = mi.check_functional_group(sdf_file, res_format='dataframe')
# print(res)

# res = mi.check_functional_group(sdf_file)
# print(res)

# ===============================
# CREATE CUSTOM FUNCTIONAL GROUP
# ===============================
# custom
custom_functional_group = [
    {'fg1': ["C1-H1", "C1-H2", "C1-O1"]},
    {'fg2': ["C1-H1", "C1-H2", "C1-C2", "C2-H3", "C2-O2"]}
]

# create custom graph
custom_g = mi.crate_custom_functional_groups(custom_functional_group)
custom_g.d("fg2")

# find custom functional groups in a compound