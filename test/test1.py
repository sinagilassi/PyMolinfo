# import packages/modules
import os
import molinfo as mi
from pprint import pprint as pp

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
sdf_file_name_1 = 'test\Structure2D_COMPOUND_CID_261.sdf'
sdf_file = os.path.join(os.getcwd(), sdf_file_name_1)

# compound
comp1 = mi.compound(sdf_file)
print(comp1)
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
