# import packages/modules
from rich import print
import os
import pyMolinfo as mi

# check version
print(mi.__version__)

# ==========================
# ! SDF FILES
# ==========================
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
# Ethinyl Estradiol
# Conformer3D_COMPOUND_CID_5991.sdf
sdf_file_name_1 = 'test\Conformer3D_COMPOUND_CID_241.sdf'
sdf_file = os.path.join(os.getcwd(), sdf_file_name_1)
sdf_name = 'Benzene'

# ===============================
# ! CREATE COMPOUND
# ===============================
# compound
comp1 = mi.compound(sdf_file)
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

# NOTE: functional groups
# print(comp1.functional_groups)

# ===============================
# ! GEOMETRY
# ===============================
# distance
# res_distance = comp1.distance_matrix(dataframe=True)
# pp(res_distance)

# _distance = comp1.distance_atoms(['O1', 'C2'])
# print(_distance)
# _distance = comp1.distance_atoms(['O1', 'H3'])
# print(_distance)
# _distance = comp1.distance_atoms(['O1', 'H4'])
# print(_distance)
# _distance = comp1.distance_atoms(['O1', 'H5'])
# print(_distance)
# _distance = comp1.distance_atoms(['O1', 'H6'])
# print(_distance)
# _distance = comp1.distance_atoms(['C2', 'H3'])
# print(_distance)

# angle
# res_angle = comp1.angle_atoms(['O1', 'C2', 'H3'])
# print(res_angle)

# res_angle = comp1.angle_atoms(['O1', 'C2', 'H4'])
# print(res_angle)
# res_angle = comp1.angle_atoms(['O1', 'C2', 'H5'])
# print(res_angle)
# res_angle = comp1.angle_atoms(['O1', 'C2', 'H6'])
# print(res_angle)
# res_angle = comp1.angle_atoms(['C2', 'H5', 'H3'])
# print(res_angle)

# dihedral angle
# res_dihedral = comp1.d_angle_atoms(['H18', 'C10', 'C7', 'H15'])
# print(res_dihedral)

# ================================
# ! CREATE GRAPH
# ================================
# NOTE: create graph
graph_1 = mi.create_graph(sdf_file, graph_name=sdf_name)
print(type(graph_1))
print(graph_1)

# view graph
mi.view_graph(graph_1)

stop=1

# * visualize compound by sdf file
# mi.g3d(sdf_file, display_bond_length=True)

# visualize compound by sdf string
# mi.g3d(sdf_string)

# visualize compound by inchi
# mi.g3d_by_inchi(
#     'InChI=1S/C14H22O6/c1-11(2)13(15)19-9-7-17-5-6-18-8-10-20-14(16)12(3)4/h1,3,5-10H2,2,4H3')

# ================================
# ! CREATE MOLECULE GRAPH
# ================================
# molecule source
# molecule_src = {
#     'MainChain': ["C1-C2", "C2-C3", "C3*{Chain1}", "C3-C4", "C4*{Chain2}", "C4-C5", "C5-C6"],
#     'Chain1': ["C1=C2", "C2-C3", "C3=*"],
#     'Chain2': ["*-C1", "C1=C2", "C2-XX3"]
# }

# molecule_src = {
#     'MainChain': ["C1-C2", "C2=C3", "C3-C4", "C3*{Chain1}", "C4=C5", "C5-C6", "C6=C1", "C6*{Chain2}"],
#     'Chain1': ["C1=C2", "C2-C3", "C3=*"],
#     'Chain2': ["*-C1", "C1=C2", "C2-XX3"]
# }

molecule_src = {
    'MainChain': ["C1-C2", "C2=C3", "C3-C4", "C3*{Chain1}", "C4=C5", "C5*{Chain1}", "C5-C6", "C6=C1", "C6*{Chain2}"],
    'Chain1': ["C1=C2", "C2-C3", "C3=*"],
    'Chain2': ["*-C1", "C1=C2", "C2-XX3"]
}

molecule_src = {
    'MainChain': ["C1*{Chain1}", "C1-C2", "C2*{Chain2}"],
    'Chain1': ["*-C1", "C1=C2", "C2-C3", "C3=C4", "C4-C5", "C5=*"],
    'Chain2': ["*-C1", "C1=C2", "C2-C3", "C3=C4", "C4-C5", "C5=*"],
}

molecule_src = {
    'MainChain': ["C1*{Chain1}", "C1-H2"],
    'Chain1': ["*=C1", "C1-C2", "C2=C3", "C3-C4", "C4=C5", "C5-*"],
}
# * create molecule graph
mol_ = mi.generate_molecule(molecule_src, molecule_name='my_molecule')
print(mol_)

# create graph
mol_graph = mol_.to_molgraph()

# display molecule graph
mol_graph.d()

stop = 1
# ================================
# ! CHECK FUNCTIONAL GROUP
# ================================
# network
# res, comp1 = mi.check_functional_group(
#     sdf_file, functional_groups=['ether'])

# dataframe format
# res, comp1 = mi.check_functional_group(sdf_file, res_format='dataframe')

# raw format
# res, comp1 = mi.check_functional_group(sdf_file)

# print(res)
# print(comp1.functional_groups)

# ================================
# ! COUNT FUNCTIONAL GROUP
# ================================
# raw format
# res, comp1 = mi.count_functional_group(sdf_file, functional_groups=[
#                                        'tertiary-alcohol'], res_format='dataframe')

# res, comp1 = mi.count_functional_group(sdf_file, res_format='dataframe')

# print(res)
# print(comp1.functional_groups)
# comp1.g3d_functional_group('hydroxyl')

# =================================
# ! CREATE CUSTOM FUNCTIONAL GROUP
# =================================
# # custom
# NOTE: list of custom functional group
custom_functional_group = [
    {'cyanide': ["C1-C2", "C2#N3"]},
    {'N#C': ["N1#C2"]},
    {'fg1': ["N1-C2", "C2-H3"]},
    {'NC=O': ["N1-C2", "C2=O3"]},
    {'HOC=C': ["H1-O2", "O2-C3", "C3=C4"]}
]

# create custom graph
# custom_g = mi.create_custom_functional_groups(custom_functional_group)
# display custom functional group
# custom_g.d("cyanide")
# print(custom_g.custom_functional_groups)

# REVIEW: create labels
# Benzene (C6H6) -> XCB

# NOTE: dict of custom functional group
custom_functional_group = {
    'cyanide': ["C1-C2", "C2#N3"],
    'N#C': ["N1#C2"],
    'fg1': ["N1-C2", "C2-H3"],
    'NC=O': ["N1-C2", "C2=O3"],
    'HOC=C': ["H1-O2", "O2-C3", "C3=C4"]
}

# custom_g = mi.create_custom_functional_groups(custom_functional_group)
# custom_g.d("fg1")
# print(custom_g.custom_functional_groups)

# NOTE: custom functional group from file
# current
custom_functional_group_file = os.path.join(
    os.getcwd(), 'test', 'custom-functional-group.yml')

# create custom graph from file
# custom_g = mi.create_custom_functional_groups(custom_functional_group_file)
# custom_g.d("METHOXY")
# custom_g.d("benzene-full")
# custom_g.d('benzene')
# custom_g.d('benzene2')
# custom_g.d('CB')

# NOTE: custom functional group from file (search a subgroup within group)
custom_functional_group_subgroup_file = os.path.join(
    os.getcwd(), 'test', 'custom-functional-group-2.yml')

# create custom graph from file
custom_g = mi.create_custom_functional_groups(
    custom_functional_group_subgroup_file)
print(custom_g)

# ==============================================
# ! FIND CUSTOM FUNCTIONAL GROUP
# ==============================================
# original
# res = mi.check_functional_group(
#     sdf_file, functional_groups=[custom_g])
# pp(res)

# # dataframe format
# res, comp1 = mi.check_functional_group(sdf_file, functional_groups=[
#                                        custom_g], res_format='dataframe')
# print(res)

# dataframe format
res, comp1 = mi.check_functional_group(sdf_file, functional_groups=[
                                       mol_graph], res_format='dataframe')
print(res)

# NOTE: display the selected functional group
# comp1.g3d_functional_group('HOC=C')
# benzene
# comp1.g3d_functional_group('benzene')

# NOTE: count the custom functional groups
res, comp1 = mi.count_functional_group(sdf_file, functional_groups=[
    mol_graph], res_format='dataframe')
print(res)

# res, comp1 = mi.count_functional_group(sdf_file, functional_groups=[
#     custom_g])
# pp(res)

# display the selected functional group
# benzene
# comp1.g3d_functional_group('benzene-full')
# comp1.g3d_functional_group('METHYLEN')
# comp1.g3d_functional_group("METHOXY")
