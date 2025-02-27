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
# toluene
# Conformer3D_COMPOUND_CID_1140.sdf
# Biphenyl
# Conformer3D_COMPOUND_CID_7095.sdf
# Ethylbenzene
# Conformer3D_COMPOUND_CID_7500.sdf
# Diphenylmethane
# Conformer3D_COMPOUND_CID_7580.sdf
# Triphenylmethane
# Conformer3D_COMPOUND_CID_10614.sdf
sdf_file_name_1 = 'test\Conformer3D_COMPOUND_CID_10614.sdf'
sdf_file = os.path.join(os.getcwd(), sdf_file_name_1)

# ================================
# ! CREATE MOLECULE GRAPH
# ================================

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

molecule_src = {
    'MainChain': ["C1*{Chain1}", "C1-C2", "C2=XX3"],
    'Chain1': ["*=C1", "C1-C2", "C2=C3", "C3-C4", "C4=C5", "C5-*"],
}

molecule_src = {
    'MainChain': ["C1*{Chain1}", "C1-C2", "C2*{Chain2}"],
    'Chain1': ["*=C1", "C1-C2", "C2=C3", "C3-C4", "C4=C5", "C5-*"],
    'Chain2': ["*=C1", "C1-C2", "C2=C3", "C3-C4", "C4=C5", "C5-*"],
}

molecule_src = {
    'MainChain': ["C1-H2", "C1-H3", "C1-C4", "C4*{Chain1}"],
    'Chain1': ["*=C1", "C1-C2", "C2=C3", "C3-C4", "C4=C5", "C5-*"],
}

molecule_src = {
    'MainChain': ["C1-H2", "C1-H3", "C1-C4", "C1-C5", "C4*{Chain1}", "C5*{Chain2}"],
    'Chain1': ["*=C1", "C1-C2", "C2=C3", "C3-C4", "C4=C5", "C5-*"],
    'Chain2': ["*=C1", "C1-C2", "C2=C3", "C3-C4", "C4=C5", "C5-*"],
}

molecule_src = {
    'MainChain': ["C1-H2", "C1-C3", "C1-C4", "C1-C5", "C3*{Chain1}", "C4*{Chain2}", "C5*{Chain3}"],
    'Chain1': ["*=C1", "C1-C2", "C2=C3", "C3-C4", "C4=C5", "C5-*"],
    'Chain2': ["*=C1", "C1-C2", "C2=C3", "C3-C4", "C4=C5", "C5-*"],
    'Chain3': ["*=C1", "C1-C2", "C2=C3", "C3-C4", "C4=C5", "C5-*"],
}

# * create molecule graph
mol_graph = mi.create_molecule_graph(molecule_src, molecule_name='my_molecule')
print(mol_graph)
# display molecule graph
mol_graph.d("my_molecule")

# =================================
# ! CREATE CUSTOM FUNCTIONAL GROUP
# =================================
# NOTE: custom functional group from file
# # current
# custom_functional_group_file = os.path.join(
#     os.getcwd(), 'test', 'custom-functional-group.yml')

# # NOTE: custom functional group from file (search a subgroup within group)
# custom_functional_group_subgroup_file = os.path.join(
#     os.getcwd(), 'test', 'custom-functional-group-2.yml')

# # create custom graph from file
# custom_g = mi.create_custom_functional_groups(
#     custom_functional_group_subgroup_file)
# print(custom_g)

# ==============================================
# ! FIND CUSTOM FUNCTIONAL GROUP
# ==============================================
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

# display the selected functional group
# benzene
# comp1.g3d_functional_group('benzene-full')
# comp1.g3d_functional_group('METHYLEN')
# comp1.g3d_functional_group("METHOXY")
