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
# Naphthalene
# Conformer3D_COMPOUND_CID_931.sdf
sdf_file_name_1 = 'test\Conformer3D_COMPOUND_CID_7095.sdf'
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

# naphthalene
molecule_src = {
    'MainChain': ["C1-C2", "C2=C3", "C3-C4", "C4=C5", "C5-C6", "C6=C1", "C1*{Chain1}", "C6*{Chain1}"],
    'Chain1': ["*-C1", "C1=C2", "C2-C3", "C3=C4", "C4-**"],
}


# NOTE: create a molecule
mol_ = mi.generate_molecule(molecule_src, molecule_name='my_molecule')
# construct molecule
constructed_molecule = mol_.constructed_molecules
print(constructed_molecule)
# main chain
chain_info = mol_.chain_info
print(chain_info)
# molecule info
molecule = mol_.molecule
print(molecule)

# NOTE: create a molecule graph (nx.Graph)
graph = mol_.to_graph()
print(graph)

# NOTE: create a molgraph (MolGraph)
mol_graph = mol_.to_molgraph()
# display molecule graph
mol_graph.d()


# SECTION: create parent & child
# CB-(CB)
molecule_src_child = {
    'MainChain': ["C1-C2", "C1*{Chain1}"],
    'Chain1': ["*=C1", "C1-C2", "C2=C3", "C3-C4", "C4=C5", "C5-*"]
}

# biphenyl
molecule_src_parent = {
    'MainChain': ["C1*{Chain1}", "C1-C2", "C2*{Chain2}"],
    'Chain1': ["*-C1", "C1=C2", "C2-C3", "C3=C4", "C4-C5", "C5=*"],
    'Chain2': ["*-C1", "C1=C2", "C2-C3", "C3=C4", "C4-C5", "C5=*"],
}

# create molecule
mol_child = mi.generate_molecule(
    molecule_src_child, molecule_name='my_molecule_child')
mol_parent = mi.generate_molecule(
    molecule_src_parent, molecule_name='my_molecule_parent')
# create graph
mol_graph_child = mol_child.to_molgraph()
mol_graph_parent = mol_parent.to_molgraph()
# display molecule graph
mol_graph_child.d()
mol_graph_parent.d()

# =================================
# ! CREATE CUSTOM FUNCTIONAL GROUP
# =================================
# NOTE: dict of custom functional group
custom_functional_group = {
    'fg_parent': mol_parent.constructed_molecule,
    'fg_child': mol_child.constructed_molecule,
    'fg_target': ['fg_parent', 'fg_child']
}

# log
print(custom_functional_group)

custom_g = mi.create_custom_functional_groups(custom_functional_group)
custom_g.d("fg_parent")
custom_g.d("fg_child")
print(custom_g.custom_functional_groups)


# ==============================================
# ! FIND CUSTOM FUNCTIONAL GROUP
# ==============================================
# dataframe format
res, comp1 = mi.check_functional_group(sdf_file, functional_groups=[
                                       custom_g], res_format='dataframe')
print(res)

# NOTE: display the selected functional group
# comp1.g3d_functional_group('HOC=C')
# benzene
# comp1.g3d_functional_group('benzene')

# NOTE: count the custom functional groups
res, comp1 = mi.count_functional_group(sdf_file, functional_groups=[
    custom_g], res_format='dataframe')
print(res)

# display the selected functional group
# benzene
# comp1.g3d_functional_group('benzene-full')
