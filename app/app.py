# import packages/modules
import streamlit as st
import molinfo as mi
# local
from graph3dKit import Graph3dKit

st.title("Molinfo")

st.write("**MolInfo** is a Python package designed for advanced molecular analysis by converting molecular structures \
         into graph representations. This package enables researchers and chemists to load various molecular file formats,\
             transform them into graphs, and extract valuable information through graph-based methods.")

mi_version = mi.__version__
# write
st.write(f"Version: {mi_version}")

# more features are available here
url = 'https://molinfo.readthedocs.io/en/latest/'
st.markdown(
    f"[More features are available]({url})")

# upload a sdf file and store in a var
uploaded_file = st.file_uploader("Upload a SDF file", type=["sdf"])

# Check if a file has been uploaded
if uploaded_file is not None:
    # Read the contents of the file
    file_contents = uploaded_file.read()
    # convert to string
    file_contents = file_contents.decode("utf-8")

    #! create new compound
    compound = mi.compound(file_contents)
    # atom elements
    __atom_elements = compound.atomElements
    # atom bonds
    __atom_bonds = compound.atomBonds
    # atom xyz
    __atom_xyz = compound.atom_xyz
    # atom xyz center
    __atom_xyz_center = compound.atom_xyz_center
    # atom bonds 1d
    __atom_bonds_1d = compound.atomBonds_1d

    # compound
    __compound = {
        'atom_elements': __atom_elements,
        'atom_block': __atom_bonds,
        'xyz_list': __atom_xyz,
        'xyz_center_list': __atom_xyz_center,
        'atom_bonds_1d': __atom_bonds_1d
    }

    #! 3dview
    # create 3d view
    Graph3dKitC = Graph3dKit(__compound)

    # 3dview
    Graph3dKitC.view3d()
else:
    st.write("No file uploaded yet.")
