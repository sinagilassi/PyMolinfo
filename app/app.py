# import packages/modules
import streamlit as st
import molinfo as mi

st.title("Molinfo")

st.write("**MolInfo** is a Python package designed for advanced molecular analysis by converting molecular structures \
         into graph representations. This package enables researchers and chemists to load various molecular file formats,\
             transform them into graphs, and extract valuable information through graph-based methods.")

# upload a sdf file and store in a var
uploaded_file = st.file_uploader("Upload a SDF file", type=["sdf"])

# Check if a file has been uploaded
if uploaded_file is not None:
    # Read the contents of the file
    file_contents = uploaded_file.read()

    # Save the file contents to a variable
    sdf_data = file_contents

    # visualize compound by sdf file
    mi.g3d(sdf_data, display_bond_length=True)

    # Do something with the SDF data (e.g., parse it, visualize it, etc.)
    st.write("SDF file uploaded successfully!")
else:
    st.write("No file uploaded yet.")
