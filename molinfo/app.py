# import packages/modules
import os
import pubchemquery as pcq
import pandas as pd
# internal
from .config import packageName
from .config import packageShortName
from .config import __version__
from .config import __description__
from .docs import MolParser
from .docs import Compound
from .docs import CustomChemGraph


def main():
    # short description and version
    _des = f"{packageShortName} {__version__}: {__description__}"
    print(_des)


def compound(f):
    '''
    Create a compound by parsing sdf file

    Parameters
    ----------
    f : str
        molecule file format (sdf) or string (sdf)

    Returns
    -------
    compound : object
        compound object
    '''
    try:
        if os.path.exists(f):
            # parse file
            MolParserC = MolParser(f)
            compound_info = MolParserC.read_file()
            # compound
            compound = Compound(compound_info)
            # res
            return compound
        else:
            # check string
            if isinstance(f, str):
                # parse string
                MolParserC = MolParser(None)
                compound_info = MolParserC.read_file(
                    sourceContent={
                        'content': f.strip(),
                        'format': 'sdf'
                    })
                # compound
                compound = Compound(compound_info)
                # res
                return compound
            else:
                raise ValueError("Invalid input file or string")
    except Exception as e:
        print(e)


def compound_by_cid(cid):
    '''
    Create a compound by cid

    Parameters
    ----------
    cid : str
        compound id

    Returns
    -------
    compound : object
        compound object
    '''
    try:
        # get sdf by cid
        sdf = pcq.get_structure_by_cid(cid)
        # check
        if sdf is None or len(sdf) == 0:
            print("sdf file not found!")
            return None
        # compound
        return compound(sdf)
    except Exception as e:
        print(e)


def compound_by_inchi(inchi):
    '''
    Create a compound by inchi

    Parameters
    ----------
    inchi : str
        compound inchi

    Returns
    -------
    compound : object
        compound object
    '''
    try:
        # get cid by inchi
        cid = pcq.get_cid_by_inchi(inchi)
        # check
        if cid:
            # get sdf by cid
            sdf = pcq.get_structure_by_cid(cid)
            # compound
            return compound(sdf)
        else:
            print('the inchi not found!')
            return None
    except Exception as e:
        print(e)


def create_graph(file):
    '''
    Convert a sdf compound file to a graph 

    Parameters
    ----------
    file : str
        molecule file format (sdf)

    Returns
    -------
    graph : object
        compound graph
    '''
    try:
        if os.path.exists(file):
            # parse file
            MolParserC = MolParser(file)
            compound_info = MolParserC.read_file()
            # compound
            compound = Compound(compound_info)
            # display 3d
            graph = compound.create_graph()
            # res
            return graph
        else:
            raise Exception("file path is not valid.")
    except Exception as e:
        print(e)


def g3d(f, fig_size=[], bg_color='#ffffff', display_legend=True, display_atom_id=True, display_bond_length=False):
    '''
    3d graph of a compound

    Parameters
    ----------
    f : str
        molecule file format (sdf) or a sdf string variable
    display_legend : bool
        display legend (default True)
    display_atom_id : bool
        display atom id (default True)
    display_bond_length : bool
        display bond length (default True)

    Returns
    -------
    None
        display 3d graph
    '''
    try:
        # create a compound
        comp = compound(f)
        # display 3d
        comp.view3d(fig_size=fig_size, display_legend=display_legend, bg_color=bg_color,
                    display_atom_id=display_atom_id, display_bond_length=display_bond_length)
    except Exception as e:
        raise Exception(f"file path/variable is not valid! {e}")


def g3d_by_inchi(inchi, fig_size=[], bg_color='#ffffff', display_legend=True, display_atom_id=True, display_bond_length=False):
    '''
    3d graph of a compound using its InChI identifier

    Parameters
    ----------
    inchi : str
        inchi code
    display_legend : bool
        display legend (default True)

    Returns
    -------
    None
        display 3d graph
    '''
    # check inchi
    if inchi is not None:
        # get cid
        cid = pcq.get_cid_by_inchi(inchi)
        # check
        if cid is not None:
            # get sdf
            sdf = pcq.get_structure_by_cid(cid)
            # check
            if sdf is not None:
                # display 3d
                # parse file
                MolParserC = MolParser(None)
                compound_info = MolParserC.read_file(
                    sourceContent={
                        'content': sdf,
                        'format': 'sdf'
                    })
                # compound
                compound = Compound(compound_info)
                # display 3d
                compound.view3d(fig_size=fig_size, display_legend=display_legend, bg_color=bg_color,
                                display_atom_id=display_atom_id, display_bond_length=display_bond_length)
            else:
                raise Exception("sdf is not found!")
        else:
            raise Exception("cid is not found!")
    else:
        raise Exception("inchi is not valid.")


def check_functional_group(file, functional_groups=[], res_format='raw'):
    '''
    Check a functional group exists in a compound

    Parameters
    ----------
    file : str
        molecule file format (sdf)
    functional_groups : list[str] or CustomChemGraph object
        functional group (default ['hydroxyl']) or CustomChemGraph object
    res_format : str
        result format (default 'raw')

    Returns
    -------
    res : dict
        a list of all count
    compound : object
        compound object (sdf file)
    '''
    # check file exists
    if os.path.exists(file):
        # parse file
        MolParserC = MolParser(file)
        compound_info = MolParserC.read_file()
        # compound
        compound = Compound(compound_info)
        # check functional group
        res = compound.check_functional_groups(functional_groups)
        # check
        if res_format == 'dataframe':
            # dataframe
            return pd.DataFrame(res), compound
        else:
            # raw
            return res, compound
    else:
        raise Exception("file path is not valid.")


def count_functional_group(file, functional_groups=[], res_format='raw'):
    '''
    Count the occurrences of functional groups within the structure of a compound.

    Parameters
    ----------
    file : str
        molecule file format (sdf)
    functional_groups : list[str] or CustomChemGraph object
        functional group (default ['hydroxyl']) or CustomChemGraph object
    res_format : str
        result format (default 'raw')

    Returns
    -------
    res : dict
        a list of all count
    compound : object
        compound object (sdf file)
    '''
    # check file exists
    if os.path.exists(file):
        # parse file
        MolParserC = MolParser(file)
        compound_info = MolParserC.read_file()
        # compound
        compound = Compound(compound_info)
        # create graph
        compound.create_graph()
        # check functional group
        res = compound.check_functional_groups(
            functional_groups, count_functional_group=True)
        # check
        if res_format == 'dataframe':
            # dataframe
            return pd.DataFrame(res), compound
        else:
            # raw
            return res, compound
    else:
        raise Exception("file path is not valid.")


def create_custom_functional_groups(functional_groups):
    '''
    create custom functional groups based on the following format:
    # CH2-O
    custom_functional_group = [
    {'fg1': ["C1-H1","C1-H2","C1-O1"]},
    {'fg2': ["C1-H1","C1-H2","C1-C2","C2-H3","C2-O2"]}
    ]

    Parameters
    ----------
    functional_groups : list[dict] or dict
        functional group

    Returns
    -------
    custom_functional_groups : list[dict]
        a list of all custom functional group
    '''
    try:
        # check format
        if isinstance(functional_groups, dict):
            # set a list
            custom_functional_groups = [functional_groups]
        elif isinstance(functional_groups, list):
            custom_functional_groups = functional_groups
        else:
            raise Exception("functional_groups is not valid.")

        # custom chem graph
        CustomChemGraphC = CustomChemGraph(custom_functional_groups)
        # res
        return CustomChemGraphC
    except Exception as e:
        raise Exception(e)


if __name__ == "__main__":
    main()
