# import packages/modules
import os
from pathlib import Path
from networkx import Graph
import pubchemquery as pcq
import pandas as pd
from typing import List, Dict, Union, Literal

# internal
from .config import packageName
from .config import packageShortName
from .config import __version__
from .config import __description__
from .config import __author__
from .docs import MolParser, Compound, CustomChemGraph, Utility


def main():
    # short description and version
    _des = f"{packageShortName} {__version__}: {__description__}"
    print(_des)


def compound(f: Union[str, Path]) -> Compound:
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
        # check str, path
        if not isinstance(f, (str, Path)):
            raise ValueError("Invalid input file path or string")

        # check file exists
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
                raise ValueError("Invalid input file path or string")
    except Exception as e:
        raise Exception(f"creating compound is failed! {e}")


def compound_by_cid(cid: Union[str, int]) -> Compound:
    '''
    Create a compound by cid

    Parameters
    ----------
    cid : str | int
        compound id

    Returns
    -------
    compound : object
        compound object
    '''
    try:
        # get sdf by cid
        sdf = pcq.get_structure_by_cid(str(cid))
        # check
        if sdf is None or len(sdf) == 0:
            # err
            raise Exception("sdf is not found!")
        # compound
        return compound(sdf)
    except Exception as e:
        raise Exception(f"creating compound is failed! {e}")


def compound_by_inchi(inchi: str) -> Compound:
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
            # check
            if sdf is None or len(sdf) == 0:
                # err
                raise Exception("sdf is not found or has no content!")
            # compound
            return compound(sdf)
        else:
            # err
            raise Exception("inchi is not valid.")
    except Exception as e:
        raise Exception(f"creating compound is failed! {e}")


def create_graph(file: Path) -> Graph:
    '''
    Converts a sdf compound file to a graph

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
        raise Exception(f"creating graph is failed! {e}")


def g3d(f: Union[str, Path], fig_size: List = [], bg_color: str = '#ffffff',
        display_legend: bool = True, display_atom_id: bool = True,
        display_bond_length: bool = False):
    '''
    3d graph of a compound

    Parameters
    ----------
    f : str
        molecule file format (sdf) or a sdf string variable
    fig_size : list
        figure size (default [])
    bg_color : str
        background color (default '#ffffff')
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


def g3d_by_inchi(inchi: str, fig_size: List = [], bg_color: str = '#ffffff', display_legend: bool = True, display_atom_id: bool = True, display_bond_length: bool = False):
    '''
    3d graph of a compound using its InChI identifier

    Parameters
    ----------
    inchi : str
        inchi code
    fig_size : list
        figure size (default [])
    bg_color : str
        background color (default '#ffffff')
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


def check_functional_group(file: Union[str, Path], functional_groups: List[Union[str, CustomChemGraph]] = [], res_format: Literal['original', 'dataframe'] = 'original'):
    '''
    Check a functional group exists in a compound

    Parameters
    ----------
    file : Path | str
        molecule file format (sdf) or string (sdf)
    functional_groups : list[str] or CustomChemGraph object
        functional group (default ['hydroxyl']) or CustomChemGraph object
    res_format : str
        result format (default 'original')

    Returns
    -------
    res : dict
        a list of all count
    compound : Compound
        compound object (sdf file)
    '''
    try:
        # create compound
        comp = compound(file)

        # check functional group
        res = comp.check_functional_groups(functional_groups)

        # check
        if res_format == 'dataframe':
            # dataframe
            df = pd.DataFrame(res)
            return df, comp
        elif res_format == 'original':
            # raw
            return res, comp
        else:
            raise Exception("res_format is not valid.")

    except Exception as e:
        raise Exception(f"checking functional group is failed! {e}")


def count_functional_group(file: Union[str, Path], functional_groups: List[Union[str, CustomChemGraph]] = [], res_format: Literal['original', 'dataframe'] = 'original'):
    '''
    Counts the occurrences of functional groups within the structure of a compound.

    Parameters
    ----------
    file : Path | str
        molecule file format (sdf) or string (sdf)
    functional_groups : list[str] or CustomChemGraph object
        functional group (default ['hydroxyl']) or CustomChemGraph object
    res_format : str
        result format (default 'original')

    Returns
    -------
    res : dict
        a list of all count
    compound : Compound
        compound object (sdf file)
    '''
    try:
        # create compound
        comp = compound(file)

        # create graph
        comp.create_graph()
        # check functional group
        res = comp.check_functional_groups(
            functional_groups, count_functional_group=True)
        # check
        if res_format == 'dataframe':
            # dataframe
            df = pd.DataFrame(res)
            return df, comp
        elif res_format == 'original':
            # raw
            return res, comp
        else:
            raise Exception("res_format is not valid.")

    except Exception as e:
        raise Exception(f"counting functional group is failed! {e}")


def create_custom_functional_groups(functional_groups: Union[Dict[str, List[str]], List[Dict[str, List[str]]], Path, str]) -> CustomChemGraph:
    '''
    Creates custom functional groups based on the following example.

    Parameters
    ----------
    functional_groups : list[dict] or dict
        functional group

    Returns
    -------
    custom_functional_groups : list[dict]
        a list of all custom functional group

    Examples
    --------
    ```python
    # fg1: CH2-O
    # fg2: CH2CHO
    # List of custom functional groups
    custom_functional_group = [
        {'fg1': ["C1-H1","C1-H2","C1-O1"]},
        {'fg2': ["C1-H1","C1-H2","C1-C2","C2-H3","C2-O2"]}
    ]

    # Dict of custom functional groups
    custom_functional_group = {
        'fg1': ["C1-H1","C1-H2","C1-O1"],
        'fg2': ["C1-H1","C1-H2","C1-C2","C2-H3","C2-O2"]
    }
    ```
    '''
    try:
        # check format
        if isinstance(functional_groups, dict):
            # set a list
            if len(functional_groups) == 1:
                custom_functional_groups = [functional_groups]
            else:
                custom_functional_groups = [
                    {k: v} for k, v in functional_groups.items()
                ]
        elif isinstance(functional_groups, list):
            # check is a list of dict
            if all(isinstance(item, dict) for item in functional_groups):
                # set
                custom_functional_groups = functional_groups
        elif isinstance(functional_groups, str):
            # check is a file yml
            custom_functional_groups = Utility.load_custom_functional_group(
                functional_groups)
        else:
            raise Exception("functional_groups format is not valid.")

        # custom chem graph
        CustomChemGraphC = CustomChemGraph(custom_functional_groups)
        # res
        return CustomChemGraphC
    except Exception as e:
        raise Exception(f'creating custom functional group is failed! {e}')
