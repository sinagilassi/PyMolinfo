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


def main():
    # short description and version
    _des = f"{packageShortName} {__version__}: {__description__}"
    print(_des)


def g3d(file, fig_size=[], display_legend=True, display_atom_id=True, display_bond_length=True):
    '''
    3d graph of a compound

    Parameters
    ----------
    file : str
        molecule file format (sdf)
    display_legend : bool
        display legend (default True)

    Returns
    -------
    None
        display 3d graph
    '''
    # check file exists
    if os.path.exists(file):
        # parse file
        MolParserC = MolParser(file)
        compound_info = MolParserC.read_file()
        # compound
        compound = Compound(compound_info)
        # display 3d
        if len(fig_size) == 0:
            compound.view3d(display_legend=display_legend,
                            display_atom_id=display_atom_id, display_bond_length=display_bond_length)
        else:
            compound.view3d(fig_size=fig_size, display_legend=display_legend,
                            display_atom_id=display_atom_id, display_bond_length=display_bond_length)
    else:
        raise Exception("file path is not valid.")


def g3d_by_inchi(inchi, display_legend=True):
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
                compound.view3d(display_legend=display_legend)
            else:
                raise Exception("sdf is not found!")
        else:
            raise Exception("cid is not found!")
    else:
        raise Exception("inchi is not valid.")


def check_functional_group(file, functional_groups=[], res_format='dict'):
    '''
    Check a functional group exists in a compound

    Parameters
    ----------
    file : str
        molecule file format (sdf)
    functional_groups : list[str]
        functional group (default ['hydroxyl'])
    res_format : str
        result format (default 'dict')

    Returns
    -------
    res : dict
        a list of all count
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
        if res_format == 'dict':
            return res
        elif res_format == 'dataframe':
            # dataframe
            return pd.DataFrame(res)
    else:
        raise Exception("file path is not valid.")


if __name__ == "__main__":
    main()
