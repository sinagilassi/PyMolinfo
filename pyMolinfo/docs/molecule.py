# GRAPHER
# ---------------
# import libs
import re
from typing import List, Dict, Union, Tuple, Optional


class Molecule():
    '''
    Build a molecule from a source
    '''

    def __init__(self, molecule_src, molecule_name):
        '''
        Initialize the molecule source

        Parameters
        ----------
        molecule_src : Dict[str,List[str]]
            The source of the molecule
        molecule_name : str
            The name of the molecule
        '''
        self.molecule_src = molecule_src
        self.molecule_name = molecule_name

    def build(self):
        '''
        Build a molecule from a source
        '''
        try:
            # SECTION: check molecule
            molecule_src_checked = self.check_molecule(self.molecule_src)

            # SECTION: construct molecule
            chain_info, molecule, constructed_molecule = self.construct_molecule(
                molecule_src_checked)

            # create molecule
            constructed_molecule_: Dict[str, List[str]] = {}
            constructed_molecule_[self.molecule_name] = constructed_molecule

            # SECTION: return
            return chain_info, molecule, constructed_molecule, constructed_molecule_
        except Exception as e:
            raise Exception(
                f"An error occurred while building the molecule: {e}")

    def extract_highest_index(self, chain):
        """
        Extracts the highest numerical index from atom labels in the given chain.

        Parameters
        ----------
        chain : list
            A list of strings representing bonds between atoms in a chain.

        Returns
        -------
        int
            The highest numerical index found in the chain.
        """
        indices = []
        for bond in chain:
            # Extract numbers from atoms
            atoms = re.findall(r'[A-Za-z]+(\d+)', bond)
            if '{' in bond or '}' in bond:
                continue  # Ignore bonds containing {Chain..}
            indices.extend(map(int, atoms))
        return max(indices) if indices else 0

    def search_for_main_chain(self, molecule_src: Dict[str, List[str]]) -> str:
        """
        Searches for the main chain in the molecule source.

        Parameters
        ----------
        molecule_src : dict
            A dictionary containing lists of strings representing bonds between atoms in a molecule.

        Returns
        -------
        str
            The name of the main chain found in the molecule source.

        Raises
        ------
        Exception
            If the main chain is not found in the molecule source.
        """
        try:
            for key, chain in molecule_src.items():
                # create a pattern to match the main chain
                pattern1 = re.compile(r'\{.*?\}\*.*')
                pattern2 = re.compile(r'.*\*\{.*?\}')

                # check if the main chain is found
                if any(pattern1.match(bond) or pattern2.match(bond) for bond in chain):
                    return str(key)

            # raise an exception if the main chain is not found
            raise Exception("Main chain not found in the molecule source.")
        except Exception as e:
            raise Exception(
                f"An error occurred while searching for the main chain: {e}")

    def check_molecule(self, molecule_src: Dict[str, List[str]]) -> Union[Dict[str, List[str]], bool]:
        """
        Checks if the molecule source is valid.

        Parameters
        ----------
        molecule_src : dict
            A dictionary containing lists of strings representing bonds between atoms in a molecule.

        Returns
        -------
        bool
            True if the molecule source is valid, False otherwise.
        """
        try:
            # create a copy of the molecule source
            molecule = {key: chain.copy()
                        for key, chain in molecule_src.items()}

            # checked molecule source
            molecule_src_checked = {}

            # check if the main chain is found
            main_chain = self.search_for_main_chain(molecule)
            if not main_chain:
                raise Exception("Main chain not found in the molecule source.")

            # add the main chain to the checked molecule source
            molecule_src_checked[main_chain] = molecule[main_chain]

            # main chain bonds
            main_chain_bonds = molecule[main_chain]

            # chain counter
            chain_counter = 0
            # looping through the molecule source
            for key, chain in molecule.items():
                if key != main_chain:
                    # check if the chain is connected to the main chain
                    # SECTION: create pattern to match the main chain
                    pattern_gate = rf"([A-Za-z]+)(\d+)\*\{{({re.escape(key)})\}}"

                    # looping through the main chain bonds
                    for i, bond in enumerate(main_chain_bonds):
                        match = re.match(pattern_gate, bond)
                        # check
                        if match:
                            # update chain counter
                            chain_counter += 1
                            # extract the matched chain
                            atom, index, key_chain = match.groups()
                            # rename key
                            key = key_chain + str(chain_counter)
                            # add the chain to the checked molecule source
                            molecule_src_checked[key] = chain

                            # update the element in the main chain
                            main_chain_bonds[i] = atom+index+"*{"+key+"}"

                    # reset chain counter
                    chain_counter = 0

            # res
            return molecule_src_checked
        except Exception as e:
            raise Exception(
                f"An error occurred while checking the molecule: {e}")

    def construct_molecule(self, molecule_src: Dict[str, List[str]]):
        """
        Constructs the molecule from the given molecule source.

        Parameters
        ----------
        molecule_src : dict
            A dictionary containing lists of strings representing bonds between atoms in a molecule.

        Returns
        -------
        dict
            A dictionary containing lists of strings representing bonds between atoms in a molecule.
        """
        try:
            # create a copy of the molecule source
            molecule = {
                key: chain.copy() for key, chain in molecule_src.items()}
            # print(f'molecule: {molecule}')

            # search for the main chain
            main_chain = self.search_for_main_chain(molecule)
            # print(f'main_chain: {main_chain}')

            chain_info: Dict[str, Dict[str, List[str]]] = {}

            # reset
            highest_index = 0

            # check if the main chain is found
            if main_chain:
                # get the highest index from the main chain
                highest_index = self.extract_highest_index(
                    molecule[main_chain])
                # print(f'highest_index: {highest_index}')

                # update index of other chains
                for key, chain in molecule.items():
                    if key != main_chain:
                        # create chain info
                        chain_info[key] = {
                            'bonds': [],
                            'gate': []
                        }

                        # update the index of the chain
                        for i, bond in enumerate(chain):

                            # SECTION: define pattern to match bonds
                            pattern = r"([A-Za-z]+)(\d+)([-=#])([A-Za-z]+)(\d+)"
                            # match the pattern
                            match_bond = re.match(pattern, bond)
                            # augment the index
                            if match_bond:
                                # extract atoms and indices
                                atom1, index1, bond_order, atom2, index2 = match_bond.groups()
                                # update the index
                                index1 = str(int(index1) + highest_index)
                                index2 = str(int(index2) + highest_index)
                                # update the bond
                                molecule[key][i] = f"{atom1}{index1}{bond_order}{atom2}{index2}"
                                # update the chain info
                                chain_info[key]['bonds'].append(
                                    molecule[key][i])

                            # SECTION: gate pattern
                            pattern_gate = r"([A-Za-z]+)(\d+)([-=#])\*"
                            # match the pattern
                            match_gate = re.match(pattern_gate, bond)
                            # augment the index
                            if match_gate:
                                # extract atoms and indices
                                atom1, index1, bond_order = match_gate.groups()
                                # update the index
                                index1 = str(int(index1) + highest_index)
                                # update the bond
                                molecule[key][i] = f"{bond_order}{atom1}{index1}"
                                # update the chain info
                                chain_info[key]['gate'].append(
                                    molecule[key][i])

                            # SECTION: gate pattern
                            pattern_gate = r"\*([-=#])([A-Za-z]+)(\d+)"
                            # match the pattern
                            match_gate = re.match(pattern_gate, bond)
                            # augment the index
                            if match_gate:
                                # extract atoms and indices
                                bond_order, atom1, index1 = match_gate.groups()
                                # update the index
                                index1 = str(int(index1) + highest_index)
                                # update the bond
                                molecule[key][i] = f"{atom1}{index1}{bond_order}"
                                # update the chain info
                                chain_info[key]['gate'].append(
                                    molecule[key][i])

                        # update highest index
                        highest_index = int(index1)
                        # print(f"highest_index: {highest_index}")

            # combine the main chain and other chains
            constructed_molecule: List[str] = []
            # find the gate atoms
            for items in molecule[main_chain]:
                # SECTION: define pattern
                pattern_gate = r"([A-Za-z]+)(\d+)\*\{([A-Za-z0-9]+)\}"
                # match the pattern
                match_gate = re.match(pattern_gate, items)
                # extract the gate atoms
                if match_gate:
                    # extract atoms and indices
                    atom1, index1, gate = match_gate.groups()
                    # find element index in molecule['main_chain']
                    element_index = molecule[main_chain].index(items)
                    # update
                    # check start with letter or number
                    if chain_info[gate]['gate'][0].startswith(('-', '=', '#')):
                        molecule[main_chain][element_index] = f"{atom1}{index1}{chain_info[gate]['gate'][0]}"
                    else:
                        molecule[main_chain][element_index] = f"{chain_info[gate]['gate'][0]}{atom1}{index1}"

                pattern_gate = r"\{([A-Za-z0-9]+)\}\*([A-Za-z]+)(\d+)"
                # match the pattern
                match_gate = re.match(pattern_gate, items)
                # extract the gate atoms
                if match_gate:
                    # extract atoms and indices
                    gate, atom1, index1 = match_gate.groups()
                    # find element index in molecule['main_chain']
                    element_index = molecule[main_chain].index(items)
                    # update
                    # check start with letter or number
                    if chain_info[gate]['gate'][0].startswith(('-', '=', '#')):
                        molecule[main_chain][element_index] = f"{atom1}{index1}{chain_info[gate]['gate'][0]}"
                    else:
                        molecule[main_chain][element_index] = f"{chain_info[gate]['gate'][0]}{atom1}{index1}"

            # combine the main chain and other chains
            for key, chain in molecule.items():
                if key == main_chain:
                    constructed_molecule.extend(chain)

            # chain info
            for key, chain in chain_info.items():
                constructed_molecule.extend(chain['bonds'])

            return chain_info, molecule, constructed_molecule
        except Exception as e:
            raise Exception(
                f"An error occurred while constructing the molecule: {e}")
