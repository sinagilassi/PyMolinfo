# MOLECULE CONSTRUCTOR
# --------------------
# import libs
import re
from typing import List, Dict, Tuple, Optional, Literal


class MoleculeConstructor:
    """
    A class for constructing complete molecular structures from component chains.

    This class provides methods to analyze molecular chains, process connections
    between chains, and construct a complete molecule from its component parts.
    """

    def __init__(self, molecule_src: Dict[str, List[str]]):
        """
        Initialize with a molecule source dictionary.

        Parameters
        ----------
        molecule_src : dict
            A dictionary containing lists of strings representing bonds between atoms in a molecule.
        """
        self.molecule_src = molecule_src
        self.molecule = {key: chain.copy()
                         for key, chain in molecule_src.items()}
        self.main_chain = self._search_for_main_chain()
        if not self.main_chain:
            raise ValueError("Main chain not found in the molecule source.")
        self.highest_index = self._extract_highest_index(
            self.molecule[self.main_chain])
        self.chain_info = {}
        self.chain_analysis = {}
        self.constructed_molecule = []

    def _search_for_main_chain(self) -> str:
        """
        Searches for the main chain in the molecule source.

        Returns
        -------
        str
            The name of the main chain found in the molecule source.
        """
        try:
            for key, chain in self.molecule.items():
                pattern1 = re.compile(r'\{.*?\}\*.*')
                pattern2 = re.compile(r'.*\*\{.*?\}')
                if any(pattern1.match(bond) or pattern2.match(bond) for bond in chain):
                    return key

            # If no chain is found, raise an error
            raise ValueError("Main chain not found in the molecule source.")
        except Exception as e:
            raise ValueError(
                "Error occurred while searching for the main chain.") from e

    def _extract_highest_index(self, chain: List[str]) -> int:
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
        try:
            indices: List[int] = []
            for bond in chain:
                # Extract numbers from atoms
                atoms = re.findall(r'[A-Za-z]+(\d+)', bond)
                if '{' in bond or '}' in bond:
                    continue  # Ignore bonds containing {Chain..}
                indices.extend(map(int, atoms))
            return max(indices) if indices else 0
        except Exception as e:
            raise Exception(
                "Error occurred while extracting the highest index.") from e

    def _analyze_chain_types(self):
        """
        Analyze the types of chains in the molecule.
        """
        chain_types = {
            "1": 'branch',
            "2": 'ring',
            "3": 'bridge',
        }

        for key, chain in self.molecule.items():
            if key != self.main_chain:
                # Count gate num
                chain_gate_num = sum(item.count('*') for item in chain)

                # Save chain type
                self.chain_analysis[key] = chain_types[str(chain_gate_num)]

                # Initialize chain info
                self.chain_info[key] = {
                    'receiver': [],
                    'bonds': [],
                    'gate': [],
                    'type': chain_types[str(chain_gate_num)],
                    'connection-port': {}
                }

    def _process_main_chain_connections(self):
        """
        Process connections between the main chain and other chains.
        """
        for chain in self.molecule[self.main_chain]:
            # Pattern for atom*{chain}
            pattern1 = r"([A-Za-z]+)(\d+)\*\{([A-Za-z0-9]+)\}"
            match = re.match(pattern1, chain)
            if match:
                atom1, index1, gate = match.groups()
                self.chain_info[gate]['receiver'].append(f"{atom1}{index1}")
                continue

            # Pattern for {chain}*atom
            pattern2 = r"\{([A-Za-z0-9]+)\}\*([A-Za-z]+)(\d+)"
            match = re.match(pattern2, chain)
            if match:
                gate, atom1, index1 = match.groups()
                self.chain_info[gate]['receiver'].append(f"{atom1}{index1}")

    def _update_chain_indices(self):
        """
        Update indices of atoms in chains and process bonds.
        """
        for key, chain in self.molecule.items():
            if key != self.main_chain:
                last_index = self.highest_index

                for i, bond in enumerate(chain):
                    # Process normal bonds
                    pattern_bond = r"([A-Za-z]+)(\d+)([-=#])([A-Za-z]+)(\d+)"
                    match = re.match(pattern_bond, bond)
                    if match:
                        atom1, index1, bond_order, atom2, index2 = match.groups()
                        index1 = str(int(index1) + self.highest_index)
                        index2 = str(int(index2) + self.highest_index)
                        self.molecule[key][i] = f"{atom1}{index1}{bond_order}{atom2}{index2}"
                        self.chain_info[key]['bonds'].append(
                            self.molecule[key][i])
                        last_index = max(int(index1), int(index2))
                        continue

                    # Process atom-* gate
                    pattern_gate1 = r"([A-Za-z]+)(\d+)([-=#])(\*+)"
                    match = re.match(pattern_gate1, bond)
                    if match:
                        atom1, index1, bond_order, gate_port = match.groups()
                        index1 = str(int(index1) + self.highest_index)
                        self.molecule[key][i] = f"{atom1}{index1}{bond_order}"
                        self.chain_info[key]['gate'].append(
                            self.molecule[key][i])
                        self.chain_info[key]['connection-port'][f"{atom1}{index1}{bond_order}"] = {
                            'port': gate_port,
                            'bond': f"{atom1}{index1}",
                            'bond-type': bond_order,
                            'bond-gate': f"{atom1}{index1}{bond_order}",
                        }
                        last_index = int(index1)
                        continue

                    # Process *-atom gate
                    pattern_gate2 = r"(\*+)([-=#])([A-Za-z]+)(\d+)"
                    match = re.match(pattern_gate2, bond)
                    if match:
                        gate_port, bond_order, atom1, index1 = match.groups()
                        index1 = str(int(index1) + self.highest_index)
                        self.molecule[key][i] = f"{bond_order}{atom1}{index1}"
                        self.chain_info[key]['gate'].append(
                            self.molecule[key][i])
                        self.chain_info[key]['connection-port'][f"{bond_order}{atom1}{index1}"] = {
                            'port': gate_port,
                            'bond': f"{atom1}{index1}",
                            'bond-type': bond_order,
                            'bond-gate': f"{bond_order}{atom1}{index1}",
                        }
                        last_index = int(index1)

                self.highest_index = last_index

    def _process_multi_gate_connection(self, element_index: int, atom1: str, index1: int | str, gate: str, chain_type_: str):
        """
        Process connections with multiple gates.

        Parameters
        ----------
        element_index : int
            Index of the element in the main chain.
        atom1 : str
            The atom symbol.
        index1 : str
            The atom index.
        gate : str
            The gate identifier.
        chain_type_ : str
            Type of the chain.
        """
        if chain_type_ == 'ring':
            for m, gate_connection in enumerate(self.chain_info[gate]['gate']):
                if gate_connection.startswith(('-', '=', '#')):
                    _connection = f"{atom1}{index1}{gate_connection}"
                else:
                    _connection = f"{gate_connection}{atom1}{index1}"

                if m == 0:
                    self.molecule[self.main_chain][element_index] = _connection
                else:
                    self.molecule[self.main_chain].append(_connection)
        elif chain_type_ == 'bridge':
            receiver_index = self.chain_info[gate]['receiver'].index(
                f"{atom1}{index1}")
            _gate_in = self.chain_info[gate]['gate'][receiver_index]

            if _gate_in.startswith(('-', '=', '#')):
                self.molecule[self.main_chain][element_index] = f"{atom1}{index1}{_gate_in}"
            else:
                self.molecule[self.main_chain][element_index] = f"{_gate_in}{atom1}{index1}"

    def _process_gate_connections(self):
        """
        Process gate connections between chains.
        """
        for i, items in enumerate(self.molecule[self.main_chain]):
            # Process atom*{chain} pattern
            pattern1 = r"([A-Za-z]+)(\d+)\*\{([A-Za-z0-9]+)\}"
            match = re.match(pattern1, items)
            if match:
                atom1, index1, gate = match.groups()
                gate_num = len(self.chain_info[gate]['gate'])
                chain_type_ = self.chain_analysis[gate]

                if gate_num == 1:
                    gate_connection = self.chain_info[gate]['gate'][0]
                    if gate_connection.startswith(('-', '=', '#')):
                        self.molecule[self.main_chain][i] = f"{atom1}{index1}{gate_connection}"
                    else:
                        self.molecule[self.main_chain][i] = f"{gate_connection}{atom1}{index1}"
                else:
                    self._process_multi_gate_connection(
                        i, atom1, index1, gate, chain_type_)
                continue

            # Process {chain}*atom pattern
            pattern2 = r"\{([A-Za-z0-9]+)\}\*([A-Za-z]+)(\d+)"
            match = re.match(pattern2, items)
            if match:
                gate, atom1, index1 = match.groups()
                gate_num = len(self.chain_info[gate]['gate'])
                chain_type_ = self.chain_analysis[gate]

                if gate_num == 1:
                    gate_connection = self.chain_info[gate]['gate'][0]
                    if gate_connection.startswith(('-', '=', '#')):
                        self.molecule[self.main_chain][i] = f"{atom1}{index1}{gate_connection}"
                    else:
                        self.molecule[self.main_chain][i] = f"{gate_connection}{atom1}{index1}"
                else:
                    self._process_multi_gate_connection(
                        i, atom1, index1, gate, chain_type_)

    def _build_constructed_molecule(self):
        """
        Build the final constructed molecule.
        """
        self.constructed_molecule = []

        # Add main chain bonds
        self.constructed_molecule.extend(self.molecule[self.main_chain])

        # Add bonds from other chains
        for key, chain in self.chain_info.items():
            self.constructed_molecule.extend(chain['bonds'])

    def construct(self):
        """
        Construct the complete molecule from its components.

        Returns
        -------
        tuple
            A tuple containing chain_info, molecule, and constructed_molecule.
        """
        # Step 1: Analyze chain types
        self._analyze_chain_types()

        # Step 2: Process main chain connections
        self._process_main_chain_connections()

        # Step 3: Update chain indices
        self._update_chain_indices()

        # Step 4: Process gate connections
        self._process_gate_connections()

        # Step 5: Build the constructed molecule
        self._build_constructed_molecule()

        return self.chain_info, self.molecule, self.constructed_molecule

    @staticmethod
    def from_source(molecule_src: Dict[str, List[str]]) -> Tuple[Dict[str, Dict[str, list]], Dict[str, List[str]], List[str]]:
        """
        Construct a molecule directly from a source dictionary.

        This is a convenience method that creates a MoleculeConstructor instance
        and immediately constructs the molecule.

        Parameters
        ----------
        molecule_src : dict
            A dictionary containing lists of strings representing bonds between atoms in a molecule.

        Returns
        -------
        tuple
            A tuple containing chain_info, molecule, and constructed_molecule.
        """
        constructor = MoleculeConstructor(molecule_src)
        return constructor.construct()
