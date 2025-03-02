# MOL GRAPH
# -----------------

# import packages/modules
import networkx as nx
import matplotlib.pyplot as plt
from typing import List, Union, Tuple, Dict, Any


class MolGraph():
    # cls variables
    # graph
    __graph: Any = None

    def __init__(self, custom_molecule: Dict[str, List[str]]):
        """_summary_

        Args:
            custom_molecule (Dict[str, List[str]]): custom molecule
        """
        self.custom_molecule = custom_molecule
        # molecule name
        self.molecule_name = list(custom_molecule.keys())[0]
        # molecule bonds
        self.molecule_bonds = list(custom_molecule.values())[0]
        # create graph
        self.__graph = self.generate_one_graph()

    @property
    def graph(self):
        return self.__graph

    def generate_one_graph(self) -> nx.Graph:
        '''
        Create a graph for a single molecule or fragment

        Parameters
        ----------
        custom_molecule : dict
            dictionary of a custom molecule

        Returns
        -------
        G : nx.Graph
            graph

        Examples
        --------
        >>> custom_molecule = {'molecule1': ["C1-H1","C1-H2","C1-O1"]}
        '''
        # group name and number of elements in the group
        custom_group_details = {list(self.custom_molecule.keys())[0]: len(
            list(self.custom_molecule.values())[0])}

        # group and subgroup collection
        g_sub_collection = {str(k).strip(): []
                            for k in custom_group_details.keys()}

        # looping through the custom functional group
        for k, v in custom_group_details.items():
            # group name
            group_name_ = str(k).strip()

            # looping through the custom functional group
            for fg in [self.custom_molecule]:
                # group name
                group_name__ = str(list(fg.keys())[0]).strip()

                # group values
                group_values__ = list(fg.values())[0]

                # looping through the group values
                for group_value___ in group_values__:
                    # check group_name_ in group_value_
                    if group_name_.lower() == group_value___.lower():
                        g_sub_collection[group_name__].append(group_name_)

        # define graph
        G = nx.Graph()

        # looping through the custom group
        for key, bonds in self.custom_molecule.items():

            # check whether sub-group exists
            if len(g_sub_collection[key]) == 0:

                for bond in bonds:
                    # check bond type
                    if '-' in bond:
                        atoms = bond.split('-')
                    elif '=' in bond:
                        atoms = bond.split('=')
                    elif '#' in bond:
                        atoms = bond.split('#')
                    else:
                        raise Exception('bond type error')

                    # check atoms
                    # Substitute empty parts with 'X'
                    # atoms = ["X0" if atom == "" else atom for atom in atoms]

                    # get letters
                    atoms_result = [{"letters": ''.join([c for c in s if c.isalpha(
                    )]), "numbers": ''.join([c for c in s if c.isdigit()])} for s in atoms]
                    # update atoms
                    atoms_name = [str(i['letters']) for i in atoms_result]
                    # atom id
                    atoms = [str(i['numbers']) for i in atoms_result]

                    # check
                    if len(atoms) == 2:
                        G.add_node(atoms[0], symbol=atoms_name[0])
                        G.add_node(atoms[1], symbol=atoms_name[1])
                        # define symbol
                        bond_id = str(atoms_name[0]).strip(
                        )+str(atoms_name[1]).strip()
                        # check bond type
                        if '=' in bond:
                            G.add_edge(atoms[0], atoms[1],
                                       symbol=bond_id, type=2)
                        elif '#' in bond:
                            G.add_edge(atoms[0], atoms[1],
                                       symbol=bond_id, type=3)
                        else:
                            G.add_edge(atoms[0], atoms[1],
                                       symbol=bond_id, type=1)

        # add a name to the graph
        G.graph['name'] = list(self.custom_molecule.keys())[0]

        # res
        return G

    def d(self):
        '''
        plots the molecule graph

        Returns
        -------
        None
        '''
        try:
            # get value from dict
            G = self.__graph
            pos = nx.spring_layout(G)
            # draw
            nx.draw_networkx(G, pos, with_labels=True, node_color='lightblue')
            labels = nx.get_edge_attributes(G, 'symbol')
            nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
            plt.show()
        except Exception as e:
            raise Exception(f"Creating graph failed: {e}")
