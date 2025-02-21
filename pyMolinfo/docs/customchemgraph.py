# CUSTOM CHEM GRAPH
# -----------------

# import packages/modules
import networkx as nx
import matplotlib.pyplot as plt
from typing import List, Union, Tuple


class CustomChemGraph():
    # cls
    # __functional_groups = []
    # custom functional groups
    _custom_functional_groups = []

    # parent group

    # child group

    def __init__(self, functional_groups):
        self.functional_groups = functional_groups

        # create graph
        __functional_groups = self.create_custom_graph(
            self.functional_groups)

        # create list
        self._custom_functional_groups = self.list_functional_groups(
            __functional_groups)

    @property
    def custom_functional_groups(self):
        return self._custom_functional_groups

    def create_custom_graph_V1(self, custom_functional_group):
        '''
        create custom graph
        bond type: - single, = double, # triple

        Parameters
        ----------
        custom_functional_group : list
            list of bonds

        Returns
        -------
        G : nx.Graph
            graph

        Examples
        --------
        >>> CH3-CH2-CH2-O
        >>> "C1-H1","C1-H2","C1-H3","C1-C2","C2-H4","C2-H5","C2-C3","C3-H6","C3-H7","C3-O"
        '''
        G = nx.Graph()
        for bond in custom_functional_group:
            atoms = bond.split('-')
            if len(atoms) == 2:
                G.add_node(atoms[0], symbol=atoms[0])
                G.add_node(atoms[1], symbol=atoms[1])
            if '=' in bond:
                G.add_edge(atoms[0], atoms[1], symbol='=', type=2)
            elif '#' in bond:
                G.add_edge(atoms[0], atoms[1], symbol='#', type=3)
            else:
                G.add_edge(atoms[0], atoms[1], symbol='-', type=1)
        return G

    def create_custom_graph(self, custom_functional_group) -> List[Tuple[str, Union[nx.Graph, dict[str, str]]]]:
        '''
        Creates a custom graph for a given functional group

        Parameters
        ----------
        custom_functional_group : list[dict]
            list of custom functional group

        Returns
        -------
        G_list : list
            list of graphs

        Examples
        --------
        >>> custom_functional_group = [
        >>> {'fg1': ["C1-H1","C1-H2","C1-O1"]},
        >>> {'fg2': ["C1-H1","C1-H2","C1-C2","C2-H3","C2-O2"]},
        >>> {'fg3': ["-C1","C1-H2"]}
        >>>    ]
        '''
        # group name and number of elements in the group
        custom_group_details = {list(fg.keys())[0]: len(
            list(fg.values())[0]) for fg in custom_functional_group}

        # group and subgroup collection
        g_sub_collection = {str(k).strip(): []
                            for k in custom_group_details.keys()}

        # looping through the custom functional group
        for k, v in custom_group_details.items():
            # group name
            group_name_ = str(k).strip()

            # looping through the custom functional group
            for fg in custom_functional_group:
                # group name
                group_name__ = str(list(fg.keys())[0]).strip()

                # group values
                group_values__ = list(fg.values())[0]

                # looping through the group values
                for group_value___ in group_values__:
                    # check group_name_ in group_value_
                    if group_name_.lower() == group_value___.lower():
                        g_sub_collection[group_name__].append(group_name_)

        # graph list
        G_list = []
        # loop for each custom functional group
        for fg in custom_functional_group:
            # custom fg -> values
            # key: fg name
            # value: list of bonds
            # define graph
            G = nx.Graph()

            # looping through the custom functional group
            for key, bonds in fg.items():

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
                    # save
                    G_list.append((key, G))
                else:
                    # save
                    G_list.append((key, g_sub_collection[key]))

        # res
        return G_list

    def list_functional_groups(self, graphs):
        '''
        make a list of dictionary

        Parameters
        ----------
        graphs : list
            list of graphs

        Returns
        -------
        functional_groups : list
            list of functional groups
        '''
        # res
        custom_functional_groups = []
        try:
            for key, graph in graphs:
                # save
                custom_functional_groups.append({key: graph})

            # res
            return custom_functional_groups
        except Exception as e:
            raise Exception(e)

    def d(self, functional_group_name):
        '''
        plots a functional group

        Parameters
        ----------
        functional_group_name : str
            name of functional group

        Returns
        -------
        None
        '''
        # check
        if len(self._custom_functional_groups) == 0:
            raise Exception('No functional groups found in the compound')

        # find graph
        functional_group_name = str(functional_group_name).strip()
        Gs = [item for item in self._custom_functional_groups if functional_group_name in item.keys()]

        # check
        if len(Gs) == 1:
            # get value from dict
            G = Gs[0][functional_group_name]
            pos = nx.spring_layout(G)
            # draw
            nx.draw_networkx(G, pos, with_labels=True, node_color='lightblue')
            labels = nx.get_edge_attributes(G, 'symbol')
            nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
            plt.show()
