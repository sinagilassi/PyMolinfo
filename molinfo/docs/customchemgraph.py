# CUSTOM CHEM GRAPH
# -----------------

# import packages/modules
import networkx as nx
import matplotlib.pyplot as plt


class CustomChemGraph():
    # cls
    # __functional_groups = []
    # custom functional groups
    _custom_functional_groups = []

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

    def create_custom_graph(self, custom_functional_group):
        '''
        create custom graph

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
        >>> {'fg2': ["C1-H1","C1-H2","C1-C2","C2-H3","C2-O2"]}
        >>>    ]
        '''
        G_list = []
        # loop for each custom functional group
        for fg in custom_functional_group:
            # custom fg -> values
            # key: fg name
            # value: list of bonds
            # define graph
            G = nx.Graph()
            for key, bonds in fg.items():
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
        plot a functional group

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
