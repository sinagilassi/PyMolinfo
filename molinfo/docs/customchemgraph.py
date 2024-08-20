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
        for fg in custom_functional_group:
            for key, bonds in fg.items():
                G = nx.Graph()
                for bond in bonds:
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
