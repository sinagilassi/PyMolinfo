# NETWORK
# -------

# import packages/modules
import pandas as pd
import networkx as nx
from networkx.algorithms import isomorphism
# local
from .chemgraphs import ChemGraphs
from .customchemgraph import CustomChemGraph


class Network(ChemGraphs):

    # functional groups
    _functional_groups = []
    _custom_functional_groups = []
    _custom_functional_group_list = {}
    # graph
    _compound_graph = None

    def __init__(self, atomElements, atomBonds, xyzList, xyzCenterList, atomBonds1d):
        self.atomElements = atomElements
        # bond block (info)
        self.atomBonds = atomBonds
        self.xyzList = xyzList
        self.xyzCenterList = xyzCenterList

        # 1d vector of atom bonds
        self.atomBonds1d = atomBonds1d

        # TODO: super
        ChemGraphs.__init__(self)

        # functional group list
        self.function_group_list = {
            'hydroxyl': [self.graph_hydroxyl()],
            'carbonyl': [self.graph_carbonyl()],
            'carboxyl': [self.graph_carboxyl()],
            'N-H': [self.graph_N_H()],
            'C-N': [self.graph_C_N_single_bond()],
            'N-O': [self.graph_N_O_single_bond()],
            'C-O': [self.graph_C_O_single_bond()],
            'C#N': [self.graph_C_N_triple_bond()],
            'methyl-group': [self.graph_methyl()],
            'methylene-group': [self.graph_methylene()],
            'methine-group': [self.graph_methine()],
            'ether': [self.graph_ether()],
            'pst-alcohols': [self.graph_primary_alcohol(), self.graph_secondary_alcohol(),
                             self.graph_tertiary_alcohol(), self.graph_secondary_alcohol_double_bond(),
                             self.graph_primary_alcohol_double_bond()],
            'primary-alcohol': [self.graph_primary_alcohol(), self.graph_primary_alcohol_double_bond()],
            'secondary-alcohol': [self.graph_secondary_alcohol(), self.graph_secondary_alcohol_double_bond()],
            'tertiary-alcohol': [self.graph_tertiary_alcohol()],
            'alkane': [self.graph_alkane()],
            'alkane CH bond (sp3)': [self.graph_alkane_CH_bond()],
            'alkene': [self.graph_alkene()],
            'alkene CH bond (sp2)': [self.graph_alkene_CH_bond()],
            'aromatic CH bond (sp2)': [self.graph_aromatic_CH_bond()],
            'aldehyde CH bond (sp2)': [self.graph_aldehyde_CH_bond()],
            'alkyne': [self.graph_alkyne()],
            'alkyne CH bond (sp)': [self.graph_alkyne_CH_bond()],
            'arene': [self.graph_arene()],
            'aldehyde': [self.graph_aldehyde()],
            'ketone': [self.graph_ketone()],
            'carboxylic-acid': [self.graph_carboxylic_acid()],
            'ester': [self.graph_ester()],
            'pst-amide': [self.graph_primary_amide(), self.graph_secondary_amide(), self.graph_tertiary_amide()],
            'primary-amide': [self.graph_primary_amide()],
            'secondary-amide': [self.graph_secondary_amide()],
            'tertiary-amide': [self.graph_tertiary_amide()],
            'pst-amine': [self.graph_primary_amine(), self.graph_secondary_amine(), self.graph_tertiary_amine()],
            'primary-amine': [self.graph_primary_amine()],
            'secondary-amine': [self.graph_secondary_amine()],
            'tertiary-amine': [self.graph_tertiary_amine()],
            'nitrile': [self.graph_nitrile()],
            'thiol': [self.graph_thiol()],
            'alkyl-halids': [self.graph_alkyl_halide('F'), self.graph_alkyl_halide('Cl'), self.graph_alkyl_halide('Br'),
                             self.graph_alkyl_halide('I'),
                             self.graph_primary_alkyl_halide(
                                 'F'), self.graph_primary_alkyl_halide('Cl'),
                             self.graph_primary_alkyl_halide('Br'), self.graph_primary_alkyl_halide('I')],
            'epoxide': [self.graph_epoxide()]
        }

        # update functional groups
        self.functional_groups = [i for i in self.function_group_list.keys()]

    # property
    @property
    def functional_groups(self):
        return self._functional_groups

    # setter
    @functional_groups.setter
    def functional_groups(self, value):
        self._functional_groups = []
        self._functional_groups = [*value]

    @property
    def custom_functional_groups(self):
        return self._custom_functional_groups

    @custom_functional_groups.setter
    def custom_functional_groups(self, value):
        self._custom_functional_groups = []
        self._custom_functional_groups = [*value]

    @property
    def custom_functional_group_list(self):
        return self._custom_functional_group_list

    @custom_functional_group_list.setter
    def custom_functional_group_list(self, value):
        self._custom_functional_group_list = value

    @property
    def compound_graph(self):
        return self._compound_graph

    @compound_graph.setter
    def compound_graph(self, value):
        self._compound_graph = value

    def check_functional_groups(self, functional_groups=[],
                                count_functional_group=False):
        '''
        Check functional groups in a compound

        Parameters
        ----------
        functional_groups : list
            list of functional groups

        Returns
        -------
        res : dict
            a list of all count
        '''
        # check functional group
        if len(functional_groups) == 0:
            functional_groups = list(self.function_group_list.keys())

        # create graph
        G = self.create_graph()

        # check
        if count_functional_group:
            res = self.count_functional_group(G, functional_groups)
        else:
            # check functional groups
            res = self.check_functional_group(G, functional_groups)
        # res
        return res

    def check_functional_group(self, G, function_groups):
        '''
        Check a functional group exists in a compound

        Parameters
        ----------
        G : graph
            graph
        function_groups : list[str]
            functional group name like hydroxyl

        Returns
        -------
        res : dict
            a list of all count
        '''
        def node_match(n1, n2):
            '''
            Define a custom node_match function to ensure element matching
            '''
            # Define dummy node pattern
            dummy_node_pattern = "XX"

            # Ignore dummy nodes
            if n1['symbol'].startswith(dummy_node_pattern) or n2['symbol'].startswith(dummy_node_pattern):
                return True

            return n1['symbol'] == n2['symbol']

        def edge_match(e1, e2):
            '''
            Define a custom edge_match function to ensure bond type matching
            '''
            return e1['type'] == e2['type']

        # res
        res_match = []

        # NOTE: for each functional group
        for item in function_groups:
            # * Check instance
            if isinstance(item, str):
                # ! check functional exists in the list
                if item in self.function_group_list:
                    # get a list of graphs
                    # REVIEW
                    function_group_graphs = self.function_group_list[item]

                    # flag to track if functional group is found
                    fg_found_any = False

                    # graph function list
                    for _fn in function_group_graphs:
                        # Create a GraphMatcher object for a functional group
                        fg_matcher = isomorphism.GraphMatcher(
                            G, _fn, node_match=node_match,
                            edge_match=edge_match)

                        # Check if a functional group is in the main graph
                        fg_found = fg_matcher.subgraph_is_isomorphic()

                        if fg_found:
                            # print(f"{item} found in the molecule!")
                            fg_found_any = True
                            # res
                            res_match.append({
                                'function_group': item,
                                'result': fg_found
                            })
                            # break
                            break
                    # check
                    if not fg_found_any:
                        # print(f"{item} not found in the molecule.")
                        # res
                        res_match.append({
                            'function_group': item,
                            'result': False
                        })
            elif isinstance(item, CustomChemGraph):
                # get a list of graphs
                for custom_functional_group in item.custom_functional_groups:
                    # check dict
                    if len(custom_functional_group) == 1:
                        # get key (custom functional group name)
                        key = list(custom_functional_group.keys())[0]
                        # get value (graph)
                        value = custom_functional_group[key]
                        print(type(value))

                        # ANCHOR: check value type
                        if isinstance(value, nx.Graph):

                            # update custom functional group
                            self.update_custom_functional_group(key, value)

                            # Create a GraphMatcher object for a functional group
                            fg_matcher = isomorphism.GraphMatcher(
                                G, value, node_match=node_match, edge_match=edge_match)

                            # Check if a functional group is in the main graph (bool)
                            fg_found = fg_matcher.subgraph_is_isomorphic()
                            # print(f"{key} found in the molecule!")

                            # res
                            res_match.append({
                                'function_group': key,
                                'result': fg_found
                            })

                        elif isinstance(value, list):
                            # SECTION: list of graphs {already created}

                            # graphs set
                            graphs_set_ = self.__look_for_group_subgroup(value)

                            # NOTE: start matching
        # REVIEW
        # get a list of functional group names
        function_group_names = [x['function_group'] for x in res_match]
        # update list
        self.functional_groups = function_group_names

        # res
        return res_match

    def count_functional_group(self, G, function_groups):
        '''
        Count the occurrences of functional groups within the structure of a compound.

        Parameters
        ----------
        G : graph
            graph
        function_groups : list[str]
            functional group name like hydroxyl

        Returns
        -------
        res : dict
            a list of all count
        '''
        def node_match(n1, n2):
            '''
            Define a custom node_match function to ensure element matching
            '''
            # Define dummy node pattern
            dummy_node_pattern = "XX"

            # Ignore dummy nodes
            if n1['symbol'].startswith(dummy_node_pattern) or n2['symbol'].startswith(dummy_node_pattern):
                return True

            return n1['symbol'] == n2['symbol']

        def edge_match(e1, e2):
            '''
            Define a custom edge_match function to ensure bond type matching
            '''
            return e1['type'] == e2['type']

        # res
        res_match = []

        # for each functional group
        for item in function_groups:
            # * Check instance
            if isinstance(item, str):
                # ! check functional exists in the list
                if item in self.function_group_list:
                    # get a list of graphs
                    # REVIEW
                    function_group_graphs = self.function_group_list[item]

                    # flag to track if functional group is found
                    fg_found_any = False
                    # count of functional group
                    fg_count = 0
                    seen_subgraphs = set()

                    # graph function list
                    for _fn in function_group_graphs:
                        # Create a GraphMatcher object for a functional group
                        fg_matcher = isomorphism.GraphMatcher(
                            G, _fn, node_match=node_match, edge_match=edge_match)

                        # Check if a functional group is in the main graph
                        # fg_found = fg_matcher.subgraph_is_isomorphic()

                        # ! Check if a functional group is in the main graph
                        for subgraph in fg_matcher.subgraph_isomorphisms_iter():
                            # Convert the subgraph to a canonical form
                            canonical_subgraph = tuple(
                                sorted(subgraph.keys()))

                            # Check if the subgraph has been seen before
                            if canonical_subgraph not in seen_subgraphs:
                                seen_subgraphs.add(canonical_subgraph)
                                fg_found_any = True
                                fg_count += 1

                    # check
                    if fg_found_any:
                        # print(f"{item} found in the molecule!")
                        res_match.append({
                            'function_group': item,
                            'result': True,
                            'count': fg_count
                        })
                    else:
                        # print(f"{item} not found in the molecule.")
                        res_match.append({
                            'function_group': item,
                            'result': False,
                            'count': 0
                        })

            elif isinstance(item, CustomChemGraph):  # ! custom functional groups
                # get a list of graphs
                for custom_functional_group in item.custom_functional_groups:
                    # check dict
                    if len(custom_functional_group) == 1:

                        # flag to track if functional group is found
                        fg_found_any = False
                        # count of functional group
                        fg_count = 0
                        seen_subgraphs = set()

                        # get key (custom functional group name)
                        key = list(custom_functional_group.keys())[0]
                        # get value (graph)
                        value = custom_functional_group[key]

                        # update custom functional group
                        self.update_custom_functional_group(key, value)

                        # Create a GraphMatcher object for a functional group
                        fg_matcher = isomorphism.GraphMatcher(
                            G, value, node_match=node_match, edge_match=edge_match)

                        # Check if a functional group is in the main graph
                        # fg_found = fg_matcher.subgraph_is_isomorphic()
                        # print(f"{key} found in the molecule!")

                        # ! Check if a functional group is in the main graph
                        for subgraph in fg_matcher.subgraph_isomorphisms_iter():
                            # Convert the subgraph to a canonical form
                            canonical_subgraph = tuple(
                                sorted(subgraph.keys()))

                            # Check if the subgraph has been seen before
                            if canonical_subgraph not in seen_subgraphs:
                                seen_subgraphs.add(canonical_subgraph)
                                fg_found_any = True
                                fg_count += 1

                        # check
                        if fg_found_any:
                            # print(f"{item} found in the molecule!")
                            res_match.append({
                                'function_group': key,
                                'result': True,
                                'count': fg_count
                            })
                        else:
                            # print(f"{item} not found in the molecule.")
                            res_match.append({
                                'function_group': key,
                                'result': False,
                                'count': 0
                            })

        # get a list of functional group names
        function_group_names = [x['function_group'] for x in res_match]
        # update list
        self.functional_groups = function_group_names

        # res
        return res_match

    def create_graph(self):
        '''
        Check functional groups in a compound

        Parameters
        ----------
        functional_groups : list
            list of functional groups

        Returns
        -------
        res : dict
            a list of all count
        '''
        # atom no
        atomNo = len(self.xyzList)
        # bond no
        bondNo = len(self.atomBonds1d)

        # Create a graph from atoms and bonds
        G = nx.Graph()

        # *** atom visualization
        for i in range(atomNo):
            # xyz
            _atom1X = self.xyzList[i, 0]
            _atom1Y = self.xyzList[i, 1]
            _atom1Z = self.xyzList[i, 2]
            _atom1XYZ = [_atom1X, _atom1Y, _atom1Z]

            # color
            # atom id
            _atomId = int(i+1)
            # symbol
            _atomSymbol = str(self.atomElements[i]).strip()
            # atom mark
            atomMark = str(_atomSymbol) + str(_atomId)

            # add node with cooridination
            # G.add_node(_atomId, symbol=_atomSymbol,
            #            x=_atom1X, y=_atom1Y, z=_atom1Z, xyz=_atom1XYZ)

            G.add_node(_atomId, symbol=_atomSymbol)

        # *** using bond block

        for i in range(bondNo):
            # atom 1
            _id1 = self.atomBonds1d[i]['id1']
            _symbol1 = self.atomBonds1d[i]['symbol1']
            # atom 2
            _id2 = self.atomBonds1d[i]['id2']
            _symbol2 = self.atomBonds1d[i]['symbol2']
            # bond type
            _bondType = self.atomBonds1d[i]['bond_type']
            # bond symbol
            _bondSymbol = self.atomBonds1d[i]['bond_symbol']

            # add edge
            G.add_edge(_id1, _id2, symbol=_bondSymbol, type=_bondType)

        # update
        self.compound_graph = G

        # res
        return G

    def node_match(self, n1, n2):
        '''
        Define a custom node_match function to ensure element matching
        '''
        # Define dummy node pattern
        dummy_node_pattern = "XX"

        # Ignore dummy nodes
        if n1['symbol'].startswith(dummy_node_pattern) or n2['symbol'].startswith(dummy_node_pattern):
            return True

        return n1['symbol'] == n2['symbol']

    def edge_match(self, e1, e2):
        '''
        Define a custom edge_match function to ensure bond type matching
        '''
        return e1['type'] == e2['type']

    def search_within_main_graph(self, functional_group):
        '''
        Search a small graph within a larger graph

        Parameters
        ----------
        functional_group : str
            functional group
        '''

        # get the graph
        G = self.compound_graph
        # sub graph
        sub_graphs = []
        sub_graphs_name = []
        # res
        res_match = []
        fg_count = 0

        # search graph in functional_groups
        for item in self.functional_groups:
            if functional_group == item:
                sub_graphs_name = [
                    i for i in self.function_group_list.keys() if i == item]
                if len(sub_graphs_name) == 1:
                    sub_graphs = self.function_group_list[str(item)]
        # search graph in custom functional_groups
        for item in self.custom_functional_groups:
            if functional_group == item:
                sub_graphs_name = [
                    i for i in self.custom_functional_group_list.keys() if i == item]
                if len(sub_graphs_name) == 1:
                    sub_graphs = self.custom_functional_group_list[str(item)]

        # check
        if len(sub_graphs) == 0:
            raise Exception('functional group not found!')

        # Create a GraphMatcher object for a functional group
        for _fn in sub_graphs:
            fg_matcher = isomorphism.GraphMatcher(
                G, _fn, node_match=self.node_match, edge_match=self.edge_match)

            # Check if a functional group is in the main graph
            for subgraph in fg_matcher.subgraph_isomorphisms_iter():
                fg_count += 1

                # Store the subgraph
                # subgraph_nodes = subgraph.items()
                # pattern subgraph
                subgraph_pattern = G.subgraph(subgraph.keys())

                # Add the subgraph and node IDs to the result
                res_match.append({
                    'function_group': _fn,
                    'result': True,
                    'count': fg_count,
                    'subgraph': subgraph,
                    'subgraph_pattern': subgraph_pattern
                })

        return res_match

    def update_custom_functional_group(self, key, value):
        '''
        Update custom functional group
        '''
        # update custom functional group
        # check key exist
        if key not in self.custom_functional_groups:
            # save
            self.custom_functional_groups.append(key)
            # dict (key: graph)
            self.custom_functional_group_list[key] = [value]

    # NOTE: search deeper within the main graph
    def __look_for_group_subgroup(self, value):
        """
        Look for a group and subgroup within the custom functional groups.

        Parameters
        ----------
        value : list
            list of graphs in a custom functional group (group and subgroup)

        Returns
        -------
        graphs_set_ : dict
            a dictionary of group and subgroup graphs
        """
        # check list (2 elements)
        if len(value) != 2:
            raise Exception(
                'custom functional group must have 2 elements')

        # elements
        graphs_ = {}

        # looping through value
        for v in value:
            # graph nodes number
            nodes_num = len(v.nodes)
            # graph edges number
            edges_num = len(v.edges)
            # graph name
            name = v.name
            # save
            graphs_[name] = {
                'nodes_num': nodes_num,
                'edges_num': edges_num,
                'graph': v
            }

        # set group and subgroup
        graphs_set_ = {
            'group': None,
            'subgroup': None
        }

        # reference
        ref_ = max([v['nodes_num']
                    for k, v in graphs_.items()])

        # looping through graphs_
        for k, v in graphs_.items():
            # check group
            if v['nodes_num'] >= ref_:
                graphs_set_['group'] = v['graph']
            else:
                graphs_set_['subgroup'] = v['graph']

        # check
        if graphs_set_['group'] is None:
            raise Exception('group not found!')

        # check
        if graphs_set_['subgroup'] is None:
            raise Exception('subgroup not found!')

        # return
        return graphs_set_
