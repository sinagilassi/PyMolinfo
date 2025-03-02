# import packages/modules
import os
from pathlib import Path
from networkx import Graph
import pubchemquery as pcq
import pandas as pd
from typing import List, Dict, Union, Literal, Optional

# internal
from .config import packageName
from .config import packageShortName
from .config import __version__
from .config import __description__
from .config import __author__
from .docs import MolParser, Compound, CustomChemGraph, Utility, Molecule


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


def create_graph(file: Path | str, graph_name: Optional[str] = 'Graph') -> Graph:
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
            graph = compound.create_graph(graph_name=graph_name)
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
        elif isinstance(functional_groups, (str, Path)):
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


def generate_molecule(molecule_src: Union[Dict[str, List[str]], Path], molecule_name: Optional[str] = 'MainChain') -> Molecule:
    '''
    Generates a `molecule object` from custom functional groups.

    Parameters
    ----------
    molecule_src : Dict[str, List[str]] | str
        molecule source
    molecule_name : str
        molecule name (default 'MainChain')

    Returns
    -------
    molecule : Molecule
        molecule object

    Examples
    --------
    ```python
    # molecule source
    molecule_src = {
        'MainChain': ["C1-C2","C2-C3","C3*{Chain1}","C3-C4","C4*{Chain2}","C4-C5","C5-C6"],
        'Chain1': ["C1=C2","C2-C3","C3=*"],
        'Chain2' : ["*-C1","C1=C2","C2-XX3"]
    }
    ```

    Notes
    -----
    - `*` is a connection point
    - `{}` is a reference to another chain
    - C3*{Chain1} means C3 is connected to Chain1 from C3=*
    - Bond index starts from 1 in all chains
    - Bond type: single bond (-), double bond (=), triple bond (#)
    '''
    try:
        # check format
        if isinstance(molecule_src, (str, Path)):
            # check is a file yml
            _mol_src = Utility.load_custom_functional_group(
                molecule_src)

        # init molecule
        MoleculeC = Molecule(molecule_src, molecule_name)

        # res
        return MoleculeC
    except Exception as e:
        raise Exception(f'creating custom functional group is failed! {e}')


def view_graph(graph: Graph, fig_size: List[int] = [10, 8], title: str = 'Graph', 
               node_size: int = 800, font_size: int = 12,
               edge_width: float = 1.5, with_labels: bool = True):
    '''
    Visualize a NetworkX Graph

    Parameters
    ----------
    graph : Graph
        NetworkX Graph object
    fig_size : List[int]
        figure size as [width, height] (default [10, 8])
    title : str
        graph title (default 'Graph')
    node_size : int
        size of nodes (default 800)
    font_size : int
        size of labels (default 12)
    edge_width : float
        width of edges (default 1.5)
    with_labels : bool
        display node labels (default True)

    Returns
    -------
    None
        displays the graph visualization
    '''
    try:
        import matplotlib.pyplot as plt
        import networkx as nx
        
        # Set up the figure
        plt.figure(figsize=tuple(fig_size))
        
        # Get node attributes if they exist
        node_labels = {}
        node_colors = []
        
        for node in graph.nodes():
            # Get symbol if available, otherwise use node id
            if 'symbol' in graph.nodes[node]:
                node_labels[node] = graph.nodes[node]['symbol']
                # Set colors based on atom type if symbol exists
                if graph.nodes[node]['symbol'] == 'C':
                    node_colors.append('#333333')  # Carbon - dark gray
                elif graph.nodes[node]['symbol'] == 'O':
                    node_colors.append('#ff0000')  # Oxygen - red
                elif graph.nodes[node]['symbol'] == 'N':
                    node_colors.append('#0000ff')  # Nitrogen - blue
                elif graph.nodes[node]['symbol'] == 'H':
                    node_colors.append('#cccccc')  # Hydrogen - light gray
                elif graph.nodes[node]['symbol'] == 'S':
                    node_colors.append('#ffff00')  # Sulfur - yellow
                else:
                    node_colors.append('#1f78b4')  # Default - blue
            else:
                node_labels[node] = str(node)
                node_colors.append('#1f78b4')  # Default - blue
        
        # Get edge attributes to style different bond types
        edge_styles = []
        edge_labels = {}
        
        for u, v, data in graph.edges(data=True):
            if 'type' in data:
                bond_type = data['type']
                if bond_type == 1:  # Single bond
                    edge_styles.append('solid')
                elif bond_type == 2:  # Double bond
                    edge_styles.append('dashed')
                elif bond_type == 3:  # Triple bond
                    edge_styles.append('dotted')
                else:
                    edge_styles.append('solid')
                
                if 'symbol' in data:
                    edge_labels[(u, v)] = data['symbol']
            else:
                edge_styles.append('solid')
        
        # Create the layout
        pos = nx.spring_layout(graph)
        
        # Draw the graph
        nx.draw(graph, pos, 
                node_color=node_colors,
                node_size=node_size,
                font_size=font_size,
                width=edge_width,
                with_labels=with_labels,
                labels=node_labels)
        
        # Draw edge labels if they exist
        if edge_labels:
            nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels)
        
        plt.title(title)
        plt.axis('off')
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        raise Exception(f"Visualizing graph failed! {e}")


def view_graph_3d(graph: Graph, fig_size: List[int] = [10, 8, 8], title: str = 'Graph 3D', 
                 node_size: int = 100, font_size: int = 10, edge_width: float = 2.0,
                 with_labels: bool = True, bg_color: str = '#ffffff'):
    '''
    Visualize a NetworkX Graph in 3D using xyz node attributes if available

    Parameters
    ----------
    graph : Graph
        NetworkX Graph object
    fig_size : List[int]
        figure size as [width, height, depth] (default [10, 8, 8])
    title : str
        graph title (default 'Graph 3D')
    node_size : int
        size of nodes (default 100)
    font_size : int
        size of labels (default 10)
    edge_width : float
        width of edges (default 2.0)
    with_labels : bool
        display node labels (default True)
    bg_color : str
        background color (default '#ffffff')

    Returns
    -------
    None
        displays the 3D graph visualization
    '''
    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        import numpy as np
        
        # Check if nodes have xyz attributes
        has_xyz = all('x' in graph.nodes[n] and 'y' in graph.nodes[n] and 'z' in graph.nodes[n] for n in graph.nodes())
        
        # Set up the 3D figure
        fig = plt.figure(figsize=tuple(fig_size[:2]))
        ax = fig.add_subplot(111, projection='3d')
        
        # Set background color
        ax.set_facecolor(bg_color)
        fig.patch.set_facecolor(bg_color)
        
        # Prepare node positions
        node_xyz = {}
        if has_xyz:
            # Use existing xyz attributes
            for node in graph.nodes():
                node_xyz[node] = (graph.nodes[node]['x'], 
                                  graph.nodes[node]['y'], 
                                  graph.nodes[node]['z'])
        else:
            # Generate 3D layout using spring layout
            import networkx as nx
            pos_2d = nx.spring_layout(graph)
            
            # Add a random z-coordinate to make it 3D
            np.random.seed(42)  # For reproducibility
            for node in graph.nodes():
                x, y = pos_2d[node]
                z = np.random.uniform(0, 1)
                node_xyz[node] = (x, y, z)
        
        # Get node attributes for visualization
        node_labels = {}
        node_colors = []
        
        for node in graph.nodes():
            # Get symbol if available, otherwise use node id
            if 'symbol' in graph.nodes[node]:
                node_labels[node] = graph.nodes[node]['symbol']
                # Set colors based on atom type if symbol exists
                if graph.nodes[node]['symbol'] == 'C':
                    node_colors.append('#333333')  # Carbon - dark gray
                elif graph.nodes[node]['symbol'] == 'O':
                    node_colors.append('#ff0000')  # Oxygen - red
                elif graph.nodes[node]['symbol'] == 'N':
                    node_colors.append('#0000ff')  # Nitrogen - blue
                elif graph.nodes[node]['symbol'] == 'H':
                    node_colors.append('#cccccc')  # Hydrogen - light gray
                elif graph.nodes[node]['symbol'] == 'S':
                    node_colors.append('#ffff00')  # Sulfur - yellow
                else:
                    node_colors.append('#1f78b4')  # Default - blue
            else:
                node_labels[node] = str(node)
                node_colors.append('#1f78b4')  # Default - blue
        
        # Plot nodes
        xs = [node_xyz[node][0] for node in graph.nodes()]
        ys = [node_xyz[node][1] for node in graph.nodes()]
        zs = [node_xyz[node][2] for node in graph.nodes()]
        
        ax.scatter(xs, ys, zs, c=node_colors, s=node_size)
        
        # Plot edges with different styles based on bond type
        for u, v, data in graph.edges(data=True):
            x = [node_xyz[u][0], node_xyz[v][0]]
            y = [node_xyz[u][1], node_xyz[v][1]]
            z = [node_xyz[u][2], node_xyz[v][2]]
            
            linestyle = '-'  # Default solid line for single bonds
            
            if 'type' in data:
                bond_type = data['type']
                if bond_type == 2:  # Double bond
                    linestyle = '--'
                elif bond_type == 3:  # Triple bond
                    linestyle = ':'
            
            # Draw the edge
            ax.plot(x, y, z, color='black', linewidth=edge_width, linestyle=linestyle)
        
        # Add node labels if requested
        if with_labels:
            for node in graph.nodes():
                x, y, z = node_xyz[node]
                ax.text(x, y, z, node_labels[node], fontsize=font_size)
        
        # Set plot title and equal aspect ratio
        ax.set_title(title)
        
        # Set equal aspect ratio for the axes
        ax.set_box_aspect([1, 1, 1])
        
        # Show the plot
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        raise Exception(f"Visualizing 3D graph failed! {e}")


def view_graph_3d_plotly(graph: Graph, fig_size: List[int] = [800, 600], title: str = 'Graph 3D',
                        node_size: int = 10, font_size: int = 10, edge_width: float = 2.0,
                        with_labels: bool = True, bg_color: str = '#ffffff'):
    '''
    Visualize a NetworkX Graph in 3D using Plotly for interactive visualization

    Parameters
    ----------
    graph : Graph
        NetworkX Graph object
    fig_size : List[int]
        figure size in pixels as [width, height] (default [800, 600])
    title : str
        graph title (default 'Graph 3D')
    node_size : int
        size of nodes (default 10)
    font_size : int
        size of labels (default 10)
    edge_width : float
        width of edges (default 2.0)
    with_labels : bool
        display node labels (default True)
    bg_color : str
        background color (default '#ffffff')

    Returns
    -------
    None
        displays the interactive 3D graph visualization
    '''
    try:
        import plotly.graph_objects as go
        import numpy as np
        import networkx as nx
        
        # Check if nodes have xyz attributes
        has_xyz = all('x' in graph.nodes[n] and 'y' in graph.nodes[n] and 'z' in graph.nodes[n] for n in graph.nodes())
        
        # Prepare node positions
        node_xyz = {}
        if has_xyz:
            # Use existing xyz attributes
            for node in graph.nodes():
                node_xyz[node] = (graph.nodes[node]['x'], 
                                 graph.nodes[node]['y'], 
                                 graph.nodes[node]['z'])
        else:
            # Generate 3D layout using spring layout
            pos_2d = nx.spring_layout(graph)
            
            # Add a random z-coordinate to make it 3D
            np.random.seed(42)  # For reproducibility
            for node in graph.nodes():
                x, y = pos_2d[node]
                z = np.random.uniform(0, 1)
                node_xyz[node] = (x, y, z)
        
        # Get node attributes for visualization
        node_labels = {}
        node_colors = []
        hover_texts = []
        
        # Color mapping for atoms
        atom_colors = {
            'C': '#333333',  # Carbon - dark gray
            'O': '#ff0000',  # Oxygen - red
            'N': '#0000ff',  # Nitrogen - blue
            'H': '#cccccc',  # Hydrogen - light gray
            'S': '#ffff00',  # Sulfur - yellow
            'F': '#90ee90',  # Fluorine - light green
            'Cl': '#a0ff20', # Chlorine - lime green
            'Br': '#a52a2a', # Bromine - brown
            'I': '#800080',  # Iodine - purple
            'P': '#ffa500'   # Phosphorus - orange
        }
        
        for node in graph.nodes():
            hover_text = f"Node: {node}"
            # Get node attributes if available
            attrs = []
            for key, value in graph.nodes[node].items():
                attrs.append(f"{key}: {value}")
                hover_text += f"<br>{key}: {value}"
            
            hover_texts.append(hover_text)
            
            # Get symbol if available, otherwise use node id
            if 'symbol' in graph.nodes[node]:
                node_labels[node] = graph.nodes[node]['symbol']
                # Set colors based on atom type if symbol exists
                symbol = graph.nodes[node]['symbol']
                if symbol in atom_colors:
                    node_colors.append(atom_colors[symbol])
                else:
                    node_colors.append('#1f78b4')  # Default - blue
            else:
                node_labels[node] = str(node)
                node_colors.append('#1f78b4')  # Default - blue
        
        # Create a trace for nodes
        node_trace = go.Scatter3d(
            x=[node_xyz[node][0] for node in graph.nodes()],
            y=[node_xyz[node][1] for node in graph.nodes()],
            z=[node_xyz[node][2] for node in graph.nodes()],
            mode='markers' + ('+text' if with_labels else ''),
            marker=dict(
                size=node_size,
                color=node_colors,
                opacity=1.0,
                line=dict(width=1, color='#000000')
            ),
            text=[node_labels[node] if with_labels else "" for node in graph.nodes()],
            textposition="middle center",
            textfont=dict(size=font_size, color='#000000'),
            hoverinfo='text',
            hovertext=hover_texts
        )
        
        # Create traces for edges with different styles based on bond type
        edge_traces = []
        
        for u, v, data in graph.edges(data=True):
            # Get edge attributes
            bond_type = data.get('type', 1)  # Default to single bond
            bond_symbol = data.get('symbol', '')
            
            # Choose line style based on bond type
            line_dash = 'solid'
            if bond_type == 2:  # Double bond
                line_dash = 'dash'
            elif bond_type == 3:  # Triple bond
                line_dash = 'dot'
            
            # Create edge coordinates
            x_edge = [node_xyz[u][0], node_xyz[v][0], None]
            y_edge = [node_xyz[u][1], node_xyz[v][1], None]
            z_edge = [node_xyz[u][2], node_xyz[v][2], None]
            
            # Create hover text
            hover_text = f"Bond: {u} - {v}<br>Type: {bond_type}"
            if bond_symbol:
                hover_text += f"<br>Symbol: {bond_symbol}"
            
            # Create edge trace
            edge_trace = go.Scatter3d(
                x=x_edge,
                y=y_edge,
                z=z_edge,
                mode='lines',
                line=dict(
                    width=edge_width,
                    color='#000000',
                    dash=line_dash
                ),
                hoverinfo='text',
                hovertext=hover_text,
                showlegend=False
            )
            
            edge_traces.append(edge_trace)
        
        # Create figure
        fig = go.Figure(data=[node_trace] + edge_traces)
        
        # Update layout
        fig.update_layout(
            title=title,
            width=fig_size[0],
            height=fig_size[1],
            scene=dict(
                xaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
                yaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
                zaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
                bgcolor=bg_color
            ),
            margin=dict(l=0, r=0, t=40, b=0),
            hovermode='closest',
            paper_bgcolor=bg_color,
            plot_bgcolor=bg_color
        )
        
        # Add camera controls
        fig.update_layout(
            scene_camera=dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=1.5, y=1.5, z=1.5)
            ),
            scene_dragmode='orbit'
        )
        
        # Show the figure
        fig.show()
        
    except Exception as e:
        raise Exception(f"Visualizing 3D graph with Plotly failed! {e}")


# Optional: Add a function to detect if a graph is likely to be a molecular graph
def is_molecular_graph(graph: Graph) -> bool:
    '''
    Check if a graph is likely to represent a molecule based on node attributes

    Parameters
    ----------
    graph : Graph
        NetworkX Graph object

    Returns
    -------
    bool
        True if the graph appears to be a molecular graph, False otherwise
    '''
    # Check if nodes have 'symbol' attributes with chemical element symbols
    chemical_symbols = {'H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I'}
    
    # Check at least some nodes for chemical symbols
    for node in graph.nodes():
        if 'symbol' in graph.nodes[node]:
            if graph.nodes[node]['symbol'] in chemical_symbols:
                return True
    
    return False
