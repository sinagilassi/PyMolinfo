# DISPLAY A CHEMICAL STRUCTURE
# -----------------------------

# import libs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import plotly as py
import plotly.express as px
import plotly.graph_objects as go
import plotly.io
import math
# internal
from .observer import Observer


class graph3d():
    '''
    3D Visualizer of a compound

    hint:
        xyzList is selected for visualization
    '''

    # properties
    _structure_type = ''
    plotScale = []

    def __init__(self, atomElements, atomBonds, xyzList, xyzCenterList, robs, tetaNo, phiNo, limits, atom_bonds_1d):
        self.atomElements = atomElements
        # bond block (info)
        self.atomBonds = atomBonds
        self.xyzList = xyzList
        self.xyzCenterList = xyzCenterList
        self.robs = robs
        self.tetaNo = tetaNo
        self.phiNo = phiNo
        self.limits = limits
        self.atomBonds_1d = atom_bonds_1d

        # set structure type
        structureType, perpendicularAxis, perpendicularVector, XYZ0 = self.StructureAnalyzer()
        self.structure_type = structureType

    @property
    def structure_type(self):
        return self._structure_type

    @structure_type.setter
    def structure_type(self, value):
        self._structure_type = value

    def StructureAnalyzer(self):
        '''
        Check geometry structure to determine 2D/3D

        Parameters
        ----------
        xyzList: list
            point coordination

        Returns
        -------
        structureType: str
            2D or 3D
        '''
        try:
            # array
            xyzList = np.array(self.xyzList)
            # set
            X = xyzList[:, 0]
            Y = xyzList[:, 1]
            Z = xyzList[:, 2]

            # check plane structure (2D, 3D)
            X0 = np.abs(X).sum()
            Y0 = np.abs(Y).sum()
            Z0 = np.abs(Z).sum()
            XYZ0 = [X0, Y0, Z0]

            # perpendicular vector [x,y,z]
            perpendicularVector = [False, False, False]
            # set
            if X0 == 0:
                perpendicularVector[0] = True

            if Y0 == 0:
                perpendicularVector[1] = True

            if Z0 == 0:
                perpendicularVector[2] = True

            # status
            if True in perpendicularVector:
                structureType = '2D'
            else:
                structureType = '3D'

            # axis selection
            perpendicularAxis = []
            if perpendicularVector[0] is True:
                perpendicularAxis.append('X')
            if perpendicularVector[1] is True:
                perpendicularAxis.append('Y')
            if perpendicularVector[2] is True:
                perpendicularAxis.append('Z')

            return structureType, perpendicularAxis, perpendicularVector, XYZ0
        except Exception as e:
            raise Exception(e)

    def set_color(self, atom_symbol):
        '''
        Set a color for each compound
        taken from https://en.wikipedia.org/wiki/CPK_coloring

        Parameters
        ----------
        atom_symbol: str
            atom symbol

        Returns
        -------
        color: str
            atom color
        '''
        colors = {
            "H": '#ffffff',
            "C": '#BCBCBC',
            "N": '#0586f6',
            "O": '#f6052a',
            "F": '#2dd930',
            "Cl": '#2dd930',
            "Br": '#950e0e',
            "I": '#360e89',
            "He": '#3dbaf1',
            "Ne": '#3dbaf1',
            "Ar": '#3dbaf1',
            "Kr": '#3dbaf1',
            "Xe": '#3dbaf1',
            "P": '#f1a03d',
            "S": '#f1ef3d',
            "B": '#efc867',
            "Li": '#6b3ccb',
            "Na": '#6b3ccb',
            "K": '#6b3ccb',
            "Rb": '#6b3ccb',
            "Cs": '#6b3ccb',
            "Fr": '#6b3ccb',
            "Be": '#1c881e',
            "Mg": '#1c881e',
            "Ca": '#1c881e',
            "Sr": '#1c881e',
            "Ba": '#1c881e',
            "Ra": '#1c881e',
            "Ti": '#3d3e40',
            "Fe": '#a48620',
            "other": '#a729ba'
        }

        # check
        _color = colors.get(str(atom_symbol))
        if _color is None:
            return colors.get(str('other'))
        else:
            return _color

    def set_size(self, symbol, _sy=300, _s=1):
        '''
        Set atom size (spherical shape)

        Parameters
        ----------
        symbol: str
            atom symbol
        _s: int
            size

        Returns
        -------
        size: int
            size
        '''
        _sy = 300
        _sx = 0.50*_sy

        sizes = {
            "H": _s*_sx,
            "C": _s*_sy,
            "N": _s*_sy,
            "O": _s*_sy,
            "F": _s*_sy,
            "Cl": _s*_sx,
            "Br": _s*_sy,
            "I": _s*_sy,
            "He": _s*_sy,
            "Ne": _s*_sy,
            "Ar": _s*_sy,
            "Kr": _s*_sy,
            "Xe": _s*_sy,
            "P": _s*_sy,
            "S": _s*_sy,
            "B": _s*_sy,
            "Li": _s*_sy,
            "Na": _s*_sy,
            "K": _s*_sy,
            "Rb": _s*_sy,
            "Cs": _s*_sy,
            "Fr": _s*_sy,
            "Be": _s*_sy,
            "Mg": _s*_sy,
            "Ca": _s*_sy,
            "Sr": _s*_sy,
            "Ba": _s*_sy,
            "Ra": _s*_sy,
            "Ti": _s*_sy,
            "Fe": _s*_sy,
            "other": _s*_sy,
        }

        _size = sizes.get(symbol)
        if _size is None:
            _sizeSet = int(sizes.get('other'))
        else:
            _sizeSet = int(_size)

        return _sizeSet

    def create_3dframe(self):
        '''
        Create 3d frame dimension
        '''
        # max length
        # x
        xMin = np.min(self.xyzList[:, 0])
        xMax = np.max(self.xyzList[:, 0])
        xLen = np.abs(np.abs(xMax) - np.abs(xMin))
        # y
        yMin = np.min(self.xyzList[:, 1])
        yMax = np.max(self.xyzList[:, 1])
        yLen = np.abs(np.abs(yMax) - np.abs(yMin))
        # z
        zMin = np.min(self.xyzList[:, 2])
        zMax = np.max(self.xyzList[:, 2])
        zLen = np.abs(np.abs(zMax) - np.abs(zMin))
        # max
        xyzLenMax = np.max([xMax, yMax, zMax])
        xyzLenMin = np.min([xMin, yMin, zMin])
        xyzR = 1/xyzLenMax

        # res
        return xyzLenMax, xyzLenMin, xyzR, xLen, yLen, zLen

    def create_bond_line(self, xyz1, xyz2, bond_type, xyzL=[1, 1, 1], xyzR=0.15):
        '''
        Create bond line (single, double, tipple)

        Parameters
        ----------
        xyz1: list
            (x,y,z) point 1
        xyz2: list
            (x,y,z) point 2
        bond_type: int
            bond type (1,2,3)
        xyzL: list
            (x,y,z) length
        xyzR: int
            (x,y,z) radius

        Returns
        -------
        bondLines: list
            list of bond lines
        '''
        # bond lines
        bondLines = []
        # bond length
        bondLength = []

        xL, yL, zL = xyzR*np.array(xyzL)

        # Calculate center-to-center vector
        center_vector = np.array(xyz2) - np.array(xyz1)

        # Calculate center-to-center line
        center_line = [[xyz1[0], xyz2[0]], [
            xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]]]

        if bond_type == 1:
            # bond line
            bondLines.append(center_line)
            # bond length
            _bond_length_res = self.calculate_distance(xyz1, xyz2)
            bondLength.append(_bond_length_res)

        elif bond_type == 2:
            # parallel line
            # y increase
            # center to center line
            # Calculate perpendicular vector to center vector
            perp_vector = np.array([center_vector[1], - center_vector[0], 0])
            perp_vector /= np.linalg.norm(perp_vector)

            # Calculate offset vectors
            offset_vector1 = 0.15 * perp_vector
            offset_vector2 = -0.15 * perp_vector

            # Calculate parallel lines
            _l1 = [[xyz1[0]+offset_vector1[0], xyz2[0]+offset_vector1[0]], [
                xyz1[1]+offset_vector1[1], xyz2[1]+offset_vector1[1]], [xyz1[2]+offset_vector1[2], xyz2[2]+offset_vector1[2]]]
            _l2 = [[xyz1[0]+offset_vector2[0], xyz2[0]+offset_vector2[0]], [
                xyz1[1]+offset_vector2[1], xyz2[1]+offset_vector2[1]], [xyz1[2]+offset_vector2[2], xyz2[2]+offset_vector2[2]]]

            # lines
            _l1_xyz1 = [xyz1[0]+offset_vector1[0], xyz1[1] +
                        offset_vector1[1], xyz1[2]+offset_vector1[2]]
            _l1_xyz2 = [xyz2[0]+offset_vector1[0], xyz2[1] +
                        offset_vector1[1], xyz2[2]+offset_vector1[2]]

            _l2_xyz1 = [xyz1[0]+offset_vector2[0], xyz1[1] +
                        offset_vector2[1], xyz1[2]+offset_vector2[2]]
            _l2_xyz2 = [xyz2[0]+offset_vector2[0], xyz2[1] +
                        offset_vector2[1], xyz2[2]+offset_vector2[2]]

            # bond lines
            bondLines.append(_l1)
            bondLines.append(_l2)
            # bond length
            _bond_length_res = self.calculate_distance(_l1_xyz1, _l1_xyz2)
            bondLength.append(_bond_length_res)
            _bond_length_res = self.calculate_distance(_l2_xyz1, _l2_xyz2)
            bondLength.append(_bond_length_res)

        elif bond_type == 3:
            # parallel line
            # y increase
            # line
            # Calculate perpendicular vector to center vector
            perp_vector = np.array([center_vector[1], -center_vector[0], 0])
            perp_vector /= np.linalg.norm(perp_vector)

            # Calculate offset vectors
            offset_vector1 = 0.125 * perp_vector
            offset_vector2 = -0.125 * perp_vector

            # Calculate parallel lines
            _l1 = [[xyz1[0]+offset_vector1[0], xyz2[0]+offset_vector1[0]], [
                xyz1[1]+offset_vector1[1], xyz2[1]+offset_vector1[1]], [xyz1[2]+offset_vector1[2], xyz2[2]+offset_vector1[2]]]
            _l2 = [[xyz1[0], xyz2[0]], [
                xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]]]

            _l3 = [[xyz1[0]+offset_vector2[0], xyz2[0]+offset_vector2[0]], [
                xyz1[1]+offset_vector2[1], xyz2[1]+offset_vector2[1]], [xyz1[2]+offset_vector2[2], xyz2[2]+offset_vector2[2]]]

            # lines
            _l1_xyz1 = [xyz1[0]+offset_vector1[0], xyz1[1] +
                        offset_vector1[1], xyz1[2]+offset_vector1[2]]
            _l1_xyz2 = [xyz2[0]+offset_vector1[0], xyz2[1] +
                        offset_vector1[1], xyz2[2]+offset_vector1[2]]
            _l2_xyz1 = xyz1
            _l2_xyz2 = xyz2
            _l3_xyz1 = [xyz1[0]+offset_vector2[0], xyz1[1] +
                        offset_vector2[1], xyz1[2]+offset_vector2[2]]
            _l3_xyz2 = [xyz2[0]+offset_vector2[0], xyz2[1] +
                        offset_vector2[1], xyz2[2]+offset_vector2[2]]

            # bond lines
            bondLines.append(_l1)
            bondLines.append(_l2)
            bondLines.append(_l3)
            # bond length
            # bond length
            _bond_length_res = self.calculate_distance(_l1_xyz1, _l1_xyz2)
            bondLength.append(_bond_length_res)
            _bond_length_res = self.calculate_distance(_l2_xyz1, _l2_xyz2)
            bondLength.append(_bond_length_res)
            _bond_length_res = self.calculate_distance(_l3_xyz1, _l3_xyz2)
            bondLength.append(_bond_length_res)

        return bondLines, bond_type, bondLength

    def create_bond_line_V2(self, xyz1, xyz2, bond_type, xyzL=[1, 1, 1], xyzR=0.15):
        '''
        Create bond line (single, double, triple)

        Parameters
        ----------
        xyz1: list
            (x,y,z) point 1
        xyz2: list
            (x,y,z) point 2
        bond_type: int
            bond type (1,2,3)
        xyzL: list
            (x,y,z) length
        xyzR: int
            (x,y,z) radius

        Returns
        -------
        bondLines: list
            list of bond lines
        '''
        bondLines = []

        xL, yL, zL = xyzR*np.array(xyzL)

        # Calculate center-to-center vector
        center_vector = np.array(xyz2) - np.array(xyz1)

        # Calculate center-to-center line
        center_line = [[xyz1[0], xyz2[0]], [
            xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]]]

        # Calculate midpoint
        midpoint = [(xyz1[0] + xyz2[0]) / 2, (xyz1[1] +
                                              xyz2[1]) / 2, (xyz1[2] + xyz2[2]) / 2]

        if bond_type == 1:
            # Break single bond into two lines
            bondLines.append(
                [[xyz1[0], midpoint[0]], [xyz1[1], midpoint[1]], [xyz1[2], midpoint[2]]])
            bondLines.append(
                [[midpoint[0], xyz2[0]], [midpoint[1], xyz2[1]], [midpoint[2], xyz2[2]]])

        elif bond_type == 2:
            # parallel line
            # y increase
            # center to center line
            # Calculate perpendicular vector to center vector
            perp_vector = np.array([center_vector[1], - center_vector[0], 0])
            perp_vector /= np.linalg.norm(perp_vector)

            # Calculate offset vectors
            offset_vector1 = 0.15 * perp_vector
            offset_vector2 = -0.15 * perp_vector

            # Calculate parallel lines
            _l1 = [[xyz1[0]+offset_vector1[0], xyz2[0]+offset_vector1[0]], [
                xyz1[1]+offset_vector1[1], xyz2[1]+offset_vector1[1]], [xyz1[2]+offset_vector1[2], xyz2[2]+offset_vector1[2]]]
            _l2 = [[xyz1[0]+offset_vector2[0], xyz2[0]+offset_vector2[0]], [
                xyz1[1]+offset_vector2[1], xyz2[1]+offset_vector2[1]], [xyz1[2]+offset_vector2[2], xyz2[2]+offset_vector2[2]]]

            # Break each parallel line into two lines
            bondLines.append([[xyz1[0]+offset_vector1[0], midpoint[0]+offset_vector1[0]], [
                xyz1[1]+offset_vector1[1], midpoint[1]+offset_vector1[1]], [xyz1[2]+offset_vector1[2], midpoint[2]+offset_vector1[2]]])
            bondLines.append([[midpoint[0]+offset_vector1[0], xyz2[0]+offset_vector1[0]], [
                midpoint[1]+offset_vector1[1], xyz2[1]+offset_vector1[1]], [midpoint[2]+offset_vector1[2], xyz2[2]+offset_vector1[2]]])
            bondLines.append([[xyz1[0]+offset_vector2[0], midpoint[0]+offset_vector2[0]], [
                xyz1[1]+offset_vector2[1], midpoint[1]+offset_vector2[1]], [xyz1[2]+offset_vector2[2], midpoint[2]+offset_vector2[2]]])
            bondLines.append([[midpoint[0]+offset_vector2[0], xyz2[0]+offset_vector2[0]], [
                midpoint[1]+offset_vector2[1], xyz2[1]+offset_vector2[1]], [midpoint[2]+offset_vector2[2], xyz2[2]+offset_vector2[2]]])

        elif bond_type == 3:
            # parallel line
            # y increase
            # line
            # Calculate perpendicular vector to center vector
            perp_vector = np.array([center_vector[1], -center_vector[0], 0])
            perp_vector /= np.linalg.norm(perp_vector)

            # Calculate offset vectors
            offset_vector1 = 0.125 * perp_vector
            offset_vector2 = -0.125 * perp_vector

            # Calculate parallel lines
            _l1 = [[xyz1[0]+offset_vector1[0], xyz2[0]+offset_vector1[0]], [
                xyz1[1]+offset_vector1[1], xyz2[1]+offset_vector1[1]], [xyz1[2]+offset_vector1[2], xyz2[2]+offset_vector1[2]]]
            _l2 = [[xyz1[0], xyz2[0]], [
                xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]]]
            _l3 = [[xyz1[0]+offset_vector2[0], xyz2[0]+offset_vector2[0]], [
                xyz1[1]+offset_vector2[1], xyz2[1]+offset_vector2[1]], [xyz1[2]+offset_vector2[2], xyz2[2]+offset_vector2[2]]]

            # Break each parallel line into two lines
            bondLines.append([[xyz1[0]+offset_vector1[0], midpoint[0]+offset_vector1[0]], [
                xyz1[1]+offset_vector1[1], midpoint[1]+offset_vector1[1]], [xyz1[2]+offset_vector1[2], midpoint[2]+offset_vector1[2]]])
            bondLines.append([[midpoint[0]+offset_vector1[0], xyz2[0]+offset_vector1[0]], [
                midpoint[1]+offset_vector1[1], xyz2[1]+offset_vector1[1]], [midpoint[2]+offset_vector1[2], xyz2[2]+offset_vector1[2]]])
            bondLines.append([[xyz1[0], midpoint[0]], [
                xyz1[1], midpoint[1]], [xyz1[2], midpoint[2]]])
            bondLines.append([[midpoint[0], xyz2[0]], [
                midpoint[1], xyz2[1]], [midpoint[2], xyz2[2]]])
            bondLines.append([[xyz1[0]+offset_vector2[0], midpoint[0]+offset_vector2[0]], [
                xyz1[1]+offset_vector2[1], midpoint[1]+offset_vector2[1]], [xyz1[2]+offset_vector2[2], midpoint[2]+offset_vector2[2]]])
            bondLines.append([[midpoint[0]+offset_vector2[0], xyz2[0]+offset_vector2[0]], [
                midpoint[1]+offset_vector2[1], xyz2[1]+offset_vector2[1]], [midpoint[2]+offset_vector2[2], xyz2[2]+offset_vector2[2]]])

        return bondLines, bond_type

    def line_mid_points(self, xyz1, xyz2):
        '''
        Divide a line in two equal parts

        Parameters
        ----------
        xyz1: list
            (x,y,z) point 1
        xyz2: list
            (x,y,z) point 2

        Returns
        -------
        midPoints: list
            list of mid points
        '''
        # Calculate midpoints
        midpoint_x = (xyz1[0] + xyz2[0]) / 2
        midpoint_y = (xyz1[1] + xyz2[1]) / 2
        midpoint_z = (xyz1[2] + xyz2[2]) / 2

        # Divide center_line into two parts
        part1 = [[xyz1[0], midpoint_x], [
            xyz1[1], midpoint_y], [xyz1[2], midpoint_z]]
        part2 = [[midpoint_x, xyz2[0]], [
            midpoint_y, xyz2[1]], [midpoint_z, xyz2[2]]]

        # Get four points
        point1 = [xyz1[0], xyz1[1], xyz1[2]]
        point2 = [midpoint_x, midpoint_y, midpoint_z]
        point3 = [midpoint_x, midpoint_y, midpoint_z]
        point4 = [xyz2[0], xyz2[1], xyz2[2]]

        # res
        return point1, point2, point3, point4, part1, part2

    def line_property(self, xyz1, xyz2):
        '''
        Check a line property with which plane is parallel
        when res contains two True, it means the False coordination contains all elements.

        Parameters
        ----------
        xyz1: list
            (x,y,z) point 1
        xyz2: list
            (x,y,z) point 2

        Returns
        -------
        res: bool
            res
        '''
        try:
            # mean value
            Xm = np.mean([xyz1[0], xyz2[0]])
            Ym = np.mean([xyz1[1], xyz2[1]])
            Zm = np.mean([xyz1[2], xyz2[2]])
            xyzMean = [Xm, Ym, Zm]

            # check plane
            isSubtractZero = np.array([False, False, False])

            # points in one line
            X = xyz1[0] - xyz2[0]
            Y = xyz1[1] - xyz2[1]
            Z = xyz1[2] - xyz2[2]
            xyzPlane = [X, Y, Z]

            # axis vector
            xyzL = np.array([0, 0, 0])

            # set
            # axis selection
            perpendicularAxis = np.array([False, False, False])
            if X == 0:
                isSubtractZero[0] = True
                perpendicularAxis[0] = True
                # xyzL = xyzL + np.array(1, 1, 0)

            if Y == 0:
                isSubtractZero[1] = True
                perpendicularAxis[1] = True
                # xyzL = xyzL + np.array(1, 1, 0)

            if Z == 0:
                isSubtractZero[2] = True
                perpendicularAxis[2] = True
                # xyzL = xyzL + np.array(1, 0, 1)

            # set xyzL
            xyzL = isSubtractZero.astype(int)

            return xyzMean, xyzPlane, isSubtractZero, xyzL, perpendicularAxis

        except Exception as e:
            raise Exception(e)

    def set_plot_scale(self):
        '''
        Set plot scale
        '''
        # atom no
        atomNo = len(self.xyzList)
        # bond no
        bondNo = len(self.atomBonds)

        # bond length list
        bondLengthList = []

        # *** using bond block
        for i in range(bondNo):
            # atom id
            _atom1Id = int(self.atomBonds[i]['id']) - 1
            # atom symbol
            _atom1Symbol = self.atomBonds[i]['symbol']
            # atom color
            _atom1Color = self.set_color(_atom1Symbol)
            # atom bond list
            _atom1BondList = self.atomBonds[i]['bonds']
            atom1BondSize = len(_atom1BondList)

            _atom1X = self.xyzList[_atom1Id, 0]
            _atom1Y = self.xyzList[_atom1Id, 1]
            _atom1Z = self.xyzList[_atom1Id, 2]
            _atom1XYZ = [_atom1X, _atom1Y, _atom1Z]

            # draw bond
            if atom1BondSize > 0:
                for j in range(atom1BondSize):
                    # atom [2] id
                    _atom2Id = int(_atom1BondList[j][0]) - 1
                    # atom [2] symbol
                    _atom2Symbol = _atom1BondList[j][1]
                    # atom color
                    _atom2Color = self.set_color(_atom2Symbol)
                    # atom [1] - atom [2] bond type
                    _bondType = int(_atom1BondList[j][3])

                    # xyz
                    _atom2X = self.xyzList[_atom2Id, 0]
                    _atom2Y = self.xyzList[_atom2Id, 1]
                    _atom2Z = self.xyzList[_atom2Id, 2]
                    _atom2XYZ = [_atom2X, _atom2Y, _atom2Z]

                    # bond connection (points)
                    _bondConnection, _bondTypeLog, _bondLengths = self.create_bond_line(
                        _atom1XYZ, _atom2XYZ, _bondType)

                    # save bond length
                    bondLengthList.append(_bondLengths)

        # check bond length list
        if len(bondLengthList) > 0:
            # flatten list
            bondLengthList_flatten = sum(bondLengthList, [])

            # max bond length
            maxBondLength = max(bondLengthList_flatten)
            # min bond length
            minBondLength = min(bondLengthList_flatten)
            # mean bond length
            meanBondLength = np.mean(bondLengthList_flatten)
            # median bond length
            medianBondLength = np.median(bondLengthList_flatten)

            # set plot scale
            self.plotScale = [minBondLength, maxBondLength,
                              meanBondLength, medianBondLength]

    def view3d(self, subgraphs=None, **kwargs):
        '''
        Draw a compound in the cartesian coordinate
        atomElements atom symbol
        atomBonds atom bonds (bond blocks)
        xyzList atom position in the cartesian coordinate
        figSize=(10, 10) plt 3d setting
        obsOption=[False, 0] display center point [0,0,0]

        Parameters
        ----------
        figSize: tuple
            figure size
        bg_color: str
            background color
        display_legend: bool
            display legend

        Returns
        -------
        fig: figure
            figure
        '''
        figSize = kwargs.get('figSize', [])
        bg_color = kwargs.get('bg_color', '#ffffff')
        display_legend = kwargs.get('display_legend', True)
        theme = kwargs.get('theme', 'light')
        display_atom_id = kwargs.get('display_atom_id', True)
        display_bond_length = kwargs.get('display_bond_length', False)
        bond_type_color = kwargs.get(
            'bond_type_color', ['#1B263B', '#EF476F', '#4361EE'])

        # plot summary
        plot_summary = []

        # Create the figure
        fig = go.Figure()

        # legend
        legend_list = []

        # marker label
        marker_labels = []

        # atom no
        atomNo = len(self.xyzList)
        # bond no
        bondNo = len(self.atomBonds_1d)

        # create 3d frame
        xyzLenMax, xyzLenMin, xyzR, xLen, yLen, zLen = self.create_3dframe()

        # *** plot scale
        plot_scale_res = self.set_plot_scale()
        # min bond length
        min_bond_length = self.plotScale[0]
        # max bond length
        max_bond_length = self.plotScale[1]
        # set marker size
        marker_size_0 = self.set_marker_size(min_bond_length, max_bond_length)

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
            # size
            _atomSize = self.set_size(_atomSymbol, _sy=marker_size_0)
            # color
            _atomColor = self.set_color(_atomSymbol)

            # atom mark
            atomMark = str(_atomSymbol) + str(_atomId)

            # atom label
            marker_labels.append(atomMark)

            if display_atom_id:  # Replace with your condition
                text_to_display = [atomMark]
            else:
                text_to_display = ['']  # Empty string to hide text

            # scatter
            fig.add_trace(go.Scatter3d(x=[_atom1X],
                                       y=[_atom1Y],
                                       z=[_atom1Z],
                                       mode='markers+text',
                                       marker=dict(
                                           color=_atomColor, size=10, sizemode='area', sizeref=1, line=dict(width=2, color='black')),
                                       hoverinfo='all',  # Display all available information
                                       hoverlabel=dict(bgcolor=bg_color),
                                       # Set custom hover text
                                       hovertext=[f'Element: {atomMark}'],
                                       text=text_to_display))

            # textfont = dict(weight='normal')

        # *** bond visualization
        # *** using bond block
        for i in range(bondNo):
            # id 1
            _atom1Id = int(self.atomBonds_1d[i]['id1']) - 1
            # id 2
            _atom2Id = int(self.atomBonds_1d[i]['id2']) - 1
            # symbol 1
            _atom1Symbol = self.atomBonds_1d[i]['symbol1']
            # symbol 2
            _atom2Symbol = self.atomBonds_1d[i]['symbol2']
            # color 1
            _atom1Color = self.set_color(_atom1Symbol)
            # color 2
            _atom2Color = self.set_color(_atom2Symbol)

            # bond symbol
            _bondSymbol = self.atomBonds_1d[i]['bond_symbol']
            # bond type
            _bondType = int(self.atomBonds_1d[i]['bond_type']-1)
            # set bond type
            _bondTypeLabel = ''
            if _bondType == 1:
                _bondTypeLabel = 'single bond'
            elif _bondType == 2:
                _bondTypeLabel = 'double bond'
            else:
                _bondTypeLabel = 'triple bond'

            # xyz atom 1
            _atom1X = self.xyzList[_atom1Id, 0]
            _atom1Y = self.xyzList[_atom1Id, 1]
            _atom1Z = self.xyzList[_atom1Id, 2]
            _atom1XYZ = [_atom1X, _atom1Y, _atom1Z]

            # xyz atom 2
            _atom2X = self.xyzList[_atom2Id, 0]
            _atom2Y = self.xyzList[_atom2Id, 1]
            _atom2Z = self.xyzList[_atom2Id, 2]
            _atom2XYZ = [_atom2X, _atom2Y, _atom2Z]

            # distance
            _distance = self.calculate_distance(_atom1XYZ, _atom2XYZ)

            # plot summary
            plot_summary.append(
                {
                    'atom1Id': _atom1Id+1,
                    'atom2Id': _atom2Id+1,
                    'atom1Symbol': str(_atom1Symbol) + str(_atom1Id+1),
                    'atom2Symbol': str(_atom2Symbol) + str(_atom2Id+1),
                    'distance': _distance
                }
            )

            # Calculate midpoint coordinates
            midX = [(_atom1X + _atom2X) / 2]
            midY = [(_atom1Y + _atom2Y) / 2]
            midZ = [(_atom1Z + _atom2Z) / 2]

            # add line
            fig.add_trace(go.Scatter3d(x=[_atom1X, _atom2X],
                                       y=[_atom1Y, _atom2Y],
                                       z=[_atom1Z, _atom2Z],
                                       mode='lines',
                                       line=dict(color=bond_type_color[_bondType], width=3), hoverinfo='none', name=_bondSymbol, showlegend=True))

            if display_bond_length is True:  # Condition to show text
                text_to_display = [f'{_distance:.3f}']
            else:
                text_to_display = ['']  # Empty string to hide text

            # Add text at midpoint
            fig.add_trace(go.Scatter3d(x=midX,
                                       y=midY,
                                       z=midZ,
                                       mode='text',
                                       # Replace with your desired text
                                       text=text_to_display,
                                       hoverinfo='text',  # Display all available information
                                       hoverlabel=dict(
                                           bgcolor=bg_color, namelength=-1),
                                       # Set custom hover text
                                       hovertext=[f'A {_bondTypeLabel} with a length of {_distance:.3f}']))

        # *** visualize subgraph
        if subgraphs is not None:
            for subgraph in subgraphs:
                # add subgraph
                subgraph_pattern = subgraph['subgraph_pattern']
                _node_list = list(subgraph_pattern.nodes())
                for i in _node_list:
                    # xyz
                    _atom1X = self.xyzList[i-1, 0]
                    _atom1Y = self.xyzList[i-1, 1]
                    _atom1Z = self.xyzList[i-1, 2]
                    _atom1XYZ = [_atom1X, _atom1Y, _atom1Z]

                    # color
                    # atom id
                    _atomId = i
                    # symbol
                    _atomSymbol = str(self.atomElements[i-1]).strip()
                    # color
                    _atomColor = self.set_color(_atomSymbol)

                    # scatter
                    fig.add_trace(go.Scatter3d(x=[_atom1X],
                                               y=[_atom1Y],
                                               z=[_atom1Z],
                                               mode='markers',
                                               marker=dict(
                                                    color=_atomColor, size=2*10, opacity=0.5,
                                                    sizemode='area', sizeref=1,
                                                    line=dict(width=2, color='black')),
                                               ))

                # edge list
                _edge_list = list(subgraph_pattern.edges())
                for i in _edge_list:
                    # xyz
                    _atom1X = self.xyzList[i[0]-1, 0]
                    _atom1Y = self.xyzList[i[0]-1, 1]
                    _atom1Z = self.xyzList[i[0]-1, 2]
                    _atom1XYZ = [_atom1X, _atom1Y, _atom1Z]

                    # xyz
                    _atom2X = self.xyzList[i[1]-1, 0]
                    _atom2Y = self.xyzList[i[1]-1, 1]
                    _atom2Z = self.xyzList[i[1]-1, 2]
                    _atom2XYZ = [_atom2X, _atom2Y, _atom2Z]

                    # add line
                    fig.add_trace(go.Scatter3d(x=[_atom1X, _atom2X],
                                               y=[_atom1Y, _atom2Y],
                                               z=[_atom1Z, _atom2Z],
                                               mode='lines',
                                               line=dict(
                                                   color='blue', width=2*5, dash='dashdot'),
                                               hoverinfo='none'))

        # Set the limits of the axes
        # set max and offset
        xyzLenMax = xyzLenMax + 1.5
        fig.update_layout(scene=dict(
            xaxis=dict(nticks=4, range=[-xyzLenMax, xyzLenMax]),
            yaxis=dict(nticks=4, range=[-xyzLenMax, xyzLenMax]),
            zaxis=dict(nticks=4, range=[-xyzLenMax, xyzLenMax])
        ))

        # Set figure size to a square
        if len(figSize) != 0:
            fig.update_layout(width=figSize[0], height=figSize[1])
        else:
            fig.update_layout(
                autosize=True
            )

        # Show legend
        # Update layout to display legend
        fig.update_layout(legend=dict(
            orientation="h",  # Horizontal legend
            yanchor="bottom",  # Anchor at bottom
            y=1.02,  # Move legend up slightly
            xanchor="right",  # Anchor at right
            x=1  # Move legend to right
        ))

        if theme == 'black':
            font_color = 'lightgray'
        else:
            font_color = 'black'

        # Update text font color
        fig.update_traces(textfont=dict(color=font_color))

        # Set background color to dark
        fig.update_layout(
            paper_bgcolor=bg_color,
            plot_bgcolor=bg_color,
            scene=dict(
                xaxis=dict(showbackground=True,
                           backgroundcolor=bg_color),
                yaxis=dict(showbackground=True,
                           backgroundcolor=bg_color),
                zaxis=dict(showbackground=True,
                           backgroundcolor=bg_color)
            )
        )

        # Remove axes and other elements
        fig.update_layout(
            scene=dict(
                xaxis=dict(showticklabels=False, showgrid=False,
                           zeroline=False, showspikes=False, title=''),
                yaxis=dict(showticklabels=False, showgrid=False,
                           zeroline=False, showspikes=False, title=''),
                zaxis=dict(showticklabels=False, showgrid=False,
                           zeroline=False, showspikes=False, title=''),
                aspectmode='manual',
                aspectratio=dict(x=1, y=1, z=1),
            ),
            showlegend=False,
            margin=dict(l=0, r=0, b=0, t=0)
        )

        # Show the plot with zoom disabled
        fig.show(config={
            'scrollZoom': True,  # Disable zoom with scroll
            'displayModeBar': True,
            'displaylogo': True,
            'modeBarButtonsToRemove': ['zoom2d', 'zoomIn2d', 'zoomOut2d', 'pan2d']
        })

        # res
        return plot_summary

    def view3dobs(self, elev=None, azim=None, figSize=(10, 10), obsOption=[True, 0]):
        '''
        Draw a compound in the cartesian coordinate with observer

        Parameters
        ----------
        elev : float
            elevation angle
        azim : float
            azimuth angle
        figSize : tuple
            figure size
        obsOption : list
            [True, 0] --> show observer, 0 --> observer radius
        '''
        # 3d plot
        fig = plt.figure(figsize=figSize)
        ax = plt.axes(projection='3d')

        # atom no
        atomNo = len(self.xyzList)
        # bond no
        bondNo = len(self.atomBonds)

        # atom visualization
        for i in range(atomNo):
            # xyz
            _atomX = self.xyzList[i, 0]
            _atomY = self.xyzList[i, 1]
            _atomZ = self.xyzList[i, 2]
            # color
            # size

            # draw atom 1
            ax.scatter3D(_atomX, _atomY, _atomZ, s=40)

        # bond visualization
        for i in range(bondNo):
            # atom id
            _atom1Id = int(self.atomBonds[i]['id']) - 1
            # atom symbol
            _atom1Symbol = self.atomBonds[i]['symbol']
            # atom bond list
            _atom1BondList = self.atomBonds[i]['bonds']
            atom1BondSize = len(_atom1BondList)

            _atom1X = self.xyzList[:, 0]
            _atom1Y = self.xyzList[:, 1]
            _atom1Z = self.xyzList[:, 2]

            # draw bond
            if atom1BondSize > 0:
                for j in range(atom1BondSize):
                    # atom [2] id
                    _atom2Id = int(_atom1BondList[j][0]) - 1
                    # atom [2] symbol
                    _atom2Symbol = _atom1BondList[j][1]
                    # atom [1] - atom [2] bond type
                    _bondType = int(_atom1BondList[j][3])

                    # set color
                    lineColor = ['k', 'b', 'c']
                    lineWidth = [1, 2, 3]

                    # xyz
                    _atom2X = self.xyzList[_atom2Id, 0]
                    _atom2Y = self.xyzList[_atom2Id, 1]
                    _atom2Z = self.xyzList[_atom2Id, 2]

                    # bond connection
                    _bondConnection = self.xyzList[[_atom1Id, _atom2Id]]
                    # draw line
                    ax.plot3D(_bondConnection[:, 0], _bondConnection[:, 1], _bondConnection[:, 2],
                              color=lineColor[_bondType-1], linewidth=lineWidth[_bondType-1])

        # obs visualization
        xyzObsList = Observer.GeneratorCircleObserver(
            self._robs, self.tetaNo, self.phiNo, self.limits['teta'])[0]
        # obs size
        xyzObsSize = len(xyzObsList)
        for i in range(xyzObsSize):
            ax.scatter3D(xyzObsList[i, :, 0],
                         xyzObsList[i, :, 1], xyzObsList[i, :, 2])
            ax.plot3D(xyzObsList[i, :, 0],
                      xyzObsList[i, :, 1], xyzObsList[i, :, 2])

        # obs show
        if obsOption[0]:
            ax.scatter3D(obsOption[1], 0, 0, s=40)

        ax.view_init(elev=elev, azim=azim)
        plt.show()

    def create_line(self, xyzList1, xyzList2, t=1):
        '''
        Create a line equation and its parallel lines

        Parameters
        ----------
        xyzList1 : list
            [x1,y1,z1]
        xyzList2 : list
            [x2,y2,z2]
        t : float
            ratio between xyzList1 and xyzList2

        Returns
        -------
        r : list
            [[x1,y1,z1], [x2,y2,z2]]
        _l1 : list
            [[x1,y1,z1], [x2,y2,z2]]
        _l2 : list
            [[x1,y1,z1], [x2,y2,z2]]
        '''
        try:
            r = np.array(xyzList2) - np.array(xyzList1)

            # matrix
            rMat = np.array([xyzList1, xyzList1])

            # x
            x = r[0]*t + xyzList1[0]
            # y
            y = r[1]*t + xyzList1[1]
            # z
            z = r[2]*t + xyzList1[2]

            xL, yL, zL = 0.1*np.ones(3)
            _l1 = [[xyzList1[0]+xL, xyzList2[0]+xL], [xyzList1[1] +
                                                      yL, xyzList2[1]+yL], [xyzList1[2]+zL, xyzList2[2]+zL]]
            _l2 = [[xyzList1[0]-xL, xyzList2[0]-xL], [xyzList1[1] -
                                                      yL, xyzList2[1]-yL], [xyzList1[2]-zL, xyzList2[2]-zL]]

            # res
            return r, _l1, _l2
        except Exception as e:
            raise Exception(e)

    def calculate_distance(self, xyz1, xyz2):
        """
        Calculate the Euclidean distance between two points.

        Parameters:
        ----------
        xyz1 : list
            Coordinates of the first point [x, y, z]
        xyz2 : list
            Coordinates of the second point [x, y, z]

        Returns:
        -------
        float
            Distance between the two points
        """
        return math.sqrt((xyz2[0] - xyz1[0])**2 +
                         (xyz2[1] - xyz1[1])**2 +
                         (xyz2[2] - xyz1[2])**2)

    def set_marker_size(self, min_distance, ref_length=1, min_marker_size=200, max_marker_size=500):
        '''
        Set the plot view based on the minimum distance.

        Parameters:
        ----------
        min_distance : float
            Minimum distance between points
        ref_length : float
            Reference length for scaling
        min_marker_size : float
            Minimum marker size
        max_marker_size : float
            Maximum

        Returns:
        -------
        float
            Marker size
        '''
        scaling_factor = min_marker_size / ref_length
        marker_size = min_distance * scaling_factor

        # Ensure the marker size is within the desired range
        marker_size = min(max_marker_size, max(min_marker_size, marker_size))

        return marker_size