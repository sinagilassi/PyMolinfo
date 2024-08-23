# MATERIAL
# ----------

# import libs
import numpy as np
import pandas as pd
from matplotlib.pyplot import xlabel
from ..config import OBSERVER_PROPERTY
from .utility import Utility
from .graph3d import graph3d
from .netwrok import Network
from .compute import Compute, CalculateMolecularMass


class Compound(graph3d, Network):
    '''
    material along with properties
    '''
    _mat_cid = ''
    _mat_name = ''
    _mat_formula = ''
    _atom_numbers = 0
    _atom_elements = []
    _atom_xyz = []
    _atom_xyz_center = []
    _atom_bond_block = []
    _atom_bond_block_1d = []
    _atom_bond_numbers = 0
    _atom_block = []
    # distance
    _distance = []

    # *** obs fixed ***
    _limits = {
        'phi': [],
        'teta': [np.deg2rad(0), np.deg2rad(180)]
    }
    _robs = 10
    _tetaNo = OBSERVER_PROPERTY.get('TETA')
    _phiNo = OBSERVER_PROPERTY.get('PHI')
    _dataNo = []
    _obsCoordinate = []

    def __init__(self, parse_prop):
        # all properties
        self.parse_prop = parse_prop
        # super class
        __atom_elements = parse_prop['atom_elements']
        __atom_block = parse_prop['atom_block']
        __atom_bonds = parse_prop['bond_block']
        __atom_xyz = parse_prop['xyz_list']
        __atom_xyz_center = parse_prop['xyz_center_list']
        # limit
        __limits = self._limits

        # TODO: convert atom_bonds to 1d vector
        __atom_bonds_1d = self.convert_atom_bonds(__atom_bonds)
        # update parse_prop
        self.parse_prop['bond_block_1d'] = __atom_bonds_1d

        # ! init parent classes
        # *** raw info (just for visualizing a structure)
        graph3d.__init__(self, __atom_elements, __atom_bonds,
                         __atom_xyz, __atom_xyz_center, self.robs, self.tetaNo, self.phiNo, __limits, __atom_bonds_1d)

        # *** network
        Network.__init__(self, __atom_elements, __atom_bonds,
                         __atom_xyz, __atom_xyz_center, __atom_bonds_1d)

        # update mat prop
        self.__update_atom_prop('mat_cid')
        self.__update_atom_prop('mat_name')
        self.__update_atom_prop('mat_formula')
        self.__update_atom_prop('atom_numbers')
        self.__update_atom_prop('atom_elements')
        self.__update_atom_prop('bond_numbers')
        self.__update_atom_prop('bond_block')
        self.__update_atom_prop('xyz_list')
        self.__update_atom_prop('xyz_center_list')
        self.__update_atom_prop('atom_block')

    def __str__(self):
        '''
        Return info about the mat
        '''
        # check
        if self.parse_prop['mat_name'] is not None:
            return str(self.parse_prop['mat_name'])
        else:
            return str(self.parse_prop['atom_elements'])

    @property
    def mat_cid(self):
        return self.parse_prop['mat_cid']

    @mat_cid.setter
    def mat_cid(self, value):
        self._mat_cid = value

    @property
    def mat_name(self):
        return self.parse_prop['mat_name']

    @mat_name.setter
    def mat_name(self, value):
        self._mat_name = value

    @property
    def mat_formula(self):
        return self.parse_prop['mat_formula']

    @mat_formula.setter
    def mat_formula(self, value):
        self._mat_formula = value

    @property
    def atom_elements(self):
        return np.array(self.parse_prop['atom_elements'])

    @atom_elements.setter
    def atom_elements(self, value):
        self._atom_elements = value

    @property
    def atom_numbers(self):
        return self.parse_prop['atom_numbers']

    @atom_numbers.setter
    def atom_numbers(self, value):
        self._atom_numbers = value

    @property
    def atom_xyz(self):
        return np.array(self.parse_prop['xyz_list'])

    @atom_xyz.setter
    def atom_xyz(self, value):
        self._atom_xyz = value

    @property
    def atom_xyz_center(self):
        return np.array(self.parse_prop['xyz_center_list'])

    @atom_xyz_center.setter
    def atom_xyz_center(self, value):
        self._atom_xyz_center = value

    @property
    def atom_bond_block(self):
        return self.parse_prop['bond_block']

    @atom_bond_block.setter
    def atom_bond_block(self, value):
        self._atom_bond_block = value

    @property
    def atom_bond_block_1d(self):
        return self.parse_prop['bond_block_1d']

    @atom_bond_block_1d.setter
    def atom_bond_block_1d(self, value):
        self._atom_bond_block_1d = value

    @property
    def atom_bond_numbers(self):
        return self.parse_prop['bond_numbers']

    @atom_bond_numbers.setter
    def atom_bond_numbers(self, value):
        self._atom_bond_numbers = value

    @property
    def atom_block(self):
        return self.parse_prop['atom_block']

    @atom_block.setter
    def atom_block(self, value):
        self._atom_block = value

    # *** observer prop
    @property
    def limits(self):
        return self._limits

    @limits.setter
    def limits(self, value):
        self._limits = []
        self._limits = value

    @property
    def robs(self):
        return self._robs

    @robs.setter
    def robs(self, value):
        self._robs = value

    @property
    def tetaNo(self):
        return self._tetaNo

    @tetaNo.setter
    def tetaNo(self, value):
        self._tetaNo = value

    @property
    def phiNo(self):
        return self._phiNo

    @phiNo.setter
    def phiNo(self, value):
        self._phiNo = value

    @property
    def dataNo(self):
        return self._dataNo

    @dataNo.setter
    def dataNo(self, value):
        self._dataNo = value

    @property
    def obsCoordinate(self):
        return self._obsCoordinate

    @obsCoordinate.setter
    def obsCoordinate(self, value):
        self._obsCoordinate = []
        self._obsCoordinate = value

    @property
    def distance(self):
        return self._distance

    @distance.setter
    def distance(self, value):
        self._distance = []
        self._distance = value

    def __update_atom_prop(self, prop_name):
        '''
        Update atom prop

        Parameters
        ----------
        prop_name: str
            atom prop name
        '''
        switchProp = {
            'header_block': 1,
            'counts_line': 1,
            'mat_cid': self.__update_mat_cid,
            'mat_name': self.__update_mat_name,
            'mat_formula': self.__update_mat_formula,
            'atom_numbers': self.__update_atom_numbers,
            'atom_elements': self.__update_atom_elements,
            'bond_numbers': self.__update_atom_bond_numbers,
            'atom_block': self.__update_atom_block,
            'bond_block': self.__update_atom_bond_block,
            'xyz_list': self.__update_atom_xyz,
            'xyz_center_list': self.__update_atom_xyz_center,
            'bond_block_1d': self.__update_atom_bond_block_1d
        }
        # select prop
        propSelection = switchProp.get(prop_name)
        # check
        if propSelection is not None:
            propSelection(self.parse_prop[prop_name])

    def __update_mat_cid(self, prop_val):
        self.mat_cid = prop_val

    def __update_mat_name(self, prop_val):
        self.mat_name = prop_val

    def __update_mat_formula(self, prop_val):
        self.mat_formula = prop_val

    def __update_atom_numbers(self, prop_val):
        self.atom_numbers = prop_val

    def __update_atom_elements(self, prop_val):
        self.atom_elements = prop_val

    def __update_atom_xyz(self, prop_val):
        self.atom_xyz = prop_val

    def __update_atom_xyz_center(self, prop_val):
        self.atom_xyz_center = prop_val

    def __update_atom_bond_block(self, prop_val):
        self.atom_bond_block = prop_val

    def __update_atom_bond_numbers(self, prop_val):
        self.atom_bond_numbers = prop_val

    def __update_atom_block(self, prop_val):
        self.atom_block = prop_val

    def __update_atom_bond_block_1d(self, prop_val):
        self.atom_bond_block_1d = prop_val

    def convert_atom_bonds(self, atom_bonds):
        '''
        Convert atom bonds to 1d list
        '''
        atom_bonds_1d = []
        for atom in atom_bonds:
            for bond in atom['bonds']:
                atom_bonds_1d.append({
                    'id1': atom['id'],
                    'symbol1': atom['symbol'],
                    'id2': bond[0],
                    'symbol2': bond[1],
                    'bond_type': bond[3],
                    'bond_symbol': bond[2]
                })
        return atom_bonds_1d

    def distance_matrix(self, dataframe=False):
        '''
        Build a matrix of atom-atom distance

        Parameters
        ----------
        dataframe: bool
            return a dataframe


        '''
        # res
        matrix, _dict = Compute.atoms_distance_matrix(
            self.xyzList, self.atom_elements)

        if dataframe:
            return pd.DataFrame(matrix, columns=self.atom_elements,
                                index=self.atom_elements)
        else:
            return matrix

    def distance_atoms(self, atoms):
        '''
        Calculate distance between two different atoms

        Parameters
        ----------
        atoms : list
            atom ids such as ['C1','H2']

        Returns
        -------
        distance: float
            distance between two atoms
        '''
        # get letters
        atoms_result = [{"letters": ''.join([c for c in s if c.isalpha(
        )]), "numbers": ''.join([c for c in s if c.isdigit()])} for s in atoms]
        # update atoms
        atom_symbols = [str(i['letters']) for i in atoms_result]
        # atom id
        atom_index = [str(int(i['numbers'])-1) for i in atoms_result]
        # res
        return Compute.atoms_distance(self.xyzList, self.atom_elements,
                                      atom_symbols, atom_index)

    def angle_atoms(self, atoms):
        '''
        Calculate angle between points p1,p2, and p3

        Parameters
        ----------
        atoms : list
            atom ids such as ['C1','H2','O3']

        Returns
        -------
        angle_degrees : float
            angle in degrees
        '''
        # check atoms
        if len(atoms) != 3:
            raise Exception('three atoms not provided!')

        # get letters
        atoms_result = [{"letters": ''.join([c for c in s if c.isalpha(
        )]), "numbers": ''.join([c for c in s if c.isdigit()])} for s in atoms]
        # update atoms
        atom_symbols = [str(i['letters']) for i in atoms_result]
        # atom id
        atom_index = [str(int(i['numbers'])-1) for i in atoms_result]

        return Compute.calculate_angle(self.xyzList, self.atom_elements,
                                       atom_symbols, atom_index)

    def d_angle_atoms(self, atoms):
        '''
        Calculate dihedral angle between points p1,p2,p3, and p4

        Parameters
        ----------
        atoms : list
            atom ids such as ['C1','H2','O3','H4']

        Returns
        -------
        angle_degrees : float
            angle in degrees
        '''
        # check atoms
        if len(atoms) != 4:
            raise Exception('4 atoms not provided!')

        # get letters
        atoms_result = [{"letters": ''.join([c for c in s if c.isalpha(
        )]), "numbers": ''.join([c for c in s if c.isdigit()])} for s in atoms]
        # update atoms
        atom_symbols = [str(i['letters']) for i in atoms_result]
        # atom id
        atom_index = [str(int(i['numbers'])-1) for i in atoms_result]

        return Compute.calculate_angle(self.xyzList, self.atom_elements,
                                       atom_symbols, atom_index)
