
# COMPUTATIONAL
# ---------------

# import libs
import numpy as np
from scipy.spatial import distance


class Compute():
    '''
    computational chemistry
    '''

    def __init__(self):
        pass

    @staticmethod
    def cal_atoms_distance(xyzAtom1, xyzAtom2):
        '''
        calculate distance between two atoms (center to center)
        '''
        return np.linalg.norm(xyzAtom1-xyzAtom2)

    @staticmethod
    def atoms_distance_matrix(xyzList, atomName):
        '''
        build a matrix containing a matrix of distance between two different atoms

        args:
            xyzList: xyz list of atoms
            atomName: atom name list such as ['C','H','H','H','H']
        '''
        # dict for atom index
        distance_res = []
        # atom no
        atomNo = len(atomName)
        atomLength = np.zeros((atomNo, atomNo))
        for i in range(atomNo):
            # atom1
            atom1XYZ = xyzList[i]
            for j in range(atomNo):
                # atom2
                atom2XYZ = xyzList[j]
                if j != i:
                    _length = Compute.cal_atoms_distance(atom1XYZ, atom2XYZ)
                else:
                    _length = 0
                # save
                atomLength[i, j] = _length

                # res
                distance_res.append(
                    {
                        'atom1': atomName[i]+str(i+1),
                        'atom2': atomName[j]+str(j+1),
                        'distance': _length
                    }
                )

        # res
        return atomLength, distance_res

    @staticmethod
    def atoms_distance(xyzList, atomName, atom_symbols, atom_index):
        '''
        calculate distance between two different atoms

        Parameters
        ----------
        xyzList : list
            xyz list of atoms
        atomName : list
            atom name list such as ['C','H','H','H','H']
        atom_symbols : list
            selected atom symbol list such as ['C','H']
        atom_index : list
            selected atom index list such as [0,1]

        Returns
        -------
        _length : float
            distance between two atoms
        '''
        try:
            # size
            atom_symbols_size = len(atom_symbols)
            atom_index_size = len(atom_index)

            # check
            if atom_symbols_size > 2 or atom_index_size > 2:
                raise Exception(
                    'atom symbol/index list must have two elements.')

            if atom_index_size == 0:
                # find the index of the first element
                atom1Index = np.where(atomName == atom_symbols[0])[0][0]
                atom2Index = np.where(atomName == atom_symbols[1])[0][0]
                # set
                atom1XYZ = xyzList[atom1Index]
                atom2XYZ = xyzList[atom2Index]
                _length = Compute.cal_atoms_distance(atom1XYZ, atom2XYZ)
            elif atom_index_size == 2:
                # set
                atom1XYZ = xyzList[int(atom_index[0])]
                atom2XYZ = xyzList[int(atom_index[1])]
                _length = Compute.cal_atoms_distance(atom1XYZ, atom2XYZ)
            elif atom_index_size > 2:
                raise Exception(
                    'atom symbol/index list must have two elements.')

            # res
            return _length
        except Exception as e:
            Exception(e)

    @staticmethod
    def calculate_angle(xyzList, atomName, atom_symbols, atom_index):
        '''
        Calculate angle/dihedral angle between points p1,p2, and p3 - p4

        Parameters
        ----------
        xyzList : list
            xyz list of atoms
        atomName : list
            atom name list such as ['C','H','H','H','H']
        atom_symbols : list
            selected atom symbol list such as ['C','H']
        atom_index : list
            selected atom index list such as [0,1]

        Returns
        -------
        angle_degrees : float
            angle in degrees
        '''
        # check
        atom_index_size = len(atom_index)

        if atom_index_size == 3:
            # angle between 3 points
            p1 = xyzList[int(atom_index[0])]
            p2 = xyzList[int(atom_index[1])]
            p3 = xyzList[int(atom_index[2])]

            # Calculate vectors
            v1 = np.array(p2) - np.array(p1)
            v2 = np.array(p2) - np.array(p3)

            # Calculate angle using cosine formula
            angle = np.arccos(
                np.dot(v1, v2) / (distance.euclidean(p1, p2) * distance.euclidean(p2, p3)))

            # Convert to degrees
            angle_degrees = np.degrees(angle)

        elif atom_index_size == 4:
            # dihedral angle 4 points
            p1 = xyzList[int(atom_index[0])]
            p2 = xyzList[int(atom_index[1])]
            p3 = xyzList[int(atom_index[2])]
            p4 = xyzList[int(atom_index[3])]

            # Calculate vectors
            v1 = p2 - p1
            v2 = p3 - p2
            v3 = p4 - p3

            # Calculate normal vectors
            n1 = np.cross(v1, v2)
            n2 = np.cross(v2, v3)

            # Normalize normal vectors
            n1 = n1 / np.linalg.norm(n1)
            n2 = n2 / np.linalg.norm(n2)

            # Calculate dihedral angle
            dihedral_angle = np.arccos(np.dot(n1, n2))

            # Calculate sign of dihedral angle
            sign = np.sign(np.dot(n1, v3))

            # Convert to degrees
            angle_degrees = np.degrees(dihedral_angle) * sign

        else:
            raise Exception(
                'atom symbol/index list must have three elements.'
            )

        # res
        return angle_degrees


# UTILITY FUNCTION
# ------------------


def CalculateMolecularMass(atom_elements, element_source):
    '''
    calculate molecular mass [g/mol]

    args:
        atom_elements: such as C, H
        element_source: periodic element table
    '''
    try:
        # # check
        # if element_source is None:
        #     # load periodic table of elements
        #     element_source = elements()
        # res
        res = 0
        # atom element size
        atomElementsSize = len(atom_elements)
        # atom property
        _atom_property = ['AtomicMass']
        # loop
        for i in range(atomElementsSize):
            _atom_symbol = atom_elements[i]
            _atom_prop_res = element_source.atom_properties(
                _atom_symbol, _atom_property)
            _atom_prop_val = _atom_prop_res.get(_atom_property[0])
            res += float(_atom_prop_val)

        return res
    except Exception as e:
        raise
