import numpy as np
import math
from atoms_dist import atoms_dist


def find_neighbors(ctrAtom, cutoff, POSCAR="POSCAR"):
    """ Return the coordinates of neighboring atoms.

    Use VASP 5 format POSCAR or CONTCAR files. Periodic boundary
    conditions are taken into accont of.

    ctrAtom : str
        Species and index of the central atom. e.g.: "O12", "Fe3".

    cutoff : float
        The cutoff radius, centered at ctrAtom (in angstrom).
        All atoms within cutoff radius of the center atom are counted.

    POSCAR : str
        Input file. Must be VASP 5 format.

    atomNames : list[str]
        Line 6 of POSCAR. Atomic species.

    atomNums : list[int]
        Line 7 of POSCAR. Atomic numbers

    atomNum_Dict : dict['str': int]
        zip atomNames and atomNums to form the dictionary.
        Each key represents one atomic species.

    atomCoor_Dict : dict['str': 2D array]
        Each key represents one atomic species. Values are 2D arrays
        of atomic coordinates. Dimension of 2D array is contingent to atomNums.

    res: dict[ 'str': List[ dict{str : list[float]} ] ]
        return value. Each key represents one atomic species. Each value is
        a list of dictionary {"atom_species + index" : coordinates}.
    """

    # center atom species and index
    #ctrAtom_name = ''.join(i for i in ctrAtom if i.isalpha())
    #ctrAtom_index = int(''.join(i for i in ctrAtom if i.isdigit()))
    if not ctrAtom[1].isalpha():
        ctrAtom_name = ctrAtom[0]
        ctrAtom_index = int(ctrAtom[1:])
    else:
        ctrAtom_name = ctrAtom[0:2]
        ctrAtom_index = int(ctrAtom[2:])

    fin = open(POSCAR, 'r')
    poscar = fin.read().splitlines()
    scaling_para = float(poscar[1])
    abc = np.array([[float(i) for i in line.split()] for line in poscar[2:5]])
    # lattice parameters in angstrom
    latt_para = abc * scaling_para
    # Lines 6 and 7 of POSCAR. atomic species and corresponding atoms numbers
    atomNames = poscar[5].split()
    atomNums = list(map(lambda x: int(x), poscar[6].split()))
    # combine atom names and numbers into a dict
    atomNum_Dict = dict(zip(atomNames, atomNums))
    # read in the coordinates of each species
    atomCoor_Dict = dict.fromkeys(atomNum_Dict, [])
    st_line = 8  # starting line number of atom coordinates
    for i in atomCoor_Dict.keys():
        end_line = st_line + atomNum_Dict[i]
        coor = np.array([[float(e) for e in line.split()[0:3]] for line in poscar[st_line: end_line]])
        st_line = end_line
        atomCoor_Dict[i] = coor
    fin.close()

    ctrAtomCoor = atomCoor_Dict[ctrAtom_name][ctrAtom_index - 1].reshape((1, 3))  # coordinates of the central atom
    length_a = np.linalg.norm(latt_para[0, :], 2)
    length_b = np.linalg.norm(latt_para[1, :], 2)
    length_c = np.linalg.norm(latt_para[2, :], 2)
    # calculate distances to the central atom
    res = dict.fromkeys(atomCoor_Dict, [])
    for i in atomCoor_Dict.keys():
        res[i] = []  # avoid name-binding problem!
        for coor in atomCoor_Dict[i]:
            currCoor = coor.reshape((1, 3))
            # dist_base = atoms_dist(ctrAtomCoor, currCoor, latt_para) # distance within the simulation cell
            index = np.where(np.all(atomCoor_Dict[i] == coor, axis=1))[0][0]
            # need to consider the periodic boundary condition
            repetition_a = math.ceil(cutoff / length_a)  # upper boundary of number of adjacent cells to search in
            repetition_b = math.ceil(cutoff / length_b)
            repetition_c = math.ceil(cutoff / length_c)
            # usually a, b and c are not big values, so this nested loop does not take much time
            for a in range(-repetition_a, repetition_a + 1):
                for b in range(-repetition_b, repetition_b + 1):
                    for c in range(-repetition_c, repetition_c + 1):
                        dist_to = currCoor + np.array([[a, 0, 0]]) + np.array([[0, b, 0]]) + np.array([[0, 0, c]])
                        if atoms_dist(ctrAtomCoor, dist_to, latt_para) <= cutoff:
                            res[i].append({i + str(int(index + 1)): list(np.squeeze(dist_to))})
    # remove the central atom from the dictionary
    res[ctrAtom_name].remove({ctrAtom: list(np.squeeze(ctrAtomCoor))})
    return(res)
