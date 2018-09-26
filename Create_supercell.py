import numpy as np
from Parse_POSCAR import parse_POSCAR


def create_supercell(POSCAR, repetitions):
    """ Create supercell from POSCAR. Supports VASP 5 only.
        Generate a outputfile named POSCAR_NxNyNz.vasp.
        Nx, Ny and Nz are the repetitions along x, y and z.

        Arguments:
        -----------------
        POSCAR : str
            Name of file with VASP POSCAR format.

        repetitions : list or tuple of int, len = 3
            Supercell size -- repetitions along x, y and z directions,
            relative to the original POSCAR or CONTCAR. Elements must
            be greater than 0.
    """

    latt_mat, _, _, atomNum_Dict, atomCoor_Dict = parse_POSCAR(POSCAR)
    # lattice matrix of the supercell, broadcasting
    output_latt_mat = latt_mat * np.array(repetitions)[:, np.newaxis]
    # compute the atom coordinates in the supercell
    output_coords = dict.fromkeys(atomCoor_Dict, [])
    for atom in atomCoor_Dict:
        coords_after_reps = repeat_N(atomCoor_Dict[atom], repetitions)
        output_coords[atom] = coords_after_reps

    # write to file
    foo = '_' + str(repetitions[0]) + str(repetitions[1]) + str(repetitions[2])
    with open("POSCAR" + foo + ".vasp", 'w') as fout:
        fout.write("supercell" + '\n'
                   "  1.000" + '\n')
        for line in output_latt_mat:
            fout.write("  " + "  ".join(("%.16f" % e) for e in line) + '\n')
        fout.write("  ".join(i for i in atomNum_Dict.keys()) + '\n' +
                   "  ".join(str(i * np.prod(repetitions)) for i in atomNum_Dict.values()) + '\n'
                   "Direct" + '\n')
        for atom in output_coords:
            for idx, line in enumerate(output_coords[atom], 1):
                fout.write("  " + "  ".join(str("%18.16f" % e) for e in line) +
                           "   " + atom + str(idx) + '\n')


def repeat_N(coords, repetitions):
    """ Compute atom coordinates for the supercell.
        Size dictated by repetitions [Nx, Ny, Nz].

        Arguments:
        -----------------
        coords : array[float], dim = (natom, 3)
            Original coordinates of atoms in the unit cell.

        repetitions : [Nx, Ny, Nz], list or tuple of int, len = 3
            Supercell size -- repetitions along x, y and z directions,
            relative to the original POSCAR or CONTCAR. Elements must
            be greater than 0.

        Returns:
        -----------------
        coords : array[float], dim = (Nx x Ny x Nz x natom, 3)
            Atom coordinates in the supercell.
    """
    # need to consider three dimensions, x, y and z axis
    for (dim, rep) in enumerate(repetitions):
        row, column = coords.shape
        # This "coords" name is local to the function.
        # Needs to be returned.
        coords = np.repeat(coords, rep, axis=0)
        added_by = np.zeros_like(coords)
        divided_by = np.ones_like(coords)
        # Add unit vector to coorespoinding coordinates
        patch = np.arange(rep)[:, np.newaxis]
        patch = np.tile(patch, (row, 1))
        added_by[:, [dim]] = patch
        divided_by[:, [dim]] = rep
        coords += added_by
        # Divide by the repetion along axis to normailize
        coords /= divided_by

    # sort the coordinate according to Z value
    coords = coords[coords[:, 2].argsort()]
    return coords
