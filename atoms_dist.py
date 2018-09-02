import numpy as np


def atoms_dist(a, b, latt_para):
    """ Return the distance between atoms a and b.
    a, b : array or list
        Two vectors of dimension (1 x 3)

    latt_para : array
        2D array with dimension (3 x 3)
    rtype : float
    """
    return np.linalg.norm(np.dot((a - b), latt_para), 2)
