def fetch_snapshots(simul_times2fetch, POTIM, input_file="XDATCAR"):
    """ Given a VASP XDATCAR file from a AIMD calculation, generate
        POSCARs at specific simulation times.

        Arguments:
        -----------------
        simul_times2fetch : list [int or float]
            List of simulation times (pico second, ps) to generate POSCARs.

        POTIM : float
            Time step (femto second, fs) used in AIMD simulation. Can be inferred from INCAR.

        input_file : str
            VASP XDATCAR-like file. By default is "XDATCAR"
    """

    with open(input_file, 'r') as fin:
        # read each line from input file
        lines = fin.read().splitlines()

    # read the ionic steps, lines starting with "Direc configurations"
    ionic_steps = [float(line.split('=')[1].strip()) for line in lines if line.startswith("Direct")]
    # convert ionic steps to simulation times (ps) by multiplying POTIM
    simul_times = [i * POTIM / 1000 for i in ionic_steps]

    for time2get in simul_times2fetch:
        # Using binary search, find the entry in simul_times that is equal or closest to simul_time
        idx = binary_search(simul_times, 0, len(simul_times) - 1, time2get)
        # Find the corresponding line number in the input file
        line2find = [line for line in lines if line.startswith("Direct") and line.endswith(' ' + str(int(ionic_steps[idx])))]
        idx_in_file = lines.index(line2find[0])
        # total number of atoms in the simulation cell
        tot_num_atoms = sum([int(i) for i in lines[idx_in_file - 1].split()])
        # starting and ending line numbers in the original XDATCAR file
        write_starts = idx_in_file - 6
        write_ends = idx_in_file + tot_num_atoms + 1
        with open("POSCAR_" + str(time2get) + ".vasp", 'w') as fout:
            fout.write("snapshot at " + str(time2get) + " ps" + '\n')
            for i in range(write_starts, write_ends):
                fout.write(lines[i] + '\n')


def binary_search(arr, left, right, x):
    """ Binary search to find index of x in arr. If x is not found,
        then return the index of value that is closest to x.

        Arguments:
        -----------------
        arr : list [float]
            List of values to search in.

        left : int
            Left limit to perform binary search.

        right : int
            Right limit to perform binary search.

        x : float
            Target value to find in arr.

        Returns: int
        -----------------
            Element index in arr.

    """

    if right < left:
        # if x is not an element in arr, then return the closest value in arr
        return left if abs(arr[left] - x) < abs(arr[right] - x) else right
    else:
        mid = left + (right - left) // 2
        if arr[mid] == x:
            return mid
        elif arr[mid] > x:
            return binary_search(arr, left, mid - 1, x)
        else:
            return binary_search(arr, mid + 1, right, x)
