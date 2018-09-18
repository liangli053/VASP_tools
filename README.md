# VASP_tools
This is a collection of python codes that process VASP input/output files. Will be updated regularly. <br/>

Parse_POSCAR.py(POSCAR="POSCAR") <br/>
    Parses files with VASP 5 POSCAR format, returns structure-related properties such as <br/>
    lattic constants, angles and atom coordinates. <br/>

Find_neighbors(ctrAtom, cutoff, POSCAR = "POSCAR") <br/>
    Finds the indices and coordinates of atoms surrounding ctrAtom. <br/>

Create_supercell.py(POSCAR, repetitions)<br/>
    Generates POSCAR file of supercell. Size determined by "repetitions". <br/>
