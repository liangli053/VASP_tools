from sys import argv, exit
import os
import shutil
import numpy as np
import random
import linecache
import subprocess
import math

""" Generate random atomic structures based on occupacies. Run Metropolis Monte Carlo
and return low-energy structures after certain iterations.

Electrostatic energyies calculated using GULP are used. May not apply to systems in which
other interactions are dominant. So far works well for layered structures. This examples
uses a four-element system. QEq parameters were fitted separately using genetic algorithm.
"""
no_Li_rmvd = int(argv[1])  # number of vacancies to generate
maxloop = 5000  # iterations of Metropolis Monte Carlo
tot_struct = 20  # number of strucutres to generate from random initializations
here = os.getcwd()
# read in pristine Li, Ir, O and Mn positions, assume Ir and O are immobile during charge
with open(here + "/Li16_pos") as textFile:
  Li16pos = np.array([[float(digit) for digit in line.split()] for line in textFile])
with open(here + "/Ir6_pos") as textFile:
  Ir6pos = np.array([[float(digit) for digit in line.split()] for line in textFile])
with open(here + "/O24_pos") as textFile:
  O24pos = np.array([[float(digit) for digit in line.split()] for line in textFile])
with open(here + "/Mn2_pos") as textFile:
  Mn2pos = np.array([[float(digit) for digit in line.split()] for line in textFile])
# read lattice constant
a_vec = linecache.getline("pristine.vasp", 3)
b_vec = linecache.getline("pristine.vasp", 4)
c_vec = linecache.getline("pristine.vasp", 5)
new_file = "gulp.in"
# Boltzmann constant (eV*K-1), room temperature
K_B = 8.6173303e-5
T = 300
run = 0
while (run < tot_struct):
  print "calculating", run, "structure"
  os.chdir(here)
  old_eng = -146 + no_Li_rmvd * 3.5
  i = 0
  print old_eng
  while (i < maxloop):
    # shuffle Li64pos, and take the first 8 as Fe positions
    # the following 48 as Li sites. Remaining are Li vacancies
    aa = Li16pos
    np.random.shuffle(aa)
    Li_pos = aa[no_Li_rmvd:16, :]
    Mn_possible_pos = np.concatenate((Mn2pos, aa[0:no_Li_rmvd, :]), axis=0)
    np.random.shuffle(Mn_possible_pos)
    Mn_pos = Mn_possible_pos[0:2, :]
##-------------------------------------------##
# write gulp input file and run QEq calculation
##-------------------------------------------##
    fout = open(new_file, 'w')
    fout.write("qeq\n"
               "title\n"
               "Mn_sub_Li2IrO3\n"
               "end\n"
               "vector\n")
    fout.write(a_vec)
    fout.write(b_vec)
    fout.write(c_vec)
    fout.write("fractional\n")
    for line in Li_pos:
      fout.write("Li   " + " ".join(str("%10.6f" % e) for e in line) + '\n')
    for line in Ir6pos:
      fout.write("Ir   " + " ".join(str("%10.6f" % e) for e in line) + '\n')
    for line in O24pos:
      fout.write("O    " + " ".join(str("%10.6f" % e) for e in line) + '\n')
    for line in Mn_pos:
      fout.write("Mn   " + " ".join(str("%10.6f" % e) for e in line) + '\n')

# QEq parameters are obtained using GA
    fout.write("space\n"
               "1\n"
               "qelectronegativity\n"
               "Li   -----   -----   -----   0.0\n"
               "Ir   -----   -----   -----   0.0\n"
               "O    -----   -----   -----   0.0\n"
               "Mn   -----   -----   -----   0.0\n"
               "qeqradius\n"
               "12.0")
    fout.close()
    os.system("/path/to/gulp <gulp.in> gulp.out")
    new_eng = (subprocess.check_output("grep 'Total lattice energy' gulp.out | head -1 \
              | awk '{print $5}'", shell=True))
    new_eng = float(new_eng)
##---------------------------------------------------------------##
# Metropolis Monte Carlo to calculate structural evolution based on
# electrostatic energies.
##---------------------------------------------------------------##
    new2old_rate = math.exp(-(new_eng - old_eng) / (K_B * T))
    # print new2old_rate
    p_rnd = random.uniform(0, 1)
    if p_rnd <= new2old_rate:
      old_eng = new_eng
      shutil.copy("gulp.in", "gulp.in_keep")
      shutil.copy("gulp.out", "gulp.out_keep")
      print old_eng
    i = i + 1

# write POSCAR based on last GULP run
  aaa = str(old_eng)
  if os.path.exists("Elec_" + aaa):
    file_name = "Elec_" + aaa + '_' + str(time.clock())
  else:
    file_name = "Elec_" + aaa
  os.mkdir(file_name)
  os.chdir(file_name)
  shutil.copy(here + "/gulp.out_keep", "gulp.out_keep")
  shutil.copy(here + "/gulp.in_keep", "gulp.in_keep")
  with open("gulp.in_keep") as textFile:
    gulp_in = textFile.read().splitlines()
  tot_atm_num = 16 + 6 + 24 + 2 - no_Li_rmvd + 9
  atm_coords = np.array([[float(digit) for digit in gulp_in[i].split()[1:4]] for i in range(9, tot_atm_num)])
  wri_POSCAR = open("POSCAR", 'w')
  wri_POSCAR.write("Mn_sub_Li2IrO3\n"
                   "1.000000\n")
  wri_POSCAR.write(a_vec)
  wri_POSCAR.write(b_vec)
  wri_POSCAR.write(c_vec)
  wri_POSCAR.write(" Li Ir O Mn\n"
                   " " + str(16 - no_Li_rmvd) + " 6 24 2\n"
                   "Direct\n")
  for line in atm_coords:
    wri_POSCAR.write("  " + " ".join(str("%10.6f" % e) for e in line) + '\n')
  wri_POSCAR.close()
  os.system("open -a VESTA.app POSCAR")
  run = run + 1
