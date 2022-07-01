#!/usr/bin/env python3

#=========================================================================================
# Peacemaker -- A Quantum Cluster Equilibrium Code.
#
# Copyright 2004-2006 Barbara Kirchner, University of Bonn
# Copyright 2007-2012 Barbara Kirchner, University of Leipzig
# Copyright 2013-2022 Barbara Kirchner, University of Bonn
#
# This file is part of Peacemaker.
#
# Peacemaker is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Peacemaker is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Peacemaker.  If not, see <http://www.gnu.org/licenses/>
#=========================================================================================

import sys
import os.path
import getopt
import numpy as np
from numpy import linalg as LA
import math


# VdW volumes taken from Bondi's compilation.
radii = {"Ac" : 2.47,
         "Al" : 1.84,
         "Am" : 2.44,
         "Sb" : 2.06,
         "Ar" : 1.88,
         "As" : 1.85,
         "At" : 2.02,
         "Ba" : 2.68,
         "Bk" : 2.44,
         "Be" : 1.53,
         "Bi" : 2.07,
         "B"  : 1.92,
         "Br" : 1.85,
         "Cd" : 2.18,
         "Ca" : 2.31,
         "Cf" : 2.45,
         "C"  : 1.7,
         "Ce" : 2.42,
         "Cs" : 3.43,
         "Cl" : 1.75,
         "Cr" : 2.06,
         "Co" : 2,
         "Cu" : 1.96,
         "Cm" : 2.45,
         "Dy" : 2.31,
         "Es" : 2.45,
         "Er" : 2.29,
         "Eu" : 2.35,
         "Fm" : 2.45,
         "F"  : 1.47,
         "Fr" : 3.48,
         "Gd" : 2.34,
         "Ga" : 1.87,
         "Ge" : 2.11,
         "Au" : 2.14,
         "Hf" : 2.23,
         "He" : 1.4,
         "Ho" : 2.3,
         "H"  : 1.1,
         "In" : 1.93,
         "I"  : 1.98,
         "Ir" : 2.13,
         "Fe" : 2.04,
         "Kr" : 2.02,
         "La" : 2.43,
         "Lr" : 2.46,
         "Pb" : 2.02,
         "Li" : 1.82,
         "Lu" : 2.24,
         "Mg" : 1.73,
         "Mn" : 2.05,
         "Md" : 2.46,
         "Hg" : 2.23,
         "Mo" : 2.17,
         "Nd" : 2.39,
         "Ne" : 1.54,
         "Np" : 2.39,
         "Ni" : 1.97,
         "Nb" : 2.18,
         "N"  : 1.55,
         "No" : 2.46,
         "Os" : 2.16,
         "O"  : 1.52,
         "Pd" : 2.1,
         "P"  : 1.8,
         "Pt" : 2.13,
         "Pu" : 2.43,
         "Po" : 1.97,
         "K"  : 2.75,
         "Pr" : 2.4,
         "Pm" : 2.38,
         "Pa" : 2.43,
         "Ra" : 2.83,
         "Rn" : 2.2,
         "Re" : 2.16,
         "Rh" : 2.1,
         "Rb" : 3.03,
         "Ru" : 2.13,
         "Sm" : 2.36,
         "Sc" : 2.15,
         "Se" : 1.9,
         "Si" : 2.1,
         "Ag" : 2.11,
         "Na" : 2.27,
         "Sr" : 2.49,
         "S"  : 1.8,
         "Ta" : 2.22,
         "Tc" : 2.16,
         "Te" : 2.06,
         "Tb" : 2.33,
         "Tl" : 1.96,
         "Th" : 2.45,
         "Tm" : 2.27,
         "Sn" : 2.17,
         "Ti" : 2.11,
         "W"  : 2.18,
         "U"  : 2.41,
         "V"  : 2.07,
         "Xe" : 2.16,
         "Yb" : 2.26,
         "Y"  : 2.32,
         "Zn" : 2.01,
         "Zr" : 2.23,
         "Xx" : 0.0
         }


# Masses taken from http://www.csudh.edu/oliver/chemdata/atmass.htm.
masses = {
        "Ac" : 227.028,
        "Al" : 26.981539,
        "Am" : 243,
        "Sb" : 121.757,
        "Ar" : 39.948,
        "As" : 74.92159,
        "At" : 210,
        "Ba" : 137.327,
        "Bk" : 247,
        "Be" : 9.012182,
        "Bi" : 208.98037,
        "Bh" : 262,
        "B"  : 10.811,
        "Br" : 79.904,
        "Cd" : 112.411,
        "Ca" : 40.078,
        "Cf" : 251,
        "C"  : 12.011,
        "Ce" : 140.115,
        "Cs" : 132.90543,
        "Cl" : 35.4527,
        "Cr" : 51.9961,
        "Co" : 58.93320,
        "Cu" : 63.546,
        "Cm" : 247,
        "Db" : 262,
        "Dy" : 162.50,
        "Es" : 252,
        "Er" : 167.26,
        "Eu" : 151.965,
        "Fm" : 257,
        "F"  : 18.9984032,
        "Fr" : 223,
        "Gd" : 157.25,
        "Ga" : 69.723,
        "Ge" : 72.61,
        "Au" : 196.96654,
        "Hf" : 178.49,
        "Hs" : 265,
        "He" : 4.002602,
        "Ho" : 164.93032,
        "H"  : 1.00794,
        "In" : 114.82,
        "I"  : 126.90447,
        "Ir" : 192.22,
        "Fe" : 55.847,
        "Kr" : 83.80,
        "La" : 138.9055,
        "Lr" : 262,
        "Pb" : 207.2,
        "Li" : 6.941,
        "Lu" : 174.967,
        "Mg" : 24.3050,
        "Mn" : 54.93805,
        "Mt" : 266,
        "Md" : 258,
        "Hg" : 200.59,
        "Mo" : 95.94,
        "Nd" : 144.24,
        "Ne" : 20.1797,
        "Np" : 237.048,
        "Ni" : 58.6934,
        "Nb" : 92.90638,
        "N"  : 14.00674,
        "No" : 259,
        "Os" : 190.2,
        "O"  : 15.9994,
        "Pd" : 106.42,
        "P"  : 30.973762,
        "Pt" : 195.08,
        "Pu" : 244,
        "Po" : 209,
        "K"  : 39.0983,
        "Pr" : 140.90765,
        "Pm" : 145,
        "Pa" : 231.0359,
        "Ra" : 226.025,
        "Rn" : 222,
        "Re" : 186.207,
        "Rh" : 102.90550,
        "Rb" : 85.4678,
        "Ru" : 101.07,
        "Rf" : 261,
        "Sm" : 150.36,
        "Sc" : 44.955910,
        "Sg" : 263,
        "Se" : 78.96,
        "Si" : 28.0855,
        "Ag" : 107.8682,
        "Na" : 22.989768,
        "Sr" : 87.62,
        "S"  : 32.066,
        "Ta" : 180.9479,
        "Tc" : 98,
        "Te" : 127.60,
        "Tb" : 158.92534,
        "Tl" : 204.3833,
        "Th" : 232.0381,
        "Tm" : 168.93421,
        "Sn" : 118.710,
        "Ti" : 47.88,
        "W"  : 183.85,
        "U"  : 238.0289,
        "V"  : 50.9415,
        "Xe" : 131.29,
        "Yb" : 173.04,
        "Y"  : 88.90585,
        "Zn" : 65.39,
        "Zr" : 91.224,
        "Xx" : 0.0
}

elements = ["C", "H", "O", "Ac", "Al", "Am", "Sb", "Ar", "As", "At", "Ba", "Bk", "Be", "Bi",
            "Bh", "B", "Br", "Cd", "Ca", "Cf", "Ce", "Cs", "Cl", "Cr", "Co", "Cu", "Cm", "Db",
            "Dy", "Es", "Er", "Eu", "Fm", "F", "Fr", "Gd", "Ga", "Ge", "Au", "Hf", "Hs", "He",
            "Ho", "In", "I", "Ir", "Fe", "Kr", "La", "Lr", "Pb", "Li", "Lu", "Mg", "Mn", "Mt",
            "Md", "Hg", "Mo", "Nd", "Ne", "Np", "Ni", "Nb", "N", "No", "Os", "Pd", "P", "Pt",
            "Pu", "Po", "K", "Pr", "Pm", "Pa", "Ra", "Rn", "Re", "Rh", "Rb", "Ru", "Rf", "Sm",
            "Sc", "Sg", "Se", "Si", "Ag", "Na", "Sr", "S", "Ta", "Tc", "Te", "Tb", "Tl", "Th",
            "Tm", "Sn", "Ti", "W", "U", "V", "Xe", "Yb", "Y", "Zn", "Zr", "Xx"]

rflag = 0
fflag = 0
xflag = 0
tflag = 0
vflag = 0
mflag = 0

try:
    opts, args = getopt.getopt(sys.argv[1:], "frxtvm")
except getopt.GetoptError:
    print("Usage:\t{} [option] <files>\n\t[-r] \t rotates molecule with inertia tensor\n\t[-x] \t creates a xyz format output\n\t[-t] \t creates a tex format output".format(sys.argv[0]))
    sys.exit(2)
for opt, arg in opts:
    if opt in ("-r"):
      rflag = 1
    if opt in ("-f"):
      fflag = 1
    if opt in ("-x"):
      xflag = 1
    if opt in ("-t"):
      tflag = 1
    if opt in ("-v"):
      vflag = 1
    if opt in ("-m"):
      mflag = 1

if len(args) == 0:
    print("Usage:\t{} [option] <files>\n\t[-r] \t rotates molecule with inertia tensor\n\t[-x] \t creates a xyz format output\n\t[-t] \t creates a tex format output".format(sys.argv[0]))
    sys.exit(1)

j = 0
nX = 0
W = np.zeros((len(args), 3))

for j in range(len(args)):
  filename = args[j]

  with open(filename, "r") as fin:
    # read xyz file
    line = fin.readline()
    natoms = int(line)
    pos = np.zeros((natoms, 3))
    symbol = []
    m = np.zeros(natoms)
    vol = np.zeros(natoms)
    line = fin.readline()
    comment = line[:-1]
    for i in range(natoms):
      line = fin.readline()
      label, x, y, z = line.split()
      label = label.lower()
      label = label.capitalize()
      if label == "Xx": nX = nX + 1
      m[i] = masses[label]
      vol[i] = 4.0/3.0*math.pi*radii[label]**3
      symbol.append(label)
      pos[i,:] = float(x), float(y), float(z)

    formula = ""
    for e in elements:
      count = symbol.count(e)
      if (count > 0):
        formula = formula + e
        if (count > 1):
          formula = formula + str(count)

    # shift COM to origin
    com = np.zeros(3)
    mtot = 0.0
    vtot = 0.0
    for i in range(natoms):
      com += m[i]*pos[i,:]
      mtot += m[i]
      vtot += vol[i]

    if (vflag):
      print(vtot)
      continue
    elif (mflag):
      print(mtot)
      continue
    elif (fflag):
      print(formula)
      continue

    com /= mtot
    for i in range(natoms):
      pos[i,:] -= com

    # calculate inertia tensor
    inertia = np.zeros((3, 3))
    for i in range(natoms):
      inertia[0, 0] += m[i]*(pos[i,1]**2 + pos[i,2]**2)
      inertia[1, 1] += m[i]*(pos[i,0]**2 + pos[i,2]**2)
      inertia[2, 2] += m[i]*(pos[i,0]**2 + pos[i,1]**2)
      inertia[0, 1] -= m[i]*pos[i,0]*pos[i,1]
      inertia[0, 2] -= m[i]*pos[i,0]*pos[i,2]
      inertia[1, 2] -= m[i]*pos[i,1]*pos[i,2]
    inertia[1,0] = inertia[0,1]
    inertia[2,0] = inertia[0,2]
    inertia[2,1] = inertia[1,2]

    # diagonalize inertia tensor and print
    w, v = LA.eig(inertia)

    p = w.argsort()
    w = w[p]
    v = v[:,p]

#    print(v)

    # if -r option is active: rotate coordinates with inertia tensor
    if rflag == 1:
      for i in range(natoms):
        pos[i,:] = v.T.dot(pos[i,:])

    # if -x option is active: print out in .xyz format
    if xflag == 1:
      print("%u" % (natoms-nX))
      print("{:s}".format(comment))
      for i in range(natoms):
        if symbol[i]=="CL": symbol[i]="Cl"
        if symbol[i]=="Xx": continue
        print("{:3s} {:11.8f} {:13.8f} {:13.8f}".format(symbol[i], pos[i,0], pos[i,1], pos[i,2]))
      print("")

    # if -t option is active: print out in tex format
    elif tflag == 1:
#     print("\documentclass[11pt]{article}")
#     print("\usepackage[a4paper,lmargin={2cm},rmargin={2cm},\ntmargin={2cm},bmargin = {2cm}]{geometry}")
#     print("\usepackage{setspace}")
#     print("\usepackage{longtable}")
#     print("\onehalfspacing")
#     print("\\renewcommand*{\\familydefault}{\sfdefault}")
#     print("\n\\begin{document}")
      print("\\begin{longtable}[l]{lrrr}")
      print("\multicolumn{{4}}{{l}}{{\\textbf{{{:s}}}}}\\\\".format(filename))
      print("%u&&&\\\\" % (natoms-nX))
      print("\multicolumn{{4}}{{l}}{{{:s}}}\\\\".format(comment))

      for i in range(natoms):
        if symbol[i]=="XX": continue
        print("{:3s} & {:11.8f} & {:13.8f} & {:13.8f} \\\\".format(symbol[i], pos[i,0], pos[i,1], pos[i,2]))
      print("\\end{longtable}\n")
#     print("\n\\end{document}")


    # default, if no option was specified
    else:
      asym = (2 * 1/w[1] - 1/w[0] - 1/w[2])/(1/w[0] - 1/w[2])

      print(formula)
      print("volume {:.4f}".format(vtot))
      print("mass {:.4f}".format(mtot))
      print("inertia {:.4f} {:.4f} {:.4f}".format(w[0], w[1], w[2]))
      print("asymmetry parameter {:.4f}\n".format(asym))
    W[j][0] = w[0]
    W[j][1] = w[1]
    W[j][2] = w[2]
  fin.close()

# if two input files where specified, calculate geometrical distance
if len(args) == 2 and xflag == 0 and tflag == 0 and vflag == 0 and fflag == 0:
  gdist = pow(pow((1/W[0][0]-1/W[1][0])/(1/W[0][0]), 2) + pow((1/W[0][1]-1/W[1][1])/(1/W[0][1]), 2) + pow((1/W[0][2]-1/W[1][2])/(1/W[0][2]), 2), 0.5)
  print("geometrical distance: {:.4f}".format(gdist))
