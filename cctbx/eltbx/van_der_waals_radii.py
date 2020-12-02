from __future__ import absolute_import, division, print_function
# Van-der-Waals radii for known elements.
#
# ******************               NOTE                *********************
# vdW radii vary for the same element depending on chemical context
# The list below may be therefore too simplistic
# A complete list of vdW radii is in chem_data/mon_lib/ener_lib.cif
# See method get_vdw_radii() of model class (mmtbx/model/model.py)
# how to get the vdw radii from ener_lib
# ****************** ******************* ******************* *******************
#
# Created: Pavel Afonine.
#
# Sources:
# 1) Bondi, J.Phys.Chem., 68, 441, 1964  for atoms:
#      Ag,Ar,As,Au,Br,Cd,Cl,Cu,F,Ga,H,He,Hg,I,In,K,Kr,
#      Li,Mg,Na,Ne,Ni,Pb,Pd,Pt,Se,Si,Sn,Te,Tl,Xe,Zn
# 2) Fokine et al, J. Appl. Cryst. (2003). 36, 352-355 for "protein" atoms:
#    C, O, N, S, P
#

class vdw(object):
  table = {
     "H":  1.20 , #wiki 1.2
     "D":  1.20 ,
     "He": 1.40 ,
     "Li": 1.82 ,
     "Be": 0.63 ,
     "B":  1.75 ,
     "C":  1.775, #wiki 1.7  #was 1.775 # CNS 2.3
     "N":  1.50 , #wiki 1.55 #was 1.5 # CNS 1.6
     "O":  1.45 , #wiki 1.52 #was 1.45 # CNS 1.6
     "F":  1.47 ,
     "Ne": 1.54 ,
     "Na": 2.27 ,
     "Mg": 1.73 ,
     "Al": 1.50 ,
     "Si": 2.10 ,
     "P":  1.90 ,
     "S":  1.80 , #wiki 1.8 # CNS 1.9
     "Cl": 1.75 ,
     "Ar": 1.88 ,
     "K":  2.75 ,
     "Ca": 1.95 ,
     "Sc": 1.32 ,
     "Ti": 1.95 ,
     "V":  1.06 ,
     "Cr": 1.13 ,
     "Mn": 1.19 ,
     "Fe": 1.26 ,
     "Co": 1.13 ,
     "Ni": 1.63 ,
     "Cu": 1.40 ,
     "Zn": 1.39 ,
     "Ga": 1.87 ,
     "Ge": 1.48 ,
     "As": 0.83 ,
     "Se": 1.90 ,
     "Br": 1.85 ,
     "Kr": 2.02 ,
     "Rb": 2.65 ,
     "Sr": 2.02 ,
     "Y":  1.61 ,
     "Zr": 1.42 ,
     "Nb": 1.33 ,
     "Mo": 1.75 ,
     "Tc": 2.00 ,
     "Ru": 1.20 ,
     "Rh": 1.22 ,
     "Pd": 1.63 ,
     "Ag": 1.72 ,
     "Cd": 1.58 ,
     "In": 1.93 ,
     "Sn": 2.17 ,
     "Sb": 1.12 ,
     "Te": 1.26 ,
     "I":  1.98 ,
     "Xe": 2.16 ,
     "Cs": 3.01 ,
     "Ba": 2.41 ,
     "La": 1.83 ,
     "Ce": 1.86 ,
     "Pr": 1.62 ,
     "Nd": 1.79 ,
     "Pm": 1.76 ,
     "Sm": 1.74 ,
     "Eu": 1.96 ,
     "Gd": 1.69 ,
     "Tb": 1.66 ,
     "Dy": 1.63 ,
     "Ho": 1.61 ,
     "Er": 1.59 ,
     "Tm": 1.57 ,
     "Yb": 1.54 ,
     "Lu": 1.53 ,
     "Hf": 1.40 ,
     "Ta": 1.22 ,
     "W":  1.26 ,
     "Re": 1.30 ,
     "Os": 1.58 ,
     "Ir": 1.22 ,
     "Pt": 1.72 ,
     "Au": 1.66 ,
     "Hg": 1.55 ,
     "Tl": 1.96 ,
     "Pb": 2.02 ,
     "Bi": 1.73 ,
     "Po": 1.21 ,
     "At": 1.12 ,
     "Rn": 2.30 ,
     "Fr": 3.24 ,
     "Ra": 2.57 ,
     "Ac": 2.12 ,
     "Th": 1.84 ,
     "Pa": 1.60 ,
     "U":  1.75 ,
     "Np": 1.71 ,
     "Pu": 1.67 ,
     "Am": 1.66 ,
     "Cm": 1.65 ,
     "Bk": 1.64 ,
     "Cf": 1.63 ,
     "Es": 1.62 ,
     "Fm": 1.61 ,
     "Md": 1.60 ,
     "No": 1.59 ,
     "Lr": 1.58 ,
      }
