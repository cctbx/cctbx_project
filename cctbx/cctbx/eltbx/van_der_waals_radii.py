# Van-der-Waals radii for known elements.
#
# Created: Pavel Afonine.
#
# Sources: 
# 1) Bondi, J.Phys.Chem., 68, 441, 1964  for atoms: 
#      Ag,Ar,As,Au,Br,Cd,Cl,Cu,F,Ga,H,He,Hg,I,In,K,Kr,
#      Li,Mg,Na,Ne,Ni,Pb,Pd,Pt,Se,Si,Sn,Te,Tl,U,Xe,Zn 
# 2) Fokine et al, J. Appl. Cryst. (2003). 36, 352-355 for "protein" atoms:
#    C, O, N, S, P 
# 3) Al, Ca -- MOPAC 2000 Manual 
# 4) Fe -- SURFNET v.1.5 Manual 
#  

class vdw:
  table = {
     "H":  1.20 ,
     "D":  2.00 ,  # default 2.00 A
     "He": 1.40 ,
     "Li": 1.82 ,
     "Be": 2.00 ,  # default 2.00 A
     "B":  2.00 ,  # default 2.00 A
     "C":  1.78 ,
     "N":  1.50 ,
     "O":  1.45 ,
     "F":  1.47 ,
     "Ne": 1.54 ,
     "Na": 2.27 ,
     "Mg": 1.73 ,
     "Al": 2.05 ,  
     "Si": 2.10 ,
     "P":  1.90 ,
     "S":  1.80 ,
     "Cl": 1.75 ,
     "Ar": 1.88 ,
     "K":  2.75 ,
     "Ca": 2.75 ,  
     "Sc": 2.00 ,  # default 2.00 A
     "Ti": 2.00 ,  # default 2.00 A
     "V":  2.00 ,  # default 2.00 A
     "Cr": 2.00 ,  # default 2.00 A
     "Mn": 2.00 ,  # default 2.00 A
     "Fe": 1.10 ,  
     "Co": 2.00 ,  # default 2.00 A
     "Ni": 1.63 ,
     "Cu": 1.40 ,
     "Zn": 1.39 ,
     "Ga": 1.87 ,
     "Ge": 2.00 ,  # default 2.00 A
     "As": 1.85 ,  # default 2.00 A
     "Se": 1.90 ,
     "Br": 1.85 ,
     "Kr": 2.02 ,
     "Rb": 2.00 ,  # default 2.00 A
     "Sr": 2.00 ,  # default 2.00 A
     "Y":  2.00 ,  # default 2.00 A
     "Zr": 2.00 ,  # default 2.00 A
     "Nb": 2.00 ,  # default 2.00 A
     "Mo": 2.00 ,  # default 2.00 A
     "Tc": 2.00 ,  # default 2.00 A
     "Ru": 2.00 ,  # default 2.00 A
     "Rh": 2.00 ,  # default 2.00 A
     "Pd": 1.63 ,
     "Ag": 1.72 ,
     "Cd": 1.58 ,
     "In": 1.93 ,
     "Sn": 2.17 ,
     "Sb": 2.00 ,  # default 2.00 A
     "Te": 2.06 ,  # default 2.00 A
     "I":  1.98 ,
     "Xe": 2.16 ,
     "Cs": 2.00 ,  # default 2.00 A
     "Ba": 2.00 ,  # default 2.00 A
     "La": 2.00 ,  # default 2.00 A
     "Ce": 2.00 ,  # default 2.00 A
     "Pr": 2.00 ,  # default 2.00 A
     "Nd": 2.00 ,  # default 2.00 A
     "Pm": 2.00 ,  # default 2.00 A
     "Sm": 2.00 ,  # default 2.00 A
     "Eu": 2.00 ,  # default 2.00 A
     "Gd": 2.00 ,  # default 2.00 A
     "Tb": 2.00 ,  # default 2.00 A
     "Dy": 2.00 ,  # default 2.00 A
     "Ho": 2.00 ,  # default 2.00 A
     "Er": 2.00 ,  # default 2.00 A
     "Tm": 2.00 ,  # default 2.00 A
     "Yb": 2.00 ,  # default 2.00 A
     "Lu": 2.00 ,  # default 2.00 A
     "Hf": 2.00 ,  # default 2.00 A
     "Ta": 2.00 ,  # default 2.00 A
     "W":  2.00 ,  # default 2.00 A
     "Re": 2.00 ,  # default 2.00 A
     "Os": 2.00 ,  # default 2.00 A
     "Ir": 2.00 ,  # default 2.00 A
     "Pt": 1.72 ,
     "Au": 1.66 ,
     "Hg": 1.55 ,
     "Tl": 1.96 ,
     "Pb": 2.02 ,
     "Bi": 2.00 ,  # default 2.00 A
     "Po": 2.00 ,  # default 2.00 A
     "At": 2.00 ,  # default 2.00 A
     "Rn": 2.00 ,  # default 2.00 A
     "Fr": 2.00 ,  # default 2.00 A
     "Ra": 2.00 ,  # default 2.00 A
     "Ac": 2.00 ,  # default 2.00 A
     "Th": 2.00 ,  # default 2.00 A
     "Pa": 2.00 ,  # default 2.00 A
     "U":  1.86 ,
     "Np": 2.00 ,  # default 2.00 A
     "Pu": 2.00 ,  # default 2.00 A
     "Am": 2.00 ,  # default 2.00 A
     "Cm": 2.00 ,  # default 2.00 A
     "Bk": 2.00 ,  # default 2.00 A
     "Cf": 2.00 ,  # default 2.00 A
     "Es": 2.00 ,  # default 2.00 A
     "Fm": 2.00 ,  # default 2.00 A
     "Md": 2.00 ,  # default 2.00 A
     "No": 2.00 ,  # default 2.00 A
     "Lr": 2.00 ,  # default 2.00 A
      }
