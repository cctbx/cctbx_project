from __future__ import division

angstrom = 'angstrom'
#angstrom = u"\u212Bngstr\u00F6m"
lower_elements =      ("c",  "n",  "o",
                                   "s",
                       )
more_lower_elements = (      "p",
                       "ge", "as", "se",
                           )
lower_elements += more_lower_elements

elements = ('X',
            'H', 'He',
            'Li','Be','B', 'C', 'N' ,'O', 'F', 'Ne',
            'Na','Mg','Al','Si','P' ,'S', 'Cl','Ar',
            'K' ,'Ca',
            'Sc','Ti','V', 'Cr','Mn','Fe','Co','Ni','Cu','Zn',
                      'Ga','Ge','As','Se','Br','Kr',
            'Rb', 'Sr',
            'Y', 'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
                      'In','Sn','Sb','Te','I', 'Xe',
            'Cs','Ba',
            'La',
            'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
                 'Hf','Ta','W', 'Re','Os','Ir','Pt','Au','Hg',
                      'Tl','Pb','Bi','Po','At','Rn',
            'Fr','Ra',
            'Ac',
            'Th','Pa','U', 'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',
                 'Rf','Ha','Sg','Ns','Hs','Mt',
            )

# Valences H-Ba
valences = ( -1,                                       # X
              1,   0,                                  # H, He
              1,   2,   3,  4,  3,  2,  1,  0,
              1,   2,   3,  4,  3,  2,  1,  0,
              1,   2,
             #-1,  -1,  -1, -1, -1, -1, -1, -1, -1, -1, # 1st transition metals
              1,   2,   3,  4,  5,  4,  3,  2,  1,  0, # 1st transition metals
                        3,  4,  3,  2,  1,  0,
              1,   2,
             #-1,  -1,  -1, -1, -1, -1, -1, -1, -1, -1, # 2nd transition metals
              1,   2,   3,  4,  5,  4,  3,  2,  1,  0, # 2nd transition metals
                        3,  4,  3,  2,  1,  0,
              1,   2,
              1,
             -1,  -1,  -1, -1, -1, -1, -1,  # rare earth 1a
             -1,  -1,  -1, -1, -1, -1, -1,  # rare earth 1b
             #     -1,  -1, -1, -1, -1, -1, -1, -1, -1, # 3rd transition metals
                   2,   3,  4,  5,  4,  3,  2,  1,  0, # 3rd transition metals
                        3,  4,  3,  2,  1,  0,
              1,   2,
              1,
             -1,  -1,   0, -1, -1, -1, -1,  # rare earth 2a
             -1,  -1,  -1, -1, -1, -1, -1,  # rare earth 2b
             #     -1,  -1, -1, -1, -1, -1, -1, -1, -1, # 4th transition metals
                   2,   3,  4,  5,  4,  3,  2,  1,  0, # 4th transition metals
                        3,  4,  3,  2,  1,  0,
            )

non_metal_indices = [0,1,2]
for i in range(5,11): non_metal_indices.append(i)
for i in range(14,19): non_metal_indices.append(i)
for i in range(32,37): non_metal_indices.append(i)
for i in range(51,55): non_metal_indices.append(i)
for i in range(84,87): non_metal_indices.append(i)

lone_pairs = ( -1,
                0,  0,
                0,  0,  0,  0,  1,  2,  3,  0,
                0,  0,  0,  0,  1,  2,  3,  0,
                0,  0,
               -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                        0,  0,  1,  2,  3,  0,
               0,   0,
               -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                        0,  0,  1,  2,  3,  0,
               0,   0,
               )

def get_valences(element, atomic_number=None, charge=0):
  assert atomic_number!=0
  assert type(charge) is type(1)
  if atomic_number is None:
    atomic_number = elements.index(element.title())
  rc = valences[atomic_number]
  if rc is None: return []
  if type(rc)==type([]):
    for i, j in enumerate(rc):
      rc[i]+=charge
    return rc
  else:
    rc += charge
    return [rc]

