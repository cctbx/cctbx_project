from __future__ import absolute_import, division, print_function
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
from mmtbx.ligands import electrons
import sys

pdbs = ['''
CRYST1   62.350   62.350   62.350  90.00  90.00  90.00 I 41 3 2
HETATM    1  C01 LIG A   1      -2.020   0.000   0.002  1.00 20.00      A    C
HETATM    2  C02 LIG A   1      -0.493   0.000   0.002  1.00 20.00      A    C
HETATM    3  C03 LIG A   1       0.223   0.000   1.188  1.00 20.00      A    C
HETATM    4  C04 LIG A   1       1.608   0.000   1.148  1.00 20.00      A    C
HETATM    5  N05 LIG A   1       2.255   0.000  -0.004  1.00 20.00      A    N
HETATM    6  C06 LIG A   1       1.608   0.000  -1.156  1.00 20.00      A    C
HETATM    7  C07 LIG A   1       0.223   0.000  -1.183  1.00 20.00      A    C
HETATM    8 H011 LIG A   1      -2.381  -0.726  -0.716  1.00 20.00      A    H
HETATM    9 H012 LIG A   1      -2.381   0.985  -0.268  1.00 20.00      A    H
HETATM   10 H013 LIG A   1      -2.381  -0.259   0.990  1.00 20.00      A    H
HETATM   11 H031 LIG A   1      -0.295   0.000   2.139  1.00 20.00      A    H
HETATM   12 H041 LIG A   1       2.169   0.000   2.075  1.00 20.00      A    H
HETATM   13 H061 LIG A   1       2.165   0.000  -2.086  1.00 20.00      A    H
HETATM   14 H071 LIG A   1      -0.300   0.000  -2.132  1.00 20.00      A    H
''',
'''
CRYST1   62.350   62.350   62.350  90.00  90.00  90.00 I 41 3 2
HETATM    1  F   014 A   1      -5.615   0.186   0.307  1.00 20.00      A    F
HETATM    2 CL   014 A   1      -5.415   0.099  -2.419  1.00 20.00      A   CL
HETATM   24  C1  014 A   1       5.194   0.337   0.293  1.00 20.00      A    C
HETATM    3  N1  014 A   1       1.933  -1.321   0.203  1.00 20.00      A    N
HETATM    4  O1  014 A   1       5.787   0.556   1.382  1.00 20.00      A    O-1
HETATM    5  C2  014 A   1       3.727  -0.086   0.294  1.00 20.00      A    C
HETATM    6  N2  014 A   1       1.549   0.011   0.256  1.00 20.00      A    N
HETATM    7  O2  014 A   1       5.807   0.481  -0.797  1.00 20.00      A    O
HETATM    8  C3  014 A   1       3.260  -1.368   0.227  1.00 20.00      A    C
HETATM    9  N3  014 A   1      -0.785  -0.001   1.194  1.00 20.00      A    N
HETATM   10  C4  014 A   1       2.649   0.754   0.311  1.00 20.00      A    C
HETATM   11  N4  014 A   1      -0.297   0.098  -0.927  1.00 20.00      A    N
HETATM   12  C5  014 A   1       0.229   0.026   0.310  1.00 20.00      A    C
HETATM   13  C6  014 A   1      -1.930   0.054   0.512  1.00 20.00      A    C
HETATM   14  C7  014 A   1      -3.277   0.057   0.926  1.00 20.00      A    C
HETATM   15  C8  014 A   1      -4.290   0.123  -0.015  1.00 20.00      A    C
HETATM   16  C10 014 A   1      -3.981   0.186  -1.361  1.00 20.00      A    C
HETATM   17  C11 014 A   1      -2.658   0.183  -1.768  1.00 20.00      A    C
HETATM   18  C12 014 A   1      -1.626   0.116  -0.810  1.00 20.00      A    C
HETATM   19  H3  014 A   1       3.867  -2.265   0.198  1.00 20.00      A    H
HETATM   20  H4  014 A   1       2.686   1.836   0.348  1.00 20.00      A    H
HETATM   21  H7  014 A   1      -3.539   0.022   1.976  1.00 20.00      A    H
HETATM   22  H11 014 A   1      -2.581  -0.031  -2.827  1.00 20.00      A    H
HETATM   23  HN3 014 A   1      -0.696  -0.053   2.190  1.00 20.00      A    H
''',
'''
CRYST1   62.350   62.350   62.350  90.00  90.00  90.00 I 41 3 2
HETATM    1  C4  03M A   1       3.552  -1.598  -2.729  1.00 20.00      A    C
HETATM    2  C6  03M A   1       2.043  -1.790   0.583  1.00 20.00      A    C
HETATM    3  C7  03M A   1       1.828  -1.480  -0.186  1.00 20.00      A    C
HETATM    4  C8  03M A   1       3.153  -1.891  -1.187  1.00 20.00      A    C
HETATM    5  C10 03M A   1       0.456  -0.918  -0.774  1.00 20.00      A    C
HETATM    6  N12 03M A   1      -0.102   1.744  -0.409  1.00 20.00      A    N
HETATM    7  C13 03M A   1      -1.175   1.811   0.798  1.00 20.00      A    C
HETATM    8  C15 03M A   1      -0.275  -0.720  -0.187  1.00 20.00      A    C
HETATM    9  C20 03M A   1      -3.871   3.297   1.952  1.00 20.00      A    C
HETATM   10  C21 03M A   1      -4.743   3.988   2.814  1.00 20.00      A    C
HETATM   11  C22 03M A   1      -6.158   3.572   2.983  1.00 20.00      A    C
HETATM   12  C24 03M A   1      -5.799   1.605   1.341  1.00 20.00      A    C
HETATM   13  C1  03M A   1       4.559  -1.965  -3.117  1.00 20.00      A    C
HETATM   14  C2  03M A   1       5.370  -2.590  -2.300  1.00 20.00      A    C
HETATM   15  C3  03M A   1       5.013  -2.795  -1.036  1.00 20.00      A    C
HETATM   16  N5  03M A   1       3.322  -2.438   0.648  1.00 20.00      A    N
HETATM   17  C9  03M A   1       3.835  -2.377  -0.591  1.00 20.00      A    C
HETATM   18  C11 03M A   1       0.328  -0.307  -0.690  1.00 20.00      A    C
HETATM   19  N14 03M A   1      -1.549   1.430   0.714  1.00 20.00      A    N
HETATM   41  O16 03M A   1      -1.480   3.068   1.213  1.00 20.00      A    O
HETATM   20  O17 03M A   1      -2.889  -0.556   0.086  1.00 20.00      A    O
HETATM   21  C18 03M A   1      -3.596   1.867   0.066  1.00 20.00      A    C
HETATM   22  C19 03M A   1      -4.353   2.092   1.174  1.00 20.00      A    C
HETATM   23  C23 03M A   1      -6.680   2.295   2.205  1.00 20.00      A    C
HETATM   24  F25 03M A   1      -7.955   2.035   2.218  1.00 20.00      A    F
HETATM   42  F26 03M A   1      -7.009   4.417   3.468  1.00 20.00      A    F
HETATM   25  C27 03M A   1       5.996  -4.017  -1.570  1.00 20.00      A    C
HETATM   26 CL   03M A   1       7.098  -2.466  -2.734  1.00 20.00      A   CL
HETATM   27  H4  03M A   1       2.664  -0.819  -3.592  1.00 20.00      A    H
HETATM   28  H6  03M A   1       1.655  -1.321   1.382  1.00 20.00      A    H
HETATM   29  H10 03M A   1      -1.103  -2.029   0.347  1.00 20.00      A    H
HETATM   30  H20 03M A   1      -2.763   3.498   2.023  1.00 20.00      A    H
HETATM   31  H21 03M A   1      -4.383   4.845   3.311  1.00 20.00      A    H
HETATM   32  H24 03M A   1      -6.063   0.910   0.937  1.00 20.00      A    H
HETATM   33  H1  03M A   1       4.934  -1.776  -4.384  1.00 20.00      A    H
HETATM   34  H18 03M A   1      -3.565   2.898  -0.557  1.00 20.00      A    H
HETATM   35 H18A 03M A   1      -4.040   1.047  -0.639  1.00 20.00      A    H
HETATM   36  H27 03M A   1       6.157  -4.645  -0.842  1.00 20.00      A    H
HETATM   37 H27A 03M A   1       7.019  -3.534  -1.908  1.00 20.00      A    H
HETATM   38 H27B 03M A   1       5.579  -4.455  -2.328  1.00 20.00      A    H
HETATM   39  HN5 03M A   1       3.703  -2.803   1.389  1.00 20.00      A    H
HETATM   40  H14 03M A   1       1.291   2.870   0.108  1.00 20.00      A    H
''',
'''
CRYST1   62.350   62.350   62.350  90.00  90.00  90.00 I 41 3 2
HETATM    1  C01 LIG A   1      -0.017   0.000  -0.006  1.00 20.00      A    C
HETATM    2  N02 LIG A   1       1.116   0.000  -0.006  1.00 20.00      A    N
HETATM    3 H011 LIG A   1      -1.100   0.000   0.013  1.00 20.00      A    H
''',
'''
CRYST1   12.216   10.000   10.019  90.00  90.00  90.00 P 1
SCALE1      0.081860  0.000000  0.000000        0.00000
SCALE2      0.000000  0.100000  0.000000        0.00000
SCALE3      0.000000  0.000000  0.099810        0.00000
HETATM    1  C01 LIG A   1       6.083   5.000   5.000  1.00 20.00      A    C
HETATM    2  N02 LIG A   1       7.216   5.000   5.000  1.00 20.00      A    N
HETATM    3 H011 LIG A   1       5.000   5.000   5.019  1.00 20.00      A    H
''',
'''
CRYST1   16.018   11.795   13.900  90.00  90.00  90.00 P 1
SCALE1      0.062430  0.000000  0.000000        0.00000
SCALE2      0.000000  0.084782  0.000000        0.00000
SCALE3      0.000000  0.000000  0.071942        0.00000
HETATM    1  C01 LIG A   1       5.361   5.910   7.544  1.00 20.00      A    C
HETATM    2  C02 LIG A   1       6.888   5.910   7.544  1.00 20.00      A    C
HETATM    3  C04 LIG A   1       9.411   5.910   7.976  1.00 20.00      A    C
HETATM    4  C06 LIG A   1      10.020   5.903   5.634  1.00 20.00      A    C
HETATM    5  C07 LIG A   1       7.744   5.910   6.455  1.00 20.00      A    C
HETATM    6  N05 LIG A   1       9.042   5.910   6.706  1.00 20.00      A    N
HETATM    7  S03 LIG A   1       7.954   5.910   8.900  1.00 20.00      A    S
HETATM    8 H011 LIG A   1       5.000   5.026   7.033  1.00 20.00      A    H
HETATM    9 H012 LIG A   1       5.000   6.795   7.033  1.00 20.00      A    H
HETATM   10 H013 LIG A   1       5.000   5.910   8.565  1.00 20.00      A    H
HETATM   11 H041 LIG A   1      10.424   5.910   8.359  1.00 20.00      A    H
HETATM   12 H061 LIG A   1      11.018   5.934   6.055  1.00 20.00      A    H
HETATM   13 H062 LIG A   1       9.905   5.000   5.046  1.00 20.00      A    H
HETATM   14 H063 LIG A   1       9.868   6.768   5.000  1.00 20.00      A    H
HETATM   15 H071 LIG A   1       7.366   5.910   5.440  1.00 20.00      A    H
''',
]
cifs = ['''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
LIG        LIG 'Unknown                  ' ligand 14 7 .
#
data_comp_LIG
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
LIG         C01    C   CH3    0    .      -2.0201    0.0000    0.0024
LIG         C02    C   CR6    0    .      -0.4930    0.0000    0.0024
LIG         C03    C   CR16   0    .       0.2233    0.0000    1.1882
LIG         C04    C   CR16   0    .       1.6080    0.0000    1.1478
LIG         N05    N   N      0    .       2.2554    0.0000   -0.0041
LIG         C06    C   CR16   0    .       1.6084    0.0000   -1.1563
LIG         C07    C   CR16   0    .       0.2233    0.0000   -1.1834
LIG        H011    H   HCH3   0    .      -2.3813   -0.7263   -0.7157
LIG        H012    H   HCH3   0    .      -2.3813    0.9850   -0.2676
LIG        H013    H   HCH3   0    .      -2.3813   -0.2587    0.9904
LIG        H031    H   HCR6   0    .      -0.2953    0.0000    2.1393
LIG        H041    H   HCR6   0    .       2.1687    0.0000    2.0747
LIG        H061    H   HCR6   0    .       2.1648    0.0000   -2.0858
LIG        H071    H   HCR6   0    .      -0.2997    0.0000   -2.1321
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
LIG   C02     C01   single        1.527 0.020
LIG   C03     C02   aromatic      1.385 0.020
LIG   C04     C03   aromatic      1.385 0.020
LIG   N05     C04   aromatic      1.321 0.020
LIG   C06     N05   aromatic      1.321 0.020
LIG   C07     C06   aromatic      1.385 0.020
LIG   C02     C07   aromatic      1.385 0.020
LIG  H011     C01   single        1.083 0.020
LIG  H012     C01   single        1.083 0.020
LIG  H013     C01   single        1.083 0.020
LIG  H031     C03   single        1.083 0.020
LIG  H041     C04   single        1.083 0.020
LIG  H061     C06   single        1.083 0.020
LIG  H071     C07   single        1.083 0.020
''',
'', # 014
'', # 03M
'''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
LIG        LIG 'Unknown                  ' ligand 3 2 .
#
data_comp_LIG
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
LIG         C01    C   CSP1   0    .      -0.0166    0.0000   -0.0063
LIG         N02    N   NS     0    .       1.1163    0.0000   -0.0063
LIG        H011    H   H      0    .      -1.0997    0.0000    0.0126
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
LIG   N02     C01   triple        1.133 0.020
LIG  H011     C01   single        1.083 0.020
''',
'''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
LIG        LIG 'Unknown                  ' ligand 3 2 .
#
data_comp_LIG
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
LIG         C01    C   CSP1   0    .      -0.0166    0.0000   -0.0063
LIG         N02    N   NS     0    .       1.1163    0.0000   -0.0063
LIG        H011    H   H      0    .      -1.0997    0.0000    0.0126
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
LIG   N02     C01   triple        1.133 0.020
LIG  H011     C01   single        1.083 0.020
''',
'''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
LIG        LIG 'Unknown                  ' ligand 15 7 .
#
data_comp_LIG
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
LIG         C01    C   CH3    0    .      -2.6390    0.0024    0.6577
LIG         C02    C   CR5    0    .      -1.1119    0.0024    0.6577
LIG         S03    S   S      0    .      -0.0464    0.0024    2.0135
LIG         C04    C   CR15   0    .       1.4106    0.0024    1.0899
LIG         N05    N   NR5    0    .       1.0418    0.0024   -0.1797
LIG         C06    C   CH3    0    .       2.0203   -0.0051   -1.2518
LIG         C07    C   CR15   0    .      -0.2560    0.0024   -0.4313
LIG        H011    H   HCH3   0    .      -3.0002   -0.8821    0.1471
LIG        H012    H   HCH3   0    .      -3.0002    0.8869    0.1471
LIG        H013    H   HCH3   0    .      -3.0001    0.0024    1.6791
LIG        H041    H   HCR5   0    .       2.4239    0.0024    1.4731
LIG        H061    H   HCH3   0    .       3.0178    0.0265   -0.8306
LIG        H062    H   HCH3   0    .       1.9050   -0.9078   -1.8397
LIG        H063    H   HCH3   0    .       1.8684    0.8603   -1.8855
LIG        H071    H   HCR5   0    .      -0.6342    0.0024   -1.4465
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
LIG   C02     C01   single        1.527 0.020
LIG   S03     C02   aromatic      1.724 0.020
LIG   C04     S03   aromatic      1.725 0.020
LIG   N05     C04   aromatic      1.322 0.020
LIG   C06     N05   single        1.452 0.020
LIG   C07     N05   aromatic      1.322 0.020
LIG   C02     C07   aromatic      1.385 0.020
LIG  H011     C01   single        1.083 0.020
LIG  H012     C01   single        1.083 0.020
LIG  H013     C01   single        1.083 0.020
LIG  H041     C04   single        1.083 0.020
LIG  H061     C06   single        1.083 0.020
LIG  H062     C06   single        1.083 0.020
LIG  H063     C06   single        1.083 0.020
LIG  H071     C07   single        1.083 0.020
''',
]
answers = [  0,
            -1, # 014
             0, # 03M
             0, # C#N asu
             0, # C#N
             1,
  ]

def main(only_i=None):
  try: only_i=int(only_i)
  except: only_i=None
  for i, (pdb_str, cif_str) in enumerate(zip(pdbs, cifs)):
    if only_i and only_i!=i+1: continue
    model_filename = 'tst_electron_01_%02d.pdb' % (i+1)
    f=open(model_filename, 'w')
    f.write(pdb_str)
    del f
    restraints_filename = 'tst_electron_01_%02d.cif' % (i+1)
    f=open(restraints_filename, 'w')
    f.write(cif_str)
    del f
    rc = run_program(program_class=electrons.Program,
                     args=[model_filename, restraints_filename],
                     logger=null_out(),
                    )
    print(rc)
    assert answers[i]==rc.total_charge
    # atom_valences = rc.atom_valences
    # for key, item in atom_valences.items():
    #   print(key,item)
  print('OK')

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
