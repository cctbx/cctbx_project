from iotbx import pdb
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from libtbx.utils import Sorry, search_for, format_cpu_times, null_out
from libtbx.test_utils import Exception_expected, block_show_diff
import libtbx.load_env
import iotbx.phil
from cStringIO import StringIO
import os

def exercise_rna_v3(mon_lib_srv, ener_lib):
  raw_records = """\
CRYST1  401.998  401.998  175.648  90.00  90.00  90.00 P 41 21 2
ATOM  25432  P     G A1206     198.794 115.083  28.238  1.00 76.49           P
ATOM  25433  OP1   G A1206     198.375 115.264  29.645  1.00 77.28           O
ATOM  25434  OP2   G A1206     200.140 114.542  27.966  1.00 76.41           O
ATOM  25435  O5'   G A1206     197.705 114.233  27.449  1.00 76.56           O
ATOM  25436  C5'   G A1206     196.418 114.761  27.170  1.00 76.79           C
ATOM  25437  C4'   G A1206     195.633 113.833  26.279  1.00 77.12           C
ATOM  25438  O4'   G A1206     196.015 114.053  24.894  1.00 77.27           O
ATOM  25439  C3'   G A1206     195.867 112.344  26.510  1.00 77.37           C
ATOM  25440  O3'   G A1206     195.077 111.810  27.556  1.00 77.97           O
ATOM  25441  C2'   G A1206     195.569 111.736  25.150  1.00 77.22           C
ATOM  25442  O2'   G A1206     194.170 111.623  24.952  1.00 76.79           O
ATOM  25443  C1'   G A1206     196.112 112.819  24.214  1.00 77.65           C
ATOM  25444  N9    G A1206     197.533 112.591  23.877  1.00 78.00           N
ATOM  25445  C8    G A1206     198.581 113.321  24.380  1.00 78.03           C
ATOM  25446  N7    G A1206     199.744 112.915  23.956  1.00 78.32           N
ATOM  25447  C5    G A1206     199.447 111.851  23.119  1.00 77.90           C
ATOM  25448  C6    G A1206     200.326 111.029  22.380  1.00 77.52           C
ATOM  25449  O6    G A1206     201.557 111.100  22.331  1.00 77.98           O
ATOM  25450  N1    G A1206     199.635 110.063  21.665  1.00 77.35           N
ATOM  25451  C2    G A1206     198.270 109.923  21.660  1.00 76.93           C
ATOM  25452  N2    G A1206     197.791 108.932  20.898  1.00 76.75           N
ATOM  25453  N3    G A1206     197.437 110.686  22.346  1.00 77.12           N
ATOM  25454  C4    G A1206     198.088 111.632  23.050  1.00 77.74           C
HETATM25455  P   2MG A1207     195.612 110.560  28.417  1.00 79.39           P
HETATM25456  OP1 2MG A1207     194.636 110.113  29.435  1.00 80.29           O
HETATM25457  OP2 2MG A1207     196.934 111.063  28.860  1.00 79.14           O
HETATM25458  O5' 2MG A1207     195.767 109.399  27.334  1.00 79.02           O
HETATM25459  C5' 2MG A1207     194.639 108.661  26.891  1.00 79.38           C
HETATM25460  C4' 2MG A1207     195.005 107.705  25.788  1.00 79.65           C
HETATM25461  O4' 2MG A1207     195.607 108.447  24.695  1.00 80.37           O
HETATM25462  C3' 2MG A1207     196.032 106.636  26.145  1.00 79.49           C
HETATM25463  O3' 2MG A1207     195.452 105.486  26.729  1.00 79.21           O
HETATM25464  C2' 2MG A1207     196.714 106.355  24.818  1.00 80.16           C
HETATM25465  O2' 2MG A1207     195.930 105.477  24.030  1.00 80.10           O
HETATM25466  C1' 2MG A1207     196.702 107.734  24.168  1.00 80.84           C
HETATM25467  N9  2MG A1207     197.946 108.483  24.449  1.00 22.09           N
HETATM25468  C8  2MG A1207     198.164 109.755  25.238  1.00 22.18           C
HETATM25469  N7  2MG A1207     199.619 110.140  25.258  1.00 22.22           N
HETATM25470  C5  2MG A1207     200.255 109.072  24.452  1.00 22.02           C
HETATM25471  C6  2MG A1207     201.664 108.956  24.144  1.00 21.94           C
HETATM25472  O6  2MG A1207     202.706 109.773  24.490  1.00 21.74           O
HETATM25473  N1  2MG A1207     201.977 107.741  23.302  1.00 22.13           N
HETATM25474  C2  2MG A1207     200.927 106.863  22.894  1.00 22.28           C
HETATM25475  N2  2MG A1207     201.271 105.673  22.035  1.00 22.47           N
HETATM25476  CM2 2MG A1207     200.196 105.202  20.911  1.00 22.37           C
HETATM25477  N3  2MG A1207     199.581 107.007  23.199  1.00 22.01           N
HETATM25478  C4  2MG A1207     199.263 108.076  23.968  1.00 22.00           C
ATOM  25479  P     C A1208     196.141 104.788  27.998  1.00 79.45           P
ATOM  25480  OP1   C A1208     195.091 104.010  28.693  1.00 79.48           O
ATOM  25481  OP2   C A1208     196.850 105.858  28.740  1.00 80.10           O
ATOM  25482  O5'   C A1208     197.196 103.794  27.347  1.00 80.15           O
ATOM  25483  C5'   C A1208     196.897 103.077  26.162  1.00 81.09           C
ATOM  25484  C4'   C A1208     198.159 102.657  25.458  1.00 82.10           C
ATOM  25485  O4'   C A1208     198.738 103.799  24.770  1.00 82.55           O
ATOM  25486  C3'   C A1208     199.268 102.155  26.367  1.00 82.68           C
ATOM  25487  O3'   C A1208     199.120 100.784  26.703  1.00 83.37           O
ATOM  25488  C2'   C A1208     200.542 102.467  25.585  1.00 82.81           C
ATOM  25489  O2'   C A1208     200.823 101.440  24.645  1.00 82.68           O
ATOM  25490  C1'   C A1208     200.147 103.738  24.821  1.00 83.10           C
ATOM  25491  N1    C A1208     200.654 104.979  25.469  1.00 83.75           N
ATOM  25492  C2    C A1208     202.030 105.208  25.539  1.00 84.04           C
ATOM  25493  O2    C A1208     202.787 104.352  25.073  1.00 85.37           O
ATOM  25494  N3    C A1208     202.504 106.336  26.117  1.00 83.20           N
ATOM  25495  C4    C A1208     201.650 107.231  26.606  1.00 83.20           C
ATOM  25496  N4    C A1208     202.152 108.331  27.170  1.00 82.99           N
ATOM  25497  C5    C A1208     200.240 107.041  26.541  1.00 83.43           C
ATOM  25498  C6    C A1208     199.794 105.918  25.966  1.00 83.75           C
ATOM  25499  P     C A1209     199.109 100.325  28.246  1.00 83.81           P
ATOM  25500  OP1   C A1209     198.251  99.114  28.342  1.00 83.30           O
ATOM  25501  OP2   C A1209     198.750 101.522  29.056  1.00 82.92           O
ATOM  25502  O5'   C A1209     200.633  99.930  28.499  1.00 84.20           O
ATOM  25503  C5'   C A1209     201.270  98.949  27.692  1.00 84.37           C
ATOM  25504  C4'   C A1209     202.741  99.244  27.501  1.00 84.82           C
ATOM  25505  O4'   C A1209     202.921 100.574  26.949  1.00 84.19           O
ATOM  25506  C3'   C A1209     203.597  99.255  28.757  1.00 85.06           C
ATOM  25507  O3'   C A1209     203.935  97.953  29.207  1.00 86.28           O
ATOM  25508  C2'   C A1209     204.805 100.090  28.332  1.00 84.19           C
ATOM  25509  O2'   C A1209     205.743  99.293  27.621  1.00 83.94           O
ATOM  25510  C1'   C A1209     204.174 101.087  27.350  1.00 83.23           C
ATOM  25511  N1    C A1209     203.980 102.422  27.942  1.00 81.76           N
ATOM  25512  C2    C A1209     205.090 103.092  28.456  1.00 81.49           C
ATOM  25513  O2    C A1209     206.191 102.525  28.418  1.00 81.76           O
ATOM  25514  N3    C A1209     204.936 104.327  28.991  1.00 81.01           N
ATOM  25515  C4    C A1209     203.734 104.903  29.006  1.00 80.86           C
ATOM  25516  N4    C A1209     203.625 106.122  29.534  1.00 79.98           N
ATOM  25517  C5    C A1209     202.587 104.247  28.470  1.00 81.59           C
ATOM  25518  C6    C A1209     202.754 103.021  27.948  1.00 81.60           C
""".splitlines()
  cif_records = """\
# electronic Ligand Builder and Optimisation Workbench (eLBOW)
#   - a module of PHENIX version dev-550 (Mon Oct  6 00:28:00 2010)
#   - file written: Wed Nov 17 14:20:56 2010
#
#   Input file: sample.pdb
#   Random seed: 3628800
#   Residue: 2MG
#
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
2MG        2MG 'Unknown                  ' ligand 40 24 .
#
data_comp_2MG
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
2MG         P      P   P     .        193.1922  110.7393   27.8108
2MG         OP1    O   OH1   .        193.2303  110.5959   29.4269
2MG         OP2    O   O     .        193.8066  112.0430   27.1100
2MG        O5'     O   O2    .        193.8090  109.3126   27.4983
2MG        C5'     C   CH2   .        195.2238  109.1932   27.1837
2MG        C4'     C   CR15  .        195.5915  108.2815   26.0928
2MG        O4'     O   O     .        196.8078  108.6487   25.5614
2MG        C3'     C   CR15  .        195.7744  106.9078   26.6204
2MG        O3'     O   OH1   .        195.9405  106.9603   27.9901
2MG        C2'     C   CR15  .        196.9455  106.4187   26.0245
2MG        O2'     O   OH1   .        196.6506  105.2541   25.2751
2MG        C1'     C   CR15  .        197.4047  107.4689   25.1175
2MG         N9     N   NR5   .        198.8548  107.5902   25.1810
2MG         C8     C   CR15  .        199.6182  107.8273   26.3017
2MG         N7     N   N     .        200.9440  107.9581   25.8975
2MG         C5     C   CR56  .        201.0035  107.8016   24.5249
2MG         C6     C   CR6   .        202.0016  107.9320   23.5649
2MG         O6     O   OH1   .        203.2939  108.1809   23.9568
2MG         N1     N   N     .        201.6911  107.8327   22.1819
2MG         C2     C   CR6   .        200.3827  107.6027   21.7653
2MG         N2     N   NH1   .        200.1198  107.3807   20.4207
2MG         CM2    C   CH3   .        201.2316  107.4071   19.4380
2MG         N3     N   N     .        199.4107  107.4740   22.6488
2MG         C4     C   CR56  .        199.6703  107.5658   24.0571
2MG        HP1     H   H     .        191.9727  110.4452   27.6374
2MG        HP11    H   HOH1  .        192.6565  111.1706   29.7890
2MG        H5'1    H   HCH2  .        195.6955  108.8897   28.0106
2MG        H5'2    H   HCH2  .        195.5673  110.1131   26.9473
2MG        H4'1    H   HCR5  .        194.9037  108.3265   25.3852
2MG        H3'1    H   HCR5  .        195.0531  106.4463   26.3998
2MG        H3'2    H   HOH1  .        195.4753  106.3269   28.3657
2MG        H2'1    H   HCR5  .        197.5835  106.2824   26.7053
2MG        H2'2    H   HOH1  .        196.7242  104.5157   25.8120
2MG        H1'1    H   HCR5  .        197.2320  107.2825   24.2174
2MG        H81     H   HCR5  .        199.2996  107.8051   27.2342
2MG        H61     H   HOH1  .        203.7851  108.3067   23.2622
2MG        H21     H   HNH1  .        199.3048  107.1092   20.1655
2MG        HM21    H   HCH3  .        201.8876  108.1500   19.6848
2MG        HM22    H   HCH3  .        200.8866  107.5682   18.5764
2MG        HM23    H   HCH3  .        201.6753  106.5837   19.4482
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
2MG   P       OP1   single        1.623 0.020
2MG   P       OP2   double        1.603 0.020
2MG   P      O5'    single        1.586 0.020
2MG  O5'     C5'    single        1.454 0.020
2MG  C5'     C4'    single        1.468 0.020
2MG  C4'     O4'    single        1.377 0.020
2MG  C4'     C3'    single        1.483 0.020
2MG  O4'     C1'    single        1.395 0.020
2MG  C3'     O3'    single        1.381 0.020
2MG  C3'     C2'    single        1.402 0.020
2MG  C2'     O2'    single        1.416 0.020
2MG  C2'     C1'    single        1.462 0.020
2MG  C1'      N9    single        1.457 0.020
2MG   N9      C8    aromatic      1.377 0.020
2MG   N9      C4    aromatic      1.389 0.020
2MG   C8      N7    aromatic      1.392 0.020
2MG   N7      C5    aromatic      1.383 0.020
2MG   C5      C6    aromatic      1.391 0.020
2MG   C5      C4    aromatic      1.432 0.020
2MG   C6      O6    single        1.373 0.020
2MG   C6      N1    aromatic      1.421 0.020
2MG   N1      C2    aromatic      1.392 0.020
2MG   C2      N2    single        1.388 0.020
2MG   C2      N3    aromatic      1.320 0.020
2MG   N3      C4    aromatic      1.435 0.020
2MG   CM2     N2    single        1.484 0.020
2MG  HP1      P     single        1.266 0.020
2MG  HP11     OP1   single        0.889 0.020
2MG  H5'1    C5'    single        0.999 0.020
2MG  H5'2    C5'    single        1.010 0.020
2MG  H4'1    C4'    single        0.988 0.020
2MG  H3'1    C3'    single        0.884 0.020
2MG  H3'2    O3'    single        0.871 0.020
2MG  H2'1    C2'    single        0.943 0.020
2MG  H2'2    O2'    single        0.916 0.020
2MG  H1'1    C1'    single        0.935 0.020
2MG  H81      C8    single        0.986 0.020
2MG  H61      O6    single        0.860 0.020
2MG  H21      N2    single        0.896 0.020
2MG  HM21     CM2   single        1.021 0.020
2MG  HM22     CM2   single        0.942 0.020
2MG  HM23     CM2   single        0.935 0.020
#
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
2MG  C5'     O5'      P           119.67 3.000
2MG   OP2     P       OP1         119.89 3.000
2MG  O5'      P       OP1          96.17 3.000
2MG  O5'      P       OP2         119.78 3.000
2MG  C4'     C5'     O5'          117.08 3.000
2MG  O4'     C4'     C5'          110.01 3.000
2MG  C3'     C4'     C5'          109.98 3.000
2MG  C1'     O4'     C4'          105.97 3.000
2MG  O3'     C3'     C4'          109.42 3.000
2MG  C2'     C3'     C4'          105.96 3.000
2MG  C3'     C4'     O4'          105.98 3.000
2MG  C2'     C1'     O4'          106.01 3.000
2MG   N9     C1'     O4'          109.98 3.000
2MG  O2'     C2'     C3'          109.75 3.000
2MG  C1'     C2'     C3'          105.99 3.000
2MG  C2'     C3'     O3'          109.54 3.000
2MG   N9     C1'     C2'          110.21 3.000
2MG  C1'     C2'     O2'          109.14 3.000
2MG   C8      N9     C1'          127.01 3.000
2MG   C4      N9     C1'          123.21 3.000
2MG   N7      C8      N9          107.94 3.000
2MG   C5      C4      N9          106.21 3.000
2MG   N3      C4      N9          133.57 3.000
2MG   C4      N9      C8          109.65 3.000
2MG   C5      N7      C8          108.57 3.000
2MG   C6      C5      N7          134.85 3.000
2MG   C4      C5      N7          107.63 3.000
2MG   O6      C6      C5          119.70 3.000
2MG   N1      C6      C5          120.55 3.000
2MG   N3      C4      C5          119.96 3.000
2MG   C4      C5      C6          117.25 3.000
2MG   C2      N1      C6          120.54 3.000
2MG   N1      C6      O6          119.74 3.000
2MG   N2      C2      N1          119.62 3.000
2MG   N3      C2      N1          120.53 3.000
2MG   CM2     N2      C2          119.78 3.000
2MG   C4      N3      C2          121.16 3.000
2MG   N3      C2      N2          119.56 3.000
2MG  HP11     OP1     P           109.46 3.000
2MG  HP1      P       OP1          97.96 3.000
2MG  HP1      P       OP2         119.89 3.000
2MG  HP1      P      O5'           97.97 3.000
2MG  H5'1    C5'     O5'          107.77 3.000
2MG  H5'2    C5'     O5'          107.86 3.000
2MG  H4'1    C4'     C5'          109.24 3.000
2MG  H5'1    C5'     C4'          107.94 3.000
2MG  H5'2    C5'     C4'          107.84 3.000
2MG  H3'1    C3'     C4'          107.10 3.000
2MG  H4'1    C4'     O4'          109.05 3.000
2MG  H1'1    C1'     O4'          113.32 3.000
2MG  H4'1    C4'     C3'          112.52 3.000
2MG  H3'2    O3'     C3'          109.63 3.000
2MG  H2'1    C2'     C3'          107.98 3.000
2MG  H3'1    C3'     O3'          111.44 3.000
2MG  H3'1    C3'     C2'          113.16 3.000
2MG  H2'2    O2'     C2'          109.64 3.000
2MG  H1'1    C1'     C2'          113.33 3.000
2MG  H2'1    C2'     O2'          113.84 3.000
2MG  H2'1    C2'     C1'          109.84 3.000
2MG  H1'1    C1'      N9          104.03 3.000
2MG  H81      C8      N9          125.95 3.000
2MG  H81      C8      N7          125.77 3.000
2MG  H61      O6      C6          109.49 3.000
2MG  H21      N2      C2          119.77 3.000
2MG  HM21     CM2     N2          109.52 3.000
2MG  HM22     CM2     N2          109.53 3.000
2MG  HM23     CM2     N2          109.42 3.000
2MG  H21      N2      CM2         119.88 3.000
2MG  H5'2    C5'     H5'1         108.04 3.000
2MG  HM22     CM2    HM21         109.38 3.000
2MG  HM23     CM2    HM21         109.45 3.000
2MG  HM23     CM2    HM22         109.53 3.000
#
loop_
_chem_comp_tor.comp_id
_chem_comp_tor.id
_chem_comp_tor.atom_id_1
_chem_comp_tor.atom_id_2
_chem_comp_tor.atom_id_3
_chem_comp_tor.atom_id_4
_chem_comp_tor.value_angle
_chem_comp_tor.value_angle_esd
_chem_comp_tor.period
2MG CONST_01       C5      N7      C8      N9           -0.02   0.0 0
2MG CONST_02       N7      C5      C4      N9           -0.02   0.0 0
2MG CONST_03       C6      C5      C4      N9         -174.95   0.0 0
2MG CONST_04       C2      N3      C4      N9          173.30   0.0 0
2MG CONST_05       C5      C4      N9      C8            0.00   0.0 0
2MG CONST_06       N3      C4      N9      C8         -173.94   0.0 0
2MG CONST_07       C6      C5      N7      C8          173.67   0.0 0
2MG CONST_08       C4      C5      N7      C8            0.03   0.0 0
2MG CONST_09       C4      N9      C8      N7            0.01   0.0 0
2MG CONST_10       N1      C6      C5      N7         -173.16   0.0 0
2MG CONST_11       N3      C4      C5      N7          174.92   0.0 0
2MG CONST_12       C2      N1      C6      C5           -0.02   0.0 0
2MG CONST_13       C2      N3      C4      C5            0.02   0.0 0
2MG CONST_14       N3      C4      C5      C6           -0.02   0.0 0
2MG CONST_15       N3      C2      N1      C6            0.02   0.0 0
2MG CONST_16       C4      C5      C6      N1            0.02   0.0 0
2MG CONST_17       C4      N3      C2      N1           -0.03   0.0 0
2MG CONST_18       N7      C8      N9     C1'         -176.01   0.0 0
2MG CONST_19       C5      C4      N9     C1'          176.21   0.0 0
2MG CONST_20       N3      C4      N9     C1'            2.27   0.0 0
2MG CONST_21       O6      C6      C5      N7            5.76   0.0 0
2MG CONST_22       N2      C2      N1      C6         -173.82   0.0 0
2MG CONST_23       C4      C5      C6      O6          178.94   0.0 0
2MG CONST_24       C2      N1      C6      O6         -178.94   0.0 0
2MG CONST_25       C4      N3      C2      N2          173.82   0.0 0
2MG CONST_26      H81      C8      N9     C1'           10.40   0.0 0
2MG CONST_27      H81      C8      N7      C5          173.58   0.0 0
2MG CONST_28      H81      C8      N9      C4         -173.58   0.0 0
2MG Var_01        C4'     C5'     O5'      P          -138.25  30.0 3
2MG Var_02        C5'     O5'      P       OP1         -98.80  30.0 3
2MG Var_03        C5'     O5'      P       OP2          30.95  30.0 3
2MG Var_04        O4'     C4'     C5'     O5'          156.12  30.0 3
2MG Var_05        C3'     C4'     C5'     O5'          -87.50  30.0 3
2MG Var_06        C1'     O4'     C4'     C5'          149.03  30.0 1
2MG Var_07        O3'     C3'     C4'     C5'          -17.76  30.0 1
2MG Var_08        C2'     C3'     C4'     C5'         -135.78  30.0 1
2MG Var_09         N9     C1'     O4'     C4'         -151.05  30.0 1
2MG Var_10        O2'     C2'     C3'     C4'         -120.19  30.0 1
2MG Var_11        O3'     C3'     C4'     O4'          101.12  30.0 1
2MG Var_12        O2'     C2'     C1'     O4'          138.94  30.0 1
2MG Var_13         C8      N9     C1'     O4'           60.18  30.0 2
2MG Var_14         C4      N9     C1'     O4'         -115.34  30.0 2
2MG Var_15         N9     C1'     C2'     C3'          139.77  30.0 1
2MG Var_16        O2'     C2'     C3'     O3'          121.88  30.0 1
2MG Var_17        C1'     C2'     C3'     O3'         -120.39  30.0 1
2MG Var_18         C8      N9     C1'     C2'          -56.35  30.0 2
2MG Var_19         C4      N9     C1'     C2'          128.13  30.0 2
2MG Var_20         N9     C1'     C2'     O2'         -102.09  30.0 1
2MG Var_21         CM2     N2      C2      N1           -0.75  30.0 2
2MG Var_22         N3      C2      N2      CM2        -174.66  30.0 2
2MG Var_23        H5'1    C5'     O5'      P            99.94  30.0 3
2MG Var_24        H5'2    C5'     O5'      P           -16.48  30.0 3
2MG Var_25        HP11     OP1     P       OP2          70.99  30.0 3
2MG Var_26        HP11     OP1     P      O5'         -159.34  30.0 3
2MG Var_27        H4'1    C4'     C5'     O5'           36.45  30.0 3
2MG Var_28        HP1      P      O5'     C5'          162.24  30.0 1
2MG Var_29        H3'1    C3'     C4'     C5'          103.17  30.0 1
2MG Var_30        H1'1    C1'     O4'     C4'           93.00  30.0 1
2MG Var_31        H3'2    O3'     C3'     C4'          136.99  30.0 3
2MG Var_32        H2'1    C2'     C3'     C4'          115.20  30.0 1
2MG Var_33        H5'1    C5'     C4'     O4'          -82.15  30.0 2
2MG Var_34        H5'2    C5'     C4'     O4'           34.35  30.0 2
2MG Var_35        H3'1    C3'     C4'     O4'         -137.95  30.0 1
2MG Var_36        H2'1    C2'     C1'     O4'          -95.61  30.0 1
2MG Var_37        H5'1    C5'     C4'     C3'           34.22  30.0 2
2MG Var_38        H5'2    C5'     C4'     C3'          150.73  30.0 2
2MG Var_39        H2'2    O2'     C2'     C3'          -89.07  30.0 3
2MG Var_40        H1'1    C1'     C2'     C3'         -104.12  30.0 1
2MG Var_41        H4'1    C4'     C3'     O3'         -139.78  30.0 1
2MG Var_42        H2'1    C2'     C3'     O3'           -2.74  30.0 1
2MG Var_43        H4'1    C4'     C3'     C2'          102.20  30.0 1
2MG Var_44        H3'2    O3'     C3'     C2'         -107.25  30.0 3
2MG Var_45        H3'1    C3'     C2'     O2'           -3.13  30.0 1
2MG Var_46        H1'1    C1'     C2'     O2'           14.02  30.0 1
2MG Var_47        H4'1    C4'     O4'     C1'          -91.18  30.0 1
2MG Var_48        H3'1    C3'     C2'     C1'          114.61  30.0 1
2MG Var_49        H2'2    O2'     C2'     C1'          155.18  30.0 3
2MG Var_50        H2'1    C2'     C1'      N9           23.37  30.0 1
2MG Var_51        H1'1    C1'      N9      C8         -178.15  30.0 1
2MG Var_52        H61      O6      C6      C5         -174.96  30.0 2
2MG Var_53        H61      O6      C6      N1            3.97  30.0 2
2MG Var_54        H21      N2      C2      N1          170.63  30.0 2
2MG Var_55        HM21     CM2     N2      C2          -38.75  30.0 3
2MG Var_56        HM22     CM2     N2      C2         -158.70  30.0 3
2MG Var_57        HM23     CM2     N2      C2           81.23  30.0 3
2MG Var_58        H21      N2      C2      N3           -3.27  30.0 2
2MG Var_59        H1'1    C1'      N9      C4            6.33  30.0 1
2MG Var_60        HP11     OP1     P      HP1          -60.37  30.0 3
2MG Var_61        H4'1    C4'     C5'     H5'1         158.18  30.0 3
2MG Var_62        H4'1    C4'     C5'     H5'2         -85.32  30.0 3
2MG Var_63        H3'1    C3'     C4'     H4'1         -18.86  30.0 1
2MG Var_64        H3'2    O3'     C3'     H3'1          18.74  30.0 3
2MG Var_65        H2'1    C2'     C3'     H3'1        -127.74  30.0 1
2MG Var_66        H2'2    O2'     C2'     H2'1          32.08  30.0 3
2MG Var_67        H1'1    C1'     C2'     H2'1         139.48  30.0 1
2MG Var_68        HM21     CM2     N2     H21          149.88  30.0 3
2MG Var_69        HM22     CM2     N2     H21           29.93  30.0 3
2MG Var_70        HM23     CM2     N2     H21          -90.14  30.0 3
#
loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
2MG chir_01   P       OP1     OP2    O5'    both
2MG chir_02  C4'     C5'     O4'     C3'    both
2MG chir_03  C3'     C4'     O3'     C2'    both
2MG chir_04  C2'     C3'     O2'     C1'    both
2MG chir_05  C1'     O4'     C2'      N9    both
#
#
loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
2MG plan-1    C1' 0.020
2MG plan-1     N9 0.020
2MG plan-1     C8 0.020
2MG plan-1     N7 0.020
2MG plan-1     C5 0.020
2MG plan-1     C6 0.020
2MG plan-1     O6 0.020
2MG plan-1     N1 0.020
2MG plan-1     C2 0.020
2MG plan-1     N2 0.020
2MG plan-1     N3 0.020
2MG plan-1     C4 0.020
2MG plan-1    H81 0.020
2MG plan-2     C2 0.020
2MG plan-2     N2 0.020
2MG plan-2    CM2 0.020
2MG plan-2    H21 0.020
"""
  import iotbx.cif
  cif_object = iotbx.cif.reader(input_string=cif_records).model()
  mon_lib_srv.process_cif_object(cif_object=cif_object)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    force_symmetry=True)
  grm = processed_pdb_file.geometry_restraints_manager()
  mod_base_5p_link_op1 = False
  mod_base_5p_link_op2 = False
  for ap in grm.angle_proxies:
    if ap.i_seqs == (8,23,24):
      mod_base_5p_link_op1 = True
    if ap.i_seqs == (8,23,25):
      mod_base_5p_link_op2 = True
  assert mod_base_5p_link_op1
  assert mod_base_5p_link_op2

def exercise_pdb_string(mon_lib_srv, ener_lib):
  raw_records = """\
CRYST1   50.066   67.126   47.862  90.00  92.41  90.00 P 1 21 1
ATOM      0  N   MET     0      18.670  12.527  40.988  1.00 52.89           N
ATOM      1  CA  MET     0      17.631  13.191  40.117  1.00 52.89           C
ATOM      2  CB  MET     0      18.110  13.132  38.643  1.00 52.89           C
ATOM      3  CG  MET     0      18.621  11.738  38.198  1.00 52.60           C
ATOM      4  SD  MET     0      19.279  11.706  36.483  1.00 52.89           S
ATOM      5  CE  MET     0      20.263  13.253  36.435  1.00 52.89           C
ATOM      6  C   MET     0      16.228  12.530  40.276  1.00 50.97           C
ATOM      7  O   MET     0      16.082  11.507  40.963  1.00 52.89           O
ATOM      8  HA  MET     0      17.568  14.136  40.375  1.00 52.89           D
ATOM      9  HB1 MET     0      18.846  13.755  38.537  1.00 45.86           D
ATOM     10  HB2 MET     0      17.375  13.382  38.060  1.00 52.35           D
ATOM     11  HG1 MET     0      17.893  11.102  38.251  1.00 47.31           D
ATOM     12  HG2 MET     0      19.339  11.465  38.791  1.00 52.89           D
ATOM     13  HE1 MET     0      20.554  13.413  35.532  1.00 52.89           D
ATOM     14  HE2 MET     0      21.028  13.143  37.008  1.00 52.89           D
ATOM     15  HE3 MET     0      19.722  13.985  36.740  1.00 52.89           D
ATOM     16  H1  MET     0      19.528  12.921  40.813  1.00 52.89           D
ATOM     17  H2  MET     0      18.706  11.586  40.792  1.00 52.89           D
ATOM     18  H3  MET     0      18.445  12.648  41.914  1.00 52.89           D
END
""".splitlines()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    force_symmetry=True)
  xray_structure = processed_pdb_file.xray_structure()
  assert xray_structure.scattering_type_registry().type_count_dict() \
      == {"S": 1, "N": 1, "C": 5, "O": 1, "D": 11}
  raw_records = [line[:66] for line in raw_records]
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    force_symmetry=True)
  xray_structure = processed_pdb_file.xray_structure()
  assert xray_structure.scattering_type_registry().type_count_dict() \
      == {"S": 1, "N": 1, "C": 5, "O": 1, "H": 11}
  #
  raw_records = """\
HETATM    1  N   NH3     1       0.000   0.000   0.000  1.00  0.00           N
HETATM    2  N   NH3     2       0.000   0.000   0.000  1.00  0.00
HETATM    3 NH3  NH3     3       0.000   0.000   0.000  1.00  0.00
HETATM    4  X   CH4     4       0.000   0.000   0.000  1.00  0.00           C
HETATM    5  C   CH4     5       0.000   0.000   0.000  1.00  0.00
HETATM    6  CH4 CH4     6       0.000   0.000   0.000  1.00  0.00
HETATM    7 U    U       7       0.000   0.000   0.000  1.00  0.00           U
END
""".splitlines()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    log=log)
  looking_for = "Ad-hoc single atom residues: "
  for line in log.getvalue().splitlines():
    line = line.strip()
    if (not line.startswith(looking_for)): continue
    counts = eval(line[len(looking_for):])
    assert counts == {"CH4": 3, "NH3": 3, "U  ": 1}
    break
  else:
    raise RuntimeError('Expected string not found in output: "%s"'
      % looking_for)
  #
  raw_records = """\
CRYST1  106.820   62.340  114.190  90.00  90.00  90.00 P 21 21 21
ATOM      7  N   SER A   4      64.059  32.579  18.554  1.00 14.95
ATOM      8  CA ASER A   4      64.561  31.418  17.802  1.00  9.04
ATOM      9  CA BSER A   4      64.804  31.635  18.406  1.00 12.07
ATOM     10  C   SER A   4      65.282  31.515  17.187  1.00 10.69
ATOM     38  N   ALA C   3     109.043  27.391  28.663  1.00 15.05
ATOM     39  CA  ALA C   3     109.073  26.531  28.433  1.00  3.14
ATOM     40  C   ALA C   3     108.930  26.867  26.637  1.00 15.22
ATOM     41  N   SER C   4     109.311  25.012  25.749  1.00  6.18
ATOM     42  CA BSER C   4     109.271  25.165  25.154  1.00 11.32
ATOM     43  CA ASER C   4     109.508  25.337  25.273  1.00  4.97
ATOM     44  C   SER C   4     109.045  24.408  24.187  1.00 15.28
END
""".splitlines()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    log=log)
  lines = search_for(
    pattern=" +bond proxies already assigned to first conformer: 3$",
    mode="re.match",
    lines=log.getvalue().splitlines())
  assert len(lines) == 1

  # test the atom_selection feature
  log = StringIO()
  pdb_inp = pdb.input(source_info = None, lines = raw_records)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    atom_selection_string="resname ser",
    pdb_inp = pdb_inp,
    crystal_symmetry=pdb_inp.crystal_symmetry(),
    log=log)
  processed_pdb_file.geometry_restraints_manager()
  lines = search_for(
    pattern=" +Nonbonded interactions: 0$",
    mode="re.match",
    lines=log.getvalue().splitlines())
  assert len(lines) == 1

  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    atom_selection_string="resname ser or resname ala",
    pdb_inp = pdb_inp,
    crystal_symmetry=pdb_inp.crystal_symmetry(),
    log=log)
  processed_pdb_file.geometry_restraints_manager()
  lines = search_for(
    pattern=" +Nonbonded interactions: 7$",
    mode="re.match",
    lines=log.getvalue().splitlines())
  assert len(lines) == 1
  #
  raw_records = """\
CRYST1   53.910   23.100   23.100  90.00 110.40  90.00 C 2
ATOM    561  N   TYR X  39      11.814  -3.041  22.005  1.00 15.52           N
ATOM    562  CA  TYR X  39      13.044  -2.727  21.259  1.00 17.08           C
ATOM    563  C   TYR X  39      13.804  -1.593  21.948  1.00 19.46           C
ATOM    564  O   TYR X  39      14.993  -1.524  21.704  1.00 28.79           O
ATOM    582  N   VAL X  40      13.093  -0.787  22.711  1.00 23.32           N
ATOM    583  CA  VAL X  40      13.567   0.356  23.491  1.00 37.02           C
ATOM    584  C   VAL X  40      13.408   0.141  24.979  1.00 47.42           C
ATOM    585  O   VAL X  40      13.919   0.926  25.809  1.00 58.31           O
ATOM    589  OXT VAL X  40      12.720  -0.842  25.372  1.00 48.85           O
HETATM  600  C1 AEOH X 200      12.974   5.558  13.017  0.54 29.75           C
HETATM  601  C1 BEOH X 200      14.446   4.322  12.408  0.46 32.58           C
HETATM  602  C2 AEOH X 200      12.259   6.535  13.707  0.54 26.57           C
HETATM  603  C2 BEOH X 200      13.905   4.879  13.576  0.46 32.29           C
END
""".splitlines()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    log=log)
  lines = search_for(
    pattern="""\
 +Inner-chain residues flagged as termini: \\['pdbres="VAL X  40 "'\\]$""",
    mode="re.match",
    lines=log.getvalue().splitlines())
  assert len(lines) == 2
  #
  raw_records = """\
CRYST1   53.910   23.100   23.100  90.00 110.40  90.00 C 2
ATOM    264  N   GLU    21       3.822  -2.789  -0.368  1.00 11.74           N
ATOM    265  CA  GLU    21       5.032  -3.595  -0.485  1.00 11.86           C
ATOM    266  C   GLU    21       6.226  -2.698  -0.821  1.00 10.70           C
ATOM    267  O   GLU    21       6.118  -1.718  -1.512  1.00 11.93           O
ATOM    268  CB AGLU    21       4.980  -4.649  -1.566  0.32 13.49           C
ATOM    269  CB BGLU    21       4.769  -4.484  -1.766  0.68 14.37           C
ATOM    270  CG AGLU    21       3.942  -5.699  -1.498  0.32 13.60           C
ATOM    271  CG BGLU    21       5.737  -5.573  -1.952  0.68 19.43           C
ATOM    272  CD AGLU    21       4.064  -6.621  -2.739  0.32 15.27           C
ATOM    273  CD BGLU    21       5.251  -6.534  -3.036  0.68 19.91           C
ATOM    274  OE1AGLU    21       5.043  -6.411  -3.478  0.32 16.46           O
ATOM    276  OE2AGLU    21       3.193  -7.474  -2.870  0.32 19.03           O
ATOM    277  OE2BGLU    21       6.037  -6.819  -3.934  0.68 23.06           O
END
""".splitlines()
  for i_pass in [0,1]:
    if (i_pass == 0):
      pattern = " +2.60 -     3.06: 2$"
      ci_counts = {0: 4, 1: 5, 2: 4}
    else:
      raw_records = raw_records[:1] + raw_records[5:]
      pattern = " +2.60 -     2.80: 1$"
      ci_counts = {1: 5, 2: 4}
    log = StringIO()
    processed_pdb_file = monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=None,
      raw_records=raw_records,
      log=log)
    grm = processed_pdb_file.geometry_restraints_manager()
    assert dict(grm.conformer_indices.counts()) == ci_counts
    lines = search_for(
      pattern=pattern, mode="re.match", lines=log.getvalue().splitlines())
    assert len(lines) == 1
  #
  raw_records = """\
REMARK    HYDROLASE                               05-SEP-07   2VB1
CRYST1   27.070   31.250   33.760  87.98 108.00 112.11 P 1           1
ATOM      1  N   LYS A   1       1.984   5.113  14.226  1.00  6.14           N
ATOM      2  CA  LYS A   1       1.811   6.069  13.092  1.00  5.23           C
ATOM      3  C   LYS A   1       2.502   7.339  13.503  1.00  4.84           C
ATOM      4  O   LYS A   1       2.414   7.746  14.659  1.00  5.74           O
ATOM      5  CB  LYS A   1       0.325   6.305  12.891  1.00  5.48           C
ATOM      6  CG  LYS A   1       0.013   7.371  11.851  1.00  5.12           C
ATOM      7  CD  LYS A   1      -1.494   7.455  11.617  1.00  5.56           C
ATOM      8  HA  LYS A   1       2.215   5.709  12.275  1.00  6.28           H
ATOM      9  HB2 LYS A   1      -0.091   5.472  12.620  1.00  6.57           H
ATOM     10  HB3 LYS A   1      -0.068   6.569  13.738  1.00  6.57           H
ATOM     11  HG2 LYS A   1       0.344   8.231  12.157  1.00  6.14           H
ATOM     12  HG3 LYS A   1       0.461   7.155  11.018  1.00  6.14           H
ATOM     13  CE ALYS A   1      -1.966   8.606  10.745  0.69  5.10           C
ATOM     14  NZ ALYS A   1      -1.548   8.473   9.287  0.69  4.56           N
ATOM     15  HD2ALYS A   1      -1.786   6.625  11.210  0.50  6.67           H
ATOM     16  HD3ALYS A   1      -1.933   7.526  12.479  0.50  6.67           H
ATOM     17  HE2ALYS A   1      -2.934   8.658  10.791  0.69  6.12           H
ATOM     18  HE3ALYS A   1      -1.610   9.436  11.100  0.69  6.12           H
ATOM     19  HZ1ALYS A   1      -1.884   7.721   8.949  0.69  6.84           H
ATOM     20  HZ2ALYS A   1      -1.854   9.169   8.825  0.69  6.84           H
ATOM     21  HZ3ALYS A   1      -0.659   8.450   9.234  0.69  6.84           H
ATOM     22  CE BLYS A   1      -1.791   8.575  10.680  0.31  5.79           C
ATOM     23  NZ BLYS A   1      -3.148   8.376  10.084  0.31  6.86           N
ATOM     24  H1 BLYS A   1       2.852   4.975  14.368  0.50  9.20           H
ATOM     25  H2 BLYS A   1       1.614   5.454  14.960  0.50  9.20           H
ATOM     26  H3 BLYS A   1       1.589   4.341  14.026  0.50  9.20           H
ATOM     27  HE2BLYS A   1      -1.763   9.419  11.157  0.31  6.95           H
ATOM     28  HE3BLYS A   1      -1.124   8.601   9.977  0.31  6.95           H
ATOM     29  HZ1BLYS A   1      -3.755   8.339  10.733  0.31 10.29           H
ATOM     30  HZ2BLYS A   1      -3.335   9.055   9.540  0.31 10.29           H
ATOM     31  HZ3BLYS A   1      -3.160   7.615   9.622  0.31 10.29           H
END
""".splitlines()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    log=log)
  lines = search_for(
    pattern="""\
  Number of resolved chirality restraint conflicts: 1""",
    mode="==",
    lines=log.getvalue().splitlines())
  assert len(lines) == 1
  #
  raw_records = """\
CRYST1  109.350  109.350  190.680  90.00  90.00 120.00 H 3 2        18
HETATM 2316 MG    MG   401      85.173  71.732  16.992  1.00 77.96          MG
HETATM 2317 MG    MG   402      79.003  69.700  14.150  1.00 73.20          MG
HETATM 2318 MG    MG   403       9.596  62.411  13.402  1.00 62.56          MG
HETATM 2319 MG    MG   404      48.026  73.068  26.732  1.00 72.59          MG
HETATM 2320 MG    MG   405      63.623  79.773  28.694  1.00 86.62          MG
HETATM 2321 MG    MG   406      25.463  68.205   2.775  1.00 67.84          MG
HETATM 2322 MG    MG   407      64.693  59.956   7.476  1.00 81.20          MG
HETATM 2323 MG    MG   408      58.001  61.900  13.126  1.00 76.59          MG
HETATM 2324 MG    MG   409      79.965  84.253  23.275  1.00 73.36          MG
HETATM 2325 MG    MG   410      22.587  66.575  19.552  1.00 62.10          MG
HETATM 2326 MG    MG   411      53.007  73.668  38.880  1.00 89.72          MG
HETATM 2327 MG    MG   412      51.568  77.963   9.128  1.00 87.46          MG
HETATM 2328  O   HOH   301      32.853  56.507   0.571  1.00 41.94           O
HETATM 2329  O   HOH   302      27.055  59.115  12.335  1.00 46.41           O
HETATM 2330  O   HOH   303      21.956  61.533  10.963  1.00 42.68           O
END
""".splitlines()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    log=log)
  lines = search_for(
    pattern="""\
          Link IDs: {None: 14}""",
    mode="==",
    lines=log.getvalue().splitlines())
  assert len(lines) == 1
  #
  raw_records = """\
HETATM 1253 CD    CD X 200      -0.849  19.743   0.211  1.00 14.82
END
""".splitlines()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    log=log)
  processed_pdb_file.xray_structure()
  lines = search_for(
    pattern="""\
     Cd      1     47.96""",
    mode="==",
    lines=log.getvalue().splitlines())
  assert len(lines) == 1
  #
  raw_records = """\
HETATM 1253 ZN    ZN X 200      -0.849  19.743   0.211  1.00 14.82
END
""".splitlines()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    log=log)
  processed_pdb_file.xray_structure()
  lines = search_for(
    pattern="""\
     Zn      1     29.99""",
    mode="==",
    lines=log.getvalue().splitlines())
  assert len(lines) == 1
  #
  raw_records = """\
CRYST1   32.100   79.470   91.830  90.00  90.00  90.00 I 2 2 2
HETATM  709 NA    NA A1398     -17.524  -7.858 -16.234  1.00 17.29          NA +
HETATM  710 NA    NA A1399     -11.813  -6.045 -17.742  1.00 25.01          NA1
HETATM  711 ZN    ZN A1400      -1.928 -11.394 -27.827  1.00 17.21          ZN+2
HETATM  712 ZN    ZN A1401      -2.733  -3.868 -16.577  1.00 20.27          ZN2+
""".splitlines()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    raw_records=raw_records)
  xray_structure = processed_pdb_file.xray_structure()
  assert xray_structure.scattering_type_registry().type_count_dict() \
      == {"Na": 1, "Na1+": 1, "Zn2+": 2}
  #
  raw_records = """\
CRYST1   14.600   26.100   29.200  90.00  90.00  90.00 P 21 21 21    4
ATOM    107  N   CYS A  16      -0.448  11.073   0.703  1.00  5.42
ATOM    108  CA  CYS A  16      -1.464  10.253   1.359  1.00  6.46
ATOM    109  C   CYS A  16      -2.718  11.036   1.746  1.00  6.39
ATOM    110  O   CYS A  16      -3.805  10.424   1.754  1.00  7.99
ATOM    111  CB  CYS A  16      -0.881   9.558   2.593  1.00  6.56
ATOM    112  SG  CYS A  16       0.273   8.212   2.183  1.00  6.50
ATOM    113  H   CYS A  16       0.296  11.370   1.301  1.00  5.42
ATOM    114  HA  CYS A  16      -1.774   9.501   0.619  1.00  6.46
ATOM    115  HB2 CYS A  16      -0.360  10.303   3.212  1.00  6.56
ATOM    116  HB3 CYS A  16      -1.704   9.154   3.201  1.00  6.56
HETATM  117  N   NH2 A  17      -2.607  12.317   2.071  1.00  6.64
HETATM  118  HN1 NH2 A  17      -1.710  12.757   2.069  1.00  6.64
HETATM  119  HN2 NH2 A  17      -3.421  12.843   2.318  1.00  6.64
""".splitlines()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    raw_records=raw_records,
    log=log)
  processed_pdb_file.xray_structure()
  lines=search_for(
    pattern="""\
          Link IDs: {'NH2_CTERM': 1}""",
    mode="==",
    lines=log.getvalue().splitlines())
  assert len(lines) == 1

def exercise_rna(
      mon_lib_srv,
      ener_lib,
      file_name,
      expected_block,
      expected_block_last_startswith=True,
      expected_modifications_used=None):
  file_path = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/"+file_name,
    test=os.path.isfile)
  if (file_path is None):
    print 'Skipping exercise_rna("%s"): input file not available' % file_name
    return
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=file_path,
    log=log)
  lines = []
  lines_modifications_used = []
  for line in log.getvalue().splitlines():
    if (line.startswith("""\
          Modifications used: {""")):
      lines_modifications_used.append(line)
    else:
      lines.append(line)
  assert not block_show_diff(
    lines, expected_block, last_startswith=expected_block_last_startswith)
  if (expected_modifications_used is None):
    assert len(lines_modifications_used) == 0
  else:
    assert len(lines_modifications_used) == len(expected_modifications_used)
    for line,expected in zip(
          lines_modifications_used, expected_modifications_used):
      modifications_used = eval(line.split(":", 1)[1])
      assert modifications_used == expected

def exercise_cns_rna(mon_lib_srv, ener_lib):
  exercise_rna(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name="cns_rna.pdb",
    expected_block= """\
  Total number of atoms: 646
  Number of models: 1
  Model: ""
    Number of chains: 1
    Chain: " "
      Number of atoms: 646
      Number of conformers: 1
      Conformer: ""
        Number of residues, atoms: 20, 646
          Classifications: {'RNA': 20}
          Link IDs: {'rna3p': 19}
  Time building chain proxies: """,
    expected_modifications_used=\
      [{'rna3p_pyr': 9, 'p5*END': 1, 'rna3p_pur': 11, '3*END': 1}])

def exercise_rna_3p_2p(mon_lib_srv, ener_lib):
  exercise_rna(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name="rna_3p_2p.pdb",
    expected_block= """\
  Total number of atoms: 63
  Number of models: 1
  Model: ""
    Number of chains: 1
    Chain: "A"
      Number of atoms: 63
      Number of conformers: 1
      Conformer: ""
        Number of residues, atoms: 3, 63
          Classifications: {'RNA': 3}
          Link IDs: {'rna3p': 1, 'rna2p': 1}
  Residues with excluded nonbonded symmetry interactions: 2
    residue:
      pdb=" P     C A 857 " occ=0.30
      ... (18 atoms not shown)
      pdb=" C6    C A 857 " occ=0.30
    residue:
      pdb=" P     U A 858 " occ=0.30
      ... (18 atoms not shown)
      pdb=" C6    U A 858 " occ=0.30
  Time building chain proxies: """,
    expected_modifications_used=[{'rna3p_pyr': 1, 'rna3p_pur': 1, 'rna2p_pyr': 1}])
  #
  exercise_rna(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name="coords_rna_gap.pdb",
    expected_block= """\
  Total number of atoms: 90
  Number of models: 1
  Model: ""
    Number of chains: 1
    Chain: "A"
      Number of atoms: 90
      Number of conformers: 1
      Conformer: ""
        Number of residues, atoms: 5, 90
          Classifications: {'RNA': 5}
          Link IDs: {'rna3p': 4}
          Chain breaks: 1
""",
    expected_block_last_startswith=False,
    expected_modifications_used=[{'rna3p_pyr': 2, 'rna3p_pur': 2}])

def exercise_rna_dna_hybrid(mon_lib_srv, ener_lib):
  exercise_rna(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name="da_a.pdb",
    expected_block= """\
  Total number of atoms: 63
  Number of models: 1
  Model: ""
    Number of chains: 2
    Chain: "A"
      Number of atoms: 30
      Number of conformers: 1
      Conformer: ""
        Number of residues, atoms: 1, 30
          Classifications: {'DNA': 1}
    Chain: "B"
      Number of atoms: 33
      Number of conformers: 1
      Conformer: ""
        Number of residues, atoms: 1, 33
          Classifications: {'RNA': 1}
  Time building chain proxies: """,
    expected_modifications_used=[{'5*END': 1}, {'rna3p_pur': 1}])

def exercise_hydrogen_deuterium_aliases():
  file_paths = []
  for file_name in ["NAD_594_HD.pdb", "NAD_594_HD.cif"]:
    file_path = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/pdb/"+file_name,
      test=os.path.isfile)
    if (file_path is None):
      print "Skipping exercise_hydrogen_deuterium_aliases():", \
        "input file not available:", file_name
      return
    file_paths.append(file_path)
  log = StringIO()
  monomer_library.pdb_interpretation.run(args=file_paths, log=log)
  assert not block_show_diff(
    log.getvalue(), """\
  Histogram of bond lengths:
        0.94 -     1.13: 56
        1.13 -     1.32: 4
        1.32 -     1.51: 56
        1.51 -     1.69: 18
        1.69 -     1.88: 1
  Bond restraints: 135
""")

def exercise_corrupt_cif_link():
  file_paths = []
  for file_name in ["hem_no.pdb", "hem_no.cif"]:
    file_path = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/misc/"+file_name,
      test=os.path.isfile)
    if (file_path is None):
      print "Skipping exercise_corrupt_cif_link():", \
        "input file not available:", file_name
      return
    file_paths.append(file_path)
  log = StringIO()
  try:
    monomer_library.pdb_interpretation.run(args=file_paths, log=log)
  except Sorry, e:
    assert str(e).startswith("Corrupt CIF link definition:")
  else: raise Exception_expected

def exercise_dna_cns_cy5_th6():
  file_paths = []
  for file_name in ["dna_cns_cy5_th6.pdb",
                    "dna_cns_cy5.cif",
                    "dna_cns_th6.cif"]:
    file_path = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/misc/"+file_name,
      test=os.path.isfile)
    if (file_path is None):
      print "Skipping exercise_dna_cns_cy5_th6():", \
        "input file not available:", file_name
      return
    file_paths.append(file_path)
  log = StringIO()
  monomer_library.pdb_interpretation.run(args=file_paths, log=log)
  assert not block_show_diff(log.getvalue(), """\
        Number of residues, atoms: 12, 244
          Classifications: %s
          Modifications used: {'5*END': 1}
          Link IDs: %s
""" % (
  str({'DNA': 12}),
  str({'rna3p': 11})))

def exercise_sym_excl_indices(mon_lib_srv, ener_lib):
  file_name = "so4_on_two_fold.pdb"
  file_path = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/"+file_name,
    test=os.path.isfile)
  if (file_path is None):
    print "Skipping exercise_sym_excl_indices():", \
      "input file not available:", file_name
    return
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=file_path,
    log=log)
  processed_pdb_file.geometry_restraints_manager()
  lines = log.getvalue().splitlines()
  assert not block_show_diff(
    lines, """\
  Residues with excluded nonbonded symmetry interactions: 1
    residue:
      pdb=" S   SO4 A  13 " occ=0.50
      ... (3 atoms not shown)
      pdb=" O4  SO4 A  13 " occ=0.50
  Time building chain proxies:
""", last_startswith=True)
  assert not block_show_diff(
    lines, """\
  Sorted by model distance:
  nonbonded pdb=" OE1 GLU A   9 "
            pdb=" O4  SO4 A  13 "
     model   vdw sym.op.
     2.134 3.040 z+1/4,y-1/4,-x+3/4
""")

def exercise_auto_alias_h_h1():
  file_paths = []
  for file_name in ["dpn.pdb", "dpn.cif"]:
    file_path = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/pdb/"+file_name,
      test=os.path.isfile)
    if (file_path is None):
      print "Skipping exercise_auto_alias_h_h1():", \
        "input file not available:", file_name
      return
    file_paths.append(file_path)
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.run(
    args=file_paths, log=log)
  assert log.getvalue().find("Modifications used: {'NH3': 1}") >= 0
  assert processed_pdb_file.all_chain_proxies.fatal_problems_message() is None

def exercise_d_aa_resnames():
  file_path = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/dpn.pdb",
    test=os.path.isfile)
  if (file_path is None):
    print "Skipping exercise_d_aa_resnames():", \
      "input file not available: dpn.pdb"
    return
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.run(
    args=[file_path], log=log)
  assert log.getvalue().find("'PEPT-D': 1") >= 0
  assert log.getvalue().find("'NH1NOTPRO': 1") >= 0
  assert processed_pdb_file.all_chain_proxies.fatal_problems_message() is None

def exercise_d_amino_acid_chain_perfect_in_box():
  file_path = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/d_amino_acid_chain_perfect_in_box.pdb",
    test=os.path.isfile)
  if (file_path is None):
    print "Skipping exercise_d_aa_resnames():", \
      "input file not available: d_amino_acid_chain_perfect_in_box.pdb"
    return
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.run(
    args=[file_path], log=log)
  grm = processed_pdb_file.geometry_restraints_manager()
  lv = log.getvalue()
  assert lv.find("Classifications: {'peptide': 16}") >= 0
  assert lv.find("'PEPT-D': 1") >= 0
  assert lv.find("'TRANS': 14") >= 0
  assert lv.find("'PCIS': 1") >= 0
  assert lv.find("""\
Simple disulfide: pdb=" SG  DCY A   4 " - pdb=" SG  DCY A  19 " distance=2.03
""") >= 0
  assert not block_show_diff(lv, """\
  Bond restraints: 122
  Sorted by residual:
  bond pdb=" CG  DPR A   7 "
       pdb=" CD  DPR A   7 "
    ideal  model  delta    sigma   weight residual
    1.503  1.507 -0.004 3.40e-02 8.65e+02 1.39e-02
""")
  assert not block_show_diff(lv, """\
  Bond angle restraints: 161
  Sorted by residual:
  angle pdb=" CA  DPR A   7 "
        pdb=" C   DPR A   7 "
        pdb=" O   DPR A   7 "
      ideal   model   delta    sigma   weight residual
     119.00  119.98   -0.98 3.00e+00 1.11e-01 1.07e-01
""")
  assert not block_show_diff(lv, """\
  Dihedral angle restraints: 47
    sinusoidal: 31
      harmonic: 16
  Sorted by residual:
  dihedral pdb=" N   DPR A   7 "
           pdb=" CG  DPR A   7 "
           pdb=" CD  DPR A   7 "
           pdb=" CB  DPR A   7 "
      ideal   model   delta sinusoidal    sigma   weight residual
      30.00   25.28    4.72     1      1.50e+01 4.44e-03 1.45e-01
""")
  assert processed_pdb_file.all_chain_proxies.fatal_problems_message() is None
  log = StringIO()
  acp = processed_pdb_file.all_chain_proxies
  grm.show_sorted(
    sites_cart=acp.sites_cart_exact(),
    site_labels=[atom.id_str() for atom in acp.pdb_atoms],
    f=log)
  lv = log.getvalue()
  assert not block_show_diff(lv, """\
chirality pdb=" CB  DIL A  10 "
          pdb=" CA  DIL A  10 "
          pdb=" CG1 DIL A  10 "
          pdb=" CG2 DIL A  10 "
  both_signs  ideal   model   delta    sigma   weight residual
    True       2.64    2.65   -0.00 2.00e-01 2.50e+01 5.99e-06
""")
  assert not block_show_diff(lv, """\
chirality pdb=" CB  DTH A  18 "
          pdb=" CA  DTH A  18 "
          pdb=" CG2 DTH A  18 "
          pdb=" OG1 DTH A  18 "
  both_signs  ideal   model   delta    sigma   weight residual
    True       2.55   -2.55   -0.00 2.00e-01 2.50e+01 4.37e-06
""")

def exercise_scale_restraints () :
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  def run_pdb_interpretation (file_name) :
    processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      params=None,
      file_name=file_name,
      strict_conflict_handling=True,
      substitute_non_crystallographic_unit_cell_if_necessary=False,
      max_atoms=None,
      log=null_out())
    return processed_pdb_file
  edits_phil = iotbx.phil.parse(
    mmtbx.monomer_library.pdb_interpretation.geometry_restraints_edits_str)
  edits = libtbx.phil.parse("""
scale_restraints {
  atom_selection = chain A and resseq 1
  scale = 2.0
}
scale_restraints {
  atom_selection = chain A and resseq 4
  scale = 0.5
  apply_to = *bond *angle *chirality
}
""")
  params_edits = edits_phil.fetch(source=edits).extract()
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/enk.pdb",
    test=os.path.isfile)
  if (not os.path.isfile(pdb_file)) :
    return
  processed_pdb_file_1 = run_pdb_interpretation(pdb_file)
  grm_1 = processed_pdb_file_1.geometry_restraints_manager()
  angle_proxies_1 = grm_1.angle_proxies
  dihedral_proxies_1 = grm_1.dihedral_proxies
  processed_pdb_file_2 = run_pdb_interpretation(pdb_file)
  grm_2 = processed_pdb_file_2.geometry_restraints_manager(
    params_edits=params_edits)
  angle_proxies_2 = grm_2.angle_proxies
  dihedral_proxies_2 = grm_2.dihedral_proxies
  assert (len(angle_proxies_2) == len(angle_proxies_1))
  assert (len(dihedral_proxies_2) == len(dihedral_proxies_1))
  for j, proxy_1 in enumerate(angle_proxies_1) :
    proxy_2 = angle_proxies_2[j]
    if (0 <= j <= 17) : # this corresponds to resseq 1
      assert (proxy_2.weight == 2*proxy_1.weight)
    elif (25 <= j <= 43) : # resseq 4 (and adjoining atoms)
      assert (proxy_2.weight == 0.5*proxy_1.weight)
    else :
      assert (proxy_2.weight == proxy_1.weight)
  for j, proxy_1 in enumerate(dihedral_proxies_1) :
    proxy_2 = dihedral_proxies_2[j]
    if (0 <= j <= 3) :
      assert (proxy_2.weight == 2*proxy_1.weight)
    else :
      assert (proxy_2.weight == proxy_1.weight)

def run(args):
  assert len(args) == 0
  exercise_scale_restraints()
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  exercise_pdb_string(mon_lib_srv, ener_lib)
  exercise_cns_rna(mon_lib_srv, ener_lib)
  exercise_rna_3p_2p(mon_lib_srv, ener_lib)
  exercise_rna_dna_hybrid(mon_lib_srv, ener_lib)
  exercise_hydrogen_deuterium_aliases()
  exercise_corrupt_cif_link()
  exercise_dna_cns_cy5_th6()
  exercise_sym_excl_indices(mon_lib_srv, ener_lib)
  exercise_auto_alias_h_h1()
  exercise_d_aa_resnames()
  exercise_d_amino_acid_chain_perfect_in_box()
  exercise_rna_v3(mon_lib_srv, ener_lib)
  print format_cpu_times()

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
