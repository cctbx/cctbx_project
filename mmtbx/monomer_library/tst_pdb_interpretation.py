from __future__ import absolute_import, division, print_function
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from libtbx.utils import Sorry, search_for, format_cpu_times, null_out
from libtbx.test_utils import Exception_expected, block_show_diff, approx_equal
import libtbx.load_env
import iotbx.phil
from libtbx import Auto
from six.moves import cStringIO as StringIO
import os
import sys
from cctbx.array_family import flex
from six.moves import zip

params = monomer_library.pdb_interpretation.master_params.extract()
params.flip_symmetric_amino_acids = False

def exercise_handle_case_insensitive(mon_lib_srv, ener_lib):
  def check(a, r, e):
    raw_records = ("""\
HETATM    1 %s    %s F   1      -7.869  17.488  18.637  1.00 30.00          %s
""" % (a,r,e)).splitlines()
    log = StringIO()
    processed_pdb_file = monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=None,
      params=params,
      raw_records=raw_records,
      log=log)
    log_lines = log.getvalue().splitlines()
    lines=search_for(
      pattern="""\
          Unusual residues: {' %s': 1}""" % r,
      mode="==",
      lines=log_lines)
    assert len(lines) == 1
    assert len(log_lines) == 14
  vars = ["ZN", "Zn"]
  for a in vars:
    for r in vars:
      for e in vars:
        check(a, r, e)

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
  log = StringIO()
  cif_object = iotbx.cif.reader(input_string=cif_records).model()
  mon_lib_srv.process_cif_object(cif_object=cif_object)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    params=params,
    force_symmetry=True,
    log=log)
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
    params=params,
    raw_records=raw_records,
    force_symmetry=True)
  xray_structure = processed_pdb_file.xray_structure()
  assert xray_structure.scattering_type_registry().type_count_dict() \
      == {"S": 1, "N": 1, "C": 5, "O": 1, "D": 11}
  raw_records = [line[:66] for line in raw_records]
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=params,
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
    params=params,
    raw_records=raw_records,
    log=log)
  looking_for = "Ad-hoc single atom residues: "
  for line in log.getvalue().splitlines():
    line = line.strip()
    if (not line.startswith(looking_for)): continue
    counts = eval(line[len(looking_for):])
    assert counts == {"U  ": 1}
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
    print('Skipping exercise_rna("%s"): input file not available' % file_name)
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
          Link IDs: {'rna2p': 1, 'rna3p': 1}
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
      print("Skipping exercise_hydrogen_deuterium_aliases():", \
        "input file not available:", file_name)
      return
    file_paths.append(file_path)
  log = StringIO()
  monomer_library.pdb_interpretation.run(args=file_paths, params=params, log=log)
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
      print("Skipping exercise_corrupt_cif_link():", \
        "input file not available:", file_name)
      return
    file_paths.append(file_path)
  log = StringIO()
  try:
    monomer_library.pdb_interpretation.run(args=file_paths, params=params, log=log)
  except Sorry as e:
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
      print("Skipping exercise_dna_cns_cy5_th6():", \
        "input file not available:", file_name)
      return
    file_paths.append(file_path)
  log = StringIO()
  monomer_library.pdb_interpretation.run(args=file_paths, params=params, log=log)
  assert not block_show_diff(log.getvalue(), """\
        Number of residues, atoms: 12, 244
          Unusual residues: %s
          Classifications: %s
          Modifications used: {'5*END': 1}
          Link IDs: %s
""" % (
  str({'CY5': 1, 'TH6': 1}),
  str({'DNA': 10, 'DNA_mixed' : 2}),
  str({'rna3p': 11})))

def exercise_sym_excl_indices(mon_lib_srv, ener_lib):
  file_name = "so4_on_two_fold.pdb"
  file_path = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/"+file_name,
    test=os.path.isfile)
  if (file_path is None):
    print("Skipping exercise_sym_excl_indices():", \
      "input file not available:", file_name)
    return
  log = StringIO()
  pdb_interpretation_params = monomer_library.pdb_interpretation.master_params.extract()
  pdb_interpretation_params.sort_atoms=False
  pdb_interpretation_params.flip_symmetric_amino_acids=False
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=file_path,
    params=pdb_interpretation_params,
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
      print("Skipping exercise_auto_alias_h_h1():", \
        "input file not available:", file_name)
      return
    file_paths.append(file_path)
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.run(
        params=params,
    args=file_paths, log=log)
  assert log.getvalue().find("Modifications used: {'NH3': 1}") >= 0
  assert processed_pdb_file.all_chain_proxies.fatal_problems_message() is None

def exercise_d_aa_resnames():
  file_path = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/dpn.pdb",
    test=os.path.isfile)
  if (file_path is None):
    print("Skipping exercise_d_aa_resnames():", \
      "input file not available: dpn.pdb")
    return
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.run(
        params=params,
    args=[file_path], log=log)
  assert log.getvalue().find("'PEPT-D': 1") >= 0
  assert log.getvalue().find("'NH1NOTPRO': 1") >= 0
  assert processed_pdb_file.all_chain_proxies.fatal_problems_message() is None

def exercise_d_amino_acid_chain_perfect_in_box_peptide_plane():
  file_path = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/d_amino_acid_chain_perfect_in_box.pdb",
    test=os.path.isfile)
  if (file_path is None):
    print("Skipping exercise_d_aa_resnames():", \
      "input file not available: d_amino_acid_chain_perfect_in_box.pdb")
    return
  log = StringIO()
  pdb_interpretation_params = monomer_library.pdb_interpretation.master_params.extract()
  pdb_interpretation_params.sort_atoms=False
  pdb_interpretation_params.flip_symmetric_amino_acids=False
  pdb_interpretation_params.peptide_link.apply_peptide_plane=True
  processed_pdb_file = monomer_library.pdb_interpretation.run(
    args=[file_path], params=pdb_interpretation_params, log=log)
  grm = processed_pdb_file.geometry_restraints_manager()
  lv = log.getvalue()
  assert lv.find("Classifications: {'peptide': 16}") >= 0
  assert lv.find("Link IDs: {'PCIS': 1, 'TRANS': 14, 'peptide plane': 15}")>=0

def exercise_d_amino_acid_chain_perfect_in_box():
  file_path = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/d_amino_acid_chain_perfect_in_box.pdb",
    test=os.path.isfile)
  if (file_path is None):
    print("Skipping exercise_d_aa_resnames():", \
      "input file not available: d_amino_acid_chain_perfect_in_box.pdb")
    return
  log = StringIO()
  pdb_interpretation_params = monomer_library.pdb_interpretation.master_params.extract()
  pdb_interpretation_params.sort_atoms=False
  pdb_interpretation_params.flip_symmetric_amino_acids=False
  processed_pdb_file = monomer_library.pdb_interpretation.run(
    args=[file_path], params=pdb_interpretation_params, log=log)
  grm = processed_pdb_file.geometry_restraints_manager()
  lv = log.getvalue()
  print(lv)
  assert lv.find("Classifications: {'peptide': 16}") >= 0
  assert lv.find("'PEPT-D': 1") >= 0
  assert lv.find("'TRANS': 14") >= 0
  assert lv.find("'PCIS': 1") >= 0
  assert lv.find("Link IDs: {'PCIS': 1, 'TRANS': 14}")>=0
  assert lv.find("""\
Simple disulfide: pdb=" SG  DCY A   4 " - pdb=" SG  DCY A  19 " distance=2.03
""") >= 0
  assert not block_show_diff(lv, """\
  Bond restraints: 121
  Sorted by residual:
  bond pdb=" CZ  DAR A  14 "
       pdb=" NH2 DAR A  14 "
    ideal  model  delta    sigma   weight residual
    1.330  1.326  0.004 1.30e-02 5.92e+03 1.15e-01
""")
  assert not block_show_diff(lv, """\
  Bond angle restraints: 161
  Sorted by residual:
  angle pdb=" NE  DAR A  14 "
        pdb=" CZ  DAR A  14 "
        pdb=" NH1 DAR A  14 "
      ideal   model   delta    sigma   weight residual
     121.50  120.07    1.43 1.00e+00 1.00e+00 2.04e+00
""")
  assert not block_show_diff(lv, """\
  Dihedral angle restraints: 50
    sinusoidal: 34
      harmonic: 16
  Sorted by residual:
  dihedral pdb=" CA  DCY A   4 "
           pdb=" CB  DCY A   4 "
           pdb=" SG  DCY A   4 "
           pdb=" SG  DCY A  19 "
      ideal   model   delta sinusoidal    sigma   weight residual
      79.00   18.71   60.29     1      2.00e+01 2.50e-03 1.21e+01
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

def exercise_asp_glu_acid():
  for resname in ["ASP", "GLU"]:
    file_name = resname.lower() + "_acid.pdb"
    file_path = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/pdb/"+file_name,
      test=os.path.isfile)
    if (file_path is None):
      print("Skipping exercise_asp_glu_acid():", \
        "input file not available:", file_name)
    else:
      log = StringIO()
      processed_pdb_file = monomer_library.pdb_interpretation.run(
        params=params,
        args=[file_path], log=log)
      if (resname == "ASP"):
        pat = "Modifications used: {'ACID-ASP': 1}"
      else:
        pat = "Modifications used: {'ACID-GLU': 1, 'COOH': 1, 'NH1NOTPRO': 1}"
      assert log.getvalue().find(pat) >= 0
      assert processed_pdb_file.all_chain_proxies \
        .fatal_problems_message() is None

# tests automatic aliasing of OP1/O1P etc. when the residue is not detected as
# RNA/DNA but the monomer entry is
def exercise_rna_dna_synonyms():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/PGP.pdb",
    test=os.path.isfile)
  cif_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/cif_files/PGP.cif",
    test=os.path.isfile)
  if (None in [pdb_file, cif_file]):
    print("Skipping exercise_rna_dna_synonyms() - input file(s) unavailable.")
  else :
    log = StringIO()
    processed_pdb_file = monomer_library.pdb_interpretation.run(
        params=params,
      args=[pdb_file, cif_file], log=log)
    msg = processed_pdb_file.all_chain_proxies.fatal_problems_message(
      ignore_unknown_scattering_types=False,
      ignore_unknown_nonbonded_energy_types=False)
    assert (msg is None)

def exercise_geostd_cif_links(mon_lib_srv, ener_lib):
  raw_records = """\
ATOM   6707  C1  MAN D9032     -82.149  64.388 127.264  1.00 20.00           C
ATOM   6708  C2  MAN D9032     -82.661  63.644 128.508  1.00 20.00           C
ATOM   6709  O2  MAN D9032     -81.755  62.603 128.885  1.00 20.00           O
ATOM   6710  C3  MAN D9032     -83.112  64.561 129.660  1.00 20.00           C
ATOM   6711  O3  MAN D9032     -83.036  63.857 130.877  1.00 20.00           O
ATOM   6712  C4  MAN D9032     -82.378  65.900 129.641  1.00 20.00           C
ATOM   6713  O4  MAN D9032     -82.870  66.786 130.611  1.00 20.00           O
ATOM   6714  C5  MAN D9032     -82.589  66.508 128.255  1.00 20.00           C
ATOM   6715  O5  MAN D9032     -81.778  65.758 127.405  1.00 20.00           O
ATOM   6716  C6  MAN D9032     -82.234  68.002 128.113  1.00 20.00           C
ATOM   6717  O6  MAN D9032     -83.296  68.629 127.421  1.00 20.00           O
ATOM   6718  O1  MAN D8032     -82.173  69.386 124.623  1.00 20.00           O
ATOM   6719  C1  MAN D8032     -82.135  70.338 125.650  1.00 20.00           C
ATOM   6720  C2  MAN D8032     -83.046  69.904 126.785  1.00 20.00           C
ATOM   6721  C3  MAN D8032     -83.124  71.033 127.800  1.00 20.00           C
ATOM   6722  O3  MAN D8032     -83.760  70.573 128.974  1.00 20.00           O
ATOM   6723  C4  MAN D8032     -81.746  71.580 128.181  1.00 20.00           C
ATOM   6724  O4  MAN D8032     -81.878  72.934 128.571  1.00 20.00           O
ATOM   6725  C5  MAN D8032     -80.638  71.462 127.125  1.00 20.00           C
ATOM   6726  O5  MAN D8032     -80.827  70.430 126.159  1.00 20.00           O
ATOM   6727  C6  MAN D8032     -79.314  71.191 127.833  1.00 20.00           C
ATOM   6728  O6  MAN D8032     -79.563  70.916 129.196  1.00 20.00           O
""".splitlines()
  link_def = """\
apply_cif_link {
  data_link = ALPHA2-6
  residue_selection_1 = chain D and resname MAN and resseq 9032
  residue_selection_2 = chain D and resname MAN and resseq 8032
}
"""
  links_phil = iotbx.phil.parse(
    mmtbx.monomer_library.pdb_interpretation.master_params_str,
    process_includes=True,
    )
  links = libtbx.phil.parse(link_def)
  params_links = links_phil.fetch(source=links)
  #params_links.show()
  params_links = params_links.extract()

  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    #force_symmetry=True,
    params=params_links,
    log=log,
    )
  #print log.getvalue()
  pat = "data_link: ALPHA2-6"
  assert log.getvalue().find(pat) >= 0

def exercise_flattened_cif_loop():
  # Previously our code would not interpret properly files where CIF loops with
  # only one row were "flattened" - i.e. listed as tag-value pairs rather than
  # in a loop construct
  pdb_string = """\
CRYST1   15.000   15.000   15.000  80.00  70.00 100.00 P 1
HETATM   13  O   HOH     1      13.866  16.009  12.098  1.00  3.00           O
HETATM   14  H1  HOH     1      13.327  16.140  12.905  1.00  3.00           H
HETATM   15  H2  HOH     1      14.117  16.901  11.777  1.00  3.00           H
"""
  restraints_cif_string = """\
#
data_comp_list
_chem_comp.id                  HOH
_chem_comp.three_letter_code   .
_chem_comp.name                'water     '
_chem_comp.group               solvent
_chem_comp.number_atoms_all    3
_chem_comp.number_atoms_nh     1
_chem_comp.desc_level          .
#
data_comp_HOH
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
 HOH           O      O    OH2      -0.408
 HOH           H1     H    HOH2      0.204
 HOH           H2     H    HOH2      0.204
loop_
_chem_comp_tree.comp_id
_chem_comp_tree.atom_id
_chem_comp_tree.atom_back
_chem_comp_tree.atom_forward
_chem_comp_tree.connect_type
 HOH      O      n/a    .      END
 HOH      H1     O      .      .
 HOH      H2     O      .      .
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 HOH      O      H1        coval       0.840    0.020    0.950
 HOH      O      H2        coval       0.840    0.020    0.950
_chem_comp_angle.comp_id           HOH
_chem_comp_angle.atom_id_1         H1
_chem_comp_angle.atom_id_2         O
_chem_comp_angle.atom_id_3         H2
_chem_comp_angle.value_angle       106.800
_chem_comp_angle.value_angle_esd   3.000
"""
  from libtbx.test_utils import open_tmp_file
  from libtbx import easy_run
  pdb_file = open_tmp_file(suffix="HOH.pdb")
  pdb_file.write(pdb_string)
  pdb_file.close()
  restraints_file = open_tmp_file(suffix="HOH.cif")
  restraints_file.write(restraints_cif_string)
  restraints_file.close()
  cmd = "phenix.pdb_interpretation \"%s\" \"%s\" write_geo_files=True" %(
    pdb_file.name, restraints_file.name)
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  geo_file = open(pdb_file.name+".geo", "r")
  geo_file_str = geo_file.read()
  geo_file.close()
  assert "Bond angle | covalent geometry | restraints: 1" in geo_file_str

def exercise_do_not_link(mon_lib_srv, ener_lib):
  """
  Do not attempt to link two residues below.
  """
  raw_records = """\
CRYST1   48.273   28.843   18.740  90.00  90.00  90.00 P 1
ATOM  17464  C1  BOG R  16      30.175 -15.554   9.755  1.00 35.95      R    C
ATOM  17465  O1  BOG R  16      31.099 -15.296   8.697  1.00 42.53      R    O
ATOM  17466  C2  BOG R  16      29.782 -17.026   9.758  1.00 39.47      R    C
ATOM  17467  O2  BOG R  16      29.087 -17.331   8.552  1.00 40.97      R    O
ATOM  17468  C3  BOG R  16      28.901 -17.357  10.950  1.00 40.73      R    C
ATOM  17469  O3  BOG R  16      28.635 -18.764  10.974  1.00 46.27      R    O
ATOM  17470  C4  BOG R  16      29.617 -16.951  12.225  1.00 40.76      R    C
ATOM  17471  O4  BOG R  16      28.755 -17.165  13.351  1.00 42.98      R    O
ATOM  17472  C5  BOG R  16      29.971 -15.472  12.137  1.00 39.41      R    C
ATOM  17473  O5  BOG R  16      30.787 -15.208  10.993  1.00 37.94      R    O
ATOM  17474  C6  BOG R  16      30.688 -14.989  13.394  1.00 40.91      R    C
ATOM  17475  O6  BOG R  16      31.999 -15.553  13.445  1.00 37.31      R    O
ATOM  17476  C1' BOG R  16      31.264 -13.898   8.462  1.00 38.79      R    C
ATOM  17477  C2' BOG R  16      32.250 -13.694   7.317  1.00 39.57      R    C
ATOM  17478  C3' BOG R  16      32.609 -12.219   7.133  1.00 44.82      R    C
ATOM  17479  C4' BOG R  16      31.515 -11.448   6.400  1.00 50.44      R    C
ATOM  17480  C5' BOG R  16      31.977 -10.037   6.047  1.00 48.33      R    C
ATOM      1  C1  BOG R  18      66.170  -7.947   8.034  1.00 55.33      R    C
ATOM      2  O1  BOG R  18      65.900  -6.739   8.753  1.00 55.00      R    O
ATOM      3  C2  BOG R  18      65.682  -9.113   8.883  1.00 55.62      R    C
ATOM      4  O2  BOG R  18      66.397  -9.125  10.123  1.00 61.39      R    O
ATOM      5  C3  BOG R  18      65.843 -10.438   8.147  1.00 60.20      R    C
ATOM      6  O3  BOG R  18      65.302 -11.514   8.924  1.00 64.46      R    O
ATOM      7  C4  BOG R  18      65.106 -10.358   6.821  1.00 55.23      R    C
ATOM      8  O4  BOG R  18      65.294 -11.575   6.088  1.00 58.17      R    O
ATOM      9  C5  BOG R  18      65.635  -9.174   6.022  1.00 55.29      R    C
ATOM     10  O5  BOG R  18      65.506  -7.953   6.761  1.00 55.00      R    O
ATOM     11  C6  BOG R  18      64.879  -9.032   4.705  1.00 50.53      R    C
ATOM     12  O6  BOG R  18      63.593  -8.452   4.947  1.00 41.70      R    O
ATOM     13  C1' BOG R  18      66.554  -5.574   8.243  1.00 52.47      R    C
ATOM     14  C2' BOG R  18      66.138  -4.389   9.113  1.00 55.51      R    C
ATOM     15  C3' BOG R  18      66.908  -3.115   8.772  1.00 56.45      R    C
ATOM     16  C4' BOG R  18      66.061  -2.210   7.826  1.00 56.42      R    C
ATOM     17  C5' BOG R  18      66.250  -0.997   7.016  1.00 58.47      R    C
ATOM     18  C6' BOG R  18      66.002   0.079   5.931  1.00 59.92      R    C
""".splitlines()
  cif_records = """
# ------------------------------------------------
#
# ---   LIST OF MONOMERS ---
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
BOG      BOG 'B-OCTYLGLUCOSIDE                    ' pyranose           48  20 .
# ------------------------------------------------------
# ------------------------------------------------------
#
# --- DESCRIPTION OF MONOMERS ---
#
data_comp_BOG
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
 BOG           C1     C    CH1       0.000      0.000    0.000    0.000
 BOG           H1     H    H         0.000      0.135   -1.086    0.103
 BOG           O1     O    O2        0.000     -1.360    0.279   -0.341
 BOG           "C1'"  C    CH2       0.000     -2.165   -0.233    0.721
 BOG           "H1'1" H    H         0.000     -1.890    0.259    1.656
 BOG           "H1'2" H    H         0.000     -2.001   -1.309    0.815
 BOG           "C2'"  C    CH2       0.000     -3.640    0.035    0.417
 BOG           "H2'1" H    H         0.000     -3.913   -0.458   -0.519
 BOG           "H2'2" H    H         0.000     -3.801    1.110    0.321
 BOG           "C3'"  C    CH2       0.000     -4.504   -0.514    1.554
 BOG           "H3'1" H    H         0.000     -4.228   -0.021    2.489
 BOG           "H3'2" H    H         0.000     -4.340   -1.589    1.649
 BOG           "C4'"  C    CH2       0.000     -5.979   -0.246    1.250
 BOG           "H4'1" H    H         0.000     -6.252   -0.737    0.314
 BOG           "H4'2" H    H         0.000     -6.141    0.830    1.155
 BOG           "C5'"  C    CH2       0.000     -6.843   -0.795    2.386
 BOG           "H5'1" H    H         0.000     -6.567   -0.303    3.322
 BOG           "H5'2" H    H         0.000     -6.679   -1.870    2.480
 BOG           "C6'"  C    CH2       0.000     -8.318   -0.527    2.083
 BOG           "H6'1" H    H         0.000     -8.591   -1.019    1.147
 BOG           "H6'2" H    H         0.000     -8.479    0.549    1.988
 BOG           "C7'"  C    CH2       0.000     -9.181   -1.075    3.219
 BOG           "H7'1" H    H         0.000     -8.906   -0.583    4.154
 BOG           "H7'2" H    H         0.000     -9.017   -2.151    3.314
 BOG           "C8'"  C    CH3       0.000    -10.656   -0.807    2.915
 BOG           "H8'3" H    H         0.000    -10.818    0.237    2.823
 BOG           "H8'2" H    H         0.000    -10.926   -1.285    2.008
 BOG           "H8'1" H    H         0.000    -11.257   -1.186    3.702
 BOG           O5     O    O2        0.000      0.320    0.635    1.235
 BOG           C5     C    CH1       0.000      1.615    0.188    1.628
 BOG           H5     H    H         0.000      1.639   -0.911    1.621
 BOG           C4     C    CH1       0.000      2.668    0.725    0.656
 BOG           H4     H    H         0.000      2.630    1.824    0.644
 BOG           O4     O    OH1       0.000      3.965    0.296    1.072
 BOG           HO4    H    H         0.000      4.629    0.637    0.457
 BOG           C3     C    CH1       0.000      2.372    0.187   -0.747
 BOG           H3     H    H         0.000      2.512   -0.903   -0.761
 BOG           O3     O    OH1       0.000      3.255    0.796   -1.691
 BOG           HO3    H    H         0.000      3.066    0.455   -2.575
 BOG           C2     C    CH1       0.000      0.919    0.523   -1.105
 BOG           H2     H    H         0.000      0.805    1.612   -1.194
 BOG           O2     O    OH1       0.000      0.577   -0.095   -2.346
 BOG           HO2    H    H         0.000     -0.339    0.119   -2.570
 BOG           C6     C    CH2       0.000      1.922    0.693    3.039
 BOG           H61    H    H         0.000      1.898    1.785    3.047
 BOG           H62    H    H         0.000      2.914    0.349    3.340
 BOG           O6     O    OH1       0.000      0.944    0.187    3.950
 BOG           HO6    H    H         0.000      1.174    0.528    4.825
loop_
_chem_comp_tree.comp_id
_chem_comp_tree.atom_id
_chem_comp_tree.atom_back
_chem_comp_tree.atom_forward
_chem_comp_tree.connect_type
 BOG      C1     n/a    O5     START
 BOG      H1     C1     .      .
 BOG      O1     C1     "C1'"  .
 BOG      "C1'"  O1     "C2'"  .
 BOG      "H1'1" "C1'"  .      .
 BOG      "H1'2" "C1'"  .      .
 BOG      "C2'"  "C1'"  "C3'"  .
 BOG      "H2'1" "C2'"  .      .
 BOG      "H2'2" "C2'"  .      .
 BOG      "C3'"  "C2'"  "C4'"  .
 BOG      "H3'1" "C3'"  .      .
 BOG      "H3'2" "C3'"  .      .
 BOG      "C4'"  "C3'"  "C5'"  .
 BOG      "H4'1" "C4'"  .      .
 BOG      "H4'2" "C4'"  .      .
 BOG      "C5'"  "C4'"  "C6'"  .
 BOG      "H5'1" "C5'"  .      .
 BOG      "H5'2" "C5'"  .      .
 BOG      "C6'"  "C5'"  "C7'"  .
 BOG      "H6'1" "C6'"  .      .
 BOG      "H6'2" "C6'"  .      .
 BOG      "C7'"  "C6'"  "C8'"  .
 BOG      "H7'1" "C7'"  .      .
 BOG      "H7'2" "C7'"  .      .
 BOG      "C8'"  "C7'"  "H8'1" .
 BOG      "H8'3" "C8'"  .      .
 BOG      "H8'2" "C8'"  .      .
 BOG      "H8'1" "C8'"  .      .
 BOG      O5     C1     .      END
 BOG      C5     O5     C6     .
 BOG      H5     C5     .      .
 BOG      C4     C5     C3     .
 BOG      H4     C4     .      .
 BOG      O4     C4     HO4    .
 BOG      HO4    O4     .      .
 BOG      C3     C4     C2     .
 BOG      H3     C3     .      .
 BOG      O3     C3     HO3    .
 BOG      HO3    O3     .      .
 BOG      C2     C3     O2     .
 BOG      H2     C2     .      .
 BOG      O2     C2     HO2    .
 BOG      HO2    O2     .      .
 BOG      C6     C5     O6     .
 BOG      H61    C6     .      .
 BOG      H62    C6     .      .
 BOG      O6     C6     .      .
 BOG      HO6    O6     .      .
 BOG      C1     C2     .    ADD
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 BOG      O1     C1        single      1.426    0.020
 BOG      C1     C2        single      1.524    0.020
 BOG      O5     C1        single      1.426    0.020
 BOG      H1     C1        single      1.099    0.020
 BOG      "C1'"  O1        single      1.426    0.020
 BOG      O2     C2        single      1.432    0.020
 BOG      C2     C3        single      1.524    0.020
 BOG      H2     C2        single      1.099    0.020
 BOG      HO2    O2        single      0.967    0.020
 BOG      O3     C3        single      1.432    0.020
 BOG      C3     C4        single      1.524    0.020
 BOG      H3     C3        single      1.099    0.020
 BOG      HO3    O3        single      0.967    0.020
 BOG      O4     C4        single      1.432    0.020
 BOG      C4     C5        single      1.524    0.020
 BOG      H4     C4        single      1.099    0.020
 BOG      HO4    O4        single      0.967    0.020
 BOG      C5     O5        single      1.426    0.020
 BOG      C6     C5        single      1.524    0.020
 BOG      H5     C5        single      1.099    0.020
 BOG      O6     C6        single      1.432    0.020
 BOG      H61    C6        single      1.092    0.020
 BOG      H62    C6        single      1.092    0.020
 BOG      HO6    O6        single      0.967    0.020
 BOG      "C2'"  "C1'"     single      1.524    0.020
 BOG      "H1'1" "C1'"     single      1.092    0.020
 BOG      "H1'2" "C1'"     single      1.092    0.020
 BOG      "C3'"  "C2'"     single      1.524    0.020
 BOG      "H2'1" "C2'"     single      1.092    0.020
 BOG      "H2'2" "C2'"     single      1.092    0.020
 BOG      "C4'"  "C3'"     single      1.524    0.020
 BOG      "H3'1" "C3'"     single      1.092    0.020
 BOG      "H3'2" "C3'"     single      1.092    0.020
 BOG      "C5'"  "C4'"     single      1.524    0.020
 BOG      "H4'1" "C4'"     single      1.092    0.020
 BOG      "H4'2" "C4'"     single      1.092    0.020
 BOG      "C6'"  "C5'"     single      1.524    0.020
 BOG      "H5'1" "C5'"     single      1.092    0.020
 BOG      "H5'2" "C5'"     single      1.092    0.020
 BOG      "C7'"  "C6'"     single      1.524    0.020
 BOG      "H6'1" "C6'"     single      1.092    0.020
 BOG      "H6'2" "C6'"     single      1.092    0.020
 BOG      "C8'"  "C7'"     single      1.513    0.020
 BOG      "H7'1" "C7'"     single      1.092    0.020
 BOG      "H7'2" "C7'"     single      1.092    0.020
 BOG      "H8'1" "C8'"     single      1.059    0.020
 BOG      "H8'2" "C8'"     single      1.059    0.020
 BOG      "H8'3" "C8'"     single      1.059    0.020
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 BOG      H1     C1     O1      109.470    3.000
 BOG      H1     C1     O5      109.470    3.000
 BOG      O1     C1     O5      109.470    3.000
 BOG      H1     C1     C2      108.340    3.000
 BOG      O1     C1     C2      109.470    3.000
 BOG      O5     C1     C2      109.470    3.000
 BOG      C1     O1     "C1'"   111.800    3.000
 BOG      O1     "C1'"  "H1'1"  109.470    3.000
 BOG      O1     "C1'"  "H1'2"  109.470    3.000
 BOG      O1     "C1'"  "C2'"   109.470    3.000
 BOG      "H1'1" "C1'"  "H1'2"  107.900    3.000
 BOG      "H1'1" "C1'"  "C2'"   109.470    3.000
 BOG      "H1'2" "C1'"  "C2'"   109.470    3.000
 BOG      "C1'"  "C2'"  "H2'1"  109.470    3.000
 BOG      "C1'"  "C2'"  "H2'2"  109.470    3.000
 BOG      "C1'"  "C2'"  "C3'"   111.000    3.000
 BOG      "H2'1" "C2'"  "H2'2"  107.900    3.000
 BOG      "H2'1" "C2'"  "C3'"   109.470    3.000
 BOG      "H2'2" "C2'"  "C3'"   109.470    3.000
 BOG      "C2'"  "C3'"  "H3'1"  109.470    3.000
 BOG      "C2'"  "C3'"  "H3'2"  109.470    3.000
 BOG      "C2'"  "C3'"  "C4'"   111.000    3.000
 BOG      "H3'1" "C3'"  "H3'2"  107.900    3.000
 BOG      "H3'1" "C3'"  "C4'"   109.470    3.000
 BOG      "H3'2" "C3'"  "C4'"   109.470    3.000
 BOG      "C3'"  "C4'"  "H4'1"  109.470    3.000
 BOG      "C3'"  "C4'"  "H4'2"  109.470    3.000
 BOG      "C3'"  "C4'"  "C5'"   111.000    3.000
 BOG      "H4'1" "C4'"  "H4'2"  107.900    3.000
 BOG      "H4'1" "C4'"  "C5'"   109.470    3.000
 BOG      "H4'2" "C4'"  "C5'"   109.470    3.000
 BOG      "C4'"  "C5'"  "H5'1"  109.470    3.000
 BOG      "C4'"  "C5'"  "H5'2"  109.470    3.000
 BOG      "C4'"  "C5'"  "C6'"   111.000    3.000
 BOG      "H5'1" "C5'"  "H5'2"  107.900    3.000
 BOG      "H5'1" "C5'"  "C6'"   109.470    3.000
 BOG      "H5'2" "C5'"  "C6'"   109.470    3.000
 BOG      "C5'"  "C6'"  "H6'1"  109.470    3.000
 BOG      "C5'"  "C6'"  "H6'2"  109.470    3.000
 BOG      "C5'"  "C6'"  "C7'"   111.000    3.000
 BOG      "H6'1" "C6'"  "H6'2"  107.900    3.000
 BOG      "H6'1" "C6'"  "C7'"   109.470    3.000
 BOG      "H6'2" "C6'"  "C7'"   109.470    3.000
 BOG      "C6'"  "C7'"  "H7'1"  109.470    3.000
 BOG      "C6'"  "C7'"  "H7'2"  109.470    3.000
 BOG      "C6'"  "C7'"  "C8'"   111.000    3.000
 BOG      "H7'1" "C7'"  "H7'2"  107.900    3.000
 BOG      "H7'1" "C7'"  "C8'"   109.470    3.000
 BOG      "H7'2" "C7'"  "C8'"   109.470    3.000
 BOG      "C7'"  "C8'"  "H8'3"  109.470    3.000
 BOG      "C7'"  "C8'"  "H8'2"  109.470    3.000
 BOG      "C7'"  "C8'"  "H8'1"  109.470    3.000
 BOG      "H8'3" "C8'"  "H8'2"  109.470    3.000
 BOG      "H8'3" "C8'"  "H8'1"  109.470    3.000
 BOG      "H8'2" "C8'"  "H8'1"  109.470    3.000
 BOG      C1     O5     C5      111.800    3.000
 BOG      O5     C5     H5      109.470    3.000
 BOG      O5     C5     C4      109.470    3.000
 BOG      O5     C5     C6      109.470    3.000
 BOG      H5     C5     C4      108.340    3.000
 BOG      H5     C5     C6      108.340    3.000
 BOG      C4     C5     C6      111.000    3.000
 BOG      C5     C4     H4      108.340    3.000
 BOG      C5     C4     O4      109.470    3.000
 BOG      C5     C4     C3      111.000    3.000
 BOG      H4     C4     O4      109.470    3.000
 BOG      H4     C4     C3      108.340    3.000
 BOG      O4     C4     C3      109.470    3.000
 BOG      C4     O4     HO4     109.470    3.000
 BOG      C4     C3     H3      108.340    3.000
 BOG      C4     C3     O3      109.470    3.000
 BOG      C4     C3     C2      111.000    3.000
 BOG      H3     C3     O3      109.470    3.000
 BOG      H3     C3     C2      108.340    3.000
 BOG      O3     C3     C2      109.470    3.000
 BOG      C3     O3     HO3     109.470    3.000
 BOG      C3     C2     H2      108.340    3.000
 BOG      C3     C2     O2      109.470    3.000
 BOG      C3     C2     C1      111.000    3.000
 BOG      H2     C2     O2      109.470    3.000
 BOG      H2     C2     C1      108.340    3.000
 BOG      O2     C2     C1      109.470    3.000
 BOG      C2     O2     HO2     109.470    3.000
 BOG      C5     C6     H61     109.470    3.000
 BOG      C5     C6     H62     109.470    3.000
 BOG      C5     C6     O6      109.470    3.000
 BOG      H61    C6     H62     107.900    3.000
 BOG      H61    C6     O6      109.470    3.000
 BOG      H62    C6     O6      109.470    3.000
 BOG      C6     O6     HO6     109.470    3.000
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
 BOG      var_1    O5     C1     O1     "C1'"    -59.799   20.000   1
 BOG      var_2    C1     O1     "C1'"  "C2'"    179.999   20.000   1
 BOG      var_3    O1     "C1'"  "C2'"  "C3'"   -179.965   20.000   3
 BOG      var_4    "C1'"  "C2'"  "C3'"  "C4'"    180.000   20.000   3
 BOG      var_5    "C2'"  "C3'"  "C4'"  "C5'"    179.946   20.000   3
 BOG      var_6    "C3'"  "C4'"  "C5'"  "C6'"    179.980   20.000   3
 BOG      var_7    "C4'"  "C5'"  "C6'"  "C7'"   -179.963   20.000   3
 BOG      var_8    "C5'"  "C6'"  "C7'"  "C8'"   -179.980   20.000   3
 BOG      var_9    "C6'"  "C7'"  "C8'"  "H8'1"  -179.975   20.000   3
 BOG      var_10   C1     O5     C5     C6       180.000   20.000   1
 BOG      var_11   O5     C5     C4     C3       -60.000   20.000   3
 BOG      var_12   C5     C4     O4     HO4     -179.949   20.000   1
 BOG      var_13   C5     C4     C3     C2        60.000   20.000   3
 BOG      var_14   C4     C3     O3     HO3     -179.973   20.000   1
 BOG      var_15   C4     C3     C2     O2       180.000   20.000   3
 BOG      var_16   C3     C2     C1     O5        60.000   20.000   3
 BOG      var_17   C3     C2     O2     HO2     -179.932   20.000   1
 BOG      var_18   O5     C5     C6     O6        59.890   20.000   3
 BOG      var_1    C5     O5     C1     C2       -55.000   20.000   1 #
loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
 BOG      chir_01  C1     O1     C2     O5        negativ
 BOG      chir_02  C2     C1     O2     C3        positiv
 BOG      chir_03  C3     C2     O3     C4        negativ
 BOG      chir_04  C4     C3     O4     C5        positiv
 BOG      chir_05  C5     C4     O5     C6        positiv
# ------------------------------------------------------
"""
  import iotbx.cif
  log = StringIO()
  cif_object = iotbx.cif.reader(input_string=cif_records).model()
  mon_lib_srv.process_cif_object(cif_object=cif_object)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    force_symmetry=False,
    log = log)
  # this is to check there is no covalent bonds between two residues
  xrs = processed_pdb_file.xray_structure()
  grm = processed_pdb_file.geometry_restraints_manager()
  rgs=list(processed_pdb_file.all_chain_proxies.pdb_hierarchy.residue_groups())
  rg1_i_seqs_1 = rgs[0].atoms().extract_i_seq()
  rg1_i_seqs_2 = rgs[1].atoms().extract_i_seq()
  bps, asu = grm.get_covalent_bond_proxies(sites_cart = xrs.sites_cart())
  for bps in bps:
    for i_seq in bps.i_seqs:
      r = []
      r.append(i_seq in rg1_i_seqs_1)
      r.append(i_seq in rg1_i_seqs_2)
      assert r.count(True)==1

def exercise_bogus_crystal_symmetry(mon_lib_srv, ener_lib):
  raw_records = """\
REMARK taken from 4arg
CRYST1    1.000    1.000    1.000   1.00   1.00   1.00 P 1           1
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      1.000000  0.000000  0.000000        0.00000
SCALE2      0.000000  1.000000  0.000000        0.00000
SCALE3      0.000000  0.000000  1.000000        0.00000
ATOM      1  CA  PRO A   1      58.034  68.065  61.446  1.00  0.00           C
ATOM      2  CA  ARG A   2      59.821  66.269  58.634  1.00  0.00           C
ATOM      3  CA  THR A   3      56.972  63.793  58.378  1.00  0.00           C
ATOM      4  CA  LEU A   4      54.332  66.493  58.625  1.00  0.00           C
ATOM      5  CA  ASN A   5      55.845  68.102  55.551  1.00  0.00           C
ATOM      6  CA  ALA A   6      55.922  64.739  53.816  1.00  0.00           C
ATOM      7  CA  TRP A   7      52.227  64.264  54.482  1.00  0.00           C
ATOM      8  CA  VAL A   8      51.590  67.865  53.504  1.00  0.00           C
ATOM      9  CA  LYS A   9      52.954  66.907  50.106  1.00  0.00           C
ATOM     10  CA  VAL A  10      50.711  63.858  50.117  1.00  0.00           C
""".splitlines()
  try:
    processed_pdb_file = monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=None,
      raw_records=raw_records,
      force_symmetry=True)
  except Sorry as e:
    assert str(e).startswith(
      "Unit cell volume is incompatible with number of atoms")
  else: raise Exception_expected

# Test automatic detection of CDL
def exercise_cdl_automatic():
  pdbstring = """\
ATOM      0  CA  GLY A   3       5.804  -2.100   7.324  1.00  1.36           C
ATOM      1  C   GLY A   3       4.651  -1.149   7.578  1.00  1.01           C
ATOM      2  O   GLY A   3       3.598  -1.553   8.071  1.00  1.38           O
ATOM      3  N   GLY A   3       6.706  -1.622   6.294  1.00  1.11           N
ATOM      4  CA  PHE A   4       3.819   1.134   7.419  1.00  0.89           C
ATOM      5  CB  PHE A   4       4.397   2.380   8.094  1.00  1.13           C
ATOM      6  C   PHE A   4       3.185   1.509   6.084  1.00  0.94           C
ATOM      7  N   PHE A   4       4.852   0.121   7.242  1.00  0.88           N
ATOM      8  O   PHE A   4       2.361   2.421   6.010  1.00  1.47           O
ATOM      9  CA  LEU A   5       3.055   1.059   3.693  1.00  0.87           C
ATOM     10  CB  LEU A   5       3.965   0.435   2.634  1.00  1.13           C
ATOM     11  C   LEU A   5       1.634   0.527   3.541  1.00  0.87           C
ATOM     12  N   LEU A   5       3.576   0.800   5.030  1.00  0.92           N
ATOM     13  O   LEU A   5       1.246  -0.440   4.196  1.00  1.23           O
"""
  pdb_1 = "REMARK   3    GEOSTD + MON.LIB. + CDL v1.2\n" + pdbstring
  pdb_2 = pdbstring
  with open("tst_cdl_auto_1.pdb", "w") as f:
    f.write(pdb_1)
  with open("tst_cdl_auto_2.pdb", "w") as f:
    f.write(pdb_2)
  cifstring = """\
data_cdl_refine
%s
loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.pdbx_PDB_model_num
  ATOM     1  N    .  GLY  A   1  ?   -9.12200   4.64529   5.82028  1.000  21.01594  N  ?  A  ?   1  1
  ATOM     2  CA   .  GLY  A   1  ?   -9.14139   4.15860   4.44999  1.000  17.74125  C  ?  A  ?   1  1
  ATOM     3  C    .  GLY  A   1  ?   -8.06917   3.11022   4.29606  1.000  17.40881  C  ?  A  ?   1  1
  ATOM     4  O    .  GLY  A   1  ?   -7.63359   2.52149   5.28628  1.000  17.89287  O  ?  A  ?   1  1
  ATOM     5  N    .  ASN  A   2  ?   -7.62610   2.88495   3.06334  1.000  16.11630  N  ?  A  ?   2  1
  ATOM     6  CA   .  ASN  A   2  ?   -6.52361   1.96527   2.82704  1.000  15.02052  C  ?  A  ?   2  1
  ATOM     7  C    .  ASN  A   2  ?   -5.24111   2.51732   3.42319  1.000  13.75642  C  ?  A  ?   2  1
  ATOM     8  O    .  ASN  A   2  ?   -5.04518   3.73206   3.51450  1.000  13.09303  O  ?  A  ?   2  1
  ATOM     9  CB   .  ASN  A   2  ?   -6.34271   1.71701   1.33679  1.000  14.87486  C  ?  A  ?   2  1
  ATOM    10  CG   .  ASN  A   2  ?   -7.61171   1.23635   0.68846  1.000  16.71339  C  ?  A  ?   2  1
  ATOM    11  OD1  .  ASN  A   2  ?   -8.11028   0.15702   1.01777  1.000  19.77837  O  ?  A  ?   2  1
  ATOM    12  ND2  .  ASN  A   2  ?   -8.16730   2.04317  -0.21593  1.000  14.75765  N  ?  A  ?   2  1
"""
  cif_1 = cifstring % """\
_refine.pdbx_stereochemistry_target_values 'GeoStd + Monomer Library + CDL v1.2'
"""
  cif_2 = cifstring % ""
  with open("tst_cdl_auto_1.cif", "w") as f:
    f.write(cif_1)
  with open("tst_cdl_auto_2.cif", "w") as f:
    f.write(cif_2)
  model_files = [
    ("tst_cdl_auto_1.pdb", "tst_cdl_auto_2.pdb"),
    ("tst_cdl_auto_1.cif", "tst_cdl_auto_2.cif"),
  ]
  for file1, file2 in model_files :
    # mmtbx.monomer_library.pdb_interpretation.run will not work with mmCIF,
    # so I'm testing via the command-line tool instead
    from mmtbx.command_line.pdb_interpretation import run
    stdout_save = sys.stdout
    sys.stdout = null_out()
    ppf_1 = run(args=[file1, "cdl=False"])[0]
    ppf_2 = run(args=[file2, "cdl=False"])[0]
    grm_1 = ppf_1.geometry_restraints_manager()
    grm_2 = ppf_2.geometry_restraints_manager()
    assert not None in [grm_1, grm_2]
    assert ppf_1.all_chain_proxies.use_cdl is None, ppf_1.all_chain_proxies.use_cdl
    assert ppf_2.all_chain_proxies.use_cdl is None, ppf_2.all_chain_proxies.use_cdl
    ppf_1 = run(args=[file1, "cdl=Auto"])[0]
    ppf_2 = run(args=[file2, "cdl=Auto"])[0]
    grm_1 = ppf_1.geometry_restraints_manager()
    grm_2 = ppf_2.geometry_restraints_manager()
    assert ppf_1.all_chain_proxies.use_cdl == True
    assert ppf_2.all_chain_proxies.use_cdl is None
    sys.stdout = stdout_save

def exercise_unk_and_cys(mon_lib_srv, ener_lib):
  raw_records = """\
CRYST1  182.247  207.944  177.997  90.00  90.00  90.00 P 21 21 21    8
ATOM  27637  N   CYS X  25     162.023 328.719 193.119  1.00 48.83           N
ATOM  27638  CA  CYS X  25     160.744 328.717 192.370  1.00 48.38           C
ATOM  27639  C   CYS X  25     160.840 329.477 191.030  1.00 48.55           C
ATOM  27640  O   CYS X  25     160.447 328.960 189.978  1.00 47.05           O
ATOM  27641  CB  CYS X  25     159.647 329.420 193.185  1.00 49.20           C
ATOM  27642  SG  CYS X  25     158.050 329.117 192.457  1.00 48.43           S
TER
HETATM28972 UNK  UNX C 262      92.003 325.919 218.014  1.00 80.36           X
END
""".splitlines()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    log=log)

def exercise_ss_bond_angles(mon_lib_srv, ener_lib):
  raw_records = """\
HEADER    HYDROLASE/IMMUNE SYSTEM                 11-MAR-07   2P45
CRYST1   73.447   73.070   42.539  90.00  90.00  90.00 P 21 21 21    4
ATOM   1119  N   CYS B  22      48.085  17.330  19.679  1.00 11.80           N
ATOM   1120  CA  CYS B  22      47.337  16.587  18.675  1.00 11.61           C
ATOM   1121  C   CYS B  22      46.034  17.288  18.406  1.00 11.47           C
ATOM   1122  O   CYS B  22      45.392  17.844  19.314  1.00 12.66           O
ATOM   1123  CB  CYS B  22      47.077  15.179  19.161  1.00 13.15           C
ATOM   1124  SG  CYS B  22      45.959  14.196  18.149  1.00 14.55           S
ATOM   1687  N   CYS B  96      48.918   9.850  18.634  1.00 12.25           N
ATOM   1688  CA  CYS B  96      47.641  10.223  19.169  1.00 11.51           C
ATOM   1689  C   CYS B  96      46.539   9.414  18.528  1.00 10.27           C
ATOM   1690  O   CYS B  96      46.654   8.958  17.371  1.00 12.04           O
ATOM   1691  CB  CYS B  96      47.460  11.706  19.105  1.00 16.73           C
ATOM   1692  SG ACYS B  96      45.815  12.414  19.205  0.50 12.31           S
ATOM   1693  SG BCYS B  96      47.259  12.563  17.544  0.50 18.06           S
END
  """.splitlines()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    log=log)
  grm = processed_pdb_file.geometry_restraints_manager()
  assert grm.angle_proxies.size() == 15
  assert grm.get_dihedral_proxies().size() == 3
  selected_dihedrals = grm.dihedral_proxies.proxy_select(
    n_seq=13, iselection=flex.size_t([4,5,10,11,12]))
  # select SS dihedrals, 2 in this case because of alternative SG atom
  assert selected_dihedrals.size() == 2
  # assert exact values from cif library...
  assert approx_equal(selected_dihedrals[0].angle_ideal, 93.0)
  assert approx_equal(selected_dihedrals[0].alt_angle_ideals[0], -86.0)
  deltas = grm.dihedral_proxies.deltas(
                  processed_pdb_file.xray_structure().sites_cart())
  assert approx_equal(deltas[3], -19.6324864704) # --> -86 degrees
  assert approx_equal(deltas[6], 23.2807269272) # --> 93 degrees

def exercise_ss_bond_angles_alt_loc(mon_lib_srv, ener_lib):
  raw_records = """\
HEADER    HYDROLASE                               19-OCT-05   2BCH
CRYST1   45.908   45.908  101.372  90.00  90.00 120.00 P 31 2 1      6
ATOM    274  N   CYS A  27      17.541   4.439  12.897  1.00 13.99           N
ATOM    275  CA  CYS A  27      16.566   5.527  12.862  1.00 14.57           C
ATOM    276  C   CYS A  27      16.236   6.026  11.467  1.00 14.53           C
ATOM    277  O   CYS A  27      15.254   6.760  11.351  1.00 16.95           O
ATOM    278  CB  CYS A  27      17.114   6.698  13.662  1.00 15.77           C
ATOM    279  SG  CYS A  27      17.230   6.332  15.443  1.00 17.57           S
ATOM   1225  CB  CYS A 123      14.607   7.591  16.260  1.00 24.16           C
ATOM   1226  SG  CYS A 123      15.316   5.939  15.946  1.00 20.05           S
ATOM   1217  N  ACYS A 123      15.023   7.279  18.624  0.58 26.40           N
ATOM   1219  CA ACYS A 123      15.266   8.190  17.491  0.58 25.69           C
ATOM   1221  C  ACYS A 123      14.764   9.599  17.776  0.58 26.33           C
ATOM   1223  O  ACYS A 123      14.197  10.238  16.886  0.58 28.70           O
ATOM   1227  OXTACYS A 123      14.975  10.081  18.878  0.58 28.31           O
ATOM   1218  N  BCYS A 123      15.023   7.288  18.685  0.42 25.68           N
ATOM   1220  CA BCYS A 123      15.108   8.205  17.548  0.42 25.86           C
ATOM   1222  C  BCYS A 123      14.270   9.460  17.813  0.42 26.42           C
ATOM   1224  O  BCYS A 123      13.915  10.125  16.837  0.42 27.75           O
ATOM   1228  OXTBCYS A 123      13.981   9.728  18.968  0.42 28.04           O
END
  """.splitlines()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    log=log)
  grm = processed_pdb_file.geometry_restraints_manager()
  assert grm.angle_proxies.size() == 21
  assert grm.get_dihedral_proxies().size() == 5, \
    "dihedrals %d" % grm.get_dihedral_proxies().size()

def exercise_bad_water(mon_lib_srv, ener_lib):
  raw_records = """\
CRYST1   34.035   17.923   13.510  90.00  90.00  90.00 P 1
SCALE1      0.029382  0.000000  0.000000        0.00000
SCALE2      0.000000  0.055794  0.000000        0.00000
SCALE3      0.000000  0.000000  0.074019        0.00000
ATOM   5961  W   HOH 21100      44.341  86.390 -10.071 1.000 50.00
ATOM   5962  W   HOH 21101      56.357  78.467 -11.839 1.000 50.00
ATOM   5963  W   HOH 21102      32.322  79.906 -13.581 1.000 50.00
END
  """.splitlines()
  e = None
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    log=log)
  msg = processed_pdb_file.all_chain_proxies.fatal_problems_message()
  assert msg == """Fatal problems interpreting model file:
  Number of atoms with unknown scattering type symbols: 3
  Number of atoms with unknown nonbonded energy type symbols: 3
    Please edit the model file to resolve the problems and/or supply a
    CIF file with matching restraint definitions, along with
    apply_cif_modification and apply_cif_link parameter definitions
    if necessary.
    It is best practice to define the element names in
    columns 77-78 of the PDB file.""", msg

def exercise_edits_parallelity(mon_lib_srv, ener_lib):
  raw_records = """\
CRYST1   47.935   37.102   30.520  90.00  90.00  90.00 P 1
ATOM     55  N   PHE A 431      31.800  18.971   9.354  1.00187.57           N
ATOM     56  CA  PHE A 431      30.664  18.262   9.908  1.00187.57           C
ATOM     57  C   PHE A 431      29.679  19.392  10.042  1.00187.57           C
ATOM     58  O   PHE A 431      29.961  20.372  10.729  1.00187.57           O
ATOM     59  CB  PHE A 431      30.976  17.709  11.288  1.00192.11           C
ATOM     60  CG  PHE A 431      32.191  16.840  11.330  1.00192.11           C
ATOM     61  CD1 PHE A 431      32.307  15.750  10.487  1.00192.11           C
ATOM     62  CD2 PHE A 431      33.209  17.101  12.221  1.00192.11           C
ATOM     63  CE1 PHE A 431      33.426  14.943  10.527  1.00192.11           C
ATOM     64  CE2 PHE A 431      34.331  16.298  12.265  1.00192.11           C
ATOM     65  CZ  PHE A 431      34.439  15.217  11.418  1.00192.11           C
ATOM    225  P   A   P  36      27.448   9.314  19.851  1.00 76.88           P
ATOM    226  OP1 A   P  36      26.648   9.676  21.048  1.00 78.59           O
ATOM    227  OP2 A   P  36      28.571   8.351  19.976  1.00 75.20           O
ATOM    228  O5' A   P  36      27.998  10.652  19.190  1.00 78.76           O
ATOM    229  C5' A   P  36      29.217  11.250  19.656  1.00 77.93           C
ATOM    230  C4' A   P  36      30.389  10.753  18.844  1.00 78.66           C
ATOM    231  O4' A   P  36      29.987  10.611  17.462  1.00 79.49           O
ATOM    232  C3' A   P  36      31.616  11.662  18.873  1.00 78.23           C
ATOM    233  O3' A   P  36      32.665  11.095  19.647  1.00 77.17           O
ATOM    234  C2' A   P  36      32.056  11.820  17.411  1.00 79.34           C
ATOM    235  O2' A   P  36      33.369  11.376  17.124  1.00 79.34           O
ATOM    236  C1' A   P  36      31.046  10.987  16.617  1.00 78.22           C
ATOM    237  N9  A   P  36      30.484  11.736  15.494  1.00 20.00           N
ATOM    238  C8  A   P  36      29.463  12.653  15.522  1.00 20.00           C
ATOM    239  N7  A   P  36      29.182  13.166  14.348  1.00 20.00           N
ATOM    240  C5  A   P  36      30.081  12.548  13.491  1.00 20.00           C
ATOM    241  C6  A   P  36      30.298  12.660  12.106  1.00 20.00           C
ATOM    242  N6  A   P  36      29.597  13.468  11.310  1.00 20.00           N
ATOM    243  N1  A   P  36      31.274  11.902  11.561  1.00 20.00           N
ATOM    244  C2  A   P  36      31.978  11.092  12.360  1.00 20.00           C
ATOM    245  N3  A   P  36      31.869  10.898  13.673  1.00 20.00           N
ATOM    246  C4  A   P  36      30.890  11.665  14.184  1.00 20.00           C
END
  """.splitlines()
  edits = """\
pdb_interpretation.geometry_restraints {
    edits {
      plane1 = chain A and resseq 431 and (name cg or name ce1 or name ce2)
      plane2 = chain P and resseq 36 and (name c2 or name c4 or name c6)
      parallelity {
        action = *add delete change
        atom_selection_1 = $plane1
        atom_selection_2 = $plane2
        sigma = 0.027
        target_angle_deg = 0
      }
    }
}"""
  gm_phil = iotbx.phil.parse(
      monomer_library.pdb_interpretation.grand_master_phil_str,
      process_includes=True)
  edits_phil = iotbx.phil.parse(edits)
  working_phil = gm_phil.fetch(edits_phil)
  params = working_phil.extract()
  # print params.geometry_restraints.edits.parallelity[0].atom_selection_1
  assert params.geometry_restraints.edits.parallelity[0].atom_selection_1 == \
      "chain A and resseq 431 and (name cg or name ce1 or name ce2)"
  processed_pdb_file = monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=None,
      raw_records=raw_records,
      params = params.pdb_interpretation,
      log=None)
  grm = processed_pdb_file.geometry_restraints_manager(
      params_edits=params.geometry_restraints.edits,
      params_remove=params.geometry_restraints.remove)
  assert grm.parallelity_proxies.size() == 1
  pp = grm.parallelity_proxies[0]
  # print dir(pp)
  assert list(pp.i_seqs) == [5,8,9]
  assert list(pp.j_seqs) == [27,30,32]
  assert approx_equal(pp.weight, 1371.74211248)
  assert approx_equal(pp.target_angle_deg, 0)

def exercise_edits_planarity(mon_lib_srv, ener_lib):
  raw_records = """\
CRYST1   47.935   37.102   30.520  90.00  90.00  90.00 P 1
ATOM     55  N   PHE A 431      31.800  18.971   9.354  1.00187.57           N
ATOM     56  CA  PHE A 431      30.664  18.262   9.908  1.00187.57           C
ATOM     57  C   PHE A 431      29.679  19.392  10.042  1.00187.57           C
ATOM     58  O   PHE A 431      29.961  20.372  10.729  1.00187.57           O
ATOM     59  CB  PHE A 431      30.976  17.709  11.288  1.00192.11           C
ATOM     60  CG  PHE A 431      32.191  16.840  11.330  1.00192.11           C
ATOM     61  CD1 PHE A 431      32.307  15.750  10.487  1.00192.11           C
ATOM     62  CD2 PHE A 431      33.209  17.101  12.221  1.00192.11           C
ATOM     63  CE1 PHE A 431      33.426  14.943  10.527  1.00192.11           C
ATOM     64  CE2 PHE A 431      34.331  16.298  12.265  1.00192.11           C
ATOM     65  CZ  PHE A 431      34.439  15.217  11.418  1.00192.11           C
ATOM    225  P   A   P  36      27.448   9.314  19.851  1.00 76.88           P
ATOM    226  OP1 A   P  36      26.648   9.676  21.048  1.00 78.59           O
ATOM    227  OP2 A   P  36      28.571   8.351  19.976  1.00 75.20           O
ATOM    228  O5' A   P  36      27.998  10.652  19.190  1.00 78.76           O
ATOM    229  C5' A   P  36      29.217  11.250  19.656  1.00 77.93           C
ATOM    230  C4' A   P  36      30.389  10.753  18.844  1.00 78.66           C
ATOM    231  O4' A   P  36      29.987  10.611  17.462  1.00 79.49           O
ATOM    232  C3' A   P  36      31.616  11.662  18.873  1.00 78.23           C
ATOM    233  O3' A   P  36      32.665  11.095  19.647  1.00 77.17           O
ATOM    234  C2' A   P  36      32.056  11.820  17.411  1.00 79.34           C
ATOM    235  O2' A   P  36      33.369  11.376  17.124  1.00 79.34           O
ATOM    236  C1' A   P  36      31.046  10.987  16.617  1.00 78.22           C
ATOM    237  N9  A   P  36      30.484  11.736  15.494  1.00 20.00           N
ATOM    238  C8  A   P  36      29.463  12.653  15.522  1.00 20.00           C
ATOM    239  N7  A   P  36      29.182  13.166  14.348  1.00 20.00           N
ATOM    240  C5  A   P  36      30.081  12.548  13.491  1.00 20.00           C
ATOM    241  C6  A   P  36      30.298  12.660  12.106  1.00 20.00           C
ATOM    242  N6  A   P  36      29.597  13.468  11.310  1.00 20.00           N
ATOM    243  N1  A   P  36      31.274  11.902  11.561  1.00 20.00           N
ATOM    244  C2  A   P  36      31.978  11.092  12.360  1.00 20.00           C
ATOM    245  N3  A   P  36      31.869  10.898  13.673  1.00 20.00           N
ATOM    246  C4  A   P  36      30.890  11.665  14.184  1.00 20.00           C
END
  """.splitlines()
  edits = """\
pdb_interpretation.geometry_restraints {
    edits {
      planarity {
        action = *add delete change
        atom_selection = resid 431 and (name N or name CA or name C)
        sigma = 0.027
      }
    }
}"""
  gm_phil = iotbx.phil.parse(
      monomer_library.pdb_interpretation.grand_master_phil_str,
      process_includes=True)
  edits_phil = iotbx.phil.parse(edits)
  working_phil = gm_phil.fetch(edits_phil)
  params = working_phil.extract()
  # print  params.geometry_restraints.edits.planarity[0].atom_selection
  assert params.geometry_restraints.edits.planarity[0].atom_selection == \
      "resid 431 and (name N or name CA or name C)"
  processed_pdb_file = monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=None,
      raw_records=raw_records,
      params = params.pdb_interpretation,
      log=None)
  grm = processed_pdb_file.geometry_restraints_manager(
      params_edits=params.geometry_restraints.edits,
      params_remove=params.geometry_restraints.remove)
  assert grm.planarity_proxies.size() == 3
  pp = grm.planarity_proxies[2]
  # print dir(pp)
  assert list(pp.i_seqs) == [0,1,2]
  assert approx_equal(list(pp.weights)[0], 1371.7421124828534)

def exercise_edits_bond(mon_lib_srv, ener_lib):
  from cctbx.geometry_restraints.linking_class import linking_class
  origin_ids = linking_class()

  raw_records = """\
CRYST1   47.935   37.102   30.520  90.00  90.00  90.00 P 1
ATOM     55  N   PHE A 431      31.800  18.971   9.354  1.00187.57           N
ATOM     56  CA  PHE A 431      30.664  18.262   9.908  1.00187.57           C
ATOM     57  C   PHE A 431      29.679  19.392  10.042  1.00187.57           C
ATOM     58  O   PHE A 431      29.961  20.372  10.729  1.00187.57           O
ATOM     59  CB  PHE A 431      30.976  17.709  11.288  1.00192.11           C
ATOM     60  CG  PHE A 431      32.191  16.840  11.330  1.00192.11           C
ATOM     61  CD1 PHE A 431      32.307  15.750  10.487  1.00192.11           C
ATOM     62  CD2 PHE A 431      33.209  17.101  12.221  1.00192.11           C
ATOM     63  CE1 PHE A 431      33.426  14.943  10.527  1.00192.11           C
ATOM     64  CE2 PHE A 431      34.331  16.298  12.265  1.00192.11           C
ATOM     65  CZ  PHE A 431      34.439  15.217  11.418  1.00192.11           C
END
  """.splitlines()
  edits = """\
pdb_interpretation.geometry_restraints {
    edits {
      bond
      {
        action = *add delete change
        atom_selection_1 = name N
        atom_selection_2 = name CB
        symmetry_operation = None
        distance_ideal = 3
        sigma = 0.1
        slack = None
        limit = 1.0
        top_out = True
      }
    }
}"""
  gm_phil = iotbx.phil.parse(
      monomer_library.pdb_interpretation.grand_master_phil_str,
      process_includes=True)
  edits_phil = iotbx.phil.parse(edits)
  working_phil = gm_phil.fetch(edits_phil)
  params = working_phil.extract()
  # print params.geometry_restraints.edits.parallelity[0].atom_selection_1
  assert params.geometry_restraints.edits.bond[0].atom_selection_1 == \
      "name N"
  processed_pdb_file = monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=None,
      raw_records=raw_records,
      params = params.pdb_interpretation,
      log=None)
  grm = processed_pdb_file.geometry_restraints_manager(
      params_edits=params.geometry_restraints.edits,
      params_remove=params.geometry_restraints.remove)
  simple = grm.pair_proxies().bond_proxies.simple
  assert simple.size() == 12
  user_defined = simple.proxy_select(origin_id=origin_ids.get_origin_id('edits'))
  assert user_defined.size() == 1
  udp = user_defined[0]
  assert list(udp.i_seqs) == [0,4]
  assert udp.limit == 1
  assert udp.top_out
  assert approx_equal(udp.distance_ideal, 3, eps=1e-4)
  assert approx_equal(udp.weight, 100, eps=1e-4)

def exercise_bad_custom_bonds(mon_lib_srv, ener_lib):
  raw_records = """
CRYST1   21.213   24.878   21.468  90.00  90.00  90.00 P 1
SCALE1      0.047141  0.000000  0.000000        0.00000
SCALE2      0.000000  0.040196  0.000000        0.00000
SCALE3      0.000000  0.000000  0.046581        0.00000
ATOM      1  N   SER A  92      10.967   6.937   5.256  1.00 19.08           N
ATOM      2  CA  SER A  92      11.206   7.027   6.703  1.00 14.61           C
ATOM      3  C   SER A  92      10.010   7.586   7.440  1.00 17.30           C
ATOM      4  O   SER A  92       9.590   7.129   8.497  1.00 21.21           O
ATOM      5  CB  SER A  92      12.423   7.902   6.959  1.00 16.36           C
ATOM      6  OG  SER A  92      12.172   9.246   6.637  1.00 23.58           O
ATOM      7  N   HIS A  93       9.395   8.633   6.883  1.00 19.20           N
ATOM      8  CA  HIS A  93       8.283   9.235   7.604  1.00 16.97           C
ATOM      9  C   HIS A  93       7.034   8.370   7.604  1.00 17.87           C
ATOM     10  O   HIS A  93       6.276   8.484   8.570  1.00 24.29           O
ATOM     11  CB  HIS A  93       7.948  10.619   7.018  1.00 17.10           C
ATOM     12  CG  HIS A  93       8.966  11.616   7.493  1.00 18.72           C
ATOM     13  ND1 HIS A  93      10.278  11.564   7.100  1.00 19.80           N
ATOM     14  CD2 HIS A  93       8.875  12.654   8.337  1.00 20.88           C
ATOM     15  CE1 HIS A  93      10.963  12.535   7.666  1.00 20.04           C
ATOM     16  NE2 HIS A  93      10.126  13.212   8.419  1.00 23.24           N
ATOM     17  N   ALA A  94       6.823   7.568   6.564  1.00 15.69           N
ATOM     18  CA  ALA A  94       5.666   6.670   6.556  1.00 18.31           C
ATOM     19  C   ALA A  94       5.873   5.431   7.418  1.00 18.40           C
ATOM     20  O   ALA A  94       5.000   5.000   8.170  1.00 19.22           O
ATOM     21  CB  ALA A  94       5.323   6.241   5.146  1.00 18.06           C
HETATM   22  C1A HEM A 154      13.653  13.778   9.638  1.00 18.41           C
HETATM   23  C1B HEM A 154      11.623  16.547   7.120  1.00 20.12           C
HETATM   24  C1C HEM A 154       8.120  16.126   9.474  1.00 13.75           C
HETATM   25  C1D HEM A 154      10.162  13.422  12.047  1.00 14.87           C
HETATM   26  C2A HEM A 154      14.872  13.734   8.859  1.00 21.34           C
HETATM   27  C2B HEM A 154      11.083  17.515   6.191  1.00 18.28           C
HETATM   28  C2C HEM A 154       6.847  16.067  10.168  1.00 20.20           C
HETATM   29  C2D HEM A 154      10.746  12.535  13.032  1.00 18.53           C
HETATM   30  C3A HEM A 154      14.709  14.571   7.805  1.00 19.89           C
HETATM   31  C3B HEM A 154       9.804  17.783   6.579  1.00 14.10           C
HETATM   32  C3C HEM A 154       7.053  15.362  11.318  1.00 20.16           C
HETATM   33  C3D HEM A 154      12.057  12.389  12.708  1.00 17.55           C
HETATM   34  C4A HEM A 154      13.376  15.131   7.909  1.00 14.01           C
HETATM   35  C4B HEM A 154       9.545  16.955   7.748  1.00 15.74           C
HETATM   36  C4C HEM A 154       8.385  14.781  11.212  1.00 11.44           C
HETATM   37  C4D HEM A 154      12.258  13.048  11.435  1.00 12.93           C
HETATM   38  CAA HEM A 154      16.084  12.838   9.160  1.00 17.09           C
HETATM   39  CAB HEM A 154       8.764  18.621   6.096  1.00 18.06           C
HETATM   40  CAC HEM A 154       6.304  15.113  12.506  1.00 17.18           C
HETATM   41  CAD HEM A 154      13.170  11.725  13.527  1.00 19.52           C
HETATM   42  CBA HEM A 154      16.213  11.573   8.272  1.00 15.28           C
HETATM   43  CBB HEM A 154       8.972  19.878   5.454  1.00 10.31           C
HETATM   44  CBC HEM A 154       5.104  15.868  12.800  1.00 23.40           C
HETATM   45  CBD HEM A 154      13.800  12.570  14.648  1.00 26.70           C
HETATM   46  CGA HEM A 154      14.900  10.799   8.148  1.00 14.87           C
HETATM   47  CGD HEM A 154      14.961  11.923  15.376  1.00 23.81           C
HETATM   48  CHA HEM A 154      13.453  12.992  10.753  1.00 13.26           C
HETATM   49  CHB HEM A 154      12.882  16.011   6.977  1.00 14.01           C
HETATM   50  CHC HEM A 154       8.314  16.918   8.363  1.00 13.33           C
HETATM   51  CHD HEM A 154       8.874  13.897  12.145  1.00 11.65           C
HETATM   52  CMA HEM A 154      15.692  14.910   6.685  1.00 14.40           C
HETATM   53  CMB HEM A 154      11.855  18.090   5.000  1.00 11.05           C
HETATM   54  CMC HEM A 154       5.571  16.713   9.646  1.00 21.38           C
HETATM   55  CMD HEM A 154      10.000  11.931  14.225  1.00 17.96           C
HETATM   56  N A HEM A 154      12.773  14.670   9.058  1.00 10.68           N
HETATM   57  N B HEM A 154      10.689  16.271   8.094  1.00 18.79           N
HETATM   58  N C HEM A 154       9.011  15.278  10.091  1.00 15.30           N
HETATM   59  N D HEM A 154      11.098  13.702  11.072  1.00 15.13           N
HETATM   60  O1A HEM A 154      14.486  10.248   9.192  1.00 22.79           O
HETATM   61  O1D HEM A 154      14.703  11.366  16.468  1.00 47.56           O
HETATM   62  O2A HEM A 154      14.357  10.789   7.028  1.00 31.78           O
HETATM   63  O2D HEM A 154      16.078  11.998  14.825  1.00 40.84           O
HETATM   64 FE   HEM A 154      10.863  14.901   9.518  1.00 17.04          Fe
HETATM   65  O   MTO A 155      11.461  16.396  10.964  1.00 20.08           O
TER
END
""".splitlines()
  edits = """\
pdb_interpretation.geometry_restraints.edits {
  bond {
    action = *add delete change
    atom_selection_1 = chain A and resname HEM and resid 154 and name FE
    atom_selection_2 = chain A and resname HIS and resid 93 and name NE2
    distance_ideal = 2.1
    sigma = 0.01
  }
  bond {
    action = *add delete change
    atom_selection_1 = chain A and resname HEM and resid 154 and name FE
    atom_selection_2 = chain A and resname MTO and resid 155 and name O
    symmetry_operation = x,y,z
    distance_ideal = 2.2
    sigma = 0.01
  }
  bond {
    atom_selection_2 = chain N and resname O and name NE
  }
  bond {
    distance_ideal = 2
  }
}"""
  bad_edits = """\
pdb_interpretation.geometry_restraints.edits {
  bond {
    distance_ideal = -1
  }
  bond {
    distance_ideal = 2
    sigma = -1
  }
}
"""
  gm_phil = iotbx.phil.parse(
      monomer_library.pdb_interpretation.grand_master_phil_str,
      process_includes=True)
  edits_phil = iotbx.phil.parse(edits+bad_edits)
  try:
    working_phil = gm_phil.fetch(edits_phil)
  except RuntimeError as e:
    assert str(e) == "geometry_restraints.edits.bond.distance_ideal element is less than the minimum allowed value: -1 < 0.001 (input line 25)", str(e)
  else: raise Exception_expected

  edits_phil = iotbx.phil.parse(edits)
  working_phil = gm_phil.fetch(edits_phil)
  params = working_phil.extract()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    params = params.pdb_interpretation,
    log = log)
  grm = processed_pdb_file.geometry_restraints_manager(
      params_edits=params.geometry_restraints.edits)
  expected = """
  Custom bonds:
    bond:
      atom 1: "HETATM   64 FE   HEM A 154 .*.    FE  "
      atom 2: "ATOM     16  NE2 HIS A  93 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   2.146
      distance_ideal:   2.100
      ideal - model:   -0.046
      slack:            0.000
      delta_slack:     -0.046
      sigma:            0.0100
    bond:
      atom 1: "HETATM   64 FE   HEM A 154 .*.    FE  "
      atom 2: "HETATM   65  O   MTO A 155 .*.     O  "
      symmetry operation: x,y,z
      distance_model:   2.164
      distance_ideal:   2.200
      ideal - model:    0.036
      slack:            0.000
      delta_slack:      0.036
      sigma:            0.0100
    Warning: Ignoring bond with distance_ideal = None:
      atom_selection_1 = None
      atom_selection_2 = "chain N and resname O and name NE"
    Warning: Ignoring bond with sigma = None:
      atom_selection_1 = None
      atom_selection_2 = None
      distance_ideal = 2
    Total number of added/changed bonds: 2"""
  log = log.getvalue().splitlines()
  found = 0
  for e in expected.splitlines():
    if(e in log): found += 1
  assert found > 29, found

def exercise_allow_polymer_cross_special_position(mon_lib_srv, ener_lib):
  raw_records = """\
CRYST1  257.160  257.160  125.289  90.00  90.00  90.00 P 42 21 2
ATOM   6125  N   GLU B 219       2.256  -0.004  52.800  1.00556.70      B    N
ATOM   6126  CA  GLU B 219       2.976  -0.009  54.101  1.00553.07      B    C
ATOM   6127  C   GLU B 219       2.020  -0.291  55.269  1.00548.53      B    C
ATOM   6128  O   GLU B 219       2.412  -1.087  56.152  1.00563.96      B    O
ATOM   6129  CB  GLU B 219       3.635   1.353  54.312  1.00549.67      B    C
ATOM   6130  N   VAL B 220       0.817   0.295  55.291  1.00579.40      B    N
ATOM   6131  CA  VAL B 220      -0.101  -0.076  56.415  1.00536.12      B    C
ATOM   6132  C   VAL B 220      -0.931  -1.312  56.040  1.00595.15      B    C
ATOM   6133  O   VAL B 220      -1.615  -1.845  56.941  1.00560.98      B    O
ATOM   6134  CB  VAL B 220      -0.969   1.109  56.872  1.00506.94      B    C
ATOM   6135  CG1 VAL B 220      -1.904   0.713  58.004  1.00516.72      B    C
ATOM   6136  CG2 VAL B 220      -0.110   2.296  57.278  1.00493.38      B    C
ATOM   6137  N   GLN B 221      -0.855  -1.759  54.781  1.00618.59      B    N
ATOM   6138  CA  GLN B 221      -1.583  -2.967  54.295  1.00643.31      B    C
ATOM   6139  C   GLN B 221      -0.837  -4.237  54.735  1.00659.80      B    C
ATOM   6140  O   GLN B 221      -1.504  -5.170  55.254  1.00664.86      B    O
ATOM   6141  CB  GLN B 221      -1.633  -2.914  52.766  1.00632.59      B    C
ATOM   6142  CG  GLN B 221      -2.226  -4.158  52.118  1.00628.14      B    C
ATOM   6143  CD  GLN B 221      -3.524  -3.878  51.400  1.00607.81      B    C
ATOM   6144  OE1 GLN B 221      -3.819  -2.746  51.025  1.00587.02      B    O
ATOM   6145  NE2 GLN B 221      -4.309  -4.921  51.188  1.00611.94      B    N
END
  """.splitlines()
  params = iotbx.phil.parse(
    monomer_library.pdb_interpretation.grand_master_phil_str,
    process_includes=True).extract()
  exp_lines = [
    "Polymer crosses special position element:",
    "ATOM      7  CA  VAL B 220       0.000  -0.000  56.415  1.00536.12      B    C",
    "Use 'allow_polymer_cross_special_position=True' to keep going."]
  cntr=0
  try:
    processed_pdb_file = monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=None,
      raw_records=raw_records,
      params = params.pdb_interpretation,
      log=None).xray_structure()
  except KeyboardInterrupt: raise
  except Exception as e:
    for line in str(e).splitlines():
      if(line.strip() in exp_lines):
        cntr+=1
  assert cntr>=2
  #
  params.pdb_interpretation.allow_polymer_cross_special_position=True
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    params = params.pdb_interpretation,
    log=None).xray_structure()

def exercise_merging_of_multiple_torsions_from_ccp4(mon_lib_srv, ener_lib):
  raw_records = '''\
CRYST1  173.131   63.222  139.747  90.00 117.44  90.00 C 1 2 1
SCALE1      0.005776  0.000000  0.002999        0.00000
SCALE2      0.000000  0.015817  0.000000        0.00000
SCALE3      0.000000  0.000000  0.008063        0.00000
ATOM     51  P   UFT D   7     -67.463  33.524   3.567  1.00112.63           P
ATOM     52  OP1 UFT D   7     -68.380  34.407   4.348  1.00110.72           O
ATOM     53  OP2 UFT D   7     -66.071  33.991   3.308  1.00125.39           O1-
ATOM     54  O5' UFT D   7     -67.227  32.153   4.352  1.00 96.88           O
ATOM     55  C5' UFT D   7     -66.275  31.194   3.850  1.00 86.30           C
ATOM     56  C4' UFT D   7     -66.867  29.809   3.998  1.00 90.89           C
ATOM     57  O4' UFT D   7     -67.018  29.186   2.691  1.00 89.89           O
ATOM     58  C3' UFT D   7     -66.068  28.802   4.846  1.00 94.41           C
ATOM     59  O3' UFT D   7     -66.587  28.769   6.175  1.00 95.64           O
ATOM     60  C2' UFT D   7     -66.247  27.476   4.085  1.00 90.12           C
ATOM     61  C1' UFT D   7     -66.300  27.971   2.642  1.00 87.60           C
ATOM     62  N1  UFT D   7     -64.948  28.191   2.070  1.00 92.31           N
ATOM     63  C2  UFT D   7     -64.148  27.069   1.895  1.00 97.52           C
ATOM     64  O2  UFT D   7     -64.517  25.929   2.181  1.00 83.90           O
ATOM     65  N3  UFT D   7     -62.899  27.333   1.376  1.00 96.36           N
ATOM     66  C4  UFT D   7     -62.378  28.569   1.029  1.00 96.57           C
ATOM     67  O4  UFT D   7     -61.235  28.641   0.578  1.00104.59           O
ATOM     68  C5  UFT D   7     -63.264  29.677   1.237  1.00 88.14           C
ATOM     69  C6  UFT D   7     -64.492  29.451   1.732  1.00 89.94           C
ATOM     70  F2' UFT D   7     -67.445  26.780   4.352  1.00 91.06           F
  '''.splitlines()
  cif_records = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
UFT UFT "2'-deoxy-2'-fluorouridine 5'-(dihydrogen phosphate)" DNA 31 21 .

data_comp_UFT
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
UFT OP3    O OP   -1 1.868  -46.867 -34.019
UFT P      P P    0  0.567  -47.432 -33.481
UFT OP1    O O    0  -0.201 -48.232 -34.516
UFT OP2    O OP   -1 -0.289 -46.386 -32.791
UFT "O5'"  O O2   0  0.984  -48.486 -32.332
UFT N1     N NR6  0  0.520  -52.299 -30.250
UFT C6     C CR16 0  0.035  -52.133 -31.532
UFT C2     C CR6  0  0.180  -53.419 -29.502
UFT O2     O O    0  0.586  -53.626 -28.364
UFT N3     N NR16 0  -0.668 -54.303 -30.131
UFT C4     C CR6  0  -1.196 -54.193 -31.405
UFT O4     O O    0  -1.946 -55.073 -31.830
UFT C5     C CR16 0  -0.794 -53.014 -32.119
UFT "F2'"  F F    0  3.749  -51.440 -29.039
UFT "C2'"  C CH1  0  2.891  -51.409 -30.143
UFT "C5'"  C CH2  0  1.715  -48.005 -31.174
UFT "C4'"  C CH1  0  2.137  -49.171 -30.313
UFT "O4'"  O O2   0  0.993  -49.989 -29.983
UFT "C1'"  C CH1  0  1.444  -51.292 -29.649
UFT "C3'"  C CH1  0  3.123  -50.141 -30.956
UFT "O3'"  O OH1  0  4.454  -49.679 -30.878
UFT H6     H H    0  0.289  -51.368 -32.013
UFT HN3    H H    0  -0.896 -55.022 -29.658
UFT H5     H H    0  -1.105 -52.862 -32.987
UFT "H2'"  H H    0  3.011  -52.211 -30.708
UFT "H5'"  H H    0  1.142  -47.394 -30.647
UFT "H5'A" H H    0  2.516  -47.504 -31.468
UFT "H4'"  H H    0  2.524  -48.816 -29.478
UFT "H1'"  H H    0  1.433  -51.381 -28.667
UFT "H3'"  H H    0  2.880  -50.294 -31.893
UFT "HO3'" H H    0  4.516  -48.968 -31.338
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
UFT C2e-chi         "O4'" "C1'" N1    C2     210.000 10.000 6
UFT C2e-nyu0        "C4'" "O4'" "C1'" "C2"   340.700 6.300  1
UFT C2e-nyu3        "C2'" "C3'" "C4'" "O4'"  22.600  4.500  1
UFT C2e-nyu4        "C3'" "C4'" "O4'" "C1'"  357.700 6.100  1
UFT C3e-chi         "O4'" "C1'" N1    C2     210.000 10.000 6
UFT C3e-nyu0        "C4'" "O4'" "C1'" "C2'"  2.8     6.100  1
UFT C3e-nyu3        "C2'" "C3'" "C4'" "O4'"  324.700 3.100  1
UFT C3e-nyu4        "C3'" "C4'" "O4'" "C1'"  20.500  5.100  1
UFT alpha           "C5'" "O5'" P     OP3    -60.000 10.00  3
UFT beta            P     "O5'" "C5'" "C4'"  180.000 10.00  3
UFT epsi            "C4'" "C3'" "O3'" "HO3'" 180.000 10.00  3
UFT gamma           "O5'" "C5'" "C4'" "C3'"  180.000 10.00  3
UFT const_11        O4    C4    C5    C6     180.000 10.0   2
UFT sp3_sp3_35      N1    "C1'" "C2'" "F2'"  180.000 10.0   3
UFT sp3_sp3_5       "F2'" "C2'" "C3'" "O3'"  60.000  10.0   3
UFT const_sp2_sp2_1 C5    C6    N1    C2     0.000   5.0    2
UFT const_23        O2    C2    N1    C6     180.000 10.0   2
UFT const_sp2_sp2_5 C4    C5    C6    N1     0.000   5.0    2
UFT const_19        O2    C2    N3    C4     180.000 10.0   2
UFT const_15        O4    C4    N3    C2     180.000 10.0   2

'''
  import iotbx.cif
  log = StringIO()
  cif_object = iotbx.cif.reader(input_string=cif_records).model()
  mon_lib_srv.process_cif_object(cif_object=cif_object,
                                 process_tor=True)
  params = iotbx.phil.parse(
    monomer_library.pdb_interpretation.grand_master_phil_str,
    process_includes=True).extract()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    params = params.pdb_interpretation,
    log=None).xray_structure()

def run(args):
  assert len(args) == 0
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  exercise_merging_of_multiple_torsions_from_ccp4(mon_lib_srv, ener_lib)
  exercise_allow_polymer_cross_special_position(mon_lib_srv, ener_lib)
  exercise_bad_custom_bonds(mon_lib_srv, ener_lib)
  exercise_bad_water(mon_lib_srv, ener_lib)
  exercise_unk_and_cys(mon_lib_srv, ener_lib)
  exercise_cdl_automatic()
  exercise_flattened_cif_loop()
  exercise_bogus_crystal_symmetry(mon_lib_srv, ener_lib)
  exercise_do_not_link(mon_lib_srv, ener_lib)
  exercise_geostd_cif_links(mon_lib_srv, ener_lib)
  exercise_handle_case_insensitive(mon_lib_srv, ener_lib)
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
  exercise_d_amino_acid_chain_perfect_in_box_peptide_plane()
  exercise_rna_v3(mon_lib_srv, ener_lib)
  exercise_asp_glu_acid()
  exercise_rna_dna_synonyms()
  exercise_ss_bond_angles(mon_lib_srv, ener_lib)
  exercise_ss_bond_angles_alt_loc(mon_lib_srv, ener_lib)
  exercise_edits_parallelity(mon_lib_srv, ener_lib)
  exercise_edits_planarity(mon_lib_srv, ener_lib)
  exercise_edits_bond(mon_lib_srv, ener_lib)
  print(format_cpu_times())

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
