from __future__ import absolute_import, division, print_function
import iotbx.pdb
from libtbx.test_utils import approx_equal
import mmtbx.refinement.rigid_body

m1_in_str = """\
CRYST1   72.240   72.010   86.990  90.00  90.00  90.00 P 21 21 21    4
ATOM      1  N   ASP L  -1      49.347 -62.804  60.380  1.00 34.60           N
ATOM      2  CA  ASP L  -1      47.975 -63.194  59.946  1.00 33.86           C
ATOM      3  C   ASP L  -1      47.122 -63.665  61.114  1.00 34.02           C
ATOM      4  O   ASP L  -1      47.573 -64.451  61.947  1.00 32.23           O
ATOM      9  N   VAL L  -2      45.889 -63.176  61.175  1.00 31.94           N
ATOM     10  CA  VAL L  -2      44.978 -63.576  62.233  1.00 29.81           C
ATOM     20  CB  GLN L   0      45.127 -68.170  63.455  1.00 29.74           C
ATOM     21  CG  GLN L   0      44.676 -69.611  63.591  1.00 39.72           C
ATOM     22  CD  GLN L   0      45.288 -70.292  64.807  1.00 51.31           C

ATOM     25  N   MET L   4      41.894 -68.026  62.256  1.00 24.27           N
ATOM     26  CA CMET L   4      40.497 -68.318  62.576  0.30 22.89           C
ATOM     27  C  CMET L   4      40.326 -69.824  62.795  0.30 21.48           C
ATOM     28  O  CMET L   4      40.633 -70.625  61.911  0.30 23.73           O

ATOM     25  N   TYR L   4      41.894 -68.026  62.256  1.00 24.27           N
ATOM     26  CA DTYR L   4      40.497 -68.318  62.576  0.70 22.89           C
ATOM     27  C  DTYR L   4      40.326 -69.824  62.795  0.70 21.48           C
ATOM     28  O  DTYR L   4      40.633 -70.625  61.911  0.70 23.73           O

ATOM    186  N   GLN L  27      43.266 -62.813  67.460  1.00 23.64           N
ATOM    187  CA  GLN L  27      43.289 -61.405  67.832  1.00 24.84           C
ATOM    197  C   SER L  27A     39.920 -58.676  67.377  1.00 22.05           C
ATOM    198  O   SER L  27A     40.712 -58.050  66.678  1.00 19.83           O

ATOM    202  CA  LEU L  27B     38.173 -58.446  65.720  1.00 20.41           C
ATOM    203  C   LEU L  27B     37.484 -57.089  65.854  1.00 21.13           C

ATOM    235  N   ASN L  28      36.659 -48.175  64.751  1.00 23.50           N
ATOM    236  CA  ASN L  28      35.418 -47.408  64.733  1.00 19.92           C
TER
ATOM   1702  N   GLY H   1      19.054 -57.059  39.111  1.00 41.20           N
ATOM   1703  CA  GLY H   1      20.002 -57.281  40.238  1.00 41.19           C
ATOM   1704  C   GLY H   1      21.452 -57.171  39.802  1.00 40.03           C

ATOM   1959  N   ASN H  35A     36.815 -57.860  44.928  1.00 15.79           N
ATOM   1960  CA  ASN H  35A     37.803 -58.882  45.265  1.00 16.70           C
ATOM   1967  N   TRP H  36      39.095 -60.772  44.593  1.00 17.23           N
ATOM   1968  CA  TRP H  36      39.440 -62.016  43.926  1.00 19.02           C
TER
ATOM   3383  O   HOH S 215      45.158 -61.466  59.209  1.00 21.71           O
ATOM   3384  O   HOH S 216      26.065 -71.512  48.872  1.00 18.56           O
ATOM   3385  O   HOH S 217      26.885 -90.624  56.632  0.92 16.47           O
ATOM   3386  O   HOH S 218      22.026 -89.723  46.422  1.04 19.34           O
TER
HETATM 3387  O   WAT L 219      33.576 -94.407  39.044  1.08 30.66
HETATM 3388  O   WAT L 220      39.377 -73.892  61.395  1.40 42.96
HETATM 3389  O  AWAT L 221      21.907 -98.365  40.731  0.40 16.45
HETATM 3390  O  BWAT L 221      28.591 -95.257  59.711  0.60 20.75
TER
ATOM   1001 AU   AU    500      14.333   3.856  26.301  1.00  7.97          Au
TER
ATOM   3342  N   VAL   226      20.817 -97.348  23.861  1.00 34.12           N
ATOM   3349  N   PRO   227      19.065-100.553  23.842  1.00 41.92           N
ATOM   3350  CA  PRO   227      19.109-101.900  24.416  1.00 46.30           C
ATOM   3351  C   PRO   227      20.202-102.722  23.753  1.00 50.65           C
ATOM   3356  N   ARG   228      20.811-103.620  24.514  1.00 56.77           N
TER
ATOM     25  N   MET X   4      41.894 -68.026  62.256  1.00 24.27           N
TER
ATOM   1001 AU   AU  K 502      11.333   3.856  26.301  1.00  7.97
ATOM   1001 AU   AU  K 100      14.333   0.856  26.301  1.00  7.97          Au
ATOM   1001 AU   AU  K 900      14.333   3.856  21.301  1.00  7.97          Au
END
"""

m1_out_str = """\
CRYST1   72.240   72.010   86.990  90.00  90.00  90.00 P 21 21 21    4
ATOM      1  N   ASP L  -1      49.347 -62.804  60.380  1.00 34.60           N
ATOM      2  CA  ASP L  -1      47.975 -63.194  59.946  1.00 33.86           C
ATOM      3  C   ASP L  -1      47.122 -63.665  61.114  1.00 34.02           C
ATOM      4  O   ASP L  -1      47.573 -64.451  61.947  1.00 32.23           O
ATOM      9  N   VAL L  -2      45.889 -63.176  61.175  1.00 31.94           N
ATOM     10  CA  VAL L  -2      44.978 -63.576  62.233  1.00 29.81           C
ATOM     20  CB  GLN L   0      45.127 -68.170  63.455  1.00 29.74           C
ATOM     21  CG  GLN L   0      44.676 -69.611  63.591  1.00 39.72           C
ATOM     22  CD  GLN L   0      45.288 -70.292  64.807  1.00 51.31           C

ATOM     25  N   MET L   4      41.894 -68.026  62.256  1.00 24.27           N
ATOM     26  CA CMET L   4      40.497 -68.318  62.576  0.30 22.89           C
ATOM     27  C  CMET L   4      40.326 -69.824  62.795  0.30 21.48           C
ATOM     28  O  CMET L   4      40.633 -70.625  61.911  0.30 23.73           O

ATOM     25  N   TYR L   4      41.894 -68.026  62.256  1.00 24.27           N
ATOM     26  CA DTYR L   4      40.497 -68.318  62.576  0.70 22.89           C
ATOM     27  C  DTYR L   4      40.326 -69.824  62.795  0.70 21.48           C
ATOM     28  O  DTYR L   4      40.633 -70.625  61.911  0.70 23.73           O

ATOM    186  N   GLN L  27      43.266 -62.813  67.460  1.00 23.64           N
ATOM    187  CA  GLN L  27      43.289 -61.405  67.832  1.00 24.84           C
ATOM    197  C   SER L  27A     39.920 -58.676  67.377  1.00 22.05           C
ATOM    198  O   SER L  27A     40.712 -58.050  66.678  1.00 19.83           O

ATOM    202  CA  LEU L  27B     38.173 -58.446  65.720  1.00 20.41           C
ATOM    203  C   LEU L  27B     37.484 -57.089  65.854  1.00 21.13           C

ATOM    235  N   ASN L  28      36.659 -48.175  64.751  1.00 23.50           N
ATOM    236  CA  ASN L  28      35.418 -47.408  64.733  1.00 19.92           C
TER
ATOM   1702  N   GLY H   1      19.054 -57.059  39.111  1.00 41.20           N
ATOM   1703  CA  GLY H   1      20.002 -57.281  40.238  1.00 41.19           C
ATOM   1704  C   GLY H   1      21.452 -57.171  39.802  1.00 40.03           C

ATOM   1959  N   ASN H  35A     36.815 -57.860  44.928  1.00 15.79           N
ATOM   1960  CA  ASN H  35A     37.803 -58.882  45.265  1.00 16.70           C
ATOM   1967  N   TRP H  36      39.095 -60.772  44.593  1.00 17.23           N
ATOM   1968  CA  TRP H  36      39.440 -62.016  43.926  1.00 19.02           C
TER
ATOM   3342  N   VAL   226      20.817 -97.348  23.861  1.00 34.12           N
ATOM   3349  N   PRO   227      19.065-100.553  23.842  1.00 41.92           N
ATOM   3350  CA  PRO   227      19.109-101.900  24.416  1.00 46.30           C
ATOM   3351  C   PRO   227      20.202-102.722  23.753  1.00 50.65           C
ATOM   3356  N   ARG   228      20.811-103.620  24.514  1.00 56.77           N
TER
END
"""

m2_in_str = """\
CRYST1   72.240   72.010   86.990  90.00  90.00  90.00 P 21 21 21    4
ATOM      1  O5'  DC A   1      29.670  19.013  -2.305  1.00 11.26           O
ATOM      2  C5'  DC A   1      30.109  19.727  -1.144  1.00  2.03           C
ATOM      3  C4'  DC A   1      30.024  19.014   0.202  1.00  4.76           C
ATOM      4  O4'  DC A   1      28.912  18.160   0.419  1.00  7.08           O
ATOM      5  C3'  DC A   1      31.204  18.275   0.705  1.00  9.25           C
ATOM      6  O3'  DC A   1      31.448  18.682   2.066  1.00 11.24           O
ATOM     20  P    DG A   2      32.339  20.004   2.406  1.00 11.95           P
ATOM     45  P    DC A   3      28.538  23.069   6.214  1.00 19.73           P
ATOM     46  OP1  DC A   3      29.232  24.056   7.091  1.00 19.15           O
ATOM     47  OP2  DC A   3      29.008  21.652   6.240  1.00 19.23           O
ATOM     48  O5'  DC A   3      26.905  23.120   6.436  1.00 16.82           O
TER
ATOM    139  O5'  DC B   7      13.405  26.346   3.336  1.00 15.19           O
ATOM    140  C5'  DC B   7      13.472  24.969   2.944  1.00  4.19           C
ATOM    158  P    DG B   8      11.130  21.728   3.702  1.00 12.13           P
ATOM    159  OP1  DG B   8      10.341  22.682   2.901  1.00 14.52           O
ATOM    186  O5'  DC B   9      16.676  16.828   1.638  1.00 13.48           O
ATOM    187  C5'  DC B   9      17.151  15.620   2.196  1.00  3.67           C
ATOM    188  C4'  DC B   9      17.855  15.745   3.510  1.00  2.48           C
TER
HETATM  277  O   HOH A  13      32.050  21.145   7.029  1.00 28.46           O
HETATM  278  O   HOH A  14      29.068  18.914   5.393  1.00 44.50           O
HETATM  279  O   HOH A  15      29.984  22.375  10.338  1.00 42.18           O
END
"""

m2_out_str = """\
CRYST1   72.240   72.010   86.990  90.00  90.00  90.00 P 21 21 21    4
ATOM      1  O5'  DC A   1      29.670  19.013  -2.305  1.00 11.26           O
ATOM      2  C5'  DC A   1      30.109  19.727  -1.144  1.00  2.03           C
ATOM      3  C4'  DC A   1      30.024  19.014   0.202  1.00  4.76           C
ATOM      4  O4'  DC A   1      28.912  18.160   0.419  1.00  7.08           O
ATOM      5  C3'  DC A   1      31.204  18.275   0.705  1.00  9.25           C
ATOM      6  O3'  DC A   1      31.448  18.682   2.066  1.00 11.24           O
ATOM     20  P    DG A   2      32.339  20.004   2.406  1.00 11.95           P
ATOM     45  P    DC A   3      28.538  23.069   6.214  1.00 19.73           P
ATOM     46  OP1  DC A   3      29.232  24.056   7.091  1.00 19.15           O
ATOM     47  OP2  DC A   3      29.008  21.652   6.240  1.00 19.23           O
ATOM     48  O5'  DC A   3      26.905  23.120   6.436  1.00 16.82           O
TER
ATOM    139  O5'  DC B   7      13.405  26.346   3.336  1.00 15.19           O
ATOM    140  C5'  DC B   7      13.472  24.969   2.944  1.00  4.19           C
ATOM    158  P    DG B   8      11.130  21.728   3.702  1.00 12.13           P
ATOM    159  OP1  DG B   8      10.341  22.682   2.901  1.00 14.52           O
ATOM    186  O5'  DC B   9      16.676  16.828   1.638  1.00 13.48           O
ATOM    187  C5'  DC B   9      17.151  15.620   2.196  1.00  3.67           C
ATOM    188  C4'  DC B   9      17.855  15.745   3.510  1.00  2.48           C
TER
END
"""

m3_in_str = """\
CRYST1   72.240   72.010   86.990  90.00  90.00  90.00 P 21 21 21    4
ATOM      1  N   ASP L  -1      49.347 -62.804  60.380  1.00 34.60           N
ATOM      2  CA  ASP L  -1      47.975 -63.194  59.946  1.00 33.86           C
ATOM      3  C   ASP L  -1      47.122 -63.665  61.114  1.00 34.02           C
ATOM      4  O   ASP L  -1      47.573 -64.451  61.947  1.00 32.23           O
ATOM      9  N   VAL L  -2      45.889 -63.176  61.175  1.00 31.94           N
ATOM     10  CA  VAL L  -2      44.978 -63.576  62.233  1.00 29.81           C
ATOM     20  CB  GLN L   0      45.127 -68.170  63.455  1.00 29.74           C
ATOM     21  CG  GLN L   0      44.676 -69.611  63.591  1.00 39.72           C
ATOM     22  CD  GLN L   0      45.288 -70.292  64.807  1.00 51.31           C
ATOM     25  N   MET L   4      41.894 -68.026  62.256  1.00 24.27           N
ATOM     26  CA CMET L   4      40.497 -68.318  62.576  0.30 22.89           C
ATOM     27  C  CMET L   4      40.326 -69.824  62.795  0.30 21.48           C
ATOM     28  O  CMET L   4      40.633 -70.625  61.911  0.30 23.73           O
TER
HETATM 3387  O   WAT L 219      33.576 -94.407  39.044  1.08 30.66
HETATM 3388  O   WAT L 220      39.377 -73.892  61.395  1.40 42.96
HETATM 3389  O  AWAT L 221      21.907 -98.365  40.731  0.40 16.45
HETATM 3390  O  BWAT L 221      28.591 -95.257  59.711  0.60 20.75
TER
ATOM   1001 AU   AU  L 500      14.333   3.856  26.301  1.00  7.97          Au
TER
END
"""

def exercise_00():
  pdb_inp = iotbx.pdb.input(source_info=None, lines=m1_in_str)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  sel_strings=mmtbx.refinement.rigid_body.rigid_groups_from_pdb_chains(
    pdb_hierarchy=pdb_hierarchy)
  answer = [
    "(chain 'L' and resid -1 through 28)",
    "(chain 'H' and resid 1 through 36)",
    "(chain ' ' and resid 226 through 228)"]
  for sel_string in sel_strings:
    assert sel_string in answer
  assert len(sel_strings) == 3
  selection = pdb_hierarchy.atom_selection_cache().selection(
    string=" or ".join(answer))
  xrs_1 = pdb_inp.xray_structure_simple().select(selection = selection)
  xrs_2 = iotbx.pdb.input(source_info=None,
    lines=m1_out_str).xray_structure_simple()
  assert approx_equal(
    xrs_1.structure_factors(d_min=2).f_calc().data(),
    xrs_2.structure_factors(d_min=2).f_calc().data())

def exercise_01():
  pdb_inp = iotbx.pdb.input(source_info=None, lines=m2_in_str)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  sel_strings=mmtbx.refinement.rigid_body.rigid_groups_from_pdb_chains(
    pdb_hierarchy=pdb_hierarchy)
  answer = [
    "(chain 'A' and resid 1 through 3)",
    "(chain 'B' and resid 7 through 9)"]
  for sel_string in sel_strings:
    assert sel_string in answer
  assert len(sel_strings) == 2
  selection = pdb_hierarchy.atom_selection_cache().selection(
    string=" or ".join(answer))
  xrs_1 = pdb_inp.xray_structure_simple().select(selection = selection)
  xrs_2 = iotbx.pdb.input(source_info=None,
    lines=m2_out_str).xray_structure_simple()
  assert approx_equal(
    xrs_1.structure_factors(d_min=2).f_calc().data(),
    xrs_2.structure_factors(d_min=2).f_calc().data())

def exercise_02():
  pdb_inp = iotbx.pdb.input(source_info=None, lines=m3_in_str)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  sel_strings=mmtbx.refinement.rigid_body.rigid_groups_from_pdb_chains(
    pdb_hierarchy=pdb_hierarchy,
    group_all_by_chain=True)
  assert sel_strings == ["(chain 'L' and not (resname HOH or resname WAT))"]

if (__name__ == "__main__"):
  exercise_00()
  exercise_01()
  exercise_02()
  print("OK")
