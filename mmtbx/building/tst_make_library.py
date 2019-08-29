
from __future__ import absolute_import, division, print_function

def exercise_misc():
  pdb_in = """\
ATOM      2  CA  LYS A   1      10.524   2.575   9.811  1.00 25.81           C
ATOM     11  CA  VAL A   2      13.845   2.559  11.678  1.00 24.92           C
ATOM     18  CA  TYR A   3      15.153  -1.021  11.773  1.00 24.62           C
ATOM     29  CA  ARG A   4      17.189  -2.536  14.547  1.00 28.84           C
ATOM     33  CA  GLY A   5      20.515  -4.049  13.441  1.00 23.05           C
ATOM     44  CA  CYS A   6      19.515  -7.712  13.919  1.00 23.43           C
ATOM     50  CA  GLU A   7      15.966  -7.099  12.737  1.00 25.93           C
ATOM     59  CA  LEU A   8      17.296  -5.750   9.429  1.00 20.70           C
ATOM     67  CA  ALA A   9      19.834  -8.568   9.089  1.00 20.77           C
ATOM     72  CA  ALA A  10      17.019 -11.103   9.422  1.00 27.53           C
ATOM     77  CA  ALA A  11      14.781  -9.269   6.951  1.00 26.18           C
ATOM     82  CA  MET A  12      17.641  -9.110   4.439  1.00 24.35           C
ATOM     90  CA  LYS A  13      18.484 -12.793   4.865  1.00 25.35           C
ATOM     99  CA  ARG A  14      14.826 -13.691   4.316  1.00 30.68           C
ATOM    110  CA  HIS A  15      14.826 -11.644   1.110  1.00 28.11           C
ATOM    120  CA  GLY A  16      17.835 -13.509  -0.317  1.00 32.23           C
ATOM    124  CA  LEU A  17      20.703 -11.051   0.228  1.00 28.04           C
ATOM    132  CA  ASP A  18      22.958 -13.359   2.253  1.00 28.33           C
ATOM    140  CA  ASN A  19      25.886 -14.090  -0.062  1.00 27.38           C
ATOM    148  CA  TYR A  20      24.031 -12.556  -3.016  1.00 24.86           C
ATOM    160  CA  ARG A  21      26.612 -12.105  -5.780  1.00 29.59           C
ATOM    171  CA  GLY A  22      29.099 -13.219  -3.147  1.00 24.97           C
ATOM    175  CA  TYR A  23      28.396 -10.439  -0.610  1.00 21.70           C
ATOM    187  CA  SER A  24      27.927 -11.740   2.927  1.00 24.53           C
ATOM    193  CA  LEU A  25      25.085 -10.427   5.048  1.00 23.54           C
ATOM    201  CA  GLY A  26      27.287  -8.119   7.099  1.00 22.24           C
ATOM    205  CA  ASN A  27      28.134  -6.159   3.940  1.00 22.17           C
"""
  from mmtbx.building import make_library
  import iotbx.pdb.hierarchy
  inp = iotbx.pdb.hierarchy.input(pdb_string=pdb_in)
  fragments = make_library.extract_peptide_fragments_by_sequence(
    pdb_hierarchy=inp.hierarchy,
    renumber_from=1,
    sequence='XYRG')
  assert (len(fragments) == 2)

if (__name__ == "__main__"):
  exercise_misc()
  print("OK")
