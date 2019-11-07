from __future__ import absolute_import, division, print_function

pdb = """HETATM    1  N   ALA A   1      -0.424   1.960   3.877  1.00 20.00      A    N+1
HETATM    2  H   ALA A   1       0.452   1.694   3.861  1.00 20.00      A    H
HETATM    3  H2  ALA A   1      -0.472   2.843   4.121  1.00 20.00      A    H
HETATM    4  H3  ALA A   1      -0.888   1.448   4.484  1.00 20.00      A    H
HETATM    5  CA  ALA A   1      -0.994   1.794   2.582  1.00 20.00      A    C
HETATM    6  HA  ALA A   1      -0.584   2.442   1.977  1.00 20.00      A    H
HETATM    7  CB  ALA A   1      -2.476   2.066   2.644  1.00 20.00      A    C
HETATM    8  HB1 ALA A   1      -2.627   2.987   2.945  1.00 20.00      A    H
HETATM    9  HB2 ALA A   1      -2.867   1.945   1.755  1.00 20.00      A    H
HETATM   10  HB3 ALA A   1      -2.896   1.444   3.272  1.00 20.00      A    H
HETATM   11  C   ALA A   1      -0.695   0.397   2.070  1.00 20.00      A    C
HETATM   12  O   ALA A   1      -1.218  -0.566   2.593  1.00 20.00      A    O
HETATM   13  N   ALA A   2       0.165   0.219   0.945  1.00 20.00      A    N
HETATM   14  H   ALA A   2       0.578   0.941   0.568  1.00 20.00      A    H
HETATM   15  CA  ALA A   2       0.460  -1.094   0.477  1.00 20.00      A    C
HETATM   16  HA  ALA A   2      -0.095  -1.682   1.027  1.00 20.00      A    H
HETATM   17  CB  ALA A   2       1.882  -1.540   0.755  1.00 20.00      A    C
HETATM   18  HB1 ALA A   2       2.509  -0.921   0.323  1.00 20.00      A    H
HETATM   19  HB2 ALA A   2       2.039  -1.544   1.721  1.00 20.00      A    H
HETATM   20  HB3 ALA A   2       2.016  -2.441   0.399  1.00 20.00      A    H
HETATM   21  C   ALA A   2      -0.001  -1.357  -0.949  1.00 20.00      A    C
HETATM   22  O   ALA A   2      -1.001  -2.020  -1.152  1.00 20.00      A    O
HETATM   23  N   ALA A   3       0.750  -0.851  -2.042  1.00 20.00      A    N
HETATM   24  H   ALA A   3       1.494  -0.350  -1.890  1.00 20.00      A    H
HETATM   25  CA  ALA A   3       0.361  -1.157  -3.372  1.00 20.00      A    C
HETATM   26  HA  ALA A   3      -0.113  -2.013  -3.369  1.00 20.00      A    H
HETATM   27  CB  ALA A   3       1.587  -1.293  -4.240  1.00 20.00      A    C
HETATM   28  HB1 ALA A   3       2.124  -0.477  -4.175  1.00 20.00      A    H
HETATM   29  HB2 ALA A   3       1.315  -1.431  -5.171  1.00 20.00      A    H
HETATM   30  HB3 ALA A   3       2.118  -2.058  -3.937  1.00 20.00      A    H
HETATM   31  C   ALA A   3      -0.576  -0.091  -3.884  1.00 20.00      A    C
HETATM   32  O   ALA A   3      -1.633  -0.418  -4.488  1.00 20.00      A    O
HETATM   33  OXT ALA A   3      -0.293   1.123  -3.725  1.00 20.00      A    O-1
"""

from six.moves import cStringIO as StringIO
from libtbx import easy_run

def run():
  f=open("tst_correct_hydrogens.pdb", "w")
  f.write(pdb)
  f.close()
  cmd = "phenix.geometry_minimization tst_correct_hydrogens.pdb"
  print(cmd)
  ero = easy_run.fully_buffered(command=cmd)
  cmd = "phenix.geometry_minimization tst_correct_hydrogens_minimized.pdb"
  print(cmd)
  ero = easy_run.fully_buffered(command=cmd)
  std = StringIO()
  ero.show_stdout(out=std)
  chiral_energy = None
  for line in std.getvalue().splitlines():
    if line.find("chirality_residual_sum")>-1:
      chiral_energy = line
  chiral_energy = chiral_energy.split()[-1]
  print('chiral energy',chiral_energy)
  assert float(chiral_energy)<0.03, 'chiral energy %0.2f greater than 0.03' % float(chiral_energy)
  print("OK")

if __name__=="__main__":
  run()
