
from __future__ import division
from __future__ import print_function
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out

def exercise () :
  import mmtbx.validation.ligands
  import iotbx.pdb.hierarchy
  pdb_1 = iotbx.pdb.hierarchy.input(pdb_string="""\
HETATM    1  C   ACT X   1       0.496   0.209   0.702  1.00 20.00      A    C
HETATM    2  O   ACT X   1       1.691   0.458   0.395  1.00 20.00      A    O
HETATM    3  OXT ACT X   1       0.100   0.388   1.883  1.00 20.00      A    O-1
HETATM    4  CH3 ACT X   1      -0.463  -0.305  -0.349  1.00 20.00      A    C
HETATM    5  H1  ACT X   1      -0.904   0.449  -0.784  1.00 20.00      A    H
HETATM    6  H2  ACT X   1      -1.135  -0.874   0.074  1.00 20.00      A    H
HETATM    7  H3  ACT X   1       0.029  -0.823  -1.014  1.00 20.00      A    H
END""")
  pdb_2 = iotbx.pdb.hierarchy.input(pdb_string="""\
HETATM    1  C   ACT X   1       0.412  -0.000   0.714  1.00 20.00      A    C
HETATM    2  O   ACT X   1       1.671  -0.000   0.714  1.00 20.00      A    O
HETATM    3  OXT ACT X   1      -0.217  -0.000   1.804  1.00 20.00      A    O-1
HETATM    4  CH3 ACT X   1      -0.344  -0.000  -0.597  1.00 20.00      A    C
HETATM    5  H1  ACT X   1      -0.507   0.920  -0.878  1.00 20.00      A    H
HETATM    6  H2  ACT X   1      -1.197  -0.460  -0.480  1.00 20.00      A    H
HETATM    7  H3  ACT X   1       0.183  -0.460  -1.277  1.00 20.00      A    H
""")
  rmsds = mmtbx.validation.ligands.compare_ligands(
    ligand_code="ACT",
    hierarchy_1=pdb_1.hierarchy,
    hierarchy_2=pdb_2.hierarchy,
    out=null_out())
  assert (len(rmsds) == 1)
  assert approx_equal(rmsds[0][0], 0.444, eps=0.0001)
  print("OK")

if (__name__ == "__main__") :
  exercise()
