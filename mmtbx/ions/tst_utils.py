
from __future__ import absolute_import, division, print_function

from libtbx.test_utils import show_diff

def exercise():
  from mmtbx.ions import utils as ion_utils
  import iotbx.pdb
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
CRYST1   20.000   60.000   50.000  90.00  90.00  90.00 P 1
HETATM 1690 ZN    ZN X   1     -14.031  10.147  -2.484  1.00 14.71      ION ZN2+
HETATM 1699 CL    CL X  10     -16.305  10.413  -3.294  0.94 15.45      ION CL1-
HETATM 1700 CL    CL X  11     -13.440   1.272 -23.660  0.64 21.31      ION CL1-
HETATM 1691 MG    MG X   2     -14.099  16.408  -0.840  1.00 14.87      ION MG2+
HETATM 1692 CD    CD X   3      -5.302  -1.809  -1.327  1.00 14.78      ION CD2+
HETATM 1693 CD    CD X   4      -8.287 -11.927 -43.776  1.00 14.70      ION CD2+
HETATM 1703 CL    CL X  14      -2.713  20.673 -12.004  0.79 20.10      ION CL1-
HETATM 1694 NI    NI X   5      -5.160  20.798 -11.755  0.93 14.92      ION NI2+
HETATM 1695 CA    CA X   6      -7.922 -11.718  -0.402  0.74 16.82      ION CA2+
HETATM 1696 CD    CD X   7     -16.886 -19.039 -34.333  0.61 15.22      ION CD2+
HETATM 1701 CL    CL X  12     -10.068 -10.650   0.239  0.53 22.83      ION CL1-
""")
  xrs = pdb_in.xray_structure_simple()
  h = pdb_in.construct_hierarchy()
  pdb_atoms = h.atoms()
  perm = ion_utils.sort_atoms_permutation(
    pdb_atoms=pdb_atoms,
    xray_structure=xrs)
  pdb_atoms = pdb_atoms.select(perm)
  xrs = xrs.select(perm)
  hierarchy = iotbx.pdb.hierarchy.root()
  model = iotbx.pdb.hierarchy.model()
  hierarchy.append_model(model)
  chain = iotbx.pdb.hierarchy.chain(id="X")
  model.append_chain(chain)
  for k, atom in enumerate(pdb_atoms):
    rg = iotbx.pdb.hierarchy.residue_group(resseq="%4d" % (k+1))
    chain.append_residue_group(rg)
    ag = iotbx.pdb.hierarchy.atom_group(resname=atom.parent().resname)
    rg.append_atom_group(ag)
    ag.append_atom(atom.detached_copy())
  assert not show_diff (hierarchy.as_pdb_string(crystal_symmetry=xrs), """\
CRYST1   20.000   60.000   50.000  90.00  90.00  90.00 P 1
SCALE1      0.050000  0.000000  0.000000        0.00000
SCALE2      0.000000  0.016667  0.000000        0.00000
SCALE3      0.000000  0.000000  0.020000        0.00000
HETATM    5 CD    CD X   1      -5.302  -1.809  -1.327  1.00 14.78      ION CD2+
HETATM    6 CD    CD X   2      -8.287 -11.927 -43.776  1.00 14.70      ION CD2+
HETATM   10 CD    CD X   3     -16.886 -19.039 -34.333  0.61 15.22      ION CD2+
HETATM    1 ZN    ZN X   4     -14.031  10.147  -2.484  1.00 14.71      ION ZN2+
HETATM    8 NI    NI X   5      -5.160  20.798 -11.755  0.93 14.92      ION NI2+
HETATM    9 CA    CA X   6      -7.922 -11.718  -0.402  0.74 16.82      ION CA2+
HETATM    2 CL    CL X   7     -16.305  10.413  -3.294  0.94 15.45      ION CL1-
HETATM    7 CL    CL X   8      -2.713  20.673 -12.004  0.79 20.10      ION CL1-
HETATM    3 CL    CL X   9     -13.440   1.272 -23.660  0.64 21.31      ION CL1-
HETATM   11 CL    CL X  10     -10.068 -10.650   0.239  0.53 22.83      ION CL1-
HETATM    4 MG    MG X  11     -14.099  16.408  -0.840  1.00 14.87      ION MG2+
""")

if (__name__ == "__main__"):
  exercise()
  print("OK")
