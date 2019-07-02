from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.monomer_library.server
import time
import mmtbx.refinement.real_space
from scitbx.array_family import flex

arg = """
ATOM      0  N   ARG A   1       8.883   8.163   7.325  1.00 11.00           N
ATOM      1  CA  ARG A   1       8.130   8.747   6.222  1.00 12.00           C
ATOM      2  C   ARG A   1       8.767   8.409   4.879  1.00 13.00           C
ATOM      3  O   ARG A   1       8.347   7.471   4.200  1.00 14.00           O
ATOM      4  CB  ARG A   1       8.027  10.265   6.388  1.00 15.00           C
ATOM      5  CG  ARG A   1       7.332  10.717   7.667  1.00 16.00           C
ATOM      6  CD  ARG A   1       5.837  10.424   7.640  1.00 17.00           C
ATOM      7  NE  ARG A   1       5.545   8.996   7.738  1.00 18.00           N
ATOM      8  CZ  ARG A   1       5.359   8.349   8.884  1.00 19.00           C
ATOM      9  NH1 ARG A   1       5.435   9.001  10.036  1.00 18.00           N
ATOM     10  NH2 ARG A   1       5.098   7.049   8.878  1.00 17.00           N
"""

arg_different_order = """
ATOM     10  NH2 ARG A   1       5.098   7.049   8.878  1.00 17.00           N
ATOM      2  C   ARG A   1       8.767   8.409   4.879  1.00 13.00           C
ATOM      4  CB  ARG A   1       8.027  10.265   6.388  1.00 15.00           C
ATOM      5  CG  ARG A   1       7.332  10.717   7.667  1.00 16.00           C
ATOM      6  CD  ARG A   1       5.837  10.424   7.640  1.00 17.00           C
ATOM      7  NE  ARG A   1       5.545   8.996   7.738  1.00 18.00           N
ATOM      0  N   ARG A   1       8.883   8.163   7.325  1.00 11.00           N
ATOM      1  CA  ARG A   1       8.130   8.747   6.222  1.00 12.00           C
ATOM      8  CZ  ARG A   1       5.359   8.349   8.884  1.00 19.00           C
ATOM      9  NH1 ARG A   1       5.435   9.001  10.036  1.00 18.00           N
ATOM      3  O   ARG A   1       8.347   7.471   4.200  1.00 14.00           O
"""

def get_object(pdb_str, backbone_sample, sort_atoms=True):
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str)
  pdb_hierarchy = pdb_inp.construct_hierarchy(sort_atoms=sort_atoms)
  pdb_atoms = pdb_hierarchy.atoms()
  pdb_atoms.reset_i_seq()
  residue = pdb_hierarchy.only_residue()
  mon_lib_srv = mmtbx.monomer_library.server.server()
  t0=time.time()
  result = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
    residue         = residue,
    mon_lib_srv     = mon_lib_srv,
    backbone_sample = backbone_sample).clusters
  print(time.time()-t0)
  return result

def exercise_00():
  #
  aa = get_object(pdb_str = arg, backbone_sample=True)
  assert aa[0].vector == [1, 4, 5, 6, 7, flex.size_t([8, 9, 10])], aa[0].vector
  #
  aa = get_object(pdb_str = arg, backbone_sample=False)
  assert aa[0].vector == [1, 4, 5, 6, 7, flex.size_t([8, 9, 10])], aa[0].vector
  #
  aa = get_object(pdb_str = arg_different_order, backbone_sample=True, sort_atoms=False)
  assert aa[0].vector == [7, 2, 3, 4, 5, flex.size_t([0, 8, 9])], aa[0].vector


if (__name__ == "__main__"):
  t0 = time.time()
  exercise_00()
  print("Time: %6.2f" % (time.time()-t0))
