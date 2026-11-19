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


asp_altloc_h = """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM    732  N   ASP A 159       9.168  28.378 -20.499  1.00 12.57           N
ATOM    733  CA  ASP A 159       8.901  27.254 -21.405  1.00 11.91           C
ATOM    734  C   ASP A 159       9.546  27.572 -22.751  1.00 11.04           C
ATOM    735  O   ASP A 159       9.135  28.491 -23.466  1.00 12.99           O
ATOM    736  CB  ASP A 159       7.373  27.114 -21.524  1.00 11.97           C
ATOM    737  H   ASP A 159       8.676  29.059 -20.683  1.00 12.57           H
ATOM    738  HA  ASP A 159       9.259  26.409 -21.092  1.00 11.91           H
ATOM    739  HB2 ASP A 159       7.012  26.903 -20.649  1.00 11.97           H
ATOM    740  HB3 ASP A 159       7.008  27.956 -21.840  1.00 11.97           H
ATOM    741  CG AASP A 159       6.916  26.027 -22.483  0.80 14.95           C
ATOM    742  OD1AASP A 159       7.698  25.447 -23.248  0.80 10.06           O
ATOM    743  OD2AASP A 159       5.687  25.772 -22.489  0.80 20.74           O
"""

def exercise_01():
  """
  Residue with a single altloc conformer (side chain altloc A, occupancy
  0.8, no B; as in PDB 2vfn ASP A 159) and H added: the conformer residue
  lists the blank-altloc atoms (incl. HB2, HB3) before the altloc A atoms, so
  the lowest-index atom of the chi1 rotating group is HB2, not CG. Such a
  residue has one conformer and is selected for real-space side-chain fitting
  in phenix.refine; with an all-H cluster selection score3 (exclude_hd=True)
  divided by zero. The cluster evaluation atom must be the heavy atom CG.
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines = asp_altloc_h)
  h = pdb_inp.construct_hierarchy()
  mon_lib_srv = mmtbx.monomer_library.server.server()
  for conformer in h.only_chain().conformers():
    residue = conformer.only_residue()
    atoms = residue.atoms()
    clusters = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
      residue         = residue,
      mon_lib_srv     = mon_lib_srv,
      backbone_sample = False).clusters
    assert len(clusters) == 2, len(clusters)
    assert conformer.altloc == "A"
    for cl in clusters:
      names = [atoms[i].name.strip() for i in cl.selection]
      assert len([i for i in cl.selection
        if not atoms[i].element_is_hydrogen()]) > 0, names
    names = [atoms[i].name.strip() for i in clusters[0].selection]
    assert names == ["CG"], (conformer.altloc, names)

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_00()
  exercise_01()
  print("Time: %6.2f" % (time.time()-t0))
