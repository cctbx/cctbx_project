from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb
import mmtbx.model
from cctbx import crystal
from libtbx.utils import null_out

pdb_str = """\
CRYST1   27.475   35.191   27.490  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   ALA E   1       6.709  31.817   5.664  1.00 20.00           N
ATOM      2  CA  ALA E   1       7.794  30.914   6.029  1.00 20.00           C
ATOM      3  C   ALA E   1       7.267  29.691   6.773  1.00 20.00           C
ATOM      4  O   ALA E   1       7.942  29.146   7.647  1.00 20.00           O
ATOM      5  CB  ALA E   1       8.831  31.643   6.870  1.00 20.00           C
ATOM      6  N   HIS E   2       6.057  29.267   6.421  1.00 20.00           N
ATOM      7  CA  HIS E   2       5.434  28.110   7.053  1.00 20.00           C
ATOM      8  C   HIS E   2       6.176  26.824   6.701  1.00 20.00           C
ATOM      9  O   HIS E   2       6.477  26.012   7.576  1.00 20.00           O
ATOM     10  CB  HIS E   2       3.965  28.000   6.640  1.00 20.00           C
ATOM     11  CG  HIS E   2       3.140  29.191   7.016  1.00 20.00           C
ATOM     12  ND1 HIS E   2       2.481  29.289   8.222  1.00 20.00           N
ATOM     13  CD2 HIS E   2       2.869  30.335   6.345  1.00 20.00           C
ATOM     14  CE1 HIS E   2       1.838  30.441   8.278  1.00 20.00           C
ATOM     15  NE2 HIS E   2       2.057  31.096   7.152  1.00 20.00           N
ATOM     16  N   CYS E   3       6.468  26.648   5.417  1.00 20.00           N
ATOM     17  CA  CYS E   3       7.193  25.473   4.949  1.00 20.00           C
ATOM     18  C   CYS E   3       8.635  25.492   5.443  1.00 20.00           C
ATOM     19  O   CYS E   3       9.223  24.444   5.712  1.00 20.00           O
ATOM     20  CB  CYS E   3       7.161  25.396   3.421  1.00 20.00           C
ATOM     21  SG  CYS E   3       5.501  25.285   2.714  1.00 20.00           S
"""

def tables(grm):
  return [(len(s.simple_edge_list()), len(s.symmetry_edge_list()))
          for s in grm.shell_sym_tables]

def exercise(max_reasonable_bond_distance):
  """
  Model is moved into a P1 box and back into its crystal cell (as real-space
  refinement does). Bonded shells of the geometry restraints manager stay the
  same, and going back does not take long even when the limit on bond length
  is large (the time used to grow as the cube of that limit).
  """
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(source_info=None, lines=pdb_str),
    log         = null_out())
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.max_reasonable_bond_distance = \
    max_reasonable_bond_distance
  model.process(pdb_interpretation_params=p, make_restraints=True)
  grm = model.get_restraints_manager().geometry
  cs = model.crystal_symmetry()
  tables_start = tables(grm)
  assert tables_start[0][0] > 15 and tables_start[0][1] == 0, tables_start
  energy_start = grm.energies_sites(sites_cart = model.get_sites_cart()).target
  shift = (10.5, 12.5, 13.5)
  model.shift_model_and_set_crystal_symmetry(
    shift_cart       = shift,
    crystal_symmetry = crystal.symmetry((60,70,60,90,90,90), "P1"))
  assert tables(grm) == tables_start
  t0 = time.process_time()
  model.shift_model_and_set_crystal_symmetry(
    shift_cart       = [-s for s in shift],
    crystal_symmetry = cs)
  t = time.process_time()-t0
  assert tables(grm) == tables_start
  energy = grm.energies_sites(sites_cart = model.get_sites_cart()).target
  assert abs(energy-energy_start) < 1.e-6*max(1, abs(energy_start))
  return t

if(__name__ == "__main__"):
  t0 = time.time()
  for d in [50, 250, None]:
    t = exercise(max_reasonable_bond_distance = d)
    print("max_reasonable_bond_distance:", d, "time: %.2f"%t)
    assert t < 0.5, [d, t]
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
