from __future__ import absolute_import, division, print_function
from cctbx import geometry_restraints
import sys

class source_info_server(object):
  def n_expected_atoms(self): return 0
  def labels(self): return "labels"

def run(args):
  assert len(args) == 0
  #
  r = geometry_restraints.angle_proxy_registry(
    strict_conflict_handling=True)
  r.initialize_table()
  r.process(
    source_info=source_info_server(),
    proxy=geometry_restraints.angle_proxy(
      i_seqs=(4,2,7), angle_ideal=0, weight=1))
  r.process(
    source_info=source_info_server(),
    proxy=geometry_restraints.angle_proxy(
      i_seqs=(5,2,3), angle_ideal=0, weight=1))
  assert r.lookup_i_proxy((4,2,7)) == 0
  assert r.lookup_i_proxy((7,2,4)) == 0
  assert r.lookup_i_proxy((2,7,4)) is None
  assert r.lookup_i_proxy((3,2,5)) == 1
  assert r.lookup_i_proxy((5,2,3)) == 1
  assert r.lookup_i_proxy((5,3,2)) is None
  #
  r = geometry_restraints.dihedral_proxy_registry(
    strict_conflict_handling=True)
  r.initialize_table()
  r.process(
    source_info=source_info_server(),
    proxy=geometry_restraints.dihedral_proxy(
      i_seqs=(4,7,2,3), angle_ideal=0, weight=1))
  r.process(
    source_info=source_info_server(),
    proxy=geometry_restraints.dihedral_proxy(
      i_seqs=(8,2,3,5), angle_ideal=0, weight=1))
  r.process(
    source_info=source_info_server(),
    proxy=geometry_restraints.dihedral_proxy(
      i_seqs=(6,3,1,9), angle_ideal=0, weight=1))
  assert r.lookup_i_proxy((3,2,7,4)) == (0, 1)
  assert r.lookup_i_proxy((4,2,7,3)) == (0, -1)
  assert r.lookup_i_proxy((4,7,2,3)) == (0, 1)
  assert r.lookup_i_proxy((3,7,2,4)) == (0, -1)
  assert r.lookup_i_proxy((5,2,3,8)) == (1, 1)
  assert r.lookup_i_proxy((6,3,1,9)) == (2, -1)
  #
  print("OK")

if (__name__ == "__main__"):
  run(sys.argv[1:])
