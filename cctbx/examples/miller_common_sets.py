from __future__ import absolute_import, division, print_function
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex

def demo():
  #
  # create toy lists of Miller indices
  #
  crystal_symmetry = crystal.symmetry(
    unit_cell=(13,14,15,90,90,90),
    space_group_symbol="P212121")
  miller_set_a = miller.set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=False,
    indices=flex.miller_index([
      (0, -1, 2),
      (-1, -2, 3),
      (2, 3, 4),
      (-3, 4, 5),
      (4, -5, 6)]))
  miller_set_b = miller_set_a.customized_copy(
    indices=flex.miller_index([
      (-5, -6, 7),
      (0, 1, 2),
      (3, -4, 5),
      (-1, -2, 3),
      (-4, -5, 6)]))
  #
  # map all indices to the asymmetric unit
  #
  asu_a = miller_set_a.map_to_asu()
  asu_b = miller_set_b.map_to_asu()
  for h in asu_a.indices():
    print("asu a:", h)
  print()
  for h in asu_b.indices():
    print("asu b:", h)
  print()
  #
  # obtain the common index sets
  #
  common_a, common_b = asu_a.common_sets(asu_b)
  for h in common_a.indices():
    print("common a:", h)
  print()
  for h in common_b.indices():
    print("common b:", h)
  print()
  #
  # obtain the "lone" index sets
  #
  lone_set_a, lone_set_b = asu_a.lone_sets(asu_b)
  for h in lone_set_a.indices():
    print("lone a:", h)
  print()
  for h in lone_set_b.indices():
    print("lone b:", h)
  print()
  #
  # now the same again, but with data (i.e. miller.array instances)
  #
  miller_array_a = miller_set_a.array(
    data=flex.random_double(size=miller_set_a.indices().size()))
  miller_array_b = miller_set_b.array(
    data=flex.random_double(size=miller_set_a.indices().size()))
  #
  # map all indices to the asymmetric unit
  #
  asu_a = miller_array_a.map_to_asu()
  asu_b = miller_array_b.map_to_asu()
  asu_a.show_array(prefix="asu a: ")
  print()
  asu_b.show_array(prefix="asu b: ")
  print()
  #
  # obtain the common index sets
  #
  common_a, common_b = asu_a.common_sets(asu_b)
  common_a.show_array(prefix="common a: ")
  print()
  common_b.show_array(prefix="common b: ")
  print()
  #
  # obtain the "lone" index sets
  #
  lone_a, lone_b = asu_a.lone_sets(asu_b)
  lone_a.show_array(prefix="lone a: ")
  print()
  lone_b.show_array(prefix="lone b: ")
  print()
  print("OK")

if (__name__ == "__main__"):
  demo()
