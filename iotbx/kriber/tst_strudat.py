from __future__ import absolute_import, division, print_function
from iotbx.kriber import strudat
from cctbx.regression import tst_direct_space_asu
from cctbx import crystal
from cctbx import uctbx
from cctbx.array_family import flex
from libtbx.test_utils import Exception_expected
from libtbx.utils import format_cpu_times
import libtbx.load_env
from six.moves import cStringIO as StringIO
import sys, os
from six.moves import zip

def exercise_basic():
  test_file = StringIO("""
*tric
Title
Reference
P 1
 11 12 13 100 110 120
Si 0.1 0.2 0.3
---------------------------------
*mono_b
Title
Reference
P 1 2 1
 11 12 13 100
Si 0.1 0.2 0.3
---------------------------------
*mono_c
Title
Reference
P 1 1 2
 11 12 13 100
Si 0.1 0.2 0.3
---------------------------------
*mono_a
Title
Reference
P 2 1 1
 11 12 13 100
Si 0.1 0.2 0.3
---------------------------------
*orth
Title
Reference
P 2 2 2
 11 12 13
Si 0.1 0.2 0.3 # remark
---------------------------------
*tetr
Title
Reference
P 4
 11 13
Si 0.1 0.2 0.3 4
---------------------------------
*trig
Title
Reference
R 3
 11 13
Si 0.1 0.2 0.3
---------------------------------
*rhom
Title
Reference
R 3 R
 11 100
Si 0.1 0.2 0.3
---------------------------------
*hexa
Title
Reference
P 6
 11 13
Si 0.1 0.2 0.3
---------------------------------
*cubi
Title
Reference
P 2 3
 11
Si 0.1 0.2 0.3
O  0.0 0.0 0.0
---------------------------------
*mixcon1
Title
Reference
P 2 3
 11
Si 0.1 0.2 0.3 4
O  0.0 0.0 0.0
---------------------------------
*mixcon2
Title
Reference
P 2 3
 11
Si 0.1 0.2 0.3 4
O  0.0 0.0 0.0
O  0.0 0.0 0.0
---------------------------------
""")
  all_entries = strudat.read_all_entries(test_file)
  for tag,cell in (("tric", (11,12,13,100,110,120)),
                   ("mono_b", (11,12,13,90,100,90)),
                   ("mono_c", (11,12,13,90,90,100)),
                   ("mono_a", (11,12,13,100,90,90)),
                   ("orth", (11,12,13,90,90,90)),
                   ("tetr", (11,11,13,90,90,90)),
                   ("trig", (11,11,13,90,90,120)),
                   ("rhom", (11,11,11,100,100,100)),
                   ("hexa", (11,11,13,90,90,120)),
                   ("cubi", (11,11,11,90,90,90))):
    assert all_entries.get(tag).unit_cell().is_similar_to(
      uctbx.unit_cell(cell))
  assert all_entries.get("orth").atoms[0].connectivity is None
  assert all_entries.get("tetr").atoms[0].connectivity == 4
  assert all_entries.get("mixcon1").connectivities() == [4, None]
  assert all_entries.get("mixcon2").connectivities() == [4, None, None]
  try:
    all_entries.get("mixcon1").connectivities(all_or_nothing=True)
  except AssertionError as e:
    assert str(e) == "Tag mixcon1: 1 atom is missing the bond count."
  else:
    raise Exception_expected
  try:
    all_entries.get("mixcon2").connectivities(all_or_nothing=True)
  except AssertionError as e:
    assert str(e) == "Tag mixcon2: 2 atoms are missing the bond count."
  else:
    raise Exception_expected
  assert all_entries.get("cubi").as_xray_structure().scatterers().size() == 2

def exercise_zeolite_atlas(distance_cutoff=3.5):
  atlas_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/misc/strudat_zeolite_atlas",
    test=os.path.isfile)
  if (atlas_file is None):
    print("Skipping exercise_zeolite_atlas(): input file not available")
    return
  with open(atlas_file) as f:
    all_entries = strudat.read_all_entries(f)
  for i,entry in enumerate(all_entries.entries):
    structure = entry.as_xray_structure()
    if ("--full" in sys.argv[1:] or i % 20 == 0):
      tst_direct_space_asu.exercise_neighbors_pair_generators(
        structure=structure,
        verbose="--Verbose" in sys.argv[1:])
    asu_mappings = structure.asu_mappings(buffer_thickness=distance_cutoff)
    pair_generator = crystal.neighbors_fast_pair_generator(
      asu_mappings=asu_mappings,
      distance_cutoff=distance_cutoff)
    bond_counts = flex.size_t(structure.scatterers().size(), 0)
    for pair in pair_generator:
      bond_counts[pair.i_seq] += 1
      if (pair.j_sym == 0):
        bond_counts[pair.j_seq] += 1
    for atom,bond_count in zip(entry.atoms, bond_counts):
      assert atom.connectivity is not None
      assert atom.connectivity == bond_count

def run():
  exercise_basic()
  exercise_zeolite_atlas()
  print(format_cpu_times())

if (__name__ == "__main__"):
  run()
