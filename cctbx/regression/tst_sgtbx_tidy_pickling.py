from __future__ import absolute_import, division, print_function

from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex
try:
  from six.moves import cPickle as pickle
except ImportError:
  import pickle

"""
Test that the spacegroup class retains order of symmetry operations when restored
from pickling. This is done in cctbx_project/cctbx/sgtbx/boost_python/space_group.cpp
by calling make_tidy() in the space_group_wrappers::setstate() function.
"""

def test_spacegroup_tidy_pickling():
  quartz_structure = xray.structure(
    crystal_symmetry=crystal.symmetry(
    unit_cell=(5.01,5.01,5.47,90,90,120),
    space_group_symbol="P6222"),
    scatterers=flex.xray_scatterer(
    [
      xray.scatterer(label="Si", site=(1/2.,1/2.,1/3.), u=0.2),
      xray.scatterer(label="O",  site=(0.197,-0.197,0.83333), u=0)
    ])
  )

  asu_mappings = quartz_structure.asu_mappings(buffer_thickness=2)
  pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  pair_asu_table.add_all_pairs(distance_cutoff=1.7)
  pair_sym_table = pair_asu_table.extract_pair_sym_table()
  new_asu_mappings = quartz_structure.asu_mappings(buffer_thickness=5)
  new_pair_asu_table = crystal.pair_asu_table(asu_mappings=new_asu_mappings)
  new_pair_asu_table.add_pair_sym_table(sym_table=pair_sym_table)
  spg = new_pair_asu_table.asu_mappings().space_group()
  pspg = pickle.loads(pickle.dumps(spg))
  mstr = ""
  pmstr = ""
  for rt in spg.all_ops():
    mstr += rt.r().as_xyz() + "\n"
  for rt in pspg.all_ops():
    pmstr += rt.r().as_xyz() + "\n"
  assert mstr == pmstr


if (__name__ == "__main__"):
  test_spacegroup_tidy_pickling()
  print("OK")
