from __future__ import absolute_import, division, print_function
from iotbx.data_manager import DataManager
import libtbx.load_env
import os
from libtbx import easy_run

def exercise():
  if not libtbx.env.has_module("phenix_regression"):
    print("phenix_regression not configured, skipping.")
    return
  fn="phenix_regression/mmtbx/extract_box_around_model_and_map/tst8.pdb"
  pdb_file_name = libtbx.env.find_in_repositories(
    relative_path=fn,
    test=os.path.isfile)
  prefix = "extract_box_around_model_and_map_tst8"
  cmd = "phenix.model_map %s output_file_name_prefix=%s grid_step=1"
  easy_run.call(cmd%(pdb_file_name,prefix))
  #
  dm = DataManager()
  mmm = dm.get_map_model_manager(
    model_file = pdb_file_name,
    map_files  = "%s.ccp4"%prefix)

  box_mmm = mmm.extract_all_maps_around_model()

  box_cs = box_mmm.crystal_symmetry()

  assert not box_mmm.crystal_symmetry().is_similar_symmetry(
     mmm.crystal_symmetry())
  assert box_mmm.unit_cell_crystal_symmetry().is_similar_symmetry(
     mmm.crystal_symmetry())

  box_mmm.remove_origin_shift_and_unit_cell_crystal_symmetry()
  assert box_mmm.crystal_symmetry().is_similar_symmetry(box_cs)
  assert box_mmm.unit_cell_crystal_symmetry().is_similar_symmetry(box_cs)
  assert not box_mmm.unit_cell_crystal_symmetry().is_similar_symmetry(
     mmm.crystal_symmetry())
  box_mmm.write_model("model.pdb")
  box_mmm.write_map("map.ccp4")

  new_mmm = dm.get_map_model_manager(
    model_file = "model.pdb",
    map_files  = "map.ccp4")

  assert new_mmm.crystal_symmetry().is_similar_symmetry(box_cs)
  assert new_mmm.unit_cell_crystal_symmetry().is_similar_symmetry(box_cs)

if (__name__ == "__main__"):
  import time
  t0 = time.time()
  exercise()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
