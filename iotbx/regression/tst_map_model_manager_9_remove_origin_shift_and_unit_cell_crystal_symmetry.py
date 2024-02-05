from __future__ import absolute_import, division, print_function
from iotbx.data_manager import DataManager
import os


def exercise():
  data_dir = os.path.dirname(os.path.abspath(__file__))
  data_ccp4 = os.path.join(data_dir, 'data',
                          'non_zero_origin_map.ccp4')
  data_model = os.path.join(data_dir, 'data',
                          'non_zero_origin_model.pdb')

  dm = DataManager()
  mmm = dm.get_map_model_manager(
    model_file = data_model,
    map_files  = data_ccp4)

  box_mmm = mmm.extract_all_maps_around_model()

  box_cs = box_mmm.crystal_symmetry()

  assert not box_mmm.crystal_symmetry().is_similar_symmetry(
     mmm.crystal_symmetry())
  assert box_mmm.unit_cell_crystal_symmetry().is_similar_symmetry(
     mmm.unit_cell_crystal_symmetry())

  box_mmm.remove_origin_shift_and_unit_cell_crystal_symmetry()
  assert box_mmm.crystal_symmetry().is_similar_symmetry(box_cs)
  assert box_mmm.unit_cell_crystal_symmetry().is_similar_symmetry(box_cs)
  assert not box_mmm.unit_cell_crystal_symmetry().is_similar_symmetry(
     mmm.crystal_symmetry())
  file_name = box_mmm.write_model("model.cif", format ='cif')
  box_mmm.write_map("map.ccp4")

  new_mmm = dm.get_map_model_manager(
    model_file = file_name,
    map_files  = "map.ccp4")

  assert new_mmm.crystal_symmetry().is_similar_symmetry(box_cs)
  assert new_mmm.unit_cell_crystal_symmetry().is_similar_symmetry(box_cs)

if (__name__ == "__main__"):
  import time
  t0 = time.time()
  exercise()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
