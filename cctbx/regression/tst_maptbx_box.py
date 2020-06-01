from __future__ import absolute_import, division, print_function
import iotbx.pdb
from libtbx.test_utils import approx_equal
from cctbx.sgtbx import space_group_info
from cctbx.development import random_structure
import cctbx.maptbx.box
from libtbx import group_args
import iotbx.pdb
from iotbx.map_manager import map_manager
import mmtbx.model

def get_random_structure_and_map():
  xrs = random_structure.xray_structure(
    space_group_info = space_group_info(19),
    volume_per_atom  = 25.,
    elements         = ('C', 'N', 'O', 'H')*10,
    min_distance     = 1.5)
  fc = xrs.structure_factors(d_min=2).f_calc()
  fft_map = fc.fft_map(resolution_factor=0.25)
  fft_map.apply_volume_scaling()
  ph = iotbx.pdb.input(
    source_info=None, lines=xrs.as_pdb_file()).construct_hierarchy()
  ph.atoms().set_xyz(xrs.sites_cart())
  map_data = fft_map.real_map_unpadded()
  mm = map_manager(
    unit_cell_grid             = map_data.accessor().all(),
    unit_cell_crystal_symmetry = fc.crystal_symmetry(),
    origin_shift_grid_units    = (0,0,0),
    map_data                   = map_data)
  model = mmtbx.model.manager(
    model_input=None, pdb_hierarchy=ph, crystal_symmetry=fc.crystal_symmetry())
  return group_args(model = model, mm = mm)

def exercise_around_model():
  mam = get_random_structure_and_map()
  map_data_orig   = mam.mm.map_data().deep_copy()
  sites_frac_orig = mam.model.get_sites_frac().deep_copy()
  sites_cart_orig = mam.model.get_sites_cart().deep_copy()
  cs_orig         = mam.model.crystal_symmetry()
  box = cctbx.maptbx.box.around_model(
    map_manager = mam.mm,
    model       = mam.model,
    cushion     = 5)
  new_mm1 = box.apply_to_map()
  new_mm2 = box.apply_to_map(map_manager=mam.mm)
  assert approx_equal(new_mm1.map_data(), new_mm2.map_data())
  new_model1 = box.apply_to_model()
  new_model2 = box.apply_to_model(model=mam.model)
  assert new_model1.crystal_symmetry().is_similar_symmetry(
         new_model2.crystal_symmetry())
  # make sure things did change
  assert new_mm2.map_data().size() != map_data_orig.size()
  # make sure things are not changed in-place
  assert approx_equal(box.map_manager.map_data(), map_data_orig)
  assert approx_equal(box.model.get_sites_frac(), sites_frac_orig)
  assert approx_equal(box.model.get_sites_cart(), sites_cart_orig)
  assert cs_orig.is_similar_symmetry(box.model.crystal_symmetry())
  assert cs_orig.is_similar_symmetry(box.map_manager.crystal_symmetry())
  #
  # IF you are about to change this - THINK TWICE!
  #
  import inspect
  r = inspect.getargspec(cctbx.maptbx.box.around_model.__init__)
  assert r.args == ['self', 'map_manager', 'model', 'cushion'], r.args

if (__name__ == "__main__"):
  exercise_around_model()
