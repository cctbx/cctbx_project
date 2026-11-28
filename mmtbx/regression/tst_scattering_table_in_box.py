from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import mmtbx.utils
from mmtbx.secondary_structure.build import ss_idealization as ssb
from cctbx import uctbx
from iotbx import map_manager, map_model_manager
from scitbx.array_family import flex
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal

"""
A model cut out in a box together with its map keeps its scattering table
(electron, neutron, ...): structure factors calculated from the model in the
box are those of the table the model came with, not of the default X-ray one.
"""

def get_inputs(table):
  ph = ssb.secondary_structure_from_sequence(
    ssb.alpha_helix_str, "ACDEFGHIKL")
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart = ph.atoms().extract_xyz(), buffer_layer = 8)
  ph.atoms().set_xyz(box.sites_cart)
  model = mmtbx.model.manager(
    model_input      = None,
    pdb_hierarchy    = ph,
    crystal_symmetry = box.crystal_symmetry(),
    log              = null_out())
  model.setup_scattering_dictionaries(scattering_table = table, d_min = 2.)
  model.set_b_iso(values = flex.double(model.size(), 40))
  fc = model.get_xray_structure().structure_factors(d_min = 2.).f_calc()
  map_data = fc.fft_map(resolution_factor = 0.25).real_map_unpadded()
  mm = map_manager.map_manager(
    map_data                   = map_data,
    unit_cell_grid             = map_data.all(),
    unit_cell_crystal_symmetry = model.crystal_symmetry(),
    wrapping                   = False)
  return model, mm

def table_of(xrs):
  return xrs.scattering_type_registry().last_table()

def exercise_extract_box(table):
  """
  mmtbx.utils.extract_box_around_model_and_map: the structure in the box.
  """
  model, mm = get_inputs(table = table)
  xrs = model.get_xray_structure()
  assert table_of(xrs) == table
  box = mmtbx.utils.extract_box_around_model_and_map(
    xray_structure = xrs,
    map_data       = mm.map_data().deep_copy(),
    box_cushion    = 3)
  xrs_box = box.xray_structure_box
  assert table_of(xrs_box) == table, table_of(xrs_box)
  # same form factors: unique gaussians of each scattering type
  g1 = xrs.scattering_type_registry().as_type_gaussian_dict()
  g2 = xrs_box.scattering_type_registry().as_type_gaussian_dict()
  assert sorted(g1.keys()) == sorted(g2.keys())
  for k in g1.keys():
    assert approx_equal(g1[k].parameters(), g2[k].parameters())
  # map in the box and the model in the box agree: amplitudes from the box
  # map against those calculated from the structure in the box
  f_obs = abs(box.box_map_coefficients(d_min = 2.5))
  f_calc = abs(f_obs.structure_factors_from_scatterers(
    xray_structure = xrs_box).f_calc())
  k = flex.sum(f_obs.data()*f_calc.data())/flex.sum(f_calc.data()*f_calc.data())
  r = flex.sum(flex.abs(f_obs.data()-k*f_calc.data()))/flex.sum(f_obs.data())
  print("%-10s extract_box: R(box map, model in box) = %6.4f"%(table, r))
  assert r < 0.05, r
  # the input structure is not modified, and the two do not share scatterers
  assert table_of(xrs) == table
  xrs_box.set_b_iso(value = 10)
  assert approx_equal(
    flex.mean(xrs.extract_u_iso_or_u_equiv())*79, 40, 0.1)

def exercise_map_model_manager(table):
  """
  iotbx.map_model_manager: model after boxing, and its deep copy.
  """
  model, mm = get_inputs(table = table)
  mmm = map_model_manager.map_model_manager(model = model, map_manager = mm)
  assert table_of(mmm.model().get_xray_structure()) == table
  mmm.box_all_maps_around_model_and_shift_origin(box_cushion = 3)
  assert table_of(mmm.model().get_xray_structure()) == table
  boxed = mmm.extract_all_maps_around_model(box_cushion = 2)
  assert table_of(boxed.model().get_xray_structure()) == table
  assert table_of(mmm.deep_copy().model().get_xray_structure()) == table
  sel = flex.bool([i%2 == 0 for i in range(mmm.model().size())])
  assert table_of(mmm.model().select(sel).get_xray_structure()) == table

if(__name__ == "__main__"):
  t0 = time.time()
  for table in ["electron", "neutron", "wk1995", "n_gaussian"]:
    exercise_map_model_manager(table = table)
    exercise_extract_box(table = table)
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
