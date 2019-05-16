from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import xray
from cctbx import crystal
from cctbx import maptbx
from cctbx.maptbx import minimization
from libtbx.test_utils import approx_equal
import random
from cctbx.development import random_structure
from cctbx import sgtbx
from six.moves import zip

if (1):
  random.seed(0)
  flex.set_random_seed(0)

def get_xrs():
  crystal_symmetry = crystal.symmetry(
    unit_cell=(10,10,10,90,90,90),
    space_group_symbol="P 1")
  return xray.structure(
    crystal_symmetry=crystal_symmetry,
    scatterers=flex.xray_scatterer([
      xray.scatterer(label="C", site=(0,0,0))]))

def get_map(xrs, d_min=1.):
  f_calc = xrs.structure_factors(d_min=d_min).f_calc()
  fft_map = f_calc.fft_map()
  fft_map.apply_sigma_scaling()
  return fft_map.real_map_unpadded(), f_calc

def exercise_00():
  """
  Exercise maptbx.target_and_gradients_diffmap .
  """
  xrs = get_xrs()
  map_data, f_calc = get_map(xrs=xrs)
  tg = maptbx.target_and_gradients_diffmap(
    unit_cell   = xrs.unit_cell(),
    map_target  = map_data,
    map_current = map_data,
    step        = 0.3,
    sites_frac  = xrs.sites_frac())
  assert approx_equal(xrs.sites_cart(), [[0,0,0]])
  assert approx_equal(tg.target(), 0)
  assert approx_equal(list(tg.gradients()), [[0,0,0]])
  xrs = xrs.translate(x=0.3, y=-0.5, z=0.7)
  assert approx_equal(xrs.sites_cart(), [[0.3,-0.5,0.7]])
  map_current, f_calc = get_map(xrs=xrs)
  tg = maptbx.target_and_gradients_diffmap(
    unit_cell   = xrs.unit_cell(),
    map_target  = map_data,
    map_current = map_current,
    step        = 0.3,
    sites_frac  = xrs.sites_frac())
  assert tg.target() > 0
  for g in tg.gradients():
    for g_ in g:
      assert abs(g_)>0.

def exercise_01(d_min=1.0):
  """
  Exercise maptbx.target_and_gradients_diffmap in action: minimization.
  """
  xrs = get_xrs()
  map_target, f_calc = get_map(xrs=xrs)
  assert approx_equal(xrs.sites_cart(), [[0,0,0]])
  for sx in [-1,0,1]:
    for sy in [-1,0,1]:
      for sz in [-1,0,1]:
        xrs_cp = xrs.deep_copy_scatterers()
        xrs_cp = xrs_cp.translate(x=0.3*sx, y=0.5*sy, z=0.7*sz)
        assert approx_equal(xrs_cp.sites_cart(), [[0.3*sx,0.5*sy,0.7*sz]],1.e-6)
        crystal_gridding = maptbx.crystal_gridding(
          unit_cell             = xrs_cp.unit_cell(),
          space_group_info      = xrs_cp.space_group_info(),
          pre_determined_n_real = map_target.accessor().all())
        o = minimization.run(
          xray_structure   = xrs_cp,
          miller_array     = f_calc,
          crystal_gridding = crystal_gridding,
          map_target       = map_target,
          step             = d_min/4,
          target_type      = "diffmap")
        assert approx_equal(xrs.sites_cart(), [[0,0,0]])


def exercise_02():
  """
  Exercise maptbx.target_and_gradients_diffmap in action: minimization
  (bigger model).
  """
  def compute_map(xray_structure, d_min=1.5, resolution_factor=1./4):
    fc = xray_structure.structure_factors(d_min = d_min).f_calc()
    fft_map = fc.fft_map(resolution_factor=resolution_factor)
    fft_map.apply_sigma_scaling()
    result = fft_map.real_map_unpadded()
    return result, fc, fft_map
  xrs = random_structure.xray_structure(
    space_group_info  = sgtbx.space_group_info("P212121"),
    elements          = ["N","C","O","S","P"]*10,
    volume_per_atom   = 50)
  map_target,tmp,tmp = compute_map(xray_structure = xrs)
  xrs_sh = xrs.deep_copy_scatterers()
  xrs_sh.shake_sites_in_place(mean_distance=0.8)
  start_error = flex.mean(xrs.distances(other = xrs_sh))
  assert start_error>0.7
  map_current, miller_array, crystal_gridding = compute_map(
    xray_structure = xrs_sh)
  for step in [miller_array.d_min()/4]*5:
    minimized = minimization.run(
      xray_structure              = xrs_sh,
      miller_array                = miller_array,
      crystal_gridding            = crystal_gridding,
      map_target                  = map_target,
      max_iterations              = 500,
      min_iterations              = 25,
      step                        = step,
      geometry_restraints_manager = None,
      target_type                 = "diffmap")
    xrs_sh = minimized.xray_structure
    map_current = minimized.map_current
    final_error = flex.mean(xrs.distances(other = minimized.xray_structure))
  assert approx_equal(start_error, 0.8, 1.e-3)
  assert final_error < 1.e-4

def exercise_03():
  """
  Exercise maptbx.target_and_gradients_simple.
  """
  def compute_map(xray_structure, d_min=1.5, resolution_factor=1./4):
    fc = xray_structure.structure_factors(d_min = d_min).f_calc()
    fft_map = fc.fft_map(resolution_factor=resolution_factor)
    fft_map.apply_sigma_scaling()
    result = fft_map.real_map_unpadded()
    return result, fc, fft_map
  xrs = random_structure.xray_structure(
    space_group_info  = sgtbx.space_group_info("P212121"),
    elements          = ["N","C","O","S","P"]*10,
    volume_per_atom   = 50)
  map_target,tmp,tmp = compute_map(xray_structure = xrs)
  xrs_sh = xrs.deep_copy_scatterers()
  xrs_sh.shake_sites_in_place(mean_distance=0.8)
  #
  t1 = maptbx.real_space_target_simple(
    unit_cell   = xrs.unit_cell(),
    density_map = map_target,
    sites_cart  = xrs_sh.sites_cart(),
    selection   = flex.bool(xrs_sh.scatterers().size(), True))
  g1 = maptbx.real_space_gradients_simple(
    unit_cell   = xrs.unit_cell(),
    density_map = map_target,
    sites_cart  = xrs_sh.sites_cart(),
    delta       = 0.25,
    selection   = flex.bool(xrs_sh.scatterers().size(), True))
  o = maptbx.target_and_gradients_simple(
    unit_cell   = xrs.unit_cell(),
    map_target  = map_target,
    sites_cart  = xrs_sh.sites_cart(),
    delta       = 0.25,
    selection   = flex.bool(xrs_sh.scatterers().size(), True))
  assert approx_equal(t1, o.target())
  for gi,gj in zip(g1, o.gradients()):
    assert approx_equal(gi, gj)

def exercise_04():
  """
  Exercise maptbx.target_and_gradients_simple in action: minimization
  (bigger model).
  """
  def compute_map(xray_structure, d_min=1., resolution_factor=1./4):
    fc = xray_structure.structure_factors(d_min = d_min).f_calc()
    fft_map = fc.fft_map(resolution_factor=resolution_factor)
    fft_map.apply_sigma_scaling()
    result = fft_map.real_map_unpadded()
    return result, fc, fft_map
  xrs = random_structure.xray_structure(
    space_group_info  = sgtbx.space_group_info("P212121"),
    elements          = ["N","C","O","S","P"]*10,
    volume_per_atom   = 150)
  map_target,tmp,tmp = compute_map(xray_structure = xrs)
  xrs_sh = xrs.deep_copy_scatterers()
  xrs_sh.shake_sites_in_place(mean_distance=0.3)
  start_error = flex.mean(xrs.distances(other = xrs_sh))
  assert start_error > 0.29
  map_current, miller_array, crystal_gridding = compute_map(
    xray_structure = xrs_sh)
  xrs_sh_ = xrs_sh.deep_copy_scatterers()
  minimized = minimization.run(
    xray_structure              = xrs_sh_,
    miller_array                = miller_array,
    crystal_gridding            = crystal_gridding,
    map_target                  = map_target,
    max_iterations              = 500,
    min_iterations              = 25,
    step                        = 0.5,
    geometry_restraints_manager = None,
    target_type                 = "simple")
  xrs_sh_ = xrs_sh_.replace_sites_cart(minimized.sites_cart)
  final_error = flex.mean(xrs.distances(other = xrs_sh_))
  assert final_error < 0.015

if (__name__ == "__main__"):
  exercise_00()
  exercise_01()
  exercise_02()
  exercise_03()
  exercise_04()
