from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from cctbx import miller
from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info
import boost_adaptbx.boost.python as bp
asu_map_ext = bp.import_ext("cctbx_asymmetric_map_ext")
from mmtbx.maps import composite_omit_map as cfom
from cctbx import maptbx
import mmtbx.f_model
import time, sys
from cctbx.development import debug_utils
from cctbx import sgtbx

def get_cc(mc1, mc2, xrs):
  crystal_gridding = mc1.crystal_gridding(
    d_min = mc1.d_min(), resolution_factor = 0.25)
  fft_map = miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = mc1)
  fft_map.apply_sigma_scaling()
  m1 = fft_map.real_map_unpadded()
  fft_map = miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = mc2)
  fft_map.apply_sigma_scaling()
  m2 = fft_map.real_map_unpadded()
  assert m1.focus()==m2.focus()
  assert m1.all()==m2.all()
  ccs = flex.double()
  for site_cart in xrs.sites_cart():
    sel = maptbx.grid_indices_around_sites(
      unit_cell  = mc1.unit_cell(),
      fft_n_real = m1.focus(),
      fft_m_real = m1.all(),
      sites_cart = flex.vec3_double([site_cart]),
      site_radii = flex.double([1.5]))
    cc = flex.linear_correlation(x=m1.select(sel), y=m2.select(sel)).coefficient()
    ccs.append(cc)
  return ccs

def run(space_group_info):
  """
  Make sure it work for all space groups and boxes with non-zero origin.
  """
  # make up data
  xrs = random_structure.xray_structure(
    space_group_info       = space_group_info,
    volume_per_atom        = 50,
    general_positions_only = False,
    u_iso                  = 0.3,
    elements               = ('C', 'N', 'O', "S")*10,
    min_distance           = 1.5)
  xrs.scattering_type_registry(table="wk1995")
  f_calc = xrs.structure_factors(d_min=2).f_calc()
  f_obs = abs(f_calc)
  # create fmodel object
  fmodel = mmtbx.f_model.manager(
    xray_structure = xrs,
    f_obs          = f_obs)
  fmodel.update_all_scales()
  mc1 = fmodel.electron_density_map().map_coefficients(
      map_type     = "2mFo-DFc",
      isotropize   = False,
      exclude_free_r_reflections=False,
      fill_missing = False)
  crystal_gridding = fmodel.f_obs().crystal_gridding(
    d_min             = fmodel.f_obs().d_min(),
    symmetry_flags    = maptbx.use_space_group_symmetry,
    resolution_factor = 1./3)
  # compute OMIT map
  r = cfom.run(
    crystal_gridding    = crystal_gridding,
    fmodel              = fmodel.deep_copy(),
    full_resolution_map = False,
    max_boxes           = 70,
    neutral_volume_box_cushion_width = 0,
    box_size_as_fraction=0.3,
    log=False)
  ccs = get_cc(mc1=mc1, mc2=r.map_coefficients(filter_noise=False), xrs=xrs)
  assert flex.mean(ccs) > 0.8
  print("  CC(min/max,mean)",ccs.min_max_mean().as_tuple())

def run_call_back(flags, space_group_info):
  run(space_group_info)

if (__name__ == "__main__"):
  t0 = time.time()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back,
    symbols_to_stdout=True, symbols_to_stderr=False)
  run(sgtbx.space_group_info("R3:R"))
  print("Time: %6.4f"%(time.time()-t0))
