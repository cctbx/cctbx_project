from __future__ import division
from scitbx.array_family import flex
from cctbx import miller
from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info
import boost.python
asu_map_ext = boost.python.import_ext("cctbx_asymmetric_map_ext")
from mmtbx.maps import composite_omit_map as cfom
from cctbx import maptbx
import mmtbx.maps
from libtbx.test_utils import approx_equal
import mmtbx.f_model
import time

class omit_p1_specific(object):
  """
  Composite full-omit map: omit entire box in real-space corresponding to Fmodel,
  which includs atomic model and non-atomic model (bulk-solvent and scales).
  ASU maps are not use, so this is limited to P1 only. Use for this test only.
  """
  def __init__(
        self,
        crystal_gridding,
        fmodel,
        map_type,
        box_size_as_fraction=None):
    sgt = fmodel.f_obs().space_group().type()
    assert sgt.number() == 1
    def get_map(fmodel, map_type, external_complete_set=None):
      f_map = fmodel.electron_density_map().map_coefficients(
          map_type     = map_type,
          isotropize   = True,
          fill_missing = False)
      fft_map = miller.fft_map(
        crystal_gridding     = crystal_gridding,
        fourier_coefficients = f_map)
      return fft_map.real_map_unpadded()
    f_model = fmodel.f_model_scaled_with_k1()
    fft_map = miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = f_model)
    f_model_map_data = fft_map.real_map_unpadded()
    zero_complex_ma = f_model.customized_copy(
      data = flex.complex_double(f_model.data().size(), 0))
    b = maptbx.boxes(
      n_real   = f_model_map_data.focus(),
      fraction = box_size_as_fraction)
    self.map_result = flex.double(flex.grid(b.n_real))
    self.r = flex.double()
    for s,e in zip(b.starts, b.ends):
      f_model_map_data_omit = maptbx.set_box_copy(
        value         = 0,
        map_data_to   = f_model_map_data,
        start         = s,
        end           = e)
      f_model_omit = f_model.structure_factors_from_map(
        map            = f_model_map_data_omit,
        use_scale      = True,
        anomalous_flag = False,
        use_sg         = False)
      fmodel_ = mmtbx.f_model.manager(
        f_obs        = fmodel.f_obs(),
        r_free_flags = fmodel.r_free_flags(),
        f_calc       = f_model_omit,
        f_mask       = zero_complex_ma)
      self.r.append(fmodel_.r_work())
      f_map_data = get_map(fmodel=fmodel_, map_type=map_type)
      etmp = [e[0]-1, e[1]-1, e[2]-1] # because .copy() includes right edge
      box = maptbx.copy(f_map_data, s, etmp)
      box.reshape(flex.grid(box.all()))
      maptbx.set_box(
        map_data_from = box,
        map_data_to   = self.map_result,
        start         = s,
        end           = e)
    sd = self.map_result.sample_standard_deviation()
    self.map_result = self.map_result/sd
    self.map_coefficients = fmodel.f_obs().structure_factors_from_map(
      map            = self.map_result,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)

def run():
  """
  This test makes sure that composite full omit maps calculated using Marat's
  ASU map code and not using ASU maps exactly match.
  """
  # make up data
  xrs = random_structure.xray_structure(
    space_group_info       = space_group_info("P 1"),
    volume_per_atom        = 250,
    general_positions_only = False,
    elements               = ('C', 'N', 'O', "S")*50,
    min_distance           = 1.0)
  xrs.scattering_type_registry(table="wk1995")
  f_obs = abs(xrs.structure_factors(d_min=2).f_calc())
  # create fmodel object
  xrs.shake_sites_in_place(mean_distance=0.3)
  sel = xrs.random_remove_sites_selection(fraction=0.1)
  xrs = xrs.select(sel)
  fmodel = mmtbx.f_model.manager(
    xray_structure = xrs,
    f_obs          = f_obs,
    r_free_flags   = f_obs.generate_r_free_flags())
  fmodel.update_all_scales(update_f_part1=False)
  crystal_gridding = fmodel.f_obs().crystal_gridding(
    d_min             = fmodel.f_obs().d_min(),
    symmetry_flags    = maptbx.use_space_group_symmetry,
    resolution_factor = 0.25)
  # compute OMIT maps
  r1 = omit_p1_specific(
    crystal_gridding = crystal_gridding,
    fmodel           = fmodel.deep_copy(),
    map_type         = "Fo",
    box_size_as_fraction=0.05)
  r2 = cfom.run(
    crystal_gridding = crystal_gridding,
    fmodel           = fmodel.deep_copy(),
    full_resolution_map = False,
    map_type         = "Fo",
    n_debias_cycles  = 1,
    neutral_volume_box_cushion_width = 0,
    box_size_as_fraction=0.05)
  assert approx_equal(r1.r, r2.r)
  def r_factor(x,y):
    x = flex.abs(abs(x).data())
    y = flex.abs(abs(y).data())
    sc = flex.sum(x*y)/flex.sum(y*y)
    return flex.sum(flex.abs(x-sc*y))/flex.sum(x+sc*y)*2
  print abs(r1.map_coefficients).data().min_max_mean().as_tuple()
  print abs(r2.map_coefficients).data().min_max_mean().as_tuple()
  assert approx_equal(flex.linear_correlation(
      x=abs(r1.map_coefficients).data(),
      y=abs(r2.map_coefficients).data()).coefficient(), 1.0)
  assert approx_equal(r_factor(
    x=r1.map_coefficients, y=r2.map_coefficients), 0.0)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print "Time: %6.4f"%(time.time()-t0)
