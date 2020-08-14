from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from cctbx import miller
from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info
import boost_adaptbx.boost.python as bp
from six.moves import zip
from six.moves import range
asu_map_ext = bp.import_ext("cctbx_asymmetric_map_ext")
from mmtbx.maps import composite_omit_map as cfom
from cctbx import maptbx
import mmtbx.maps
from libtbx.test_utils import approx_equal
import mmtbx.f_model
import time
import sys
import cctbx

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
        max_boxes,
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
      fraction = 0.3,
      n_real   = f_model_map_data.focus(),
      max_boxes=max_boxes,
      log=sys.stdout)
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

class omit_general_obsolete(object):
  """
  Composite full-omit map: omit entire box in real-space corresponding to
  Fmodel, which includs atomic model and non-atomic model (bulk-solvent and
  scales).  This is much faster than the refinement-based method, and is used
  by default in phenix.composite_omit_map.
  """
  def __init__(
        self,
        crystal_gridding,
        fmodel,
        map_type,
        box_size_as_fraction=0.03,
        max_boxes=100,
        n_debias_cycles=2,
        neutral_volume_box_cushion_width=1,
        full_resolution_map=True,
        log=sys.stdout):
    self.crystal_gridding = crystal_gridding
    # assert compatibility of symops with griding
    assert self.crystal_gridding._symmetry_flags is not None
    self.sgt = fmodel.f_obs().space_group().type()
    self.zero_cmpl_ma = fmodel.f_calc().customized_copy(
      data = flex.complex_double(fmodel.f_calc().size(), 0))
    # embedded utility functions
    def get_map(fmodel, map_type, crystal_gridding, asu=True):
      f_map = fmodel.electron_density_map().map_coefficients(
        map_type                   = map_type,
        isotropize                 = True,
        exclude_free_r_reflections = True,
        fill_missing               = False)
      fft_map = cctbx.miller.fft_map(
        crystal_gridding     = crystal_gridding,
        fourier_coefficients = f_map)
      if(asu): return asu_map_ext.asymmetric_map(self.sgt,
        fft_map.real_map_unpadded()).data()
      else:
        return fft_map.real_map_unpadded()
    # f_model map
    f_model_map_data = fmodel.f_model_scaled_with_k1().fft_map(
      symmetry_flags   = maptbx.use_space_group_symmetry,
      crystal_gridding = self.crystal_gridding).real_map_unpadded()
    self.n_real = f_model_map_data.focus()
    # extract asu map from full P1
    f_model_map_data_asu=asu_map_ext.asymmetric_map(
      self.sgt, f_model_map_data).data()
    self.acc = f_model_map_data_asu.accessor()
    f_model_map_data_asu = f_model_map_data_asu.shift_origin()
    # set up boxes
    b = maptbx.boxes(
      n_real   = f_model_map_data_asu.focus(),
      fraction = box_size_as_fraction,
      max_boxes= max_boxes,
      log      = log)
    self.map_result_asu = flex.double(flex.grid(b.n_real))
    assert f_model_map_data_asu.focus()==b.n_real
    assert b.n_real==self.map_result_asu.focus()
    n_real_asu = b.n_real
    self.r = flex.double() # for regression test only
    n_boxes = len(b.starts)
    i_box = 0
    for s,e in zip(b.starts, b.ends):
      i_box+=1
      # define wide box: neutral + phased volumes
      if(neutral_volume_box_cushion_width>0):
        sh = neutral_volume_box_cushion_width
        ss = [max(s[i]-sh,0) for i in [0,1,2]]
        ee = [min(e[i]+sh,n_real_asu[i]) for i in [0,1,2]]
      else: ss,ee = s,e
      # omit wide box from f_model map, repeat n_debias_cycles times
      f_model_map_data_asu_ = f_model_map_data_asu.deep_copy()
      for i in range(n_debias_cycles):
        f_model_omit, f_model_map_data_asu_ = self.omit_box(s=ss, e=ee,
          md_asu=f_model_map_data_asu_)
      # get fmodel for omit map calculation
      fmodel_ = mmtbx.f_model.manager(
        f_obs        = fmodel.f_obs(),
        r_free_flags = fmodel.r_free_flags(),
        f_calc       = f_model_omit,
        f_mask       = self.zero_cmpl_ma)
      rw = fmodel_.r_work()
      self.r.append(rw) # for regression test only
      f_map_data_asu = get_map(fmodel=fmodel_, map_type=map_type,
        crystal_gridding=self.crystal_gridding)
      f_map_data_asu = f_map_data_asu.shift_origin()
      if(log):
        print("box %2d of %2d:"%(i_box, n_boxes), s, e, "%6.4f"%rw, file=log)
      assert f_map_data_asu.focus() == self.map_result_asu.focus()
      maptbx.copy_box(
        map_data_from = f_map_data_asu,
        map_data_to   = self.map_result_asu,
        start         = s,
        end           = e)
    # result
    self.map_result_asu.reshape(self.acc)
    self.asu_map_omit = asu_map_ext.asymmetric_map(
      self.sgt, self.map_result_asu, self.n_real)
    self.map_coefficients = self.zero_cmpl_ma.customized_copy(
      indices = self.zero_cmpl_ma.indices(),
      data    = self.asu_map_omit.structure_factors(self.zero_cmpl_ma.indices()))
    # full resolution map (reflections in sphere, not in box!)
    if(full_resolution_map):
      cs = self.zero_cmpl_ma.complete_set(d_min=self.zero_cmpl_ma.d_min())
      asu_map_omit = asu_map_ext.asymmetric_map(
        self.sgt,self.map_result_asu,self.n_real)
      fill = self.zero_cmpl_ma.customized_copy(
        indices = cs.indices(),
        data    = asu_map_omit.structure_factors(cs.indices()))
      self.map_coefficients = self.map_coefficients.complete_with(
        other=fill, scale=True)

  def omit_box(self, s, e, md_asu):
    md_asu_omit = maptbx.set_box_copy(value = 0, map_data_to = md_asu,
      start = s, end = e)
    md_asu_omit.reshape(self.acc)
    asu_map_omit = asu_map_ext.asymmetric_map(self.sgt, md_asu_omit, self.n_real)
    ma_omit = self.zero_cmpl_ma.customized_copy(
      indices = self.zero_cmpl_ma.indices(),
      data    = asu_map_omit.structure_factors(self.zero_cmpl_ma.indices()))
    fft_map = cctbx.miller.fft_map(
      crystal_gridding = self.crystal_gridding, fourier_coefficients = ma_omit)
    md = fft_map.real_map_unpadded()
    asu_map = asu_map_ext.asymmetric_map(self.sgt, md)
    md_asu_omit = asu_map.data()
    md_asu_omit = md_asu_omit.shift_origin()
    return ma_omit, md_asu_omit

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
    u_iso                  = 0.1,
    min_distance           = 1.0)
  xrs.scattering_type_registry(table="wk1995")
  f_obs = abs(xrs.structure_factors(d_min=2).f_calc())
  # create fmodel object
  xrs.shake_sites_in_place(mean_distance=0.3)
  sel = xrs.random_remove_sites_selection(fraction=0.1)
  xrs = xrs.select(sel)
  fmodel = mmtbx.f_model.manager(
    xray_structure = xrs,
    f_obs          = f_obs)
  fmodel.update_all_scales(update_f_part1=False)
  crystal_gridding = fmodel.f_obs().crystal_gridding(
    d_min             = fmodel.f_obs().d_min(),
    symmetry_flags    = maptbx.use_space_group_symmetry,
    resolution_factor = 0.25)
  # compute OMIT maps
  r1 = omit_p1_specific(
    crystal_gridding = crystal_gridding,
    fmodel           = fmodel.deep_copy(),
    max_boxes=70,
    map_type         = "Fo")
  r2 = omit_general_obsolete(
    crystal_gridding = crystal_gridding,
    fmodel           = fmodel.deep_copy(),
    full_resolution_map = False,
    map_type         = "Fo",
    n_debias_cycles  = 1,
    neutral_volume_box_cushion_width = 0,
    box_size_as_fraction=0.3,
    max_boxes=70,
    log=sys.stdout)
  r3 = cfom.run(
    crystal_gridding = crystal_gridding,
    fmodel           = fmodel.deep_copy(),
    full_resolution_map = False,
    neutral_volume_box_cushion_width = 0,
    box_size_as_fraction=0.3,
    max_boxes=70,
    log=sys.stdout)
  assert approx_equal(r1.r, r2.r)
  def r_factor(x,y):
    x = flex.abs(abs(x).data())
    y = flex.abs(abs(y).data())
    sc = flex.sum(x*y)/flex.sum(y*y)
    return flex.sum(flex.abs(x-sc*y))/flex.sum(x+sc*y)*2
  print(abs(r1.map_coefficients).data().min_max_mean().as_tuple())
  print(abs(r2.map_coefficients).data().min_max_mean().as_tuple())
  cc1=flex.linear_correlation(
      x=abs(r1.map_coefficients).data(),
      y=abs(r2.map_coefficients).data()).coefficient()
  assert approx_equal(cc1, 1.0)
  cc2=flex.linear_correlation(
      x=abs(r1.map_coefficients).data(),
      y=abs(r3.map_coefficients(filter_noise=False)).data()).coefficient()
  assert cc2 > 0.8, cc2
  assert approx_equal(r_factor(
    x=r1.map_coefficients, y=r2.map_coefficients), 0.0)
  cc3=flex.linear_correlation(x=r1.r, y=r3.r).coefficient()
  assert cc3>0.95
  for cutoff in [0.5,0.6,0.7,0.8,0.9,0.95,0.99]:
    print(maptbx.cc_peak(
      cutoff       = cutoff,
      map_coeffs_1 = r1.map_coefficients,
      map_coeffs_2 = r3.map_coefficients(filter_noise=False)), "CCpeak", cutoff)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.4f"%(time.time()-t0))
