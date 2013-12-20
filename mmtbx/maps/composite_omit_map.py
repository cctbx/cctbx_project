from __future__ import division
from scitbx.array_family import flex
from cctbx import miller
from cctbx import maptbx
import boost.python
asu_map_ext = boost.python.import_ext("cctbx_asymmetric_map_ext")
import mmtbx.f_model
import cctbx.miller
import sys

class run(object):
  """
  Composite full-omit map: omit entire box in real-space corresponding to Fmodel,
  which includs atomic model and non-atomic model (bulk-solvent and scales).
  """
  def __init__(
        self,
        crystal_gridding,
        fmodel,
        map_type,
        box_size_as_fraction=0.03,
        n_debias_cycles=2,
        resolution_factor=0.25,
        neutral_colume_box_cushion_width=1,
        full_resolution_map=False,
        log=sys.stdout):
    sgt = fmodel.f_obs().space_group().type()
    def get_map(fmodel, map_type, fft_map_ref):
      f_map = fmodel.electron_density_map(
        update_f_part1=False).map_coefficients(
          map_type     = map_type,
          isotropize   = True,
          fill_missing = False)
      fft_map = cctbx.miller.fft_map(
        crystal_gridding     = fft_map_ref,
        fourier_coefficients = f_map)
      asu_map_fmd = asu_map_ext.asymmetric_map(sgt, fft_map.real_map_unpadded())
      return asu_map_fmd.data()
    f_model = fmodel.f_model_scaled_with_k1()
    fft_map_ref = f_model.fft_map(
      symmetry_flags    = maptbx.use_space_group_symmetry,
      resolution_factor = resolution_factor)
    f_model_map_data = fft_map_ref.real_map_unpadded()
    zero_complex_ma = f_model.customized_copy(
      data = flex.complex_double(f_model.data().size(), 0))
    n_real = f_model_map_data.focus()
    # extract asu map from full P1
    asu_map = asu_map_ext.asymmetric_map(sgt, f_model_map_data)
    f_model_map_data_asu = asu_map.data()
    b = maptbx.boxes(
      n_real   = f_model_map_data_asu.focus(),
      fraction = box_size_as_fraction,
      log      = log)
    self.map_result_asu = flex.double(flex.grid(b.n_real))
    assert f_model_map_data_asu.focus()==b.n_real
    assert b.n_real==self.map_result_asu.focus()
    n_real_asu = b.n_real
    self.r = flex.double() # for regression test only
    n_boxes = len(b.starts)
    i_box = 0
    for s,e in zip(b.starts, b.ends):
      if(log): print "box %2d of %2d:"%(i_box, n_boxes), s, e
      i_box+=1
      # define wide box: neutral + phased volumes
      if(neutral_colume_box_cushion_width>0):
        sh = neutral_colume_box_cushion_width
        ss = [max(s[i]-sh,0) for i in [0,1,2]]
        ee = [min(e[i]+sh,n_real_asu[i]) for i in [0,1,2]]
      else: ss,ee = s,e
      # omit wide box from f_model map, repeat n_debias_cycles times
      f_model_map_data_asu_ = f_model_map_data_asu.deep_copy()
      for i in xrange(n_debias_cycles):
        f_model_omit, f_model_map_data_asu_ = self.omit_box(s=ss, e=ee, sgt=sgt,
          ma=zero_complex_ma, cg=fft_map_ref, md_asu=f_model_map_data_asu_,
          n_real=n_real)
      # get fmodel for omit map calculation
      fmodel_ = mmtbx.f_model.manager(
        f_obs        = fmodel.f_obs(),
        r_free_flags = fmodel.r_free_flags(),
        f_calc       = f_model_omit,
        f_mask       = zero_complex_ma)
      self.r.append(fmodel_.r_work()) # for regression test only
      f_map_data_asu = get_map(fmodel=fmodel_, map_type=map_type,
        fft_map_ref=fft_map_ref)
      assert f_map_data_asu.focus() == self.map_result_asu.focus()
      maptbx.copy_box(
        map_data_from = f_map_data_asu,
        map_data_to   = self.map_result_asu,
        start         = s,
        end           = e)
    # result
    sd = self.map_result_asu.sample_standard_deviation()
    self.map_result_asu = self.map_result_asu/sd
    asu_map_omit = asu_map_ext.asymmetric_map(sgt, self.map_result_asu, n_real)
    self.map_coefficients = f_model.customized_copy(
      indices = f_model.indices(),
      data    = asu_map_omit.structure_factors(f_model.indices()))
    # filling
    if(full_resolution_map):
      cs = f_model.complete_set(d_min=f_model.d_min())
      map_result_asu_trunc = self.map_result_asu.set_selected(
        self.map_result_asu<1, 0)
      asu_map_omit = asu_map_ext.asymmetric_map(sgt,map_result_asu_trunc,n_real)
      fill = f_model.customized_copy(
        indices = cs.indices(),
        data    = asu_map_omit.structure_factors(cs.indices()))
      self.map_coefficients = self.map_coefficients.complete_with(
        other=fill, scale=True)

  def omit_box(self, s, e, sgt, ma, cg, md_asu, n_real):
    md_asu_omit = maptbx.set_box_copy(value = 0, map_data_to = md_asu,
      start = s, end = e)
    asu_map_omit = asu_map_ext.asymmetric_map(sgt, md_asu_omit, n_real)
    ma_omit = ma.customized_copy(
      indices = ma.indices(),
      data    = asu_map_omit.structure_factors(ma.indices()))
    fft_map = cctbx.miller.fft_map(
      crystal_gridding = cg, fourier_coefficients = ma_omit)
    md = fft_map.real_map_unpadded()
    asu_map = asu_map_ext.asymmetric_map(sgt, md)
    md_asu_omit = asu_map.data()
    return ma_omit, md_asu_omit
