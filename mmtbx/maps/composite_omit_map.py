from __future__ import division
from scitbx.array_family import flex
from cctbx import miller
from cctbx import maptbx
import boost.python
asu_map_ext = boost.python.import_ext("cctbx_asymmetric_map_ext")
import mmtbx.f_model
import cctbx.miller

class run(object):
  """
  Composite full-omit map: omit entire box in real-space corresponding to Fmodel,
  which includs atomic model and non-atomic model (bulk-solvent and scales).
  """
  def __init__(
        self,
        crystal_gridding,
        fmodel,
        map_type="2mFo-DFc",
        box_size_as_unit_cell_fraction=0.03,
        n_debias_cycles=10,
        reset_below_sigma=None,
        equalize=False,
        external_complete_set=None):
    sgt = fmodel.f_obs().space_group().type()
    def get_map(fmodel, map_type, external_complete_set=None):
      f_map = fmodel.electron_density_map(
        update_f_part1=False).map_coefficients(
          map_type     = map_type,
          isotropize   = True,
          fill_missing = False)
      if(external_complete_set is not None):
        f_map = f_map.complete_with(external_complete_set, scale=True)
      fft_map = cctbx.miller.fft_map(
        crystal_gridding     = crystal_gridding,
        fourier_coefficients = f_map)
      return fft_map.real_map_unpadded(), f_map
    f_model = fmodel.f_model_scaled_with_k1()
    fft_map = cctbx.miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = f_model)
    f_model_map_data = fft_map.real_map_unpadded()
    zero_complex_ma = f_model.customized_copy(
      data = flex.complex_double(f_model.data().size(), 0))
    n_real = f_model_map_data.focus()
    # extract asu map from full P1
    asu_map = asu_map_ext.asymmetric_map(sgt, f_model_map_data)
    f_model_map_data_asu = asu_map.data()
    b = maptbx.boxes(
      n_real     = f_model_map_data_asu.focus(),
      unit_cell  = fmodel.f_obs().crystal_symmetry().unit_cell(),
      box_size_as_unit_cell_fraction = box_size_as_unit_cell_fraction)
    self.map_result = flex.double(flex.grid(b.n_real))
    self.r = flex.double()
    for s,e in zip(b.starts, b.ends):
      print "box:", s,e
      # de-bias map outside the box start
      f_model_map_data_asu_ = f_model_map_data_asu.deep_copy()
      for i in xrange(n_debias_cycles): # 1 is equivalnet to usual scenario
        f_model_map_data_asu_omit = maptbx.set_box_copy(
          value         = 0,
          map_data_to   = f_model_map_data_asu_,
          start         = s,
          end           = e)
        asu_map_omit = asu_map_ext.asymmetric_map(sgt,f_model_map_data_asu_omit,
          n_real)
        f_model_omit = f_model.customized_copy(
          indices = f_model.indices(),
          data    = asu_map_omit.structure_factors(f_model.indices()))
        #
        fft_map = cctbx.miller.fft_map(
          crystal_gridding     = crystal_gridding,
          fourier_coefficients = f_model_omit)
        f_model_map_data_ = fft_map.real_map_unpadded()
        asu_map = asu_map_ext.asymmetric_map(sgt, f_model_map_data_)
        f_model_map_data_asu_ = asu_map.data()
      # de-bias map outside the box end
      fmodel_ = mmtbx.f_model.manager(
        f_obs        = fmodel.f_obs(),
        r_free_flags = fmodel.r_free_flags(),
        f_calc       = f_model_omit,
        f_mask       = zero_complex_ma)
      self.r.append(fmodel_.r_work())
      f_map_data, tmp_mc = get_map(fmodel=fmodel_, map_type=map_type)
      etmp = [e[0]-1, e[1]-1, e[2]-1] # because .copy() includes right edge
      box = maptbx.copy(f_map_data, s, etmp)
      box.reshape(flex.grid(box.all()))
      if(reset_below_sigma is not None):
        sd = box.sample_standard_deviation()
        box = box/sd
        maptbx.reset(
          data                   = box,
          substitute_value       = 0,
          less_than_threshold    = reset_below_sigma,
          greater_than_threshold = flex.min(box)-1.e-3,
          use_and                = True)
      if(equalize):
        box = maptbx.volume_scale(map = box, n_bins = 1000).map_data()
      maptbx.set_box(
        map_data_from = box,
        map_data_to   = self.map_result,
        start         = s,
        end           = e)
    sd = self.map_result.sample_standard_deviation()
    self.map_result = self.map_result/sd
    asu_map_omit = asu_map_ext.asymmetric_map(sgt, self.map_result, n_real)
    self.map_coefficients = f_model.customized_copy(
      indices = f_model.indices(),
      data    = asu_map_omit.structure_factors(f_model.indices()))
