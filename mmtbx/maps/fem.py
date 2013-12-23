from __future__ import division
from cctbx import miller
from cctbx import maptbx
import mmtbx.maps.kick
import mmtbx.maps

class run(object):

  def __init__(self, fmodel, resolution_factor=0.25):
    crystal_gridding = fmodel.f_obs().crystal_gridding(
      d_min              = fmodel.f_obs().d_min(),
      resolution_factor  = resolution_factor)
    ko = mmtbx.maps.kick.run(fmodel = fmodel, crystal_gridding=crystal_gridding)
    # cut obvious noise
    m = ko.map_data_result
    maptbx.reset(
      data=m,
      substitute_value=0.0,
      less_than_threshold=0.5,
      greater_than_threshold=-9999,
      use_and=True)
    # smooth binary edges
    for i in xrange(3):
      maptbx.map_box_average(
        map_data   = m,
        cutoff     = 0.6,
        index_span = 1)
    # deblur using unsharp mask
    maptbx.sharpen(map_data=m, index_span=1, n_averages=2)
    # equalize histogram
    m  = maptbx.volume_scale(map = m,  n_bins = 10000).map_data()
    # compute map coefficients from map
    mc = ko.mc_result.structure_factors_from_map(
      map            = m,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
    # b-factor sharpen
    self.mc_result=mmtbx.maps.b_factor_sharpening_by_map_kurtosis_maximization(
      map_coeffs=mc)
    fft_map = miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = self.mc_result)
    fft_map.apply_sigma_scaling()
    m_tmp = fft_map.real_map_unpadded()
    # remove new noise due to sharpening
    s = m<0.9
    mask = m.set_selected(s, 0)
    mask = mask.set_selected(~s, 1)
    m_tmp = m_tmp * mask
    self.mc_result = self.mc_result.structure_factors_from_map(
      map            = m_tmp,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
