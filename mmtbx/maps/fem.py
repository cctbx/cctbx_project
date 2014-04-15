from __future__ import division
from cctbx import miller
from cctbx import maptbx
import mmtbx.maps
import mmtbx.map_tools
import random
from scitbx.array_family import flex
from cctbx import adptbx
from libtbx.test_utils import approx_equal
from mmtbx.maps import kick
import sys
from libtbx import adopt_init_args

def get_map(mc, cg=None):
  if(cg is None):
    fft_map = mc.fft_map(resolution_factor=0.25)
  else:
    fft_map = miller.fft_map(
      crystal_gridding     = cg,
      fourier_coefficients = mc)
  fft_map.apply_sigma_scaling()
  return fft_map.real_map_unpadded()

class counter(object):
  def __init__(self, n1, n2, log):
    adopt_init_args(self, locals())
    self.n=0
    assert self.n1*self.n2>0
    self.progress_scale = 100./(self.n1*self.n2)
  def show(self):
    self.n += 1
    self.log.write(
      "\r%s %d%%" %("FEM loop: done so far:",
      int(self.n*self.progress_scale)))
    self.log.flush()

class run(object):

  def __init__(self, fmodel, signal_threshold, use_omit=False, sharp=True):
    # make sure common B moved to k_total
    if(sharp):
      b_iso_min = flex.min(
        fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1))
      assert approx_equal(b_iso_min, 0, 1.e-3)
    # usual 2mFo-DFc mc
    self.mc = mmtbx.map_tools.electron_density_map(
      fmodel=fmodel).map_coefficients(
        map_type     = "2mFo-DFc",
        isotropize   = False,
        fill_missing = False)
    self.mc_iso = mmtbx.map_tools.electron_density_map(
      fmodel=fmodel).map_coefficients(
        map_type         = "2mFo-DFc",
        isotropize       = True,
        fill_missing     = False)
    # define gridding
    crystal_gridding = fmodel.f_obs().crystal_gridding(
      d_min              = fmodel.f_obs().d_min(),
      resolution_factor  = 0.25)
    # Fill missing
    missing_reflections_manager = mmtbx.map_tools.model_missing_reflections(
      coeffs=self.mc_iso, fmodel=fmodel)
    missing = missing_reflections_manager.get_missing(deterministic=True)
    self.mc_filled        = self.mc_iso.complete_with(other=missing, scale=True)
    self.mc_filled_no_iso = self.mc.complete_with(other=missing, scale=True)
    # composite OMIT map
    m_omit = None
    if(use_omit):
      como = mmtbx.maps.composite_omit_map.run(
        map_type                  = "mFo-DFc",
        crystal_gridding          = crystal_gridding,
        fmodel                    = fmodel.deep_copy(), # XXX
        reset_to_zero_below_sigma = None,
        use_shelx_weight          = False)
      mc_omit = como.map_coefficients
      mc_omit = mc_omit.complete_with(other=self.mc_filled, scale=True)
      m_omit = get_map(mc_omit, crystal_gridding)
      maptbx.reset(
        data=m_omit,
        substitute_value=0.0,
        less_than_threshold=0.5,
        greater_than_threshold=-9999,
        use_and=True)
      m_omit = m_omit.set_selected(m_omit< 0.5, 0)
      m_omit = m_omit.set_selected(m_omit>=0.5, 1)
    # Extra filter
    m_filter = None
    for mc in [self.mc, self.mc_iso, self.mc_filled, self.mc_filled_no_iso]:
      m_filter_ = get_map(mc = mc, cg=crystal_gridding)
      maptbx.binarize(
        map_data               = m_filter_,
        threshold              = 0.5,
        substitute_value_below = 0,
        substitute_value_above = 1)
      if(m_filter is None): m_filter = m_filter_
      else:                 m_filter = m_filter * m_filter_
    # FEM loop
    progress_counter = counter(n1=10, n2=16, log=sys.stdout)
    self.map_result = fem_loop(
      fmodel           = fmodel,
      mc               = self.mc_iso,
      n                = 10,
      n_cycles         = 50,
      n_macro_cycles   = 16,
      signal_threshold = signal_threshold,
      crystal_gridding = crystal_gridding,
      m_filter         = m_filter,
      progress_counter = progress_counter)
    print
    #
    if(m_omit is not None):
      self.map_result = self.map_result*m_omit
    self.mc_result = self.mc_filled.structure_factors_from_map(
      map            = self.map_result,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)

def ccp4_map(map_data, unit_cell, space_group, n_real, file_name):
  from iotbx import ccp4_map
  ccp4_map.write_ccp4_map(
      file_name=file_name,
      unit_cell=unit_cell,
      space_group=space_group,
      #gridding_first=(0,0,0),# This causes a bug (map gets shifted)
      #gridding_last=n_real,  # This causes a bug (map gets shifted)
      map_data=map_data,
      labels=flex.std_string([""]))

def fem_loop(fmodel, mc, n, crystal_gridding, signal_threshold, n_cycles,
             n_macro_cycles, m_filter, progress_counter=None):
  mac = maptbx.map_accumulator(n_real = crystal_gridding.n_real())
  for i in xrange(n_macro_cycles):
    m = fem_loop_(fmodel=fmodel, mc=mc, n=n, crystal_gridding=crystal_gridding,
      n_cycles=n_cycles, signal_threshold=signal_threshold,
      progress_counter=progress_counter)
    m = m*m_filter
    mac.add(map_data=m)
#    ccp4_map(map_data=m, unit_cell=mc.unit_cell(),
#      space_group=mc.space_group(), n_real=m.all(), file_name="%s.ccp4"%str(i))
#  ccp4_map(map_data=mm, unit_cell=mc.unit_cell(),
#      space_group=mc.space_group(), n_real=m.all(), file_name="ave.ccp4")
  return mac.as_median_map()

def fem_loop_(fmodel, mc, n, crystal_gridding, n_cycles, signal_threshold,
              progress_counter=None):
  m = None
  average_peak_volume = None
  for it in xrange(n):
    if(progress_counter is not None): progress_counter.show()
    #
    missing_reflections_manager = mmtbx.map_tools.model_missing_reflections(
      coeffs=mc, fmodel=fmodel)
    missing = missing_reflections_manager.get_missing(
      deterministic=random.choice([True, False]))
    wam=kick.weighted_average(fmodel=fmodel,map_coefficients=mc,missing=missing)
    #
    mc_wa = wam.random_weight_averaged_map_coefficients(
      random_scale  = 5.,
      random_seed   = int(random.random()*10000000),
      n_cycles      = n_cycles,
      fraction_keep = 0.95)
    m_ = get_map(mc=mc_wa, cg=crystal_gridding)
    maptbx.reset(
      data                   = m_,
      substitute_value       = 0.0,
      less_than_threshold    = signal_threshold,
      greater_than_threshold = -9999,
      use_and                = True)
    if(m is None):
      m = m_
      average_peak_volume = int(maptbx.peak_volume_estimate(
        map_data         = m_,
        sites_cart       = fmodel.xray_structure.sites_cart(),
        crystal_symmetry = fmodel.xray_structure.crystal_symmetry(),
        cutoff           = signal_threshold)*0.7)*+1
    else:
      m = m + m_
  m = m / n
  co = maptbx.connectivity(map_data=m, threshold=signal_threshold)
  m = m*co.volume_cutoff_mask(volume_cutoff=average_peak_volume).as_double()
  m = maptbx.volume_scale(map = m,  n_bins = 10000).map_data()
  return m
