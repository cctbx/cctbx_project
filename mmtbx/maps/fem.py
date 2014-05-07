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
      msk = truncate_with_roots2(m=m_filter_,fmodel=fmodel,c1=0.5,c2=0.25,
        cutoff=0.5,scale=0.7, as_int=True)
      maptbx.truncate_special(mask=msk, map_data=m_filter_)
      maptbx.binarize(
        map_data               = m_filter_,
        threshold              = 0.0,
        substitute_value_below = 0,
        substitute_value_above = 1)
      if(m_filter is None): m_filter = m_filter_
      else:                 m_filter = m_filter * m_filter_
    # FEM loop
    progress_counter = counter(n1=10, n2=16, log=sys.stdout)
    self.map_result = fem_loop(
      fmodel           = fmodel,
      complete_set = self.mc_filled,
      mc               = self.mc_iso,
      n                = 10,
      n_cycles         = 50,
      n_macro_cycles   = 16,
      signal_threshold = signal_threshold,
      crystal_gridding = crystal_gridding,
      m_filter         = m_filter,
      progress_counter = progress_counter)
#
    maptbx.sharpen(map_data=self.map_result, index_span=2, n_averages=1)
#
    mc = self.mc_filled.structure_factors_from_map(
      map            = self.map_result,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
    m = get_map(mc=mc, cg = crystal_gridding)

    average_peak_volume = int(maptbx.peak_volume_estimate(
      map_data         = m,
      sites_cart       = fmodel.xray_structure.sites_cart(),
      crystal_symmetry = fmodel.xray_structure.crystal_symmetry(),
      cutoff           = 1)*0.7)
    co = maptbx.connectivity(map_data=m, threshold=1.0)

    F = co.volume_cutoff_mask(volume_cutoff=average_peak_volume).as_double()
    self.map_result = self.map_result*F
    self.map_result = maptbx.volume_scale(map = self.map_result,  n_bins = 10000).map_data()
####
    keep_interblob=True
    F = None
    for c1 in [0.95, 0.9, 0.85, 0.8,0.7,0.6,0.5,0.4,0.3,0.2]:
      if(c1<0.9): keep_interblob=False
      f = truncate_with_roots2(m=self.map_result,fmodel=fmodel,c1=c1,c2=c1-0.1,
        cutoff=c1,scale=1, keep_interblob=keep_interblob)
      if(F is None): F = f
      else: F = F * f
#
    mc = self.mc_filled.structure_factors_from_map(
      map            = self.map_result,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
    m = get_map(mc=mc, cg = crystal_gridding)
    keep_interblob=True
    for c1 in [1.0,0.9,0.8,0.7,0.6,0.5]:
      if(c1<0.9): keep_interblob=False
      f = truncate_with_roots2(m=m,fmodel=fmodel,c1=c1,c2=c1-0.1,cutoff=c1,scale=1.0,
        keep_interblob=keep_interblob)
      if(F is None): F = f
      else: F = F * f
#
    self.map_result = self.map_result*F
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
             n_macro_cycles, m_filter, complete_set, progress_counter=None):
  mac = maptbx.map_accumulator(n_real = crystal_gridding.n_real())
  missing_reflections_manager = mmtbx.map_tools.model_missing_reflections(
    coeffs=mc, fmodel=fmodel)
  weighted_average_manager = kick.weighted_average(
    fmodel = fmodel, map_coefficients = mc)
  #
  b_sharp = get_b_sharp(fmodel=fmodel, mc=mc, n_cycles=n_cycles,
      crystal_gridding=crystal_gridding, signal_threshold=signal_threshold,
      weighted_average_manager=weighted_average_manager,
      missing_reflections_manager=missing_reflections_manager)
  print "b_sharp:", b_sharp
  #
  for i in xrange(n_macro_cycles):
    m = fem_loop_(
      fmodel                      = fmodel,
      complete_set = complete_set,
      mc                          = mc,
      n                           = n,
      crystal_gridding            = crystal_gridding,
      n_cycles                    = n_cycles,
      signal_threshold            = signal_threshold,
      progress_counter            = progress_counter,
      missing_reflections_manager = missing_reflections_manager,
      b_sharp                     = b_sharp,
      weighted_average_manager    = weighted_average_manager)
    m = m*m_filter
    mac.add(map_data=m)
#    ccp4_map(map_data=m, unit_cell=mc.unit_cell(),
#      space_group=mc.space_group(), n_real=m.all(), file_name="%s.ccp4"%str(i))
#  ccp4_map(map_data=mm, unit_cell=mc.unit_cell(),
#      space_group=mc.space_group(), n_real=m.all(), file_name="ave.ccp4")
  m = mac.as_median_map()
  m = maptbx.volume_scale(map = m,  n_bins = 10000).map_data()
  return m

def one_iteration(
      fmodel,
      mc,
      n_cycles,
      crystal_gridding,
      signal_threshold,
      missing_reflections_manager,
      b_sharp,
      weighted_average_manager):
  mc_wa = weighted_average_manager.random_weight_averaged_map_coefficients(
    random_scale  = 5.,
    random_seed   = int(random.random()*10000000),
    n_cycles      = n_cycles,
    missing       = missing_reflections_manager.get_missing_fast(),
    fraction_keep = 0.95)
  if(b_sharp is not None): # and random.choice([True, False])):
  #if(b_sharp is not None and random.choice([True, False])):
    mc_wa = mmtbx.maps.b_factor_sharpening_by_map_kurtosis_maximization(
      map_coeffs = mc_wa, show = False, b_sharp_best = b_sharp)
  m = get_map(mc=mc_wa, cg=crystal_gridding)
  ####
  #ccp4_map(map_data=m, unit_cell=mc.unit_cell(),
  #    space_group=mc.space_group(),
  #    n_real=m.all(),
  #    file_name="%s.ccp4"%random.choice(xrange(10000)))
  ####
  if(signal_threshold is not None):
    maptbx.reset(
      data                   = m,
      substitute_value       = 0.0,
      less_than_threshold    = signal_threshold,
      greater_than_threshold = -9999,
      use_and                = True)
  return m

def fem_loop_(
      fmodel,
      mc,
      n,
      crystal_gridding,
      n_cycles,
      signal_threshold,
      weighted_average_manager,
      missing_reflections_manager,
      b_sharp,
      complete_set,
      progress_counter=None,
      cut_by_connectivity=True):
  m = None
  average_peak_volume = None
  average_peak_volume_he = None
  for it in xrange(n):
    if(progress_counter is not None): progress_counter.show()
    m_ = one_iteration(
      fmodel           = fmodel,
      mc               = mc,
      n_cycles         = n_cycles,
      crystal_gridding = crystal_gridding,
      signal_threshold = signal_threshold,
      b_sharp=b_sharp,
      missing_reflections_manager=missing_reflections_manager,
      weighted_average_manager = weighted_average_manager)
    if(m is None): m = m_
    else:          m = m + m_
  m = m / n
  #
  if(0):
    average_peak_volume = int(maptbx.peak_volume_estimate(
      map_data         = m,
      sites_cart       = fmodel.xray_structure.sites_cart(),
      crystal_symmetry = fmodel.xray_structure.crystal_symmetry(),
      cutoff           = signal_threshold)*0.7)
    co = maptbx.connectivity(map_data=m, threshold=signal_threshold)
    m = m*co.volume_cutoff_mask(volume_cutoff=average_peak_volume).as_double()
    m = maptbx.volume_scale(map = m,  n_bins = 10000).map_data()
  else:
  #
    if(cut_by_connectivity):
      msk = truncate_with_roots2(
        m=m,fmodel=fmodel,c1=0.5,c2=0.25,cutoff=0.5,scale=0.7, as_int=True)
      maptbx.truncate_special(mask=msk, map_data=m)
      #
      F = None
      for c1 in [0.7,0.6,]:
        f = truncate_with_roots2(m=m,fmodel=fmodel,c1=c1,c2=c1-0.1,
          cutoff=c1,scale=0.7, keep_interblob=True)
        if(F is None): F = f
        else: F = F * f
      m = m * F
      #
    m = maptbx.volume_scale(map = m,  n_bins = 10000).map_data()
    #if(random.choice([True, False, False])):
    if(random.choice([True, False])):
      F = None
      keep_interblob=True
      for c1 in [0.9, 0.85, 0.8, 0.7]:
        if(c1<0.9): keep_interblob=False
        f = truncate_with_roots2(m=m,fmodel=fmodel,c1=c1,c2=c1-0.1,
          cutoff=c1,scale=1, keep_interblob=keep_interblob)
        if(F is None): F = f
        else: F = F * f
      m = m * F
  return m

def truncate_with_roots(m, fmodel, c1, c2, cutoff, scale, keep_interblob=False):
  assert c1>=c2
  average_peak_volume = int(maptbx.peak_volume_estimate(
    map_data         = m,
    sites_cart       = fmodel.xray_structure.sites_cart(),
    crystal_symmetry = fmodel.xray_structure.crystal_symmetry(),
    cutoff           = cutoff)*scale)
  co1 = maptbx.connectivity(map_data=m, threshold=c1)
  co2 = maptbx.connectivity(map_data=m, threshold=c2)
  res_mask = co2.noise_elimination_two_cutoffs(connectivity_t1=co1,
    volume_threshold_t1=average_peak_volume, keep_interblob=keep_interblob).as_double()
  return m*res_mask

def truncate_with_roots2(m, fmodel, c1, c2, cutoff, scale, keep_interblob=False, as_int=False):
  assert c1>=c2
  average_peak_volume = int(maptbx.peak_volume_estimate(
    map_data         = m,
    sites_cart       = fmodel.xray_structure.sites_cart(),
    crystal_symmetry = fmodel.xray_structure.crystal_symmetry(),
    cutoff           = cutoff)*scale)
  co1 = maptbx.connectivity(map_data=m, threshold=c1)
  co2 = maptbx.connectivity(map_data=m, threshold=c2)
  result = co2.noise_elimination_two_cutoffs(connectivity_t1=co1,
    volume_threshold_t1=average_peak_volume, keep_interblob=keep_interblob)
  if(as_int): return result
  else:       return result.as_double()

def get_b_sharp(fmodel, mc, n_cycles, crystal_gridding, signal_threshold,
      weighted_average_manager, missing_reflections_manager):
  result = flex.double()
  #return None
  for i in xrange(1):
    mc_wa = weighted_average_manager.random_weight_averaged_map_coefficients(
      random_scale  = 5.,
      random_seed   = int(random.random()*10000000),
      n_cycles      = n_cycles,
      missing       = missing_reflections_manager.get_missing_fast(),
      fraction_keep = 0.95)
    m = get_map(mc=mc_wa, cg=crystal_gridding)
    maptbx.reset(
      data                   = m,
      substitute_value       = 0.0,
      less_than_threshold    = signal_threshold,
      greater_than_threshold = -9999,
      use_and                = True)
    tmp = mc_wa.structure_factors_from_map(
      map            = m,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
    b =  mmtbx.maps.b_factor_sharpening_by_map_kurtosis_maximization(
      map_coeffs=tmp, show=False, b_only=True)
    result.append(b)
  return flex.min(result)
