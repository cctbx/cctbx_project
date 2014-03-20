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

def get_map(mc, cg=None):
  if(cg is None):
    fft_map = mc.fft_map(resolution_factor=0.25)
  else:
    fft_map = miller.fft_map(
      crystal_gridding     = cg,
      fourier_coefficients = mc)
  fft_map.apply_sigma_scaling()
  return fft_map.real_map_unpadded()

def scale(x,y):
  n = flex.sum( flex.abs(x)*flex.abs(y) )
  d = flex.sum( flex.abs(y)*flex.abs(y) )
  return n/d

def compute_map_and_combine(
      map_coeffs,
      map_data,
      cg,
      intersection_thresholds):
  m = get_map(mc=map_coeffs, cg=cg)
  #maptbx.reset(
  #    data=m,
  #    substitute_value=0.0,
  #    less_than_threshold=0.0,
  #    greater_than_threshold=-9999,
  #    use_and=True)
  #maptbx.map_box_average(
  #  map_data   = m,
  #  cutoff     = 0.0,
  #  index_span = 1)
  if(map_data is None): map_data = m
  else:
    if(intersection_thresholds is not None):
      map_data=intersection(m1=map_data,m2=m,thresholds=intersection_thresholds)
    else:
      map_data = (m+map_data)/2
  return map_data

def intersection(m1,m2,thresholds=flex.double([i/10. for i in xrange(9)]+[1.]),
      average=True):
  assert [m1,m2].count(None) != 2
  if(  m1 is None): return m2
  elif(m2 is None): return m1
  else:
    maptbx.intersection(
      map_data_1 = m1,
      map_data_2 = m2,
      thresholds = thresholds,
      average    = average)
  return m1

class run(object):

  def __init__(self, fmodel, use_omit=False, use_resolve=False, sharp=True):
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
    # compute bs mask and define gridding
    crystal_gridding = fmodel.f_obs().crystal_gridding(
      d_min              = fmodel.f_obs().d_min(),
      resolution_factor  = 0.25)
    bs_mask = mmtbx.masks.mask_from_xray_structure(
      xray_structure           = fmodel.xray_structure,
      p1                       = True,
      solvent_radius           = 1.,#4,
      shrink_truncation_radius = 1.,#2,
      for_structure_factors    = False,
      n_real                   = crystal_gridding.n_real()).mask_data
    crystal_gridding = maptbx.crystal_gridding(
      unit_cell             = fmodel.f_obs().unit_cell(),
      space_group_info      = fmodel.f_obs().space_group_info(),
      pre_determined_n_real = bs_mask.all())
    # Fill missing
    if(use_resolve):
      self.mc_filled = mmtbx.map_tools.resolve_dm_map(
        fmodel       = fmodel,
        map_coeffs   = self.mc_iso,
        mask_cycles  = 5,
        minor_cycles = 5,
        pdb_inp      = None,
        use_model_hl = True,
        fill         = True)
    else:
      missing_reflections_manager = mmtbx.map_tools.model_missing_reflections(
        coeffs=self.mc_iso, fmodel=fmodel)
      missing = missing_reflections_manager.get_missing(deterministic=True)
      self.mc_filled = self.mc_iso.complete_with(other=missing, scale=True)
    m_filter = kick.randomize_completeness(
      map_coeffs       = self.mc_filled,
      crystal_gridding = crystal_gridding)
    # Filter masks
    m_filter = maptbx.conditional_solvent_region_filter(
      bulk_solvent_mask = bs_mask,
      map_data          = m_filter,
      threshold         = 0.5)
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
      m_filter = m_filter * m_omit
    #
    # FEM loop
    wam = kick.weighted_average(fmodel=fmodel, map_coefficients=self.mc,
      missing=self.mc_filled)#missing)
    self.map_result=None
    cntr = 0
    n_loop_1 = 10
    n_loop_2 = 10
    progress_scale = 100./(n_loop_1*n_loop_2)
    for it2 in xrange(n_loop_1):
      m,cntr= fem_loop_2(wam=wam, n=n_loop_2, crystal_gridding=crystal_gridding,
        cntr=cntr, progress_scale=progress_scale)
      maptbx.reset(
        data=m,
        substitute_value=0.0,
        less_than_threshold=0.25,
        greater_than_threshold=-9999,
        use_and=True)
      if(sharp):
        maptbx.sharpen(map_data=m, index_span=1, n_averages=2)
      m  = maptbx.volume_scale(map = m,  n_bins = 10000).map_data()
      self.map_result = intersection(m1 = self.map_result, m2 = m)
      #
      mc_ = self.mc_filled.structure_factors_from_map(
        map            = self.map_result,
        use_scale      = True,
        anomalous_flag = False,
        use_sg         = False)
      m_=get_map(mc_, crystal_gridding)
      m_ = m_.set_selected(m_< 1, 0)
      m_ = m_.set_selected(m_>=1, 1)
      m_filter = m_filter * m_
      #
    sys.stdout.write("\n")
    sys.stdout.flush()
    #
    self.map_result = self.map_result*m_filter
    self.mc_result = self.mc_filled.structure_factors_from_map(
      map            = self.map_result,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)

def fem_loop_2(wam, n, crystal_gridding, n_cycles=100, progress_scale=None,
               cntr=None, intersection_thresholds=None):
  m = None
  for it in xrange(n):
    if(cntr is not None):
      cntr+=1
      sys.stdout.write(
        "\r%s %d%%" %("FEM loop: done so far:", int(cntr*progress_scale)))
      sys.stdout.flush()
    mc_wa = wam.random_weight_averaged_map_coefficients(
      random_scale  = 5.,
      random_seed   = int(random.random()*10000000),
      n_cycles      = n_cycles,
      fraction_keep = 0.95)
    m = compute_map_and_combine(
      map_coeffs = mc_wa,
      map_data   = m,
      cg         = crystal_gridding,
      intersection_thresholds = intersection_thresholds)
  return m, cntr
