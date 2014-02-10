from __future__ import division
from cctbx import miller
from cctbx import maptbx
import mmtbx.maps
import mmtbx.map_tools
import random
from scitbx.array_family import flex
from cctbx import adptbx
from libtbx.test_utils import approx_equal

def get_map(mc, cg=None):
  if(cg is None):
    fft_map = mc.fft_map(resolution_factor=0.25)
  else:
    fft_map = miller.fft_map(
      crystal_gridding     = cg,
      fourier_coefficients = mc)
  fft_map.apply_sigma_scaling()
  return fft_map.real_map_unpadded()

def compute_map_and_combine(
      map_coeffs,
      map_data):
  m = get_map(mc=map_coeffs)
  maptbx.reset(
      data=m,
      substitute_value=0.0,
      less_than_threshold=0.0,
      greater_than_threshold=-9999,
      use_and=True)
  maptbx.map_box_average(
    map_data   = m,
    cutoff     = 0.0,
    index_span = 1)
  if(map_data is None): map_data = m
  else:
    map_data = (m+map_data)/2
  return map_data

def intersection(m1,m2,thresholds=flex.double([i/10. for i in xrange(9)]+[1.]),
      average=True):
  assert [m1,m2].count(None) != 2
  if(  m1 is None): return m2
  elif(m2 is None):  return m1
  else:
    maptbx.intersection(
      map_data_1 = m1,
      map_data_2 = m2,
      thresholds = thresholds,
      average    = average)
  return m1

class run(object):

  def __init__(self, fmodel, use_omit=True):
    # make sure common B moved to k_total
    b_iso_min = flex.min(
      fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1))
    assert approx_equal(b_iso_min, 0, 1.e-3)
    # usual 2mFo-DFc mc
    self.mc = mmtbx.map_tools.electron_density_map(
      fmodel=fmodel).map_coefficients(
        map_type     = "2mFo-DFc",
        isotropize   = False,
        fill_missing = False)
    self.mc_all_scales = mmtbx.map_tools.electron_density_map(
      fmodel=fmodel).map_coefficients(
        map_type         = "2mFo-DFc",
        isotropize       = True,
        use_shelx_weight = True,
        fill_missing     = False)
    # compute bs mask and define gridding
    crystal_gridding = fmodel.f_obs().crystal_gridding(
      d_min              = fmodel.f_obs().d_min(),
      resolution_factor  = 0.25)
    bs_mask = mmtbx.masks.mask_from_xray_structure(
      xray_structure           = fmodel.xray_structure,
      p1                       = True,
      solvent_radius           = 1.4,
      shrink_truncation_radius = 1.2,
      for_structure_factors    = False,
      n_real                   = crystal_gridding.n_real()).mask_data
    crystal_gridding = maptbx.crystal_gridding(
      unit_cell             = fmodel.f_obs().unit_cell(),
      space_group_info      = fmodel.f_obs().space_group_info(),
      pre_determined_n_real = bs_mask.all())
    # Resolve DM map
    self.mc_resolve = mmtbx.map_tools.resolve_dm_map(
      fmodel       = fmodel,
      map_coeffs   = self.mc_all_scales,
      pdb_inp      = None,
      use_model_hl = True,
      fill         = True)
    #XXX assert 100% complete XXX
    mc_missing = self.mc_resolve.lone_set(self.mc)
    m_resolve = get_map(self.mc_resolve, crystal_gridding)
    # Resolve filter mask
    m_resolve = m_resolve.set_selected(m_resolve< 0.25, 0)
    m_resolve = m_resolve.set_selected(m_resolve>=0.25, 1)
    s1 = bs_mask == 1
    s0 = bs_mask == 0
    assert s1.count(True)+s0.count(True) == bs_mask.size()
    bs_mask = bs_mask.set_selected(s1, 0)
    bs_mask = bs_mask.set_selected(s0, 1)
    m_resolve = m_resolve + bs_mask
    s2 = m_resolve > 1
    m_resolve = m_resolve.set_selected(s2, 1)
    # composite OMIT map
    m_omit = None
    if(use_omit):
      como = mmtbx.maps.composite_omit_map.run(
        map_type                  = "mFo-DFc",
        crystal_gridding          = crystal_gridding,
        fmodel                    = fmodel.deep_copy(), # XXX
        reset_to_zero_below_sigma = 0,
        use_shelx_weight          = True)
      mc_omit = como.map_coefficients
      mc_omit = mc_omit.complete_with(other=self.mc_resolve, scale=True)
      m_omit = get_map(mc_omit, crystal_gridding)
      maptbx.reset(
        data=m_omit,
        substitute_value=0.0,
        less_than_threshold=0.5,
        greater_than_threshold=-9999,
        use_and=True)
      m_omit = m_omit.set_selected(m_omit< 0.5, 0)
      m_omit = m_omit.set_selected(m_omit>=0.5, 1)
    #
    # FEM loop 3
    self.map_result = None
    for i in [1,2,3]:
      mc_w = loop_1_2(fmodel=fmodel, mc=self.mc, mc_fill=self.mc_resolve,
        he_cycles=5, random_weigh_cycles=5)
      m=get_map(mc_w, crystal_gridding)

#       ################
#      # investigate more / run conditionally
#      mc_result=mmtbx.maps.b_factor_sharpening_by_map_kurtosis_maximization(
#        map_coeffs=mc_w)
#      m_tmp=get_map(mc_result, crystal_gridding)
#      # remove new noise due to sharpening
#      s = m_tmp<0.9
#      m = m_tmp.deep_copy()
#      mask = m.set_selected(s, 0)
#      mask = mask.set_selected(~s, 1)
#      m = m_tmp * mask
#      ################

      maptbx.reset(
        data=m,
        substitute_value=0.0,
        less_than_threshold=0.5,
        greater_than_threshold=-9999,
        use_and=True)
      m  = maptbx.volume_scale(map = m,  n_bins = 10000).map_data()
      m = m*m_resolve
      if(m_omit is not None): m = m*m_omit
      self.map_result = intersection(m1 = self.map_result, m2 = m)
    self.mc_result = self.mc_resolve.structure_factors_from_map(
      map            = self.map_result,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)

def loop_1_2(fmodel, mc, mc_fill, he_cycles=10, random_weigh_cycles=10):
  mr = None
  for i in xrange(he_cycles):
    print i
    m = None
    for scale in xrange(random_weigh_cycles):
      print "  ", scale
      weight = mmtbx.map_tools.shelx_weight(
        f_obs   = fmodel.f_obs(),
        f_model = fmodel.f_model_scaled_with_k1(),
        weight_parameter = None)
      data = mc.data()*weight
      mc_w = mc.customized_copy(data = data)
      if(mc_fill is not None):
        mc_w = mc_w.complete_with(other=mc_fill, scale=True)
      if(random.choice([True, False])):
        mc_w = mc_w.select(flex.random_bool(mc_w.size(), 0.95))
      m = compute_map_and_combine(map_coeffs=mc_w, map_data=m)
    maptbx.reset(
      data=m,
      substitute_value=0.0,
      less_than_threshold=0.5,
      greater_than_threshold=-9999,
      use_and=True)
    # run conditionally: low-res/high B
    maptbx.sharpen(map_data=m, index_span=1, n_averages=2)
    m  = maptbx.volume_scale(map = m,  n_bins = 10000).map_data()
    mr = intersection(m1 = mr, m2 = m)
  return mc_w.structure_factors_from_map(
    map            = mr,
    use_scale      = True,
    anomalous_flag = False,
    use_sg         = False)
