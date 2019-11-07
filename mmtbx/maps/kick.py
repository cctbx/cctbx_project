from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import random
from cctbx import miller
from cctbx import maptbx
from libtbx.test_utils import approx_equal
import mmtbx.f_model
from cctbx import maptbx
import mmtbx.maps.composite_omit_map
from six.moves import range

def get_map(map_coeffs, crystal_gridding):
  fft_map = miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = map_coeffs)
  fft_map.apply_sigma_scaling()
  return fft_map.real_map_unpadded()

def compute_map_and_combine(
      map_coeffs,
      crystal_gridding,
      map_data,
      thresholds = [0,0.1,0.2,0.3,0.4,0.5]):
  m = get_map(map_coeffs=map_coeffs, crystal_gridding=crystal_gridding)
  if(map_data is None): map_data = m
  else:
    maptbx.intersection(
      map_data_1 = m,
      map_data_2 = map_data,
      thresholds = flex.double(thresholds),
      average    = True)
  return map_data

class weighted_average(object):
  def __init__(self, fmodel, map_coefficients):
    self.map_coefficients_data = map_coefficients.data()
    self.map_coefficients = map_coefficients
    f_o = fmodel.f_obs().data()
    i_o = f_o * f_o
    i_c = abs(fmodel.f_model_scaled_with_k1()).data() * \
          abs(fmodel.f_model_scaled_with_k1()).data()
    sig = fmodel.f_obs().sigmas()
    self.r = abs(i_o-i_c)/abs(i_o+i_c)
    if(sig is None):
      sig = flex.double(f_o.size(), 1)
    self.so = sig/f_o

  def random_weight_averaged_map_coefficients(self,
        random_scale, random_seed, n_cycles, fraction_keep, missing=None):
    assert self.map_coefficients_data.size() == \
           self.r.size() == \
           self.so.size()
    mc_data = maptbx.fem_averaging_loop(
      map_coefficients = self.map_coefficients_data,
      r_factors        = self.r,
      sigma_over_f_obs = self.so,
      random_scale     = random_scale,
      random_seed      = random_seed,
      n_cycles         = n_cycles)
    mc_wa = self.map_coefficients.customized_copy(data = mc_data)
    if(missing is not None):
      mc_wa = mc_wa.complete_with(other=missing, scale=True)
    return mc_wa.select(flex.random_bool(mc_wa.size(), fraction_keep))

def randomize_completeness2(map_coeffs, crystal_gridding, n_cycles=10):
  map_filter = None
  from mmtbx.maps import fem
  for cycle in range(n_cycles):
    if(cycle>0):
      mc = map_coeffs.select(flex.random_bool(map_coeffs.size(), 0.99))
    else:
      mc = map_coeffs.deep_copy()
    m = get_map(map_coeffs=mc, crystal_gridding=crystal_gridding)
    maptbx.reset(
      data=m,
      substitute_value=0.0,
      less_than_threshold=0.5,
      greater_than_threshold=-9999,
      use_and=True)
    m  = maptbx.volume_scale(map = m,  n_bins = 10000).map_data()
    if(map_filter is not None):
      map_filter  = maptbx.volume_scale(map = map_filter,  n_bins = 10000).map_data()
    map_filter = fem.intersection(m1=map_filter, m2=m,
      thresholds=flex.double([i/100. for i in range(0,55,5)]),
      average=True)
  #map_filter = map_filter.set_selected(map_filter< 0.5, 0)
  #map_filter = map_filter.set_selected(map_filter>=0.5, 1)
  return map_filter

def randomize_completeness(map_coeffs, crystal_gridding, n_cycles=10,
      thresholds=[0,0.1,0.2,0.3,0.4,0.5]):
  map_data = None
  for cycle in range(n_cycles):
    if(cycle>0):
      mc = map_coeffs.select(flex.random_bool(map_coeffs.size(), 0.95))
    else:
      mc = map_coeffs.deep_copy()
    map_data = compute_map_and_combine(
      map_coeffs        = mc,
      crystal_gridding  = crystal_gridding,
      map_data          = map_data,
      thresholds        = thresholds)
  return map_data

def randomize_struture_factors(map_coeffs, number_of_kicks, phases_only=False):
  map_coeff_data = None
  for kick in range(number_of_kicks):
    rc, ar, pr = random.choice([(0.1, 0.10, 40),
                                (0.2, 0.09, 40),
                                (0.3, 0.08, 30),
                                (0.4, 0.07, 30),
                                (0.5, 0.06, 20),
                                (0.6, 0.05, 20),
                                (0.7, 0.04, 20),
                                (0.8, 0.03, 10),
                                (0.9, 0.02, 10),
                                (1.0, 0.01, 10)
                               ])
    if(phases_only): ar = 0
    sel = flex.random_bool(map_coeffs.size(), rc)
    mc = map_coeffs.randomize_amplitude_and_phase(
      amplitude_error=ar, phase_error_deg=pr, selection=sel)
    if(map_coeff_data is None): map_coeff_data = mc.data()
    else:                       map_coeff_data = map_coeff_data + mc.data()
  map_coeff_data = map_coeff_data/number_of_kicks
  return miller.set(
    crystal_symmetry = map_coeffs.crystal_symmetry(),
    indices          = map_coeffs.indices(),
    anomalous_flag   = False).array(data = map_coeff_data)

def kick_map_coeffs(
      map_coeffs,
      crystal_gridding,
      number_of_kicks,
      macro_cycles,
      missing           = None,
      kick_completeness = 0.95,
      phases_only       = False):
  map_data = None
  if(macro_cycles==0):
    map_data = compute_map_and_combine(
      map_coeffs       = map_coeffs,
      crystal_gridding = crystal_gridding,
      map_data         = map_data)
  for it in range(macro_cycles):
    print("  %d"%it)
    if(number_of_kicks>0):
      mc = randomize_struture_factors(map_coeffs=map_coeffs,
        number_of_kicks=number_of_kicks, phases_only=phases_only)
    else:
      mc = map_coeffs.deep_copy()
    if(missing is not None):
      mc = mc.complete_with(other=missing, scale=True)
    if(kick_completeness):
      mc = mc.select(flex.random_bool(mc.size(), kick_completeness))
    map_data = compute_map_and_combine(
      map_coeffs       = mc,
      crystal_gridding = crystal_gridding,
      map_data         = map_data)
  return map_data

def kick_fmodel(
      fmodel,
      map_type,
      crystal_gridding,
      number_of_kicks,
      macro_cycles,
      missing           = None,
      kick_completeness = 0.95):
  f_model = fmodel.f_model_no_scales()
  zero = fmodel.f_calc().customized_copy(data =
    flex.complex_double(fmodel.f_calc().data().size(), 0))
  fmodel_dc  = mmtbx.f_model.manager(
    f_obs         = fmodel.f_obs(),
    r_free_flags  = fmodel.r_free_flags(),
    k_isotropic   = fmodel.k_isotropic(),
    k_anisotropic = fmodel.k_anisotropic(),
    f_calc        = fmodel.f_model_no_scales(),
    f_part1       = fmodel.f_part1(),
    f_part2       = fmodel.f_part2(),
    f_mask        = zero)
  r1 = fmodel.r_work()
  r2 = fmodel_dc.r_work()
  assert approx_equal(r1, r2, 1.e-4), [r1, r2]
  def get_mc(fm):
   return fm.electron_density_map().map_coefficients(
       map_type     = map_type,
       isotropize   = True,
       fill_missing = False)
  def recreate_r_free_flags(fmodel):
    rc = random.choice([0.05, 0.9])
    r_free_flags = flex.random_bool(fmodel.f_obs().indices().size(), rc)
    fmodel._r_free_flags._data = r_free_flags
    return fmodel
  map_data = None
  for it in range(macro_cycles):
    print("  %d"%it)
    f_model_kick = randomize_struture_factors(map_coeffs=f_model,
      number_of_kicks=number_of_kicks)
    fmodel_dc = recreate_r_free_flags(fmodel = fmodel_dc)
    fmodel_dc.update(f_calc = f_model_kick)
    mc = get_mc(fm=fmodel_dc)
    if(missing is not None):
      mc = mc.complete_with(missing, scale=True)
    if(kick_completeness):
      mc = mc.select(flex.random_bool(mc.size(), kick_completeness))
    map_data = compute_map_and_combine(
      map_coeffs       = mc,
      crystal_gridding = crystal_gridding,
      map_data         = map_data)
  return map_data

class run(object):
  """
  Note 1: Assume Fmodel has correct scaling already (all f_parts etc).
  Note 2: More macro_cycles can substantially clean map in extremely bad cases.
  """

  def __init__(
      self,
      fmodel,
      map_type           = "2mFo-DFc",
      mask_data          = None,
      crystal_gridding   = None,
      number_of_kicks    = 100,
      macro_cycles       = 10,
      kick_completeness  = 0.95,
      omit               = True):
    fmodel = self.convert_to_non_anomalous(fmodel=fmodel)
    from mmtbx import map_tools
    self.mc_orig = map_tools.electron_density_map(
      fmodel=fmodel).map_coefficients(
        map_type     = "2mFo-DFc",
        isotropize   = True,
        fill_missing = False)
    md_orig = get_map(map_coeffs=self.mc_orig, crystal_gridding=crystal_gridding)
    # model missing
    self.complete_set = map_tools.resolve_dm_map(
      fmodel       = fmodel,
      map_coeffs   = self.mc_orig,
      pdb_inp      = None,
      use_model_hl = True,
      fill         = True)
    md_complete_set = get_map(map_coeffs=self.complete_set, crystal_gridding=crystal_gridding)
    self.missing = self.complete_set.lone_set(self.mc_orig)
    # Kick map coefficients
    md_kick = kick_map_coeffs(
      map_coeffs        = self.mc_orig,
      crystal_gridding  = crystal_gridding,
      number_of_kicks   = number_of_kicks,
      macro_cycles      = macro_cycles,
      missing           = self.missing,
      kick_completeness = kick_completeness)
    self.mc_kick = self.map_coeffs_from_map(map_data=md_kick)
    # Kick fmodel
    md_fm = kick_fmodel(
      fmodel            = fmodel,
      map_type          = map_type,
      crystal_gridding  = crystal_gridding,
      number_of_kicks   = number_of_kicks,
      macro_cycles      = macro_cycles,
      missing           = self.missing,
      kick_completeness = kick_completeness)
    self.mc_fm = self.map_coeffs_from_map(map_data=md_fm)
    if(omit):
    # Kick OMIT map
      com1 = mmtbx.maps.composite_omit_map.run(
        crystal_gridding     = crystal_gridding,
        fmodel               = fmodel.deep_copy(), # XXX
        full_resolution_map  = True,
        box_size_as_fraction = 0.03)
      md_com1 = kick_map_coeffs(
        map_coeffs        = com1.map_coefficients(),
        crystal_gridding  = crystal_gridding,
        number_of_kicks   = number_of_kicks,
        macro_cycles      = macro_cycles,
        phases_only       = True,
        missing           = self.missing,
        kick_completeness = kick_completeness)
      self.mc_com1 = self.map_coeffs_from_map(map_data=md_com1)
      # Kick OMIT map 2
      com2 = mmtbx.maps.composite_omit_map.run(
        crystal_gridding     = crystal_gridding,
        fmodel               = fmodel.deep_copy(), # XXX
        full_resolution_map  = True,
        box_size_as_fraction = 0.03)
      md_com2 = kick_map_coeffs(
        map_coeffs        = com2.map_coefficients(),
        crystal_gridding  = crystal_gridding,
        number_of_kicks   = number_of_kicks,
        macro_cycles      = macro_cycles,
        phases_only       = True,
        missing           = self.missing,
        kick_completeness = kick_completeness)
      self.mc_com2 = self.map_coeffs_from_map(map_data=md_com2)
    # combine maps
    def intersect(m1,m2, use_average):
      maptbx.intersection(
        map_data_1 = m1,
        map_data_2 = m2,
        thresholds = flex.double([0,0.1,0.2,0.3,0.4,0.5]),
        average    = False)
      if(use_average): return (m1+m2)/2
      else:            return m1
    m = (md_kick + md_fm)/2
    m = intersect(m, md_kick,         use_average=True)
    m = intersect(m, md_fm,           use_average=True)
    if(omit):
      m = intersect(m, md_com1,         use_average=False)
      m = intersect(m, md_com2,         use_average=False)
    m = intersect(m, md_orig,         use_average=False)
    m = intersect(m, md_complete_set, use_average=False)
    self.map_data_result = m
    self.mc_result  = self.map_coeffs_from_map(map_data=self.map_data_result)

  def write_mc(self, file_name="mc.mtz"):
    mtz_dataset = self.mc_orig.as_mtz_dataset(column_root_label="mc_orig")
    mtz_dataset.add_miller_array(
      miller_array=self.mc_kick,
      column_root_label="mc_kick")
    mtz_dataset.add_miller_array(
      miller_array=self.mc_fm,
      column_root_label="mc_fm")
    if(omit):
      mtz_dataset.add_miller_array(
        miller_array=self.mc_com1,
        column_root_label="mc_com1")
      mtz_dataset.add_miller_array(
        miller_array=self.mc_com2,
        column_root_label="mc_com2")
    mtz_dataset.add_miller_array(
      miller_array=self.mc_result,
      column_root_label="mc_result")
    mtz_dataset.add_miller_array(
      miller_array=self.complete_set,
      column_root_label="complete_set")
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = file_name)

  def map_coeffs_from_map(self, map_data):
    return self.mc_orig.structure_factors_from_map(
      map            = map_data,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)

  def convert_to_non_anomalous(self, fmodel):
    if(fmodel.f_obs().anomalous_flag()):
      f_obs        = fmodel.f_obs().average_bijvoet_mates()
      r_free_flags = fmodel.r_free_flags().average_bijvoet_mates()
      fmodel = mmtbx.f_model.manager(
        f_obs = f_obs,
        r_free_flags = r_free_flags,
        xray_structure = fmodel.xray_structure)
      fmodel.update_all_scales()
    return fmodel
