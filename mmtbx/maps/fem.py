from __future__ import absolute_import, division, print_function
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
from cctbx import adptbx
from mmtbx import map_tools
from cctbx import miller
from libtbx import Auto
import mmtbx.maps.composite_omit_map
from six.moves import range

class run(object):

  def __init__(
        self,
        f_obs,
        r_free_flags,
        xray_structure,
        use_resolve,
        use_omit,
        use_max_map,
        sharp,
        use_unsharp_masking,
        resolution_factor,
        n_inner_loop = 10,
        n_outer_loop = 10,
        log          = None):
    adopt_init_args(self, locals())
    if(self.log is None): self.log = sys.stdout
    print("Start FEM...", file=self.log)
    self.prepare_f_obs_and_flags()
    self.mc_orig = self.compute_original_map()
    self.b_overall = None
    print("Create un-sharpened fmodel...", file=self.log)
    self.fmodel_nosharp = self.create_fmodel(show=True).deep_copy()
    if(self.sharp): self.remove_common_isotropic_adp()
    print("Create fmodel...", file=self.log)
    self.fmodel = self.create_fmodel(show=True)
    self.crystal_gridding = self.f_obs.crystal_gridding(
      d_min             = self.fmodel.f_obs().d_min(),
      symmetry_flags    = maptbx.use_space_group_symmetry,
      resolution_factor = self.resolution_factor)
    # initial maps
    print("Compute initial maps...", file=self.log)
    self.mc = map_tools.electron_density_map(
      fmodel=self.fmodel).map_coefficients(
        map_type     = "2mFo-DFc",
        isotropize   = True,
        fill_missing = False)
    self.mc_fill = map_tools.electron_density_map(
      fmodel=self.fmodel).map_coefficients(
        map_type     = "2mFo-DFc",
        isotropize   = True,
        fill_missing = True)
    print("Finding typical atom volumes...", file=self.log)
    self.selection = good_atoms_selection(
      crystal_gridding = self.crystal_gridding,
      map_coeffs       = self.mc_fill,
      xray_structure   = self.fmodel.xray_structure)
    # main FEM loop
    m = self.outer_loop()
    m = low_volume_density_elimination(m=m, fmodel=self.fmodel,
      selection=self.selection)
    self.zero_below_threshold(m = m)
    # HE
    m = maptbx.volume_scale(map = m,  n_bins = 10000).map_data()
    # create Fourier ripples filter
    sel = m<0.5
    msk = m.set_selected(sel, 0)
    sel = msk>=0.5
    msk = msk.set_selected(sel, 1)

#XXX    maptbx.sharpen(map_data=m, index_span=2, n_averages=1,
#XXX          allow_negatives=False)


    # Use Resolve filter
    m_resolve = self.resolve_filter_map()
    if(m_resolve is not None): m = m * m_resolve
    # Use OMIT
    if(self.use_omit):
      comit = mmtbx.maps.composite_omit_map.run(
        fmodel                           = self.fmodel_nosharp,
        crystal_gridding                 = self.crystal_gridding,
        box_size_as_fraction             = 0.2,
        max_boxes                        = 2000,
        neutral_volume_box_cushion_width = 2,
        full_resolution_map              = True,
        log                              = self.log)
      omit_map = comit.as_p1_map()
      #ccp4_map(cg=self.crystal_gridding, file_name="omit1.ccp4", map_data=omit_map)
      omit_map = low_volume_density_elimination(m=omit_map, fmodel=self.fmodel,
        selection=self.selection,end=16)
      #ccp4_map(cg=self.crystal_gridding, file_name="omit2.ccp4", map_data=omit_map)
      sel      = omit_map<1.5
      omit_map = omit_map.set_selected(sel, 0)
      sel      = omit_map>=1.5
      omit_map = omit_map.set_selected(sel, 1)
      m = m * omit_map
      #
      # Extra filter: seems to ne redundant, TBD.
      #
      #omit_mc = comit.result_as_sf()
      #omit_map = get_map(mc=omit_mc, cg=self.crystal_gridding)
      #omit_map = low_volume_density_elimination(m=omit_map, fmodel=self.fmodel, selection=self.selection)
      #ccp4_map(cg=self.crystal_gridding, file_name="omit2.ccp4", map_data=omit_map)
      #sel      = omit_map<0.5
      #omit_map = omit_map.set_selected(sel, 0)
      #sel      = omit_map>=0.5
      #omit_map = omit_map.set_selected(sel, 1)
      #m = m * omit_map

    #
    self.mc_result = self.mc_fill.structure_factors_from_map(
      map            = m,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
    if(self.sharp):
      self.mc_result=mmtbx.maps.b_factor_sharpening_by_map_kurtosis_maximization(
        map_coeffs=self.mc_result, show=True, b_only=False)
    self.map_result = get_map(mc=self.mc_result, cg=self.crystal_gridding)*msk

  def resolve_filter_map(self):
    m_resolve = None
    cmpl = self.fmodel.f_obs().resolution_filter(d_min=6).completeness()
    if((self.use_resolve is Auto and
       (self.fmodel.r_work()>0.2 or self.b_overall>30.) or cmpl<0.7) or
       self.use_resolve is True):
      print("Running Resolve density modificaiton", file=self.log)
      mc_resolve = self.fmodel.resolve_dm_map_coefficients()
      m_resolve = get_map(mc=mc_resolve, cg=self.crystal_gridding)
      m_resolve = low_volume_density_elimination(m=m_resolve, fmodel=self.fmodel,
        selection=self.selection)
      m_resolve = m_resolve.set_selected(m_resolve < 0.25, 0)
      m_resolve = m_resolve.set_selected(m_resolve >=0.25, 1)
      print("Obtained Resolve filter", file=self.log)
    return m_resolve

  def write_output_files(self, mtz_file_name, ccp4_map_file_name, fem_label,
                         orig_label):
    mtz_dataset = self.mc_result.as_mtz_dataset(column_root_label=fem_label)
    mtz_dataset.add_miller_array(
      miller_array      = self.mc_orig,
      column_root_label = orig_label)
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = mtz_file_name)
    ccp4_map(cg=self.crystal_gridding, file_name=ccp4_map_file_name,
      map_data=self.map_result)

  def compute_original_map(self):
    return map_tools.electron_density_map(
      fmodel = self.create_fmodel(update_f_part1=False)).map_coefficients(
        map_type     = "2mFo-DFc",
        isotropize   = False,
        fill_missing = False)

  def outer_loop(self):
    missing_reflections_manager = mmtbx.map_tools.model_missing_reflections(
      coeffs=self.mc, fmodel=self.fmodel)
    missing = missing_reflections_manager.get_missing(deterministic=True)
    wam = kick.weighted_average(fmodel=self.fmodel, map_coefficients=self.mc)
    progress_counter = counter(
      n1=self.n_inner_loop, n2=self.n_outer_loop, log=self.log)
    map_accumulator = maptbx.map_accumulator(
      n_real = self.crystal_gridding.n_real(), smearing_b=1, max_peak_scale=100,
      smearing_span=5, use_max_map=self.use_max_map)
    for i in range(self.n_outer_loop):
      m = inner_loop(
        fmodel           = self.fmodel,
        wam              = wam,
        missing          = missing,
        crystal_gridding = self.crystal_gridding,
        n                = self.n_inner_loop,
        progress_counter = progress_counter,
        use_max_map      = self.use_max_map)
      m = low_volume_density_elimination(m=m, fmodel=self.fmodel,
        selection=self.selection)
      if(self.sharp and self.use_unsharp_masking):
        maptbx.sharpen(map_data=m, index_span=1, n_averages=2,
          allow_negatives=False)
        maptbx.gamma_compression(map_data=m, gamma=0.1)
      self.zero_below_threshold(m = m)
      m = m/flex.max(m)
      map_accumulator.add(map_data=m)
    m = map_accumulator.as_median_map()
    sd = m.sample_standard_deviation()
    print(file=self.log)
    return m/sd

  def remove_common_isotropic_adp(self):
    xrs = self.xray_structure
    b_iso_min = flex.min(xrs.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1))
    self.b_overall = b_iso_min
    print("Max B subtracted from atoms and used to sharpen map:", b_iso_min, file=self.log)
    xrs.shift_us(b_shift=-b_iso_min)
    b_iso_min = flex.min(xrs.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1))
    assert approx_equal(b_iso_min, 0, 1.e-3)

  def prepare_f_obs_and_flags(self):
    if not self.f_obs:
      from libtbx.utils import Sorry
      raise Sorry("No FOBS available?")
    if not self.r_free_flags:
      from libtbx.utils import Sorry
      raise Sorry("No r_free_flags (test set) available?")
    sel = self.f_obs.data()>0
    self.f_obs = self.f_obs.select(sel)
    self.r_free_flags = self.r_free_flags.select(sel)
    #
    merged = self.f_obs.as_non_anomalous_array().merge_equivalents()
    self.f_obs = merged.array().set_observation_type(self.f_obs)
    #
    merged = self.r_free_flags.as_non_anomalous_array().merge_equivalents()
    self.r_free_flags = merged.array().set_observation_type(self.r_free_flags)
    #
    self.f_obs, self.r_free_flags = self.f_obs.common_sets(self.r_free_flags)

  def create_fmodel(self, update_f_part1=True, show=False):
    fmodel = mmtbx.f_model.manager(
      f_obs          = self.f_obs,
      r_free_flags   = self.r_free_flags,
      xray_structure = self.xray_structure)
    fmodel.update_all_scales(update_f_part1 = update_f_part1)
    if(show):
      fmodel.show(show_header=False, show_approx=False)
      print("r_work=%6.4f r_free=%6.4f" % (fmodel.r_work(),
        fmodel.r_free()), file=self.log)
    return fmodel

  def zero_below_threshold(self, m):
    maptbx.reset(
      data                   = m,
      substitute_value       = 0.0,
      less_than_threshold    = 0.24,
      greater_than_threshold = -9999,
      use_and                = True)

def low_volume_density_elimination(m, fmodel, selection, end=11):
  rr = [i/10. for i in range(5,end)]#+[1.5, 1.75, 2.]
  rr.reverse()
  for c in rr:
    sh=0.25
    zero_all_interblob_region = False
    if(c<=0.5): zero_all_interblob_region=True
    msk = truncate_with_roots(
      m=m, fmodel=fmodel, c1=c, c2=c-sh, cutoff=c, scale=1, as_int=False,
      zero_all_interblob_region=zero_all_interblob_region,
      selection = selection)
    if(msk is not None):
      m = m * msk
  return m

def truncate_with_roots(
      m, fmodel, c1, c2, cutoff, scale, zero_all_interblob_region=True,
      as_int=False, average_peak_volume=None, selection=None):
  assert c1>=c2
  if(average_peak_volume is None):
    sites_cart = fmodel.xray_structure.sites_cart()
    if(selection is not None):
      sites_cart = sites_cart.select(selection)
    average_peak_volume = maptbx.peak_volume_estimate(
      map_data         = m,
      sites_cart       = sites_cart,
      crystal_symmetry = fmodel.xray_structure.crystal_symmetry(),
      cutoff           = cutoff)
  if(average_peak_volume is None or int(average_peak_volume*scale)-1==0):
    return None
  average_peak_volume = int(average_peak_volume*scale/2)-1 # XXX "/2" is ad hoc and I don't know why!
  co1 = maptbx.connectivity(map_data=m, threshold=c1)
  co2 = maptbx.connectivity(map_data=m, threshold=c2)
  result = co2.noise_elimination_two_cutoffs(
    connectivity_object_at_t1=co1,
    elimination_volume_threshold_at_t1=average_peak_volume,
    zero_all_interblob_region=zero_all_interblob_region)
  if(as_int): return result
  else:       return result.as_double()

def good_atoms_selection(
      crystal_gridding,
      map_coeffs,
      xray_structure):
  #XXX copy from model_missing_reflections map_tools.py, consolidate later
  #XXX Also look for similar crap in f_model.py
  fft_map = miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = map_coeffs)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  rho_atoms = flex.double()
  for site_frac in xray_structure.sites_frac():
    rho_atoms.append(map_data.eight_point_interpolation(site_frac))
  #rho_mean = flex.mean_default(rho_atoms.select(rho_atoms>1.0), 1.0)
  sel_exclude = rho_atoms < 1.0 # XXX ??? TRY 0.5!
  sites_cart = xray_structure.sites_cart()
  #
  f_calc = map_coeffs.structure_factors_from_scatterers(
    xray_structure = xray_structure).f_calc()
  fft_map = miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = f_calc)
  fft_map.apply_sigma_scaling()
  map_data2 = fft_map.real_map_unpadded()
  #
  hd_sel = xray_structure.hd_selection()
  for i_seq, site_cart in enumerate(sites_cart):
    selection = maptbx.grid_indices_around_sites(
      unit_cell  = map_coeffs.unit_cell(),
      fft_n_real = map_data.focus(),
      fft_m_real = map_data.all(),
      sites_cart = flex.vec3_double([site_cart]),
      site_radii = flex.double([1.5]))
    cc = flex.linear_correlation(x=map_data.select(selection),
      y=map_data2.select(selection)).coefficient()
    if(cc<0.7 or hd_sel[i_seq]): sel_exclude[i_seq] = True
  return ~sel_exclude

def inner_loop(fmodel, wam, missing, crystal_gridding, n, progress_counter,
               use_max_map):
  mac = maptbx.map_accumulator(n_real = crystal_gridding.n_real(),
    smearing_b=1, max_peak_scale=100, smearing_span=5, use_max_map=use_max_map)
  for j in range(n):
    mc_w = wam.random_weight_averaged_map_coefficients(
      missing       = missing,
      random_scale  = random.choice([0,1,2,3,4,5]),
      random_seed   = random.choice(range(1, 9754365, 10000)),
      n_cycles      = 100,
      fraction_keep = random.choice([0.9, 0.95, 1.0]))
    if(random.choice([True,False])):
      mc_w = kick.randomize_struture_factors(map_coeffs=mc_w,
        number_of_kicks=100)
    m = get_map(mc=mc_w, cg=crystal_gridding)
    m = m/flex.max(m)
    sel = m<0
    m = m.set_selected(sel, 0)
    mac.add(map_data=m)
    progress_counter.show()
  mm = mac.as_median_map()
  sd = mm.sample_standard_deviation()
  r = mm/sd
  return r

def get_map(mc, cg):
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

def ccp4_map(cg, file_name, mc=None, map_data=None):
  assert [mc, map_data].count(None)==1
  if(map_data is None):
    map_data = get_map(mc=mc, cg=cg)
  from iotbx import mrcfile
  mrcfile.write_ccp4_map(
      file_name=file_name,
      unit_cell=cg.unit_cell(),
      space_group=cg.space_group(),
      #gridding_first=(0,0,0),# This causes a bug (map gets shifted)
      #gridding_last=n_real,  # This causes a bug (map gets shifted)
      map_data=map_data,
      labels=flex.std_string([""]))
