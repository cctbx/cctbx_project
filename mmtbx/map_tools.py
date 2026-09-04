"""
This module contains a variety of functions related to map calculation, many
of which are access via the mmtbx.f_model.manager API.  It has some overlap
with the separate mmtbx.maps module.
"""

from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import miller
from cctbx import maptbx
from libtbx.utils import null_out
import mmtbx
import libtbx
import random
import boost_adaptbx.boost.python as bp
from six.moves import range
mmtbx_f_model_ext = bp.import_ext("mmtbx_f_model_ext")
import mmtbx.masks

def shelx_weight(
      f_obs,
      f_model,
      weight_parameter=None):
  if(weight_parameter is not None):
    sc = weight_parameter
  else:
    sc = flex.double()
    for sc_ in range(f_obs.data().size()):
      sc.append(random.choice([1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2]))
  f_m = abs(f_model).data()
  f_o = abs(f_obs).data()
  i_c = f_m * f_m
  i_f = f_o * f_o
  sig_f = f_obs.sigmas()
  assert sig_f.size() == f_o.size()
  w = 1./( sig_f*sig_f*sig_f*sig_f + (sc * i_f)*(sc * i_f) )
  sigma_i = 1/flex.sqrt(w)
  r1 = i_c*i_c/(i_c*i_c + sigma_i*sigma_i)
  r2 = 1 / (1 + sig_f*sig_f*sig_f*sig_f/i_c/i_c + sc*sc*i_f*i_f/i_c/i_c)
  return r1

class fo_fc_scales(object):
  def __init__(self,
               fmodel,
               map_type_str,
               acentrics_scale = 2.0,
               centrics_scale = 1.0):
    """
    Compute x and y for x*Fobs-y*Fmodel given map type string
    """
    self.fmodel = fmodel
    mnm = mmtbx.map_names(map_name_string = map_type_str)
    # R.Read, SIGMAA: 2mFo-DFc (acentrics) & mFo (centrics)
    centric_flags  = self.fmodel.f_obs().centric_flags().data()
    acentric_flags = ~centric_flags
    if(mnm.k != 0):
      self.fo_scale = flex.double(self.fmodel.f_obs().size(), 1.0)
    else: self.fo_scale = flex.double(self.fmodel.f_obs().size(), 0.0)
    if(mnm.n != 0):
      self.fc_scale = flex.double(self.fmodel.f_obs().size(), 1.0)
    else: self.fc_scale = flex.double(self.fmodel.f_obs().size(), 0.0)
    abs_kn = abs(mnm.k*mnm.n)
    if(mnm.k != abs(mnm.n) and abs_kn > 1.e-6):
      self.fo_scale.set_selected(acentric_flags, mnm.k)
      self.fo_scale.set_selected( centric_flags, max(mnm.k-centrics_scale,0.))
      self.fc_scale.set_selected(acentric_flags, mnm.n)
      self.fc_scale.set_selected( centric_flags, max(mnm.n-centrics_scale,0.))
    elif(mnm.k == abs(mnm.n) and abs_kn > 1.e-6):
      fo_scale_k = self.fo_scale*mnm.k
      self.fc_scale_n = self.fc_scale*mnm.n
      self.fo_scale.set_selected(acentric_flags, fo_scale_k*acentrics_scale)
      self.fo_scale.set_selected( centric_flags, fo_scale_k)
      self.fc_scale.set_selected(acentric_flags, self.fc_scale_n*acentrics_scale)
      self.fc_scale.set_selected( centric_flags, self.fc_scale_n)
    else:
      self.fo_scale *= mnm.k
      self.fc_scale *= mnm.n

class combine(object):
  def __init__(self,
               fmodel,
               map_type_str,
               fo_scale,
               fc_scale,
               use_shelx_weight,
               shelx_weight_parameter,
               map_calculation_helper=None):
    self.mch = map_calculation_helper
    self.mnm = mmtbx.map_names(map_name_string = map_type_str)
    self.fmodel = fmodel
    self.fc_scale = fc_scale
    self.fo_scale = fo_scale
    self.f_obs = None
    self.f_model = None
    self.use_shelx_weight = use_shelx_weight
    self.shelx_weight_parameter = shelx_weight_parameter
    if(not self.mnm.ml_map):
      self.f_obs = self.fmodel.f_obs().data()*fo_scale
    else:
      if(self.mch is None): self.mch = self.fmodel.map_calculation_helper()
      if(self.fmodel.hl_coeffs() is None):
        self.f_obs = self.mch.f_obs.data()*fo_scale*self.mch.fom
      else:
        cp = fmodel.combine_phases(map_calculation_helper = self.mch)
        self.f_obs = self.mch.f_obs.data()*fo_scale*cp.f_obs_phase_and_fom_source()

  def map_coefficients(self, f_model=None):
    def compute(fo,fc,miller_set):
      if(type(fo) == flex.double):
        fo = miller.array(miller_set=miller_set, data=fo).phase_transfer(
          phase_source = self.f_model).data()
      return miller.array(miller_set=miller_set, data=fo+fc)
    if(f_model is None):
      if(not self.mnm.ml_map):
        self.f_model = self.fmodel.f_model_scaled_with_k1()
        f_model_data = self.f_model.data()*self.fc_scale
      else:
        self.f_model = self.mch.f_model
        f_model_data = self.f_model.data()*self.fc_scale*self.mch.alpha.data()
    else:
      self.f_model = f_model
      if(not self.mnm.ml_map):
        f_model_data = self.f_model.data()*self.fc_scale
      else:
        f_model_data = self.f_model.data()*self.fc_scale*self.mch.alpha.data()
    # f_model_data may be multiplied by scales like "-1", so it cannot be
    # phase source !
    result = compute(fo=self.f_obs, fc=f_model_data,
      miller_set=self.fmodel.f_obs())
    if(self.use_shelx_weight):
      sw = shelx_weight(
        f_obs   = self.fmodel.f_obs(),
        f_model = self.fmodel.f_model_scaled_with_k1(),
        weight_parameter = self.shelx_weight_parameter)
      result = result.customized_copy(data = result.data()*sw)
    return result

class electron_density_map_llgi(object):
  """ LLGI-native counterpart to electron_density_map -- builds
  2mFo-DFc/mFo-DFc-style map coefficients using Feff and the LLGI
  target's own D/fom (see mmtbx.f_model.manager.
  map_calculation_helper_llgi's docstring for the full derivation and
  correspondence to the ordinary ML alpha/beta/fom), instead of Fobs and
  a fresh generic ML alpha/beta/fom fit.

  Deliberately a SEPARATE class from electron_density_map, not a
  target-aware branch inside it -- mirrors mmtbx.f_model.manager.
  r_work_llgi()'s own separation from r_work() (see
  llgi_r_factors_available()'s docstring): electron_density_map is
  relied on by many non-LLGI callers (ml, mlhl, twin targets) that
  assume it is always Fobs-based, so it is left completely untouched;
  this class is reached only by callers that explicitly ask for the
  LLGI-native map (see mmtbx.f_model.manager.map_coefficients_llgi()).

  Reuses fo_fc_scales/combine UNCHANGED (both already only depend on
  centric_flags/.size(), which are identical whether sourced from f_obs
  or Feff, since Feff is always index-matched to f_obs -- see
  llgi_r_factors_available()'s own docstring) -- only the
  map_calculation_helper passed in differs (map_calculation_helper_llgi()
  instead of map_calculation_helper()).

  Special-case maps (anomalous, anomalous_residual, phaser_sad_llg,
  Fcalc-only) and use_shelx_weight are NOT supported here -- none of
  those are meaningful extensions of the LLGI target's own likelihood
  (anomalous/SAD analysis and SHELX weighting are independent features
  layered on the classic ML machinery, not yet ported to LLGI); raises
  NotImplementedError if requested, rather than silently falling back to
  an F-obs-based answer.

  fill_missing/fill_missing_method: like electron_density_map's own
  fill_missing/fill_missing_method, but only "f_model" (the default for
  BOTH classes) is implemented here, via fill_missing_f_obs_llgi/
  model_missing_reflections_llgi (Dobs*sigmaA*Emodel*sqrt(TEPS)*RESN,
  the LLGI-native DFc-analog, in place of D*Fc -- see that class's own
  docstring). "dsf"/"resolve_dm" raise NotImplementedError here (they
  are density-based, not model/target-based, so they would work
  unchanged, but are not yet wired through this class -- reachable via
  the ordinary electron_density_map path instead for now).
  """

  def __init__(self, fmodel):
    self.fmodel = fmodel
    self.mch = None

  def map_coefficients(self,
                       map_type,
                       acentrics_scale = 2.0,
                       centrics_pre_scale = 1.0,
                       exclude_free_r_reflections=False,
                       fill_missing=False,
                       fill_missing_method="f_model",
                       isotropize=True,
                       sharp=False):
    map_name_manager = mmtbx.map_names(map_name_string = map_type)
    if(map_name_manager.anomalous or map_name_manager.anomalous_residual
       or map_name_manager.phaser_sad_llg):
      raise NotImplementedError(
        "electron_density_map_llgi does not support anomalous/"
        "anomalous_residual/phaser_sad_llg map types; these features "
        "have not been ported to the LLGI target.")
    if(fill_missing and fill_missing_method not in ("f_model", None, False)):
      raise NotImplementedError(
        "electron_density_map_llgi only supports fill_missing_method="
        "'f_model' (Dobs*sigmaA*Emodel*sqrt(TEPS)*RESN); '%s' is not "
        "yet wired through the LLGI path." % fill_missing_method)
    mnm = mmtbx.map_names(map_name_string = map_type)
    if(mnm.k==0 and abs(mnm.n)==1):
      # Fcalc-only map: no observed-amplitude dependence at all, so the
      # F-obs vs Feff distinction is moot -- reuse electron_density_map
      # unchanged rather than duplicating this special case.
      return electron_density_map(fmodel=self.fmodel).map_coefficients(
        map_type=map_type, fill_missing=fill_missing)
    if(self.mch is None):
      self.mch = self.fmodel.map_calculation_helper_llgi()
    ffs = fo_fc_scales(
      fmodel          = self.fmodel,
      map_type_str    = map_type,
      acentrics_scale = acentrics_scale,
      centrics_scale  = centrics_pre_scale)
    fo_scale, fc_scale = ffs.fo_scale, ffs.fc_scale
    coeffs = combine(
      fmodel                 = self.fmodel,
      map_type_str           = map_type,
      fo_scale               = fo_scale,
      fc_scale               = fc_scale,
      map_calculation_helper = self.mch,
      use_shelx_weight       = False,
      shelx_weight_parameter = None).map_coefficients()
    r_free_flags = None
    scale_default = 1. / (self.fmodel.k_isotropic()*self.fmodel.k_anisotropic())
    scale_array = coeffs.customized_copy(data=scale_default)
    if (exclude_free_r_reflections):
      if (coeffs.anomalous_flag()):
        coeffs = coeffs.average_bijvoet_mates()
      r_free_flags = self.fmodel.r_free_flags()
      if (r_free_flags.anomalous_flag()):
        r_free_flags = r_free_flags.average_bijvoet_mates()
        scale_array = scale_array.average_bijvoet_mates()
      coeffs = coeffs.select(~r_free_flags.data())
      scale_array = scale_array.select(~r_free_flags.data())
    if(isotropize):
      if (scale_array.anomalous_flag()) and (not coeffs.anomalous_flag()):
        scale_array = scale_array.average_bijvoet_mates()
      scale = scale_array.data()
      coeffs = coeffs.customized_copy(data = coeffs.data()*scale)
    if(fill_missing):
      if(coeffs.anomalous_flag()):
        coeffs = coeffs.average_bijvoet_mates()
      coeffs = fill_missing_f_obs_llgi(coeffs=coeffs, fmodel=self.fmodel)
    if(sharp):
      ss = 1./flex.pow2(coeffs.d_spacings().data()) / 4.
      from cctbx import adptbx
      b = flex.mean(self.fmodel.xray_structure.extract_u_iso_or_u_equiv() *
        adptbx.u_as_b(1))/2
      k_sharp = 1./flex.exp(-ss * b)
      coeffs = coeffs.customized_copy(data = coeffs.data()*k_sharp)
    return coeffs

  def fft_map(self,
              resolution_factor = 1/3.,
              symmetry_flags = None,
              map_coefficients = None,
              map_type = None,
              acentrics_scale = 2.0,
              centrics_pre_scale = 1.0):
    if(map_coefficients is None):
      map_coefficients = self.map_coefficients(
        map_type           = map_type,
        acentrics_scale    = acentrics_scale,
        centrics_pre_scale = centrics_pre_scale)
    return map_coefficients.fft_map(
      resolution_factor = resolution_factor,
      symmetry_flags    = symmetry_flags)

class electron_density_map(object):

  def __init__(self,
               fmodel,
               map_calculation_helper = None):
    self.fmodel = fmodel
    self.anom_diff = None
    self.mch = map_calculation_helper
    if(self.fmodel.f_obs().anomalous_flag()):
      self.anom_diff = self.fmodel.f_obs().anomalous_differences()
      f_model = self.fmodel.f_model().as_non_anomalous_array().\
        merge_equivalents().array()
      fmodel_match_anom_diff, anom_diff_common = \
        f_model.common_sets(other =  self.anom_diff)
      assert anom_diff_common.indices().size()==self.anom_diff.indices().size()
      self.anom_diff = anom_diff_common.phase_transfer(
        phase_source = fmodel_match_anom_diff)

  def map_coefficients(self,
                       map_type,
                       acentrics_scale = 2.0,
                       centrics_pre_scale = 1.0,
                       exclude_free_r_reflections=False,
                       fill_missing=False,
                       fill_missing_method="f_model",
                       isotropize=True,
                       sharp=False,
                       pdb_hierarchy=None, # XXX required for map_type=llg
                       merge_anomalous=None,
                       use_shelx_weight=False,
                       shelx_weight_parameter=1.5):
    map_name_manager = mmtbx.map_names(map_name_string = map_type)
    # Special case #1: anomalous map
    if(map_name_manager.anomalous):
      if(self.anom_diff is not None):
        # Formula from page 141 in "The Bijvoet-Difference Fourier Synthesis",
        # Jeffrey Roach, METHODS IN ENZYMOLOGY, VOL. 374
        return miller.array(miller_set = self.anom_diff,
                            data       = self.anom_diff.data()/(2j))
      else: return None
    # Special case #2: anomalous residual map
    elif (map_name_manager.anomalous_residual):
      if (self.anom_diff is not None):
        return anomalous_residual_map_coefficients(
          fmodel=self.fmodel,
          exclude_free_r_reflections=exclude_free_r_reflections)
      else : return None
    # Special case #3: Phaser SAD LLG map
    elif (map_name_manager.phaser_sad_llg):
      if (pdb_hierarchy is None):
        raise RuntimeError("pdb_hierarchy must not be None when a Phaser SAD "+
          "LLG map is requested.")
      if (self.anom_diff is not None):
        return get_phaser_sad_llg_map_coefficients(
          fmodel=self.fmodel,
          pdb_hierarchy=pdb_hierarchy)
      else :
        return None
    # Special case #4: Fcalc map
    mnm = mmtbx.map_names(map_name_string = map_type)
    if(mnm.k==0 and abs(mnm.n)==1):
      if(fill_missing):
        return self.fmodel.xray_structure.structure_factors(
          d_min = self.fmodel.f_obs().d_min()).f_calc()
      else:
        return self.fmodel.f_obs().structure_factors_from_scatterers(
          xray_structure = self.fmodel.xray_structure).f_calc()
    #
    if(self.mch is None):
      self.mch = self.fmodel.map_calculation_helper()
    ffs = fo_fc_scales(
      fmodel          = self.fmodel,
      map_type_str    = map_type,
      acentrics_scale = acentrics_scale,
      centrics_scale  = centrics_pre_scale)
    fo_scale, fc_scale = ffs.fo_scale, ffs.fc_scale
    coeffs = combine(
      fmodel                 = self.fmodel,
      map_type_str           = map_type,
      fo_scale               = fo_scale,
      fc_scale               = fc_scale,
      map_calculation_helper = self.mch,
      use_shelx_weight       = use_shelx_weight,
      shelx_weight_parameter = shelx_weight_parameter).map_coefficients()
    r_free_flags = None
    # XXX the default scale array (used for the isotropize option) needs to be
    # calculated and processed now to avoid array size errors
    scale_default = 1. / (self.fmodel.k_isotropic()*self.fmodel.k_anisotropic())
    scale_array = coeffs.customized_copy(data=scale_default)
    if (exclude_free_r_reflections):
      if (coeffs.anomalous_flag()):
        coeffs = coeffs.average_bijvoet_mates()
      r_free_flags = self.fmodel.r_free_flags()
      if (r_free_flags.anomalous_flag()):
        r_free_flags = r_free_flags.average_bijvoet_mates()
        scale_array = scale_array.average_bijvoet_mates()
      coeffs = coeffs.select(~r_free_flags.data())
      scale_array = scale_array.select(~r_free_flags.data())
    scale=None
    if(isotropize):
      if(scale is None):
        if (scale_array.anomalous_flag()) and (not coeffs.anomalous_flag()):
          scale_array = scale_array.average_bijvoet_mates()
        scale = scale_array.data()
      coeffs = coeffs.customized_copy(data = coeffs.data()*scale)
    if(fill_missing):
      if(coeffs.anomalous_flag()):
        coeffs = coeffs.average_bijvoet_mates()
      coeffs = fill_missing_f_obs(
        coeffs = coeffs,
        fmodel = self.fmodel,
        method = fill_missing_method)
    if(sharp):
      ss = 1./flex.pow2(coeffs.d_spacings().data()) / 4.
      from cctbx import adptbx
      b = flex.mean(self.fmodel.xray_structure.extract_u_iso_or_u_equiv() *
        adptbx.u_as_b(1))/2
      k_sharp = 1./flex.exp(-ss * b)
      coeffs = coeffs.customized_copy(data = coeffs.data()*k_sharp)
    if (merge_anomalous) and (coeffs.anomalous_flag()):
      return coeffs.average_bijvoet_mates()
    return coeffs

  def fft_map(self,
              resolution_factor = 1/3.,
              symmetry_flags = None,
              map_coefficients = None,
              other_fft_map = None,
              map_type = None,
              force_anomalous_flag_false = None,
              acentrics_scale = 2.0,
              centrics_pre_scale = 1.0,
              use_all_data = True):
    if(map_coefficients is None):
      map_coefficients = self.map_coefficients(
        map_type           = map_type,
        acentrics_scale    = acentrics_scale,
        centrics_pre_scale = centrics_pre_scale)
      if(force_anomalous_flag_false):
        map_coefficients = map_coefficients.average_bijvoet_mates()
    if(force_anomalous_flag_false):
      map_coefficients = map_coefficients.average_bijvoet_mates()
    if(not use_all_data):
      map_coefficients = map_coefficients.select(self.fmodel.arrays.work_sel)
    if(other_fft_map is None):
      return map_coefficients.fft_map(
        resolution_factor = resolution_factor,
        symmetry_flags    = symmetry_flags)
    else:
      return miller.fft_map(
        crystal_gridding     = other_fft_map,
        fourier_coefficients = map_coefficients)

def resolve_dm_map(
      fmodel,
      map_coeffs,
      xrs,
      use_model_hl,
      fill,
      solvent_content=None,
      solvent_content_attenuator=0.1,
      mask_cycles  = 2,
      minor_cycles = 2,
      mask_from_model = None,  # use xrs as model_mask in density modification
      add_mask=None,  # alternative to mask_from_model: use xrs as added mask
      denmod_with_model = True, # use xrs in density modification
      input_text   = None):
  """
  Compute Resolve DM map
  """
  input_text="""
keep_missing
"""
  from solve_resolve.resolve_python import density_modify_in_memory
  if(solvent_content is None):
    import mmtbx.utils
    solvent_content = mmtbx.utils.f_000(xray_structure =
      fmodel.xray_structure).solvent_fraction-solvent_content_attenuator
  else:
    solvent_content = solvent_content-solvent_content_attenuator
  hl_model = None
  if(use_model_hl):
    hl_model = miller.set(crystal_symmetry=fmodel.f_obs().crystal_symmetry(),
      indices = fmodel.f_obs().indices(),
      anomalous_flag=False).array(
        data=fmodel.f_model_phases_as_hl_coefficients(
          map_calculation_helper=None))
  f_obs = fmodel.f_obs()
  f_obs = f_obs.array(
      data = f_obs.data(),
      sigmas = flex.double(f_obs.indices().size(), 1.0))  # must be present
  if(fill):
    data = f_obs.data()
    sigmas = f_obs.sigmas()
    complete_set = f_obs.complete_set()
    n_missing =  complete_set.indices().size() - f_obs.indices().size()
    complete_set = complete_set.array(
      data = flex.double(complete_set.indices().size(), 0.01), # must be > 0.0
      sigmas = flex.double(complete_set.indices().size(), -1.0))  # must be -1.0
    f_obs = f_obs.complete_with(other=complete_set)
  if(add_mask or mask_from_model):
    assert not (add_mask and mask_from_model)
    model_xrs = xrs # add_mask using xrs to define it. Could instead supply
                    #  a different xrs to define add_mask or model_mask
  else:
    model_xrs = None
  cmn=density_modify_in_memory.run(
    fp_sigfp            = f_obs.deep_copy(),
    hendrickson_lattman = hl_model,
    rad_mask            = max(2.5, f_obs.d_min()),
    map_coeffs_start    = map_coeffs,
    solvent_content     = solvent_content,
    mask_cycles         = mask_cycles,
    minor_cycles        = minor_cycles,
    xrs                 = xrs,
    model_xrs           = model_xrs,
    denmod_with_model   = denmod_with_model,
    add_mask            = add_mask,
    mask_from_model     = mask_from_model,
    verbose             = False,
    input_text          = input_text,
    out                 = null_out())
  result = cmn.map_coeffs_out_as_miller_array
  assert f_obs.indices().size()==result.indices().size()
  return result

class model_missing_reflections(object):
  def __init__(
        self,
        fmodel,
        coeffs):
    # XXX see f_model.py: duplication! Consolidate.
    self.fmodel = fmodel
    self.coeffs = coeffs
    crystal_gridding = fmodel.f_obs().crystal_gridding(
      d_min              = self.fmodel.f_obs().d_min(),
      resolution_factor  = 1./3)
    fft_map = miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = self.coeffs)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()
    rho_atoms = flex.double()
    for site_frac in self.fmodel.xray_structure.sites_frac():
      rho_atoms.append(map_data.eight_point_interpolation(site_frac))
    rho_mean = flex.mean_default(rho_atoms.select(rho_atoms>0.5), 0.5)
    sel_exclude = rho_atoms > min(rho_mean/2., 1)
    sites_cart = fmodel.xray_structure.sites_cart()
    #
    fft_map = miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = self.fmodel.f_model())
    fft_map.apply_sigma_scaling()
    map_data2 = fft_map.real_map_unpadded()
    #
    for i_seq, site_cart in enumerate(sites_cart):
      selection = maptbx.grid_indices_around_sites(
        unit_cell  = self.coeffs.unit_cell(),
        fft_n_real = map_data.focus(),
        fft_m_real = map_data.all(),
        sites_cart = flex.vec3_double([site_cart]),
        site_radii = flex.double([1.5]))
      cc = flex.linear_correlation(x=map_data.select(selection),
        y=map_data2.select(selection)).coefficient()
      if(cc<0.7): sel_exclude[i_seq] = False
    #
    del map_data, fft_map, rho_atoms
    self.d_min = fmodel.f_obs().d_min()
    cs = self.fmodel.f_obs().average_bijvoet_mates().complete_set(d_min=self.d_min)
    self.complete_set = cs.array(data = flex.double(cs.indices().size(), 0))
    self.xray_structure_cut = self.fmodel.xray_structure.select(sel_exclude)
    self.missing_set = self.complete_set.common_set(self.coeffs)
    #
    self.f_calc_missing = self.complete_set.structure_factors_from_scatterers(
      xray_structure = self.xray_structure_cut).f_calc()
    self.ss_missing = 1./flex.pow2(self.f_calc_missing.d_spacings().data()) / 4.
    mask_manager = mmtbx.masks.manager(
      miller_array      = self.f_calc_missing,
      miller_array_twin = None,
      mask_params       = None)
    self.f_mask_missing = mask_manager.shell_f_masks(
      xray_structure = self.xray_structure_cut,
      force_update   = True)
    self.zero_data = flex.complex_double(self.f_calc_missing.data().size(), 0)

  def get_missing_fast(self):
    k_sol = random.choice([i/100. for i in range(20,41)])
    b_sol = random.choice([i for i in range(20,85, 5)])
    f_calc = self.kick.randomize_struture_factors(
      map_coeffs=self.f_calc_missing, number_of_kicks=10)
    data = mmtbx_f_model_ext.core(
      f_calc        = f_calc.data(),
      shell_f_masks = [self.f_mask_missing[0].data(),],
      k_sols        = [k_sol,],
      b_sol         = b_sol,
      f_part1       = self.zero_data,
      f_part2       = self.zero_data,
      u_star        = [0,0,0,0,0,0],
      hkl           = self.f_calc_missing.indices(),
      uc            = self.f_calc_missing.unit_cell(),
      ss            = self.ss_missing)
    return miller.array(miller_set=self.f_calc_missing, data=data.f_model)

  def get_missing(self, deterministic=False):
    if(deterministic):
      xrs = self.xray_structure_cut
    else:
      sel = flex.random_bool(self.xray_structure_cut.scatterers().size(), 0.9)
      xrs = self.xray_structure_cut.select(sel)
    if(deterministic):
      # XXX This may result in very bad maps at low resolution
      #fm = mmtbx.f_model.manager(
      #  f_obs          = self.complete_set,
      #  xray_structure = xrs,
      #  k_sol          = 0.35,
      #  b_sol          = 46.0)
      #result = fm.f_model_no_scales()
      # This is slower but better
      # Never change this unless checked with CK case
      fm = self.fmodel.deep_copy()
      r = fm.update_all_scales(fast=False, refine_hd_scattering=False)
      fm = mmtbx.f_model.manager(
        f_obs          = self.complete_set,
        xray_structure = xrs,
        k_sol          = r.k_sol[0],
        b_sol          = r.b_sol[0],
        b_cart         = r.b_cart)
      result = fm.f_model_no_scales()
    else:
      if(random.choice([True, False])):
        k_sol = random.choice([i/100. for i in range(20,41)])
        b_sol = random.choice([i for i in range(20,85, 5)])
        fm = mmtbx.f_model.manager(
          f_obs          = self.complete_set,
          xray_structure = xrs,
          k_sol          = k_sol,
          b_sol          = b_sol)
        result = fm.f_model_no_scales()
      else:
        result = xrs.structure_factors(d_min = self.d_min).f_calc()
    return result

def fill_missing_f_obs_1(coeffs, fmodel):
  mro = model_missing_reflections(coeffs=coeffs, fmodel=fmodel)
  missing = mro.get_missing(deterministic=True)
  return coeffs.complete_with(other = missing, scale=True)

def fill_missing_f_obs_2(coeffs, fmodel):
  scale_to = fmodel.f_obs().average_bijvoet_mates()
  dsf = coeffs.double_step_filtration(
    vol_cutoff_plus_percent=5.0,
    vol_cutoff_minus_percent=5.0,
    scale_to=scale_to)
  return coeffs.complete_with(other = dsf, scale=True)

def fill_missing_f_obs_3(coeffs, fmodel):
  if(not libtbx.env.has_module("solve_resolve")): return coeffs
  mc_dm_filled = resolve_dm_map(
    fmodel       = fmodel,
    map_coeffs   = coeffs,
    pdb_inp      = None,
    use_model_hl = True,
    fill         = True)
  return coeffs.complete_with(other = mc_dm_filled, scale=True)

def fill_missing_f_obs(coeffs, fmodel, method):
  assert method in ["resolve_dm", "dsf", "f_model", None, False]
  if(method == "f_model"):
    return fill_missing_f_obs_1(coeffs=coeffs, fmodel=fmodel)
  elif(method == "dsf"):
    return fill_missing_f_obs_2(coeffs=coeffs, fmodel=fmodel)
  elif(method == "resolve_dm"):
    return fill_missing_f_obs_3(coeffs=coeffs, fmodel=fmodel)
  elif(method in [None, False]):
    return coeffs
  else:
    raise RuntimeError("Invalid arg of fill_missing_f_obs: method:"%(method))

class model_missing_reflections_llgi(object):
  """ LLGI-native counterpart to model_missing_reflections -- fills
  missing/unmeasured reflections with a DFc-analog term built from
  Dobs*sigmaA*Emodel*sqrt(TEPS)*RESN (E-scale, anisotropy-free), rather
  than model_missing_reflections' own D*Fc (ML alpha-weighted, F-scale).
  See mmtbx.f_model.manager.map_calculation_helper_llgi's docstring for
  the same D/Emodel/RESN quantities used elsewhere in the LLGI map-
  coefficient/sigmaA machinery, reused UNCHANGED here.

  Reuses model_missing_reflections' own __init__ entirely (the atom-
  correlation filtering, complete_set/xray_structure_cut/f_calc_missing/
  f_mask_missing construction, is genuinely target-agnostic -- it only
  needs an fmodel and the CURRENT map coefficients being completed, not
  anything ML-specific) via composition rather than duplicating it,
  since model_missing_reflections doesn't currently expose a hook for a
  subclass to override just the "how do we weight/value a missing
  reflection" step.

  A missing reflection has, by definition, no actual measurement, so
  three of the quantities the LLGI target normally needs at each
  reflection (Dobs, RESN, and -- less obviously -- k_isotropic itself,
  since fmodel.k_isotropic() is a plain per-reflection array from bss's
  own fit, not a closed-form function of resolution) have no natural
  value there either. Each is handled differently, per what's actually
  knowable:
    Dobs   -- treated as 1 (full reliability): there is no actual
              measurement to distrust, so sigmaA(d) alone (evaluated at
              the missing reflection's own resolution, via
              e_sigmaa_target_evaluator.evaluate_at) stands in for the
              usual Dobs*sigmaA product -- mirrors how the ML fill
              itself applies NO alpha/fom weighting at all to its D*Fc
              term (see model_missing_reflections.get_missing), just a
              milder version of the same idea (weighted by model
              confidence at that resolution, not left completely
              unweighted).
    RESN   -- nacelle never computed this for missing reflections (it's
              inherently tied to the experimental data collection), and
              this module has no independent model of its resolution
              shape to extrapolate with -- each missing reflection
              borrows the RESN value of its NEAREST observed reflection
              in resolution (1D nearest-neighbour in d*^2), rather than
              fitting/assuming a smooth curve we have no basis for.
    k_iso  -- fmodel.k_isotropic() has no value at an unobserved
              resolution, but (unlike RESN) IS expected to be a smooth,
              well-behaved function of resolution alone (it is exactly
              this module's own SigmaP/Emodel machinery's premise) --
              fit a plain B-spline curve to the OBSERVED k_isotropic(ss)
              array (reusing mmtbx.refinement.llgi_e_bulk_solvent's own
              _b_spline_design_matrix, log-space to keep it positive,
              same convention as sigmaA's own z=log-space-ish
              parameterisation) and evaluate that fitted curve at each
              missing reflection's own ss, clamped to the observed
              range exactly as e_sigmaa_target_evaluator.evaluate_at
              clamps sigmaA.
  k_mask/k_sol/b_sol themselves ARE well-defined at any resolution
  (k_mask(ss) = k_sol*exp(-b_sol*ss) is a closed-form function by
  construction -- see mmtbx.refinement.llgi_e_bulk_solvent.k_mask_and_
  gradients) using fmodel.k_sol_b_sol_from_k_mask()'s current fit (the
  LIVE fmodel's own already-current bulk-solvent state -- NOT re-fit
  here, unlike the ML "careful" deterministic=True path's fresh, slow
  update_all_scales() re-fit: the LLGI macrocycle pipeline already keeps
  k_sol/b_sol current every cycle via run_inner_loop/estimate_e_sigmaa_
  fixed_bulk_solvent, so re-fitting here would be redundant work solving
  an already-solved problem).

  Requires fmodel.llgi_data() (FEFF/DOBS/TEPS/RESN) AND llgi_data.sigmaa
  already attached (same precondition as map_calculation_helper_llgi) --
  raises AttributeError otherwise.
  """

  def __init__(self, fmodel, coeffs):
    self.fmodel = fmodel
    self.coeffs = coeffs
    self._base = model_missing_reflections(fmodel=fmodel, coeffs=coeffs)
    llgi_data = fmodel.llgi_data()
    if(llgi_data is None):
      raise AttributeError(
        "model_missing_reflections_llgi requires llgi_data (FEFF/DOBS/"
        "TEPS/RESN) -- none attached.")
    sigmaa = getattr(llgi_data, "sigmaa", None)
    if(sigmaa is None):
      raise AttributeError(
        "model_missing_reflections_llgi requires sigmaa already "
        "attached to llgi_data (see update_llgi_sigmaa_scatfrac()) -- "
        "not yet available.")
    self.llgi_data = llgi_data

  def _fit_k_isotropic_curve(self):
    """ Fit a smooth B-spline curve to the OBSERVED fmodel.k_isotropic()
    array vs ss (log-space, so the fitted curve stays positive), for
    later evaluation at the missing reflections' own ss -- see this
    class's own docstring. Returns (coeffs, ss_range, n_coeffs, degree)
    ready for _eval_k_isotropic_curve.
    """
    import numpy as np
    from mmtbx.refinement.llgi_e_bulk_solvent import (
      ss_from_f_obs, _b_spline_design_matrix)
    f_obs = self.fmodel.f_obs()
    ss = ss_from_f_obs(f_obs).as_numpy_array()
    k_iso = np.asarray(self.fmodel.k_isotropic(), dtype=float)
    ss_range = (float(ss.min()), float(ss.max()))
    n_coeffs, degree = 8, 3
    design = _b_spline_design_matrix(ss, n_coeffs, degree, x_range=ss_range)
    log_k_iso = np.log(np.clip(k_iso, 1.e-12, None))
    coeffs, _res, _rank, _sv = np.linalg.lstsq(design, log_k_iso, rcond=None)
    return coeffs, ss_range, n_coeffs, degree

  def _eval_k_isotropic_curve(self, ss_missing, fit):
    import numpy as np
    from mmtbx.refinement.llgi_e_bulk_solvent import _b_spline_design_matrix
    coeffs, ss_range, n_coeffs, degree = fit
    ss_low, ss_high = ss_range
    ss_np = np.clip(np.asarray(ss_missing, dtype=float), ss_low, ss_high)
    design = _b_spline_design_matrix(
      ss_np, n_coeffs, degree, x_range=ss_range)
    log_k_iso = design.dot(coeffs)
    return flex.double(np.exp(log_k_iso).tolist())

  def _nearest_resn(self, d_star_sq_missing):
    """ RESN for each missing reflection, borrowed from its nearest
    OBSERVED reflection in resolution (1D nearest-neighbour in d*^2) --
    see this class's own docstring on why RESN is not fit/extrapolated
    the way sigmaA/k_isotropic are.
    """
    import numpy as np
    f_obs = self.fmodel.f_obs()
    d_star_sq_obs = np.asarray(f_obs.d_star_sq().data(), dtype=float)
    resn_obs = np.asarray(self.llgi_data.resn.data(), dtype=float)
    order = np.argsort(d_star_sq_obs)
    d_sorted = d_star_sq_obs[order]
    resn_sorted = resn_obs[order]
    d_missing_np = np.asarray(d_star_sq_missing, dtype=float)
    idx = np.searchsorted(d_sorted, d_missing_np)
    idx = np.clip(idx, 0, len(d_sorted) - 1)
    idx_prev = np.clip(idx - 1, 0, len(d_sorted) - 1)
    use_prev = (
      (idx == len(d_sorted)) |
      (np.abs(d_sorted[idx_prev] - d_missing_np) <
       np.abs(d_sorted[idx] - d_missing_np)))
    nearest_idx = np.where(use_prev, idx_prev, idx)
    return flex.double(resn_sorted[nearest_idx].tolist())

  def get_missing(self):
    """ Dobs*sigmaA*Emodel*sqrt(TEPS)*RESN (Dobs=1, see this class's own
    docstring) for the missing-reflection set, on the SAME index set/
    order as self._base.f_calc_missing (mirrors model_missing_
    reflections.get_missing's own return shape, so fill_missing_f_obs_
    llgi can use it identically to fill_missing_f_obs_1).
    """
    import mmtbx.refinement.llgi_e_bulk_solvent as llgi_e_bs
    base = self._base
    f_calc_missing = base.f_calc_missing
    f_mask_missing = base.f_mask_missing[0]
    ss_missing = flex.double(base.ss_missing)
    d_star_sq_missing = f_calc_missing.d_star_sq().data()
    epsilons_missing = f_calc_missing.epsilons().data().as_double()

    k_sol, b_sol = self.fmodel.k_sol_b_sol_from_k_mask()
    k_mask_result = llgi_e_bs.k_mask_and_gradients(ss_missing, k_sol, b_sol)
    k_iso_fit = self._fit_k_isotropic_curve()
    k_iso_missing = self._eval_k_isotropic_curve(ss_missing, k_iso_fit)

    fmnas_missing = k_iso_missing * (
      f_calc_missing.data() + k_mask_result.k_mask * f_mask_missing.data())

    sigma_p_missing = llgi_e_bs.build_sigma_p(
      llgi_e_bs.f_model_no_aniso_scale(self.fmodel).data(),
      self.fmodel.f_obs().d_star_sq().data(),
      d_star_sq_eval=d_star_sq_missing)
    e_model_missing = fmnas_missing * (
      1.0 / flex.sqrt(epsilons_missing * sigma_p_missing))
    e_model_abs_missing = flex.abs(e_model_missing)

    sigmaa_result = self.llgi_data.sigmaa
    # sigmaa_result here is the miller.array attached to llgi_data (an
    # already-evaluated curve, not the fitted evaluator itself) -- the
    # fitted evaluator/x_range needed for evaluate_at() come from a
    # FRESH sigmaA fit against the CURRENT Emodel instead (see below),
    # since llgi_data.sigmaa's own fit is not guaranteed still current
    # against this exact fmodel state, and evaluate_at is only exposed
    # on estimate_e_sigmaa's own OWN return value, not persisted
    # anywhere on llgi_data itself.
    f_obs = self.fmodel.f_obs()
    epsilons_obs = f_obs.epsilons().data().as_double()
    d_star_sq_obs = f_obs.d_star_sq().data()
    fmnas_obs = llgi_e_bs.f_model_no_aniso_scale(self.fmodel).data()
    e_model_obs = llgi_e_bs.build_e_model(fmnas_obs, epsilons_obs, d_star_sq_obs)
    e_eff_obs = llgi_e_bs.build_e_eff(
      self.llgi_data.feff.data(), self.llgi_data.resn.data())
    sigmaa_refit = llgi_e_bs.estimate_e_sigmaa(
      e_eff=e_eff_obs, r_free_flags=self.fmodel.r_free_flags().data(),
      e_model=flex.abs(e_model_obs.e_model), dobs=self.llgi_data.dobs.data(),
      centric_flags=f_obs.centric_flags().data(), d_star_sq=d_star_sq_obs)
    sigmaa_missing = sigmaa_refit.evaluate_at(
      d_star_sq_missing, x_range=sigmaa_refit.x_range)

    resn_missing = self._nearest_resn(d_star_sq_missing)
    teps_missing = flex.double(d_star_sq_missing.size(), 1.0)  # TEPS==1 only
    sqrt_teps_resn = flex.sqrt(teps_missing) * resn_missing

    fill_data = sigmaa_missing * e_model_abs_missing * sqrt_teps_resn
    # Phase from fmnas_missing (real-valued magnitude above needs a
    # phase source -- matches map_calculation_helper_llgi's own
    # .f_model = f_model_no_aniso_scale convention, whose phase equals
    # fmodel.f_model()'s own phase since k_anisotropic is real/positive).
    fill_complex = miller.array(
      miller_set=f_calc_missing, data=fill_data).phase_transfer(
        phase_source=fmnas_missing).data()
    return f_calc_missing.customized_copy(data=fill_complex)

def fill_missing_f_obs_llgi(coeffs, fmodel):
  """ LLGI-native counterpart to fill_missing_f_obs_1 -- fills missing
  reflections with Dobs*sigmaA*Emodel*sqrt(TEPS)*RESN (see
  model_missing_reflections_llgi) instead of D*Fc. The ONLY LLGI fill-
  missing method implemented (mirrors fill_missing_f_obs_1 being the
  default/"f_model" ML method) -- the "dsf"/"resolve_dm" alternatives
  are density-based, not model/target-based, and already work unchanged
  on any map coefficients regardless of ML/LLGI provenance, so they stay
  reachable only through the ordinary (non-LLGI) fill_missing_f_obs
  dispatch for now.
  """
  mro = model_missing_reflections_llgi(coeffs=coeffs, fmodel=fmodel)
  missing = mro.get_missing()
  return coeffs.complete_with(other=missing, scale=True)

def sharp_evaluation_target(sites_frac, map_coeffs, resolution_factor = 0.25):
  fft_map = map_coeffs.fft_map(resolution_factor=resolution_factor)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  target = 0
  cntr = 0
  for site_frac in sites_frac:
    mv = map_data.eight_point_interpolation(site_frac)
    if(mv >0):
      target += mv
      cntr += 1
  t_mean = target/cntr
  return t_mean

def sharp_map(sites_frac, map_coeffs, ss = None, b_sharp=None, b_min = -150,
              b_max = 150, step = 10):
  if(ss is None):
    ss = 1./flex.pow2(map_coeffs.d_spacings().data()) / 4.
  from cctbx import miller
  if(b_sharp is None):
    t=-1
    map_coeffs_best = None
    b_sharp_best = None
    for b_sharp in range(b_min,b_max,step):
      map_coeffs_ = map_coeffs.deep_copy()
      sc2 = flex.exp(b_sharp*ss)
      map_coeffs_ = map_coeffs_.customized_copy(data = map_coeffs_.data()*sc2)
      t_=sharp_evaluation_target(sites_frac=sites_frac, map_coeffs=map_coeffs_)
      if(t_>t):
        t=t_
        b_sharp_best = b_sharp
        map_coeffs_best = map_coeffs_.deep_copy()
    print("b_sharp:", b_sharp_best, t)
  else:
    scale = flex.exp(b_sharp*ss)
    map_coeffs_best = map_coeffs.customized_copy(data=map_coeffs.data()*scale)
    b_sharp_best = b_sharp
  return map_coeffs_best, b_sharp_best

def get_phaser_sad_llg_map_coefficients(
    fmodel,
    pdb_hierarchy,
    log=None,
    verbose=False):
  """
  Calculates an anomalous log-likelihood gradient (LLG) map using the SAD
  target in Phaser.  This is essentially similar to an anomalous difference-
  difference map, but more sensitive.
  """
  if (not libtbx.env.has_module("phaser")):
    raise Sorry("Phaser not available - required for SAD LLG maps.")
  from phaser.phenix_adaptors import sad_target
  assert (fmodel.f_model().anomalous_flag())
  f_obs = fmodel.f_obs().select(fmodel.f_obs().data()>0)
  r_free_flags = fmodel.r_free_flags().common_set(other=f_obs)
  data = sad_target.data_adaptor(
    f_obs=f_obs,
    r_free_flags=r_free_flags,
    verbose=True)
  if (verbose) and (log is not None):
    data.output.setPackagePhenix(log)
  else :
    data.output.setPackagePhenix(null_out())
  t = data.target(
    xray_structure=fmodel.xray_structure,
    pdb_hierarchy=pdb_hierarchy,
    log=null_out())
  t.set_f_calc(fmodel.f_model())
  map_coeffs = t.llg_map_coeffs()
  return map_coeffs

def anomalous_residual_map_coefficients(fmodel, weighted=False,
    exclude_free_r_reflections=True):
  """
  EXPERIMENTAL

  Calculates map coefficients showing the difference in anomalous scattering
  between F-obs and F-model.  Similar to the Phaser SAD LLG map, but appears
  to be less sensitive.
  """
  assert (fmodel.f_obs().anomalous_flag())
  f_obs_anom = fmodel.f_obs().anomalous_differences()
  f_model_anom = abs(fmodel.f_model()).anomalous_differences()
  if (weighted):
    mch = fmodel.map_calculation_helper()
    fom = fmodel.f_obs().customized_copy(
      data=mch.fom, sigmas=None).average_bijvoet_mates()
    alpha = mch.alpha.average_bijvoet_mates()
    fom = fom.common_set(other=f_obs_anom)
    alpha = alpha.common_set(other=f_obs_anom)
    f_obs_anom = f_obs_anom.customized_copy(data=f_obs_anom.data()*fom.data())
    f_model_anom = f_model_anom.customized_copy(
      data=f_model_anom.data()*alpha.data())
  anom_diff_diff = f_obs_anom.customized_copy(
    data=f_obs_anom.data() - f_model_anom.data())
  f_model = fmodel.f_model().as_non_anomalous_array().\
    merge_equivalents().array()
  fmodel_match_anom_diff, anom_diff_diff_common = \
    f_model.common_sets(other =  anom_diff_diff)
  assert (anom_diff_diff_common.indices().size() ==
          anom_diff_diff.indices().size())
  phases_tmp = miller.array(
    miller_set=anom_diff_diff_common,
    data=flex.double(anom_diff_diff_common.indices().size(), 1)
    ).phase_transfer(phase_source=fmodel_match_anom_diff)
  map_coeffs = miller.array(
    miller_set=anom_diff_diff_common,
    data = anom_diff_diff_common.data() * phases_tmp.data())
  if (exclude_free_r_reflections):
    r_free_flags = fmodel.r_free_flags().average_bijvoet_mates()
    r_free_flags, map_coeffs = r_free_flags.common_sets(map_coeffs)
    map_coeffs = map_coeffs.select(~(r_free_flags.data()))
  return miller.array(
    miller_set=map_coeffs,
    data=map_coeffs.data()/(2j))
