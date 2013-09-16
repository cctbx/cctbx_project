from __future__ import division
from cctbx.array_family import flex
from cctbx import miller
from cctbx import maptbx
import libtbx.phil
from libtbx.utils import null_out
from libtbx import adopt_init_args
import mmtbx

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
               map_calculation_helper=None):
    self.mch = map_calculation_helper
    self.mnm = mmtbx.map_names(map_name_string = map_type_str)
    self.fmodel = fmodel
    self.fc_scale = fc_scale
    self.fo_scale = fo_scale
    self.f_obs = None
    self.f_model = None
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
    return compute(fo=self.f_obs, fc=f_model_data, miller_set=self.fmodel.f_obs())

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
                       ncs_average=False,
                       isotropize=True,
                       sharp=False,
                       post_processing_callback=None,
                       pdb_hierarchy=None, # XXX required for map_type=llg
                       merge_anomalous=None) :
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
    elif (map_name_manager.anomalous_residual) :
      if (self.anom_diff is not None) :
        return anomalous_residual_map_coefficients(
          fmodel=self.fmodel,
          exclude_free_r_reflections=exclude_free_r_reflections)
      else : return None
    # Special case #3: Phaser SAD LLG map
    elif (map_name_manager.phaser_sad_llg) :
      if (pdb_hierarchy is None) :
        raise RuntimeError("pdb_hierarchy must not be None when a Phaser SAD "+
          "LLG map is requested.")
      if (self.anom_diff is not None) :
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
      map_calculation_helper = self.mch).map_coefficients()
    r_free_flags = None
    # XXX the default scale array (used for the isotropize option) needs to be
    # calculated and processed now to avoid array size errors
    scale_default = 1. / (self.fmodel.k_isotropic()*self.fmodel.k_anisotropic())
    scale_array = coeffs.customized_copy(data=scale_default)
    if (exclude_free_r_reflections) :
      if (coeffs.anomalous_flag()) :
        coeffs = coeffs.average_bijvoet_mates()
      r_free_flags = self.fmodel.r_free_flags()
      if (r_free_flags.anomalous_flag()) :
        r_free_flags = r_free_flags.average_bijvoet_mates()
        scale_array = scale_array.average_bijvoet_mates()
      coeffs = coeffs.select(~r_free_flags.data())
      scale_array = scale_array.select(~r_free_flags.data())
    scale=None
    if (ncs_average) and (post_processing_callback is not None) :
      # XXX NCS averaging done here
      assert hasattr(post_processing_callback, "__call__")
      coeffs = post_processing_callback(
        map_coeffs=coeffs,
        fmodel=self.fmodel,
        map_type=map_type)
    if(isotropize):
      if(scale is None):
        if (scale_array.anomalous_flag()) and (not coeffs.anomalous_flag()) :
          scale_array = scale_array.average_bijvoet_mates()
        scale = scale_array.data()
      coeffs = coeffs.customized_copy(data = coeffs.data()*scale)
    if(fill_missing):
      if(coeffs.anomalous_flag()):
        coeffs = coeffs.average_bijvoet_mates()
      coeffs = fill_missing_f_obs(coeffs, self.fmodel)
    if(sharp):
      ss = 1./flex.pow2(coeffs.d_spacings().data()) / 4.
      from cctbx import adptbx
      b = flex.mean(self.fmodel.xray_structure.extract_u_iso_or_u_equiv() *
        adptbx.u_as_b(1))/2
      k_sharp = 1./flex.exp(-ss * b)
      coeffs = coeffs.customized_copy(data = coeffs.data()*k_sharp)
    if (merge_anomalous) and (coeffs.anomalous_flag()) :
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

def fill_missing_f_obs_1(coeffs, fmodel):
  cs = fmodel.f_obs().average_bijvoet_mates().complete_set(
    d_min = fmodel.f_obs().d_min())
  ca = cs.array(data = flex.double(cs.indices().size(), 0))
  fm = mmtbx.f_model.manager(
    f_obs = ca,
    xray_structure = fmodel.xray_structure,
    k_sol = 0.35,
    b_sol = 46.0)
  f_map = fm.f_model_no_scales()
  return coeffs.complete_with(other = f_map, scale=True)

def fill_missing_f_obs_2(coeffs, fmodel):
  scale_to = fmodel.f_obs().average_bijvoet_mates()
  dsf = coeffs.double_step_filtration(
    vol_cutoff_plus_percent=10.0,
    vol_cutoff_minus_percent=10.0,
    scale_to=scale_to)
  return coeffs.complete_with(other = dsf, scale=True)

def fill_missing_f_obs(coeffs, fmodel):
  if(fmodel.xray_structure is not None):
    return fill_missing_f_obs_1(coeffs=coeffs, fmodel=fmodel)
  else:
    return fill_missing_f_obs_2(coeffs=coeffs, fmodel=fmodel)

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
    print "b_sharp:", b_sharp_best, t
  else:
    scale = flex.exp(b_sharp*ss)
    map_coeffs_best = map_coeffs.customized_copy(data=map_coeffs.data()*scale)
    b_sharp_best = b_sharp
  return map_coeffs_best, b_sharp_best

def get_phaser_sad_llg_map_coefficients (
    fmodel,
    pdb_hierarchy,
    log=None,
    verbose=False) :
  """
  Calculates an anomalous log-likelihood gradient (LLG) map using the SAD
  target in Phaser.  This is essentially similar to an anomalous difference-
  difference map, but more sensitive.
  """
  if (not libtbx.env.has_module("phaser")) :
    raise Sorry("Phaser not available - required for SAD LLG maps.")
  from phaser.phenix_adaptors import sad_target
  assert (fmodel.f_model().anomalous_flag())
  f_obs = fmodel.f_obs().select(fmodel.f_obs().data()>0)
  r_free_flags = fmodel.r_free_flags().common_set(other=f_obs)
  data = sad_target.data_adaptor(
    f_obs=f_obs,
    r_free_flags=r_free_flags,
    verbose=True)
  if (verbose) and (log is not None) :
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

def anomalous_residual_map_coefficients (fmodel, weighted=False,
    exclude_free_r_reflections=True) :
  """
  EXPERIMENTAL

  Calculates map coefficients showing the difference in anomalous scattering
  between F-obs and F-model.  Similar to the Phaser SAD LLG map, but appears
  to be less sensitive.
  """
  assert (fmodel.f_obs().anomalous_flag())
  f_obs_anom = fmodel.f_obs().anomalous_differences()
  f_model_anom = abs(fmodel.f_model()).anomalous_differences()
  if (weighted) :
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
  if (exclude_free_r_reflections) :
    r_free_flags = fmodel.r_free_flags().average_bijvoet_mates()
    r_free_flags, map_coeffs = r_free_flags.common_sets(map_coeffs)
    map_coeffs = map_coeffs.select(~(r_free_flags.data()))
  return miller.array(
    miller_set=map_coeffs,
    data=map_coeffs.data()/(2j))

ncs_averaging_params = """
resolution_factor = 0.25
  .type = float
use_molecule_mask = False
  .type = bool
averaging_radius = 5.0
  .type = float
solvent_content = 0.5
  .type = float
exclude_hd = True
  .type = bool
skip_difference_map = Auto
  .type = bool
"""

# XXX it would be more useful to have this integrated with the rest of the
# code, instead of making map averaging an afterthought.  however, the
# external overhead is currently substantial.
class ncs_averager (object) :
  def __init__ (self, ncs_object, params=None, log=None, verbose=False) :
    if (params is None) :
      params = libtbx.phil.parse(ncs_averaging_params).extract()
    if (log is None) :
      log = null_out()
    adopt_init_args(self, locals())
    self.mask = None

  def __call__ (self,
                map_coeffs,
                fmodel,
                generate_new_mask=False,
                map_type=None) :
    # XXX probably not a good idea to average anomalous maps
    #if (map_type is not None) and (map_type.lower().startswith("anom")) :
    #  return map_coeffs
    from solve_resolve.resolve_python.resolve_utils import get_map_mask_sg_cell
    from solve_resolve.resolve_python.ncs_average import ncs_average
    from cctbx import maptbx
    from scitbx.array_family import flex
    if (map_coeffs.anomalous_flag()) :
      map_coeffs = map_coeffs.average_bijvoet_mates()
    fft_map = map_coeffs.fft_map(
      symmetry_flags=maptbx.use_space_group_symmetry,
      resolution_factor=self.params.resolution_factor)
    map = fft_map.apply_volume_scaling().real_map_unpadded().as_float()
    if (self.verbose) :
      out = self.log
    else :
      out = null_out()
    if (self.mask is None) or (generate_new_mask) :
      if (self.params.use_molecule_mask) :
        self.mask = flex.float(real_map.size(), 0)
        sites_cart = fmodel.xray_structure.sites_cart()
        if (self.params.exclude_hd) :
          sites_cart = sites_cart.select(~fmodel.xray_structure.hd_selection())
        indices = maptbx.grid_indices_around_sites(
          unit_cell=map_coeffs.unit_cell(),
          fft_n_real=real_map.focus(),
          fft_m_real=real_map.all(),
          sites_cart=sites_cart,
          site_radii=flex.double(sites_cart.size(),
            self.params.averaging_radius))
        mask.set_selected(indices, 1)
        mask.reshape(real_map.accessor())
      else :
        mask_map_coeffs = fmodel.electron_density_map(
          update_f_part1=True).map_coefficients(
            map_type="2mFo-DFc")
        mask_fft_map = mask_map_coeffs.fft_map(
          symmetry_flags=maptbx.use_space_group_symmetry,
          resolution_factor=self.params.resolution_factor)
        mask_map = mask_fft_map.apply_volume_scaling().real_map_unpadded().as_float()
        map_db,mask_map_db,space_group_object,unit_cell_object=\
          get_map_mask_sg_cell(
            map_coeffs=mask_map_coeffs,
            map=mask_map,
            space_group=map_coeffs.space_group(),
            unit_cell=map_coeffs.unit_cell(),
            solvent_content=self.params.solvent_content,
            wang_radius=self.params.averaging_radius,
            resolution=map_coeffs.d_min(),
            out=out,
            resolve_command_list=None)
        #map = map_db.map
        self.mask = mask_map_db.map
    averaged = ncs_average(
      map=map,
      mask=self.mask,
      ncs_object=self.ncs_object,
      space_group=map_coeffs.space_group(),
      unit_cell=map_coeffs.unit_cell(),
      resolution=map_coeffs.d_min(),
      out=out)
    new_map_coeffs = map_coeffs.structure_factors_from_map(
      map=averaged.average_map.as_double(),
      use_sg=True)
    return new_map_coeffs
