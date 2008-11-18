from cctbx import miller
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import xray
from cctbx import eltbx
import cctbx.xray.structure_factors
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, host_and_user, multi_out
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from iotbx.pdb import xray_structure
from iotbx import pdb
import libtbx.phil.command_line
from cStringIO import StringIO
from scitbx.python_utils import easy_pickle
from scitbx.math import matrix
from cctbx import adptbx
import sys, os, math
from mmtbx import monomer_library
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.monomer_library.server
import iotbx.phil
import libtbx.phil.command_line
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx.option_parser import iotbx_option_parser
from iotbx import crystal_symmetry_from_any
from iotbx import pdb
from iotbx.pdb import crystal_symmetry_from_pdb
from iotbx import mtz
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
from scitbx.python_utils.misc import user_plus_sys_time, show_total_time
from libtbx.str_utils import show_string
from libtbx.utils import Sorry, date_and_time, host_and_user, multi_out
from libtbx import adopt_init_args
import random, sys, os
from libtbx.test_utils import approx_equal
from mmtbx.refinement import print_statistics
import libtbx.load_env
from mmtbx import max_lik
from mmtbx.max_lik import maxlik
from scitbx import fftpack
from scitbx import matrix
from cctbx import maptbx
from mmtbx import masks
import boost.python
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss

ext = boost.python.import_ext("mmtbx_f_model_ext")

class map_helper(object):
  def __init__(self, fmodel, free_reflections_per_bin = 100,
               interpolation = True):
    ss = 1./flex.pow2(fmodel.f_obs.d_spacings().data())/4.
    fb_cart  = fmodel.fb_cart()
    scale_k1 = fmodel.scale_k1()
    f_obs_scale   = 1.0 / (fb_cart * scale_k1)
    f_model_scale = 1.0 / fb_cart
    f_obs_data_scaled = fmodel.f_obs.data() * f_obs_scale
    f_model_data_scaled = fmodel.f_model().data() * f_model_scale
    self.f_obs_scaled = fmodel.f_obs.array(data = f_obs_data_scaled)
    self.f_model_scaled = fmodel.f_obs.array(data = f_model_data_scaled)
    self.alpha, beta = maxlik.alpha_beta_est_manager(
      f_obs                    = self.f_obs_scaled,
      f_calc                   = self.f_model_scaled,
      free_reflections_per_bin = free_reflections_per_bin,
      flags                    = fmodel.r_free_flags.data(),
      interpolation            = interpolation).alpha_beta()
    self.fom = max_lik.fom_and_phase_error(
      f_obs          = f_obs_data_scaled,
      f_model        = flex.abs(f_model_data_scaled),
      alpha          = self.alpha.data(),
      beta           = beta.data(),
      space_group    = fmodel.r_free_flags.space_group(),
      miller_indices = fmodel.r_free_flags.indices()).fom()

class kick_map(object):

  def __init__(self, fmodel,
                     map_type,
                     number_of_kicks = 50,
                     update_bulk_solvent_and_scale = False,
                     resolution_factor = 1./4,
                     symmetry_flags = None,
                     other_fft_map = None,
                     real_map_unpadded = True,
                     real_map = False,
                     average_maps = False):
    assert [real_map_unpadded, real_map].count(True) == 1
    fmodel_tmp = fmodel.deep_copy()
    bss_params = bss.master_params.extract()
    map_helper_obj = map_helper(fmodel = fmodel)
    ss = 1./flex.pow2(fmodel.f_obs.d_spacings().data())/4.
    self.map_data = None
    self.map_coeffs = None
    map_coeff_data = None
    assert number_of_kicks > 0
    counter = 0
    for kick_size in [0,0.1,0.3,0.5,0.7]:
      b_sharp = None # 8 * math.pi**2 * kick_size**2 #XXX needs testing
      for trial in xrange(number_of_kicks):
        xray_structure = fmodel.xray_structure.deep_copy_scatterers()
        xray_structure.shake_sites_in_place(mean_distance = kick_size)
        max_kick = fmodel.xray_structure.max_distance(other = xray_structure)
        print kick_size, max_kick
        model_to_map_obj = model_to_map(
          fmodel_tmp                    = fmodel_tmp,
          xray_structure                = xray_structure,
          fmodel                        = fmodel,
          map_type                      = map_type,
          update_bulk_solvent_and_scale = update_bulk_solvent_and_scale,
          resolution_factor             = resolution_factor,
          other_fft_map                 = other_fft_map,
          symmetry_flags                = symmetry_flags,
          bss_params                    = bss_params,
          alpha_fom_source              = map_helper_obj)
        if(average_maps):
          self.fft_map = model_to_map_obj.fft_map()
          if(real_map):
            tmp_result = self.fft_map.real_map()
          elif(real_map_unpadded):
            tmp_result = self.fft_map.real_map_unpadded()
          if(self.map_data is None): self.map_data = tmp_result
          else: self.map_data += tmp_result
        else:
          map_coeff_ = model_to_map_obj.map_coefficients()
          map_coeff_data_ = map_coeff_.data()
          if(b_sharp is not None):
            map_coeff_data_ = map_coeff_data_ * flex.exp(b_sharp*ss)
          if(map_coeff_data is None):
            map_coeff_data = map_coeff_data_
          else:
            map_coeff_data = map_coeff_data + map_coeff_data_
        counter += 1
    if(self.map_data is not None):
      self.map_data = self.map_data/counter
      # produce sigma scaled map: copied from miller.py
      from cctbx import maptbx
      statistics = maptbx.statistics(self.map_data)
      self.average = statistics.mean()
      self.standard_deviation = statistics.sigma()
      self.map_data /= self.standard_deviation
    if(map_coeff_data is not None):
      self.map_coeffs = map_coeff_.customized_copy(
        data = map_coeff_data/counter)

class model_to_map(object):

  def __init__(self, fmodel_tmp,
                     xray_structure,
                     fmodel,
                     map_type,
                     update_bulk_solvent_and_scale,
                     resolution_factor,
                     other_fft_map,
                     symmetry_flags,
                     bss_params,
                     alpha_fom_source):
    self.map_type = map_type
    self.resolution_factor = resolution_factor
    self.other_fft_map = other_fft_map
    self.alpha_fom_source = alpha_fom_source
    self.symmetry_flags = symmetry_flags
    fmodel_tmp.update_xray_structure(xray_structure = xray_structure,
                                     update_f_calc  = True,
                                     update_f_mask  = True)
    if(update_bulk_solvent_and_scale):
      fmodel_tmp.update(
        b_cart            = fmodel.b_cart(),
        k_sol             = fmodel.k_sol(),
        b_sol             = fmodel.b_sol(),
        alpha_beta_params = fmodel.alpha_beta_params)
      bss_params.k_sol_b_sol_grid_search = False
      bss_params.number_of_macro_cycles = 1
      fmodel_tmp.update_solvent_and_scale(params = bss_params)
    self.electron_density_map_manager = fmodel_tmp.electron_density_map()

  def fft_map(self):
    result = self.electron_density_map_manager.\
      fft_map(map_type                   = self.map_type,
              other_fft_map              = self.other_fft_map,
              symmetry_flags             = self.symmetry_flags,
              resolution_factor          = self.resolution_factor,
              alpha_fom_source           = self.alpha_fom_source,
              force_anomalous_flag_false = True)
    result.apply_sigma_scaling()
    return result

  def map_coefficients(self):
    return self.electron_density_map_manager.map_coefficients(
      map_type = self.map_type, alpha_fom_source = self.alpha_fom_source)

class electron_density_map(object):

  def __init__(self, fmodel, fill_missing_f_obs = False):
    self.fmodel = fmodel.deep_copy()
    self.fill_missing_f_obs = fill_missing_f_obs
    self.anom_diff = None
    if(self.fmodel.f_obs.anomalous_flag()):
      self.anom_diff = self.fmodel.f_obs.anomalous_differences()
      f_model = self.fmodel.f_model().as_non_anomalous_array().\
        merge_equivalents().array()
      fmodel_match_anom_diff, anom_diff_common = \
        f_model.common_sets(other =  self.anom_diff)
      assert anom_diff_common.indices().size()==self.anom_diff.indices().size()
      self.anom_diff = self._phase_transfer(miller_array = anom_diff_common,
        phase_source = fmodel_match_anom_diff)
    if(self.fill_missing_f_obs): self.fmodel = self._fill_f_obs()
    self.map_helper_obj = map_helper(fmodel  = self.fmodel)
    #del self.fmodel # XXX

  def map_coefficients(self, map_type, alpha_fom_source = None):
    map_name_manager = mmtbx.map_names(map_name_string = map_type)
    if(map_name_manager.anomalous):
      if(self.anom_diff is not None):
        # Formula from page 141 in "The Bijvoet-Difference Fourier Synthesis",
        # Jeffrey Roach, METHODS IN ENZYMOLOGY, VOL. 374
        return miller.array(miller_set = self.anom_diff,
                            data       = self.anom_diff.data()/(2j))
      else: return None
    if(not map_name_manager.ml_map):
       return self._map_coeff(
         f_obs         = self.map_helper_obj.f_obs_scaled,
         f_model       = self.map_helper_obj.f_model_scaled,
         f_obs_scale   = map_name_manager.k,
         f_model_scale = map_name_manager.n)
    if(map_name_manager.ml_map):
      if(alpha_fom_source is not None):
        alpha = alpha_fom_source.alpha.data()
        fom = alpha_fom_source.alpha.data()
      else:
        alpha = self.map_helper_obj.alpha.data()
        fom = self.map_helper_obj.alpha.data()
      return self._map_coeff(
        f_obs         = self.map_helper_obj.f_obs_scaled,
        f_model       = self.map_helper_obj.f_model_scaled,
        f_obs_scale   = map_name_manager.k*fom,
        f_model_scale = map_name_manager.n*alpha)

  def _fill_f_obs(self):
    f_model = self.fmodel.f_model()
    n_refl_orig = f_model.data().size()
    complete_set = f_model.complete_set(d_min = f_model.d_min(), d_max=None)
    f_calc_atoms = complete_set.structure_factors_from_scatterers(
      xray_structure = self.fmodel.xray_structure).f_calc()
    f_calc_atoms_lone = f_calc_atoms.lone_set(other = f_model)
    n_refl_lone = f_calc_atoms_lone.data().size()
    f_mask_lone = masks.manager(
      miller_array = f_calc_atoms_lone,
      xray_structure = self.fmodel.xray_structure).f_mask()
    ss = 1./flex.pow2(f_mask_lone.d_spacings().data())/4.
    r_free_flags_lone = f_mask_lone.array(
      data = flex.bool(f_mask_lone.size(), False))
    f_model_core = ext.core(f_calc = f_calc_atoms_lone.data(),
      f_mask = f_mask_lone.data(),
      b_cart = self.fmodel.b_cart(),
      k_sol  = self.fmodel.k_sol(),
      b_sol  = self.fmodel.b_sol(),
      hkl    = f_calc_atoms_lone.indices(),
      uc     = f_mask_lone.unit_cell(),
      ss     = ss)
    f_obs_orig = self.fmodel.f_obs.deep_copy()
    r_free_flags_orig = self.fmodel.r_free_flags
    # compose new fileld fmodel
    f_model_lone = abs(miller.array(
      miller_set = f_mask_lone,
      data       = f_model_core.f_model * self.fmodel.scale_k1()))
    new_f_obs = self.fmodel.f_obs.concatenate(other = f_model_lone)
    new_r_free_flags = self.fmodel.r_free_flags.concatenate(
      other = r_free_flags_lone)
    self.fmodel = mmtbx.f_model.manager(
      xray_structure = self.fmodel.xray_structure,
      r_free_flags   = new_r_free_flags,
      target_name    = "ls_wunit_k1",
      f_obs          = new_f_obs)
    self.fmodel.update_solvent_and_scale()
    # replace 'F_obs' -> alpha * 'F_obs' for filled F_obs
    alpha, beta = maxlik.alpha_beta_est_manager(
      f_obs                    = self.fmodel.f_obs,
      f_calc                   = self.fmodel.f_model_scaled_with_k1(),
      free_reflections_per_bin = 100,
      flags                    = self.fmodel.r_free_flags.data(),
      interpolation            = True).alpha_beta()
    apply_alpha_sel = flex.bool(n_refl_orig, False).concatenate(
      flex.bool(n_refl_lone, True)) # assume order did not change
    assert apply_alpha_sel.size() == self.fmodel.f_obs.data().size()
    alpha = alpha.select(apply_alpha_sel)
    # compose new fileld fmodel
    f_model_lone = abs(miller.array(
      miller_set = f_mask_lone,
      data       = f_model_core.f_model * self.fmodel.scale_k1()*alpha.data()))
    new_f_obs = f_obs_orig.concatenate(other = f_model_lone)
    new_r_free_flags = r_free_flags_orig.concatenate(
      other = r_free_flags_lone)
    self.fmodel = mmtbx.f_model.manager(
      xray_structure = self.fmodel.xray_structure,
      r_free_flags   = new_r_free_flags,
      target_name    = "ls_wunit_k1",
      f_obs          = new_f_obs)
    self.fmodel.update_solvent_and_scale()
    #
    return self.fmodel

  def _phase_transfer(self, miller_array, phase_source):
    # XXX could be a method in miller.py under a better name in future
    tmp = miller.array(miller_set = miller_array,
      data = flex.double(miller_array.indices().size(), 1)
      ).phase_transfer(phase_source = phase_source)
    return miller.array(miller_set = miller_array,
      data = miller_array.data() * tmp.data())

  def _map_coeff(self, f_obs, f_model, f_obs_scale, f_model_scale):
    obs = miller.array(miller_set = f_model,
                       data       = f_obs.data()*f_obs_scale)
    obs_phi_calc = self._phase_transfer(miller_array = obs,
      phase_source = f_model)
    result = miller.array(
      miller_set = obs_phi_calc,
      data       = obs_phi_calc.data()-f_model.data()*f_model_scale)
    # I don't understand why, but this really improves the maps.
    # CNS does the same.
    centrics  = result.select_centric()
    acentrics = result.select_acentric()
    acentrics_data = acentrics.data() * 2.0
    centrics_data  = centrics.data()
    result = acentrics.customized_copy(
      indices = acentrics.indices().concatenate(centrics.indices()),
      data    = acentrics_data.concatenate(centrics_data))
    return result

  def fft_map(self,
              resolution_factor = 1/3.,
              symmetry_flags = None,
              map_coefficients = None,
              other_fft_map = None,
              map_type = None,
              alpha_fom_source = None,
              force_anomalous_flag_false = None):
    if(map_coefficients is None):
      map_coefficients = self.map_coefficients(
        map_type         = map_type,
        alpha_fom_source = alpha_fom_source)
      if(force_anomalous_flag_false):
        map_coefficients = map_coefficients.average_bijvoet_mates()
    if(force_anomalous_flag_false):
      map_coefficients = map_coefficients.average_bijvoet_mates()
    if(other_fft_map is None):
      return map_coefficients.fft_map(
        resolution_factor = resolution_factor,
        symmetry_flags    = symmetry_flags)
    else:
      return miller.fft_map(
        crystal_gridding     = other_fft_map,
        fourier_coefficients = map_coefficients)
