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

ext = boost.python.import_ext("mmtbx_f_model_ext")


class kick_map(object):

  def __init__(self, fmodel,
                     map_type,
                     kick_size,
                     number_of_kicks,
                     update_bulk_solvent_and_scale,
                     resolution_factor,
                     symmetry_flags,
                     other_fft_map = None,
                     real_map_unpadded = True,
                     real_map = False):
    assert [real_map_unpadded, real_map].count(True) == 1
    self.map_data = None
    assert number_of_kicks > 0
    counter = 0
    for kick_size in [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7]:
      b_sharp = 8 * math.pi**2 * kick_size**2
      for trial in xrange(number_of_kicks):
        xray_structure = fmodel.xray_structure.deep_copy_scatterers()
        xray_structure.shake_sites_in_place(mean_distance = kick_size)
        self.fft_map = model_to_map(
          xray_structure                = xray_structure,
          fmodel                        = fmodel,
          map_type                      = map_type,
          update_bulk_solvent_and_scale = update_bulk_solvent_and_scale,
          resolution_factor             = resolution_factor,
          other_fft_map                 = other_fft_map,
          symmetry_flags                = symmetry_flags,
          b_sharp                       = b_sharp).fft_map
        if(real_map):
          tmp_result = self.fft_map.real_map()
        elif(real_map_unpadded):
          tmp_result = self.fft_map.real_map_unpadded()
        if(self.map_data is None): self.map_data = tmp_result
        else: self.map_data += tmp_result
        counter += 1
    self.map_data = self.map_data/counter
    # produce sigma scaled map: copied from miller.py
    from cctbx import maptbx
    statistics = maptbx.statistics(self.map_data)
    self.map_data /= statistics.sigma()


class model_to_map(object):

  def __init__(self, xray_structure,
                     fmodel,
                     map_type,
                     update_bulk_solvent_and_scale,
                     resolution_factor,
                     other_fft_map,
                     symmetry_flags,
                     b_sharp):
    fmodel_result = mmtbx.f_model.manager(
      xray_structure = xray_structure,
      r_free_flags   = fmodel.r_free_flags,
      target_name    = fmodel.target_name,
      f_obs          = fmodel.f_obs)
    if(update_bulk_solvent_and_scale):
      fmodel_result.update_solvent_and_scale()
    else:
      fmodel_result.update(
        f_mask            = fmodel.f_mask(),
        b_cart            = fmodel.b_cart(),
        k_sol             = fmodel.k_sol(),
        b_sol             = fmodel.b_sol(),
        alpha_beta_params = fmodel.alpha_beta_params)
    self.fft_map = fmodel_result.electron_density_map(
      map_type          = map_type,
      resolution_factor = resolution_factor,
      other_fft_map     = other_fft_map,
      symmetry_flags    = symmetry_flags,
      b_sharp           = b_sharp)
    self.fft_map.apply_sigma_scaling()

###############################################################################

class electron_density_map(object):

  def __init__(self, fmodel, b_sharp = None, kick_map = False,
                     fill_missing_f_obs = False):
    self.fmodel = fmodel.deep_copy()
    self.b_sharp = b_sharp
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
    if(self.fill_missing_f_obs): self._fill_f_obs()
    ss = 1./flex.pow2(self.fmodel.f_obs.d_spacings().data())/4.
    fb_cart  = self.fmodel.fb_cart()
    scale_k1 = self.fmodel.scale_k1()
    f_obs_scale   = 1.0 / (fb_cart * scale_k1)
    f_model_scale = 1.0 / fb_cart
    self.f_obs_data_scaled = self.fmodel.f_obs.data() * f_obs_scale
    self.f_model_data_scaled = self.fmodel.f_model().data() * f_model_scale
    if(b_sharp is not None): # XXX determine automatically as suggested by Axel.
      self.f_obs_data_scaled *= flex.exp(b_sharp*ss)
      self.f_model_data_scaled *= flex.exp(b_sharp*ss)
    self.f_obs_scaled = self.fmodel.f_obs.array(data = self.f_obs_data_scaled)
    self.f_model_scaled = self.fmodel.f_obs.array(data = self.f_model_data_scaled)
    self.r_free_flags = self.fmodel.r_free_flags
    del self.fmodel

  def map_coefficients(self, map_type):
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
         f_obs         = self.f_obs_scaled,
         f_model       = self.f_model_scaled,
         f_obs_scale   = map_name_manager.k,
         f_model_scale = map_name_manager.n)
    if(map_name_manager.ml_map):
      alpha, beta = maxlik.alpha_beta_est_manager(
        f_obs                    = self.f_obs_scaled,
        f_calc                   = self.f_model_scaled,
        free_reflections_per_bin = 100,
        flags                    = self.r_free_flags.data(),
        interpolation            = True).alpha_beta()
      fom = max_lik.fom_and_phase_error(
        f_obs          = self.f_obs_data_scaled,
        f_model        = flex.abs(self.f_model_data_scaled),
        alpha          = alpha.data(),
        beta           = beta.data(),
        space_group    = self.r_free_flags.space_group(),
        miller_indices = self.r_free_flags.indices()).fom()
      return self._map_coeff(
        f_obs         = self.f_obs_scaled,
        f_model       = self.f_model_scaled,
        f_obs_scale   = map_name_manager.k*fom,
        f_model_scale = map_name_manager.n*alpha.data())

  def _fill_f_obs(self):
    f_model = self.fmodel.f_model()
    complete_set = f_model.complete_set(d_min = f_model.d_min(), d_max=None)
    f_calc_atoms = complete_set.structure_factors_from_scatterers(
      xray_structure = self.fmodel.xray_structure).f_calc()
    f_calc_atoms_lone = f_calc_atoms.lone_set(other = f_model)
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
    f_model_lone = abs(miller.array(
      miller_set = f_mask_lone,
      data       = f_model_core.f_model * self.fmodel.scale_k1()))
    new_f_obs = self.fmodel.f_obs.concatenate(other = f_model_lone)
    new_r_free_flags = self.fmodel.r_free_flags.concatenate(
      other = r_free_flags_lone)
    self.fmodel = mmtbx.f_model.manager(
      xray_structure = self.fmodel.xray_structure,
      r_free_flags   = new_r_free_flags,
      target_name    = "ml",
      f_obs          = new_f_obs)
    self.fmodel.update_solvent_and_scale()

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
              map_type = None):
    if(map_coefficients is None):
      map_coefficients = self.map_coefficients(map_type = map_type)
    if(other_fft_map is None):
      return map_coefficients.fft_map(
        resolution_factor = resolution_factor,
        symmetry_flags    = symmetry_flags)
    else:
      return miller.fft_map(
        crystal_gridding     = other_fft_map,
        fourier_coefficients = map_coefficients)
