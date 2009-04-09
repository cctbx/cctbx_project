from cctbx import miller
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import xray
from cctbx import eltbx
from cctbx import maptbx
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
from libtbx.utils import user_plus_sys_time, show_total_time
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
    map_helper_obj = fmodel.map_calculation_helper()
    ss = fmodel.ss
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
        if 0:#(average_maps): # XXX disabled till further investigation (does not work)
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
    if(map_coeff_data is not None):
      self.map_coeffs = map_coeff_.customized_copy(
        data = map_coeff_data/counter)
    ###
    self.fft_map = self.map_coeffs.fft_map(resolution_factor = resolution_factor)
    if(real_map):
      self.map_data = self.fft_map.real_map()
    else:
      self.map_data = self.fft_map.real_map_unpadded()

    statistics = maptbx.statistics(self.map_data)
    self.average = statistics.mean()
    self.standard_deviation = statistics.sigma()
    self.map_data /= self.standard_deviation

    if 0: # XXX works but no much impact... Find out WHY?
      for i in xrange(500):
        print i
        sel = flex.random_bool(size = self.map_coeffs.data().size(), threshold = 0.5)
        map_coeffs = self.map_coeffs.select(sel)

        self.fft_map = miller.fft_map(crystal_gridding = self.fft_map,
          fourier_coefficients = map_coeffs)

        #
        if(real_map):
          map_data = self.fft_map.real_map()
        else:
          map_data = self.fft_map.real_map_unpadded()
        if(map_data is not None):
          # produce sigma scaled map: copied from miller.py
          statistics = maptbx.statistics(map_data)
          self.average = statistics.mean()
          self.standard_deviation = statistics.sigma()
          map_data /= self.standard_deviation
        self.map_data += map_data
      self.map_data = self.map_data / (i+1)

      statistics = maptbx.statistics(self.map_data)
      self.average = statistics.mean()
      self.standard_deviation = statistics.sigma()
      self.map_data /= self.standard_deviation
      #
      self.map_coeffs = self.map_coeffs.structure_factors_from_map(
        map = self.map_data,
        use_scale = True)
      ###

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

  def __init__(self, fmodel, fill_missing_f_obs = False, filled_f_obs_file_name = None, fill_mode = None):
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
    if(self.fill_missing_f_obs):
      self.fmodel = self.fmodel.fill_missing_f_obs(fill_mode = fill_mode)
      assert filled_f_obs_file_name is not None
      if 0: # XXX make it an option
        self.fmodel.export_filled_f_obs(file_name = filled_f_obs_file_name)
    self.map_helper_obj = self.fmodel.map_calculation_helper()
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
    centric_flags = self.fmodel.f_obs.centric_flags().data()
    cf_scale = flex.double(self.fmodel.f_obs.size(), 1.0)
    acf_scale = flex.double(self.fmodel.f_obs.size(), 1.0)
    fo_scale = flex.double(self.fmodel.f_obs.size(), 1.0)
    if(map_name_manager.k != map_name_manager.n and
       abs(map_name_manager.k*map_name_manager.n) > 1.e-6):
      cf_scale = (~centric_flags).as_double()
      fo_scale.set_selected(~centric_flags, map_name_manager.k)
      fo_scale.set_selected(centric_flags, max(map_name_manager.k-1.,0.))
    else:
      acf_scale.set_selected(~centric_flags, 2.0)
    if(not map_name_manager.ml_map):
       return self._map_coeff(
         f_obs         = self.map_helper_obj.f_obs_scaled,
         f_model       = self.map_helper_obj.f_model_scaled,
         f_obs_scale   = fo_scale*acf_scale,
         f_model_scale = map_name_manager.n*cf_scale*acf_scale)
    if(map_name_manager.ml_map):
      if(alpha_fom_source is not None):
        alpha = alpha_fom_source.alpha.data()
        fom = alpha_fom_source.alpha.data()
      else:
        alpha = self.map_helper_obj.alpha.data()
        fom = self.map_helper_obj.alpha.data()
      if(self.fmodel.abcd is None):
        return self._map_coeff(
          f_obs         = self.map_helper_obj.f_obs_scaled,
          f_model       = self.map_helper_obj.f_model_scaled,
          f_obs_scale   = fo_scale*fom*acf_scale,
          f_model_scale = map_name_manager.n*alpha*cf_scale*acf_scale)
      else:
        comb_p = self.fmodel.combine_phases()
        fo_all_scales = self.map_helper_obj.f_obs_scaled.data()*fo_scale*\
          acf_scale*comb_p.f_obs_phase_and_fom_source()
        fc_all_scales = self.map_helper_obj.f_model_scaled.data()*\
          map_name_manager.n*alpha*cf_scale*acf_scale
        return miller.array(
          miller_set = self.fmodel.f_calc(),
          data       = fo_all_scales - fc_all_scales)

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
