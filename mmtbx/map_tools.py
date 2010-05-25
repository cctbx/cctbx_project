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
                     kick_sizes         = [0,0.2,0.4,0.6,0.8],
                     number_of_kicks    = 30,
                     acentrics_scale    = 2.0,
                     centrics_pre_scale = 1.0):
    fmodel_tmp = fmodel.deep_copy()
    map_helper_obj = fmodel.map_calculation_helper()
    map_coeff_data = None
    self.map_coeffs = None
    counter = 0
    for kick_size in kick_sizes:
      for kick in xrange(number_of_kicks):
        counter += 1
        print kick_size, kick
        xray_structure = fmodel.xray_structure.deep_copy_scatterers()
        xray_structure.shake_sites_in_place(mean_distance = kick_size)
        fmodel_tmp.update_xray_structure(xray_structure = xray_structure,
                                         update_f_calc  = True,
                                         update_f_mask  = True,
                                         force_update_f_mask = True)
        self.map_coeffs = fmodel_tmp.electron_density_map().map_coefficients(
          map_type         = map_type,
          alpha_fom_source = map_helper_obj)
        if(map_coeff_data is None):
          map_coeff_data = self.map_coeffs.data()
        else:
          map_coeff_data = map_coeff_data + self.map_coeffs.data()
    if(self.map_coeffs is not None):
      self.map_coeffs = self.map_coeffs.customized_copy(
        data = map_coeff_data/counter)

def b_sharp_map(map_coefficients,
                resolution_factor = 1/4.,
                b_sharp_min  = 0,
                b_sharp_max  = 260,
                b_sharp_step = 25,
                b_sharp      = None):
  ss = 1./flex.pow2(map_coefficients.d_spacings().data()) / 4.
  fft_map_start = map_coefficients.fft_map(resolution_factor=resolution_factor)
  def sigma_scale(map_data):
    statistics = maptbx.statistics(map_data)
    average = statistics.mean()
    standard_deviation = statistics.sigma()
    map_data /= standard_deviation
    return map_data
  if(b_sharp is not None):
    scale = flex.exp(b_sharp*ss)
    return map_coefficients.customized_copy(data=map_coefficients.data()*scale)
  else:
    fft_map_start.apply_sigma_scaling()
    map_data_start = fft_map_start.real_map_unpadded()
    counter = 0
    mcd = None
    map_data = map_data_start.deep_copy()
    for b_sharp in range(b_sharp_min, b_sharp_max, b_sharp_step):
      counter += 1
      scale = flex.exp(b_sharp*ss)
      mc = map_coefficients.customized_copy(data=map_coefficients.data()*scale)
      fft_map = miller.fft_map(crystal_gridding = fft_map_start,
        fourier_coefficients = mc)
      fft_map.apply_sigma_scaling()
      map_data_ = fft_map.real_map_unpadded()
      map_data = maptbx.combine_and_maximize_maps(
        map_data_1 = map_data_,
        map_data_2 = map_data,
        n_real     = fft_map.n_real())
      map_data = sigma_scale(map_data)
    map_data = maptbx.combine_and_maximize_maps(
      map_data_1 = map_data,
      map_data_2 = map_data_start,
      n_real     = fft_map.n_real())
    result = map_coefficients.structure_factors_from_map(map = map_data,
      use_scale = True)
    return result


class electron_density_map(object):

  def __init__(self,
               fmodel,
               fill_missing_f_obs = False,
               filled_f_obs_file_name = None,
               fill_mode = None,
               reverse_scale = True):
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
      if 0: # XXX make it an option
        self.fmodel.export_filled_f_obs(file_name = filled_f_obs_file_name)
    self.map_helper_obj = self.fmodel.map_calculation_helper(
      reverse_scale = reverse_scale)
    #del self.fmodel # XXX

  def map_coefficients(self, map_type, alpha_fom_source = None,
                       acentrics_scale = 2.0, centrics_pre_scale = 1.0):
    map_name_manager = mmtbx.map_names(map_name_string = map_type)
    if(map_name_manager.anomalous):
      if(self.anom_diff is not None):
        # Formula from page 141 in "The Bijvoet-Difference Fourier Synthesis",
        # Jeffrey Roach, METHODS IN ENZYMOLOGY, VOL. 374
        return miller.array(miller_set = self.anom_diff,
                            data       = self.anom_diff.data()/(2j))
      else: return None
    #
    # R.Read, SIGMAA: 2mFo-DFc (acentrics) & mFo (centrics)
    #
    centric_flags = self.fmodel.f_obs.centric_flags().data()
    if(map_name_manager.k != 0):
      fo_scale = flex.double(self.fmodel.f_obs.size(), 1.0)
    else: fo_scale = flex.double(self.fmodel.f_obs.size(), 0.0)
    if(map_name_manager.n != 0):
      fc_scale = flex.double(self.fmodel.f_obs.size(), 1.0)
    else: fc_scale = flex.double(self.fmodel.f_obs.size(), 0.0)
    if(map_name_manager.k != abs(map_name_manager.n) and
       abs(map_name_manager.k*map_name_manager.n) > 1.e-6):
      fo_scale.set_selected(~centric_flags, map_name_manager.k)
      fo_scale.set_selected(centric_flags, max(map_name_manager.k-centrics_pre_scale,0.))
      fc_scale.set_selected(~centric_flags, map_name_manager.n)
      fc_scale.set_selected(centric_flags, max(map_name_manager.n-centrics_pre_scale,0.))
    elif(map_name_manager.k == abs(map_name_manager.n) and
       abs(map_name_manager.k*map_name_manager.n) > 1.e-6):
      fo_scale.set_selected(~centric_flags, fo_scale*map_name_manager.k*acentrics_scale)
      fo_scale.set_selected( centric_flags, fo_scale*map_name_manager.k)
      fc_scale.set_selected(~centric_flags, fc_scale*map_name_manager.n*acentrics_scale)
      fc_scale.set_selected( centric_flags, fc_scale*map_name_manager.n)
    else:
      fo_scale *= map_name_manager.k
      fc_scale *= map_name_manager.n
    if(not map_name_manager.ml_map):
       return self._map_coeff(
         f_obs         = self.map_helper_obj.f_obs_scaled,
         f_model       = self.map_helper_obj.f_model_scaled,
         f_obs_scale   = fo_scale,
         f_model_scale = fc_scale)
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
          f_obs_scale   = fo_scale*fom,
          f_model_scale = fc_scale*alpha)
      else:
        comb_p = self.fmodel.combine_phases()
        fo_all_scales = self.map_helper_obj.f_obs_scaled.data()*fo_scale*\
          comb_p.f_obs_phase_and_fom_source()
        fc_all_scales = self.map_helper_obj.f_model_scaled.data()*\
          fc_scale*alpha
        return miller.array(
          miller_set = self.fmodel.f_calc(),
          data       = fo_all_scales + fc_all_scales)

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
      data       = obs_phi_calc.data()+f_model.data()*f_model_scale)
    return result

  def fft_map(self,
              resolution_factor = 1/3.,
              symmetry_flags = None,
              map_coefficients = None,
              other_fft_map = None,
              map_type = None,
              alpha_fom_source = None,
              force_anomalous_flag_false = None,
              acentrics_scale = 2.0,
              centrics_pre_scale = 1.0,
              use_all_data = True):
    if(map_coefficients is None):
      map_coefficients = self.map_coefficients(
        map_type           = map_type,
        acentrics_scale    = acentrics_scale,
        centrics_pre_scale = centrics_pre_scale,
        alpha_fom_source   = alpha_fom_source)
      if(force_anomalous_flag_false):
        map_coefficients = map_coefficients.average_bijvoet_mates()
    if(force_anomalous_flag_false):
      map_coefficients = map_coefficients.average_bijvoet_mates()
    if(not use_all_data):
      map_coefficients = map_coefficients.select(self.fmodel.work)
    if(other_fft_map is None):
      return map_coefficients.fft_map(
        resolution_factor = resolution_factor,
        symmetry_flags    = symmetry_flags)
    else:
      return miller.fft_map(
        crystal_gridding     = other_fft_map,
        fourier_coefficients = map_coefficients)
