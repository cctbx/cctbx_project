from cctbx.array_family import flex
from cctbx import miller
from cctbx import maptbx
import boost.python
import mmtbx

ext = boost.python.import_ext("mmtbx_f_model_ext")

class kick_map(object):

  def __init__(self, fmodel,
                     map_type,
                     kick_sizes         = [0.5],
                     number_of_kicks    = 50,
                     acentrics_scale    = 2.0,
                     centrics_pre_scale = 1.0,
                     isotropize         = False,
                     sharp              = False,
                     exclude_free_r_reflections = False):
    self.map_coeffs = None
    fmodel_tmp = fmodel.deep_copy()
    isotropize_helper = fmodel.isotropize_helper()
    map_helper_obj = fmodel.map_calculation_helper()
    map_coeff_data = None
    counter = 0
    b_sh_cntr = flex.double()
    ss = 1./flex.pow2(
      fmodel_tmp.f_obs().average_bijvoet_mates().d_spacings().data()) / 4.
    sites_frac = fmodel.xray_structure.sites_frac()
    for kick_size in kick_sizes:
      b_ext = None
      for kick in xrange(number_of_kicks):
        print "kick", kick
        counter += 1
        xray_structure = fmodel.xray_structure.deep_copy_scatterers()
        xray_structure.shake_sites_in_place(mean_distance = kick_size)
        fmodel_tmp.update_xray_structure(xray_structure = xray_structure,
                                         update_f_calc  = True,
                                         update_f_mask  = False)
        # XXX need to call average_bijvoet_mates() regardless of whether the
        # data are anomalous - otherwise the assertion below will fail with
        # some datasets.  there may be a quicker way to do this, but in the
        # context of kicked map computation it does not appear to matter.
        self.map_coeffs = fmodel_tmp.electron_density_map().map_coefficients(
            map_type = map_type,
            external_alpha_fom_source = map_helper_obj).average_bijvoet_mates()
        assert isotropize_helper.iso_scale.indices().all_eq(
          self.map_coeffs.indices())
        #
        if(isotropize):
          self.map_coeffs = self.map_coeffs.array(
            data = self.map_coeffs.data()*isotropize_helper.iso_scale.data())
        ##############################################
        if(0): #XXX experimental
          fft_map = self.map_coeffs.fft_map(resolution_factor=1./3)
          fft_map.apply_sigma_scaling()
          map_data = fft_map.real_map_unpadded()
          map_data = maptbx.denmod_simple(
            map_data = map_data,
            n_real   = fft_map.n_real())
          self.map_coeffs = self.map_coeffs.structure_factors_from_map(map=map_data,
            use_scale = True, anomalous_flag = False, use_sg = True)
        ##############################################
        if(sharp):
          self.map_coeffs, b_ext = mmtbx.maps.sharp_map(
            sites_frac = sites_frac,
            map_coeffs = self.map_coeffs,
            ss         = ss,
            b_sharp    = b_ext)
          b_sh_cntr.append(b_ext)
          if(b_sh_cntr.size()<2): b_ext = None
          else: b_ext = flex.mean(b_sh_cntr)
        if(map_coeff_data is None):
          map_coeff_data = self.map_coeffs.data()
        else:
          map_coeff_data = map_coeff_data + self.map_coeffs.data()
    if(sharp):
      self.map_coeffs, b_ext = mmtbx.maps.sharp_map(
        sites_frac = sites_frac,
        map_coeffs = self.map_coeffs,
        ss         = ss)
    if(self.map_coeffs is not None):
      self.map_coeffs = miller.set(
        crystal_symmetry = self.map_coeffs.crystal_symmetry(),
        indices          = self.map_coeffs.indices(),
        anomalous_flag   = self.map_coeffs.anomalous_flag()).array(
          data = map_coeff_data/counter)
      if (exclude_free_r_reflections) :
        self.map_coeffs = self.map_coeffs.select(fmodel.arrays.work_sel)

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
      self.anom_diff = self._phase_transfer(miller_array = anom_diff_common,
        phase_source = fmodel_match_anom_diff)

  def map_coefficients(self,
                       map_type,
                       acentrics_scale = 2.0,
                       centrics_pre_scale = 1.0,
                       external_alpha_fom_source = None):
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
    centric_flags = self.fmodel.f_obs().centric_flags().data()
    if(map_name_manager.k != 0):
      fo_scale = flex.double(self.fmodel.f_obs().size(), 1.0)
    else: fo_scale = flex.double(self.fmodel.f_obs().size(), 0.0)
    if(map_name_manager.n != 0):
      fc_scale = flex.double(self.fmodel.f_obs().size(), 1.0)
    else: fc_scale = flex.double(self.fmodel.f_obs().size(), 0.0)
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
       self.mch = self.fmodel.map_calculation_helper()
       return self._map_coeff(
         f_obs         = self.fmodel.f_obs(),
         f_model       = self.fmodel.f_model_scaled_with_k1(),
         f_obs_scale   = fo_scale,
         f_model_scale = fc_scale)
    if(map_name_manager.ml_map):
      if(self.mch is None):
        self.mch = self.fmodel.map_calculation_helper()
      alpha = self.mch.alpha
      fom = self.mch.fom
      if(external_alpha_fom_source is not None):
        alpha = external_alpha_fom_source.alpha
        fom = external_alpha_fom_source.fom
      if(self.fmodel.hl_coeffs() is None):
        return self._map_coeff(
          f_obs         = self.mch.f_obs,
          f_model       = self.mch.f_model,
          f_obs_scale   = fo_scale*fom,
          f_model_scale = fc_scale*alpha.data())
      else:
        comb_p = self.fmodel.combine_phases(map_calculation_helper=self.mch)
        fo_all_scales = self.mch.f_obs.data()*fo_scale*\
          comb_p.f_obs_phase_and_fom_source()
        fc_all_scales = self.mch.f_model.data()*\
          fc_scale*alpha.data()
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
