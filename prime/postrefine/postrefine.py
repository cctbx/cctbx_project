from __future__ import division
from __future__ import print_function
from cctbx.array_family import flex
from cctbx import miller, crystal
from six.moves import cPickle as pickle
from .mod_leastsqr import leastsqr_handler
from .mod_results import postref_results
from cctbx.crystal import symmetry
from scitbx.matrix import sqr
from cctbx import statistics
from .mod_partiality import partiality_handler
from .mod_lbfgs_partiality import lbfgs_partiality_handler
from .mod_mx import mx_handler
import math, os
from .mod_input import read_frame

class postref_handler(object):
  """
  handle post-refinement
  - read-in and store input in input_handler object
  - generate a mean-intensity-scaled mtz file as a reference set
  - perform post-refinement
  """
  def __init__(self):
    """
    Constructor
    """

  def organize_input(self, observations_pickle, iparams, avg_mode, pickle_filename=None):
    """Given the pickle file, extract and prepare observations object and
    the alpha angle (meridional to equatorial).
    """
    #get general parameters
    if iparams.isoform_name is not None:
      if "identified_isoform" not in observations_pickle:
        return None, "No identified isoform"
      if observations_pickle["identified_isoform"] != iparams.isoform_name:
        return None, "Identified isoform(%s) is not the requested isoform (%s)"%(observations_pickle["identified_isoform"], iparams.isoform_name)
    if iparams.flag_weak_anomalous:
      if avg_mode == 'final':
        target_anomalous_flag = iparams.target_anomalous_flag
      else:
        target_anomalous_flag = False
    else:
      target_anomalous_flag = iparams.target_anomalous_flag
    img_filename_only = ''
    if pickle_filename: img_filename_only = os.path.basename(pickle_filename)
    txt_exception = ' {0:40} ==> '.format(img_filename_only)
    #for dials integration pickles - also look for experimentxxx.json
    if "miller_index" in observations_pickle:
      from dxtbx.model.experiment_list import ExperimentListFactory
      exp_json_file = os.path.join(os.path.dirname(pickle_filename),img_filename_only.split('_')[0]+'_refined_experiments.json')
      if os.path.isfile(exp_json_file):
        experiments = ExperimentListFactory.from_json_file(exp_json_file)
        dials_crystal = experiments[0].crystal
        detector = experiments[0].detector
        beam = experiments[0].beam
        crystal_symmetry = crystal.symmetry(
            unit_cell=dials_crystal.get_unit_cell().parameters(),
            space_group_symbol=iparams.target_space_group)
        miller_set_all=miller.set(
                    crystal_symmetry=crystal_symmetry,
                    indices=observations_pickle['miller_index'],
                    anomalous_flag=target_anomalous_flag)
        observations = miller_set_all.array(
                  data=observations_pickle['intensity.sum.value'],
                  sigmas=flex.sqrt(observations_pickle['intensity.sum.variance'])).set_observation_type_xray_intensity()
        detector_distance_mm = detector[0].get_distance()
        alpha_angle_obs = flex.double([0]*len(observations.data()))
        wavelength = beam.get_wavelength()
        spot_pred_x_mm = observations_pickle['s1'] #a disguise of s1
        spot_pred_y_mm = flex.double([0]*len(observations.data()))
        #calculate the crystal orientation
        O = sqr(dials_crystal.get_unit_cell().orthogonalization_matrix()).transpose()
        R = sqr(dials_crystal.get_U()).transpose()
        from cctbx.crystal_orientation import crystal_orientation, basis_type
        crystal_init_orientation = crystal_orientation(O*R, basis_type.direct)
      else:
        txt_exception += exp_json_file+' not found'
        print(txt_exception)
        return None, txt_exception
    else:
      #for cctbx.xfel proceed as usual
      observations = observations_pickle["observations"][0]
      detector_distance_mm = observations_pickle['distance']
      mm_predictions = iparams.pixel_size_mm*(observations_pickle['mapped_predictions'][0])
      xbeam = observations_pickle["xbeam"]
      ybeam = observations_pickle["ybeam"]
      alpha_angle_obs = flex.double([math.atan(abs(pred[0]-xbeam)/abs(pred[1]-ybeam)) \
                                     for pred in mm_predictions])
      spot_pred_x_mm = flex.double([pred[0]-xbeam for pred in mm_predictions])
      spot_pred_y_mm = flex.double([pred[1]-ybeam for pred in mm_predictions])
      #Polarization correction
      wavelength = observations_pickle["wavelength"]
      crystal_init_orientation = observations_pickle["current_orientation"][0]
    #continue reading...
    if iparams.flag_LP_correction and "observations" in observations_pickle:
      fx = 1 - iparams.polarization_horizontal_fraction
      fy = 1 - fx
      if fx > 1.0 or fx < 0:
        print('Horizontal polarization fraction is not correct. The value must be >= 0 and <= 1')
        print('No polarization correction. Continue with post-refinement')
      else:
        phi_angle_obs = flex.double([math.atan2(pred[1]-ybeam, pred[0]-xbeam) \
                                         for pred in mm_predictions])
        bragg_angle_obs = observations.two_theta(wavelength).data()
        P = ((fx*((flex.sin(phi_angle_obs)**2)+((flex.cos(phi_angle_obs)**2)*flex.cos(bragg_angle_obs)**2)))+\
          (fy*((flex.cos(phi_angle_obs)**2)+((flex.sin(phi_angle_obs)**2)*flex.cos(bragg_angle_obs)**2))))
        I_prime = observations.data()/P
        sigI_prime =observations.sigmas()/P
        observations = observations.customized_copy(data=flex.double(I_prime),
                                                    sigmas=flex.double(sigI_prime))
    #set observations with target space group - !!! required for correct
    #merging due to map_to_asu command.
    if iparams.target_crystal_system is not None:
      target_crystal_system = iparams.target_crystal_system
    else:
      target_crystal_system = observations.crystal_symmetry().space_group().crystal_system()
    lph = lbfgs_partiality_handler()
    if iparams.flag_override_unit_cell:
      uc_constrained_inp = lph.prep_input(iparams.target_unit_cell.parameters(), target_crystal_system)
    else:
      uc_constrained_inp = lph.prep_input(observations.unit_cell().parameters(), target_crystal_system)
    uc_constrained = list(lph.prep_output(uc_constrained_inp, target_crystal_system))
    try:
      #apply constrain using the crystal system
      miller_set = symmetry(
          unit_cell=uc_constrained,
          space_group_symbol=iparams.target_space_group
        ).build_miller_set(
          anomalous_flag=target_anomalous_flag,
          d_min=iparams.merge.d_min)
      observations = observations.customized_copy(anomalous_flag=target_anomalous_flag,
                      crystal_symmetry=miller_set.crystal_symmetry())
    except Exception:
      a,b,c,alpha,beta,gamma = uc_constrained
      txt_exception += 'Mismatch spacegroup (%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f)'%(a,b,c,alpha,beta,gamma)
      print(txt_exception)
      return None, txt_exception
    #reset systematic absence
    sys_absent_negate_flags = flex.bool([sys_absent_flag[1]==False for sys_absent_flag in observations.sys_absent_flags()])
    observations = observations.select(sys_absent_negate_flags)
    alpha_angle_obs = alpha_angle_obs.select(sys_absent_negate_flags)
    spot_pred_x_mm = spot_pred_x_mm.select(sys_absent_negate_flags)
    spot_pred_y_mm = spot_pred_y_mm.select(sys_absent_negate_flags)

    #remove observations from rejection list
    if iparams.rejections:
      if pickle_filename in iparams.rejections:
        miller_indices_ori_rejected = iparams.rejections[pickle_filename]
        i_sel_flag = flex.bool([True]*len(observations.data()))
        cnrej = 0
        for miller_index_ori_rejected in miller_indices_ori_rejected:
          for i_index_ori, miller_index_ori in enumerate(observations.indices()):
            if miller_index_ori_rejected == miller_index_ori:
              i_sel_flag[i_index_ori] = False
              cnrej += 1
        observations = observations.customized_copy(indices=observations.indices().select(i_sel_flag),
            data=observations.data().select(i_sel_flag),
            sigmas=observations.sigmas().select(i_sel_flag))
        alpha_angle_obs = alpha_angle_obs.select(i_sel_flag)
        spot_pred_x_mm = spot_pred_x_mm.select(i_sel_flag)
        spot_pred_y_mm = spot_pred_y_mm.select(i_sel_flag)

    #filter resolution
    i_sel_res = observations.resolution_filter_selection(d_max=iparams.merge.d_max, d_min=iparams.merge.d_min)
    observations = observations.select(i_sel_res)
    alpha_angle_obs = alpha_angle_obs.select(i_sel_res)
    spot_pred_x_mm = spot_pred_x_mm.select(i_sel_res)
    spot_pred_y_mm = spot_pred_y_mm.select(i_sel_res)

    #Filter weak
    i_sel = (observations.data()/observations.sigmas()) > iparams.merge.sigma_min
    observations = observations.select(i_sel)
    alpha_angle_obs = alpha_angle_obs.select(i_sel)
    spot_pred_x_mm = spot_pred_x_mm.select(i_sel)
    spot_pred_y_mm = spot_pred_y_mm.select(i_sel)

    #filter icering (if on)
    if iparams.icering.flag_on:
      miller_indices = flex.miller_index()
      I_set = flex.double()
      sigI_set = flex.double()
      alpha_angle_obs_set = flex.double()
      spot_pred_x_mm_set = flex.double()
      spot_pred_y_mm_set = flex.double()
      for miller_index, d, I, sigI, alpha, spot_x, spot_y in zip(observations.indices(), observations.d_spacings().data(),
                                        observations.data(), observations.sigmas(), alpha_angle_obs,
                                        spot_pred_x_mm, spot_pred_y_mm):
        if d > iparams.icering.d_upper or d < iparams.icering.d_lower:
          miller_indices.append(miller_index)
          I_set.append(I)
          sigI_set.append(sigI)
          alpha_angle_obs_set.append(alpha)
          spot_pred_x_mm_set.append(spot_x)
          spot_pred_y_mm_set.append(spot_y)
      observations = observations.customized_copy(indices=miller_indices,
          data=I_set, sigmas=sigI_set)
      alpha_angle_obs = alpha_angle_obs_set[:]
      spot_pred_x_mm = spot_pred_x_mm_set[:]
      spot_pred_y_mm = spot_pred_y_mm_set[:]
    #replacing sigI (if set)
    if iparams.flag_replace_sigI:
      observations = observations.customized_copy(sigmas=flex.sqrt(observations.data()))
    inputs = observations, alpha_angle_obs, spot_pred_x_mm, spot_pred_y_mm, detector_distance_mm, wavelength, crystal_init_orientation
    return inputs, 'OK'

  def get_observations_non_polar(self, observations_original, pickle_filename, iparams):
    #return observations with correct polarity
    if iparams.indexing_ambiguity.index_basis_in is None:
      return observations_original.map_to_asu(), 'h,k,l'
    ind_pickle = iparams.indexing_ambiguity.index_basis_in
    if pickle_filename not in ind_pickle:
      return observations_original.map_to_asu(), 'Not Found'
    from cctbx import sgtbx
    cb_op = sgtbx.change_of_basis_op(ind_pickle[pickle_filename])
    observations_alt = observations_original.map_to_asu().change_basis(cb_op).map_to_asu()
    return observations_alt, ind_pickle[pickle_filename]

  def postrefine_by_frame(self, frame_no, pickle_filename, iparams, miller_array_ref, pres_in, avg_mode):
    #1. Prepare data
    observations_pickle = read_frame(pickle_filename)
    pickle_filepaths = pickle_filename.split('/')
    img_filename_only = pickle_filepaths[len(pickle_filepaths)-1]
    txt_exception = ' {0:40} ==> '.format(img_filename_only)
    if observations_pickle is None:
      txt_exception += 'empty or bad input file\n'
      return None, txt_exception
    inputs, txt_organize_input = self.organize_input(observations_pickle, iparams, avg_mode, pickle_filename=pickle_filename)
    if inputs is not None:
      observations_original, alpha_angle, spot_pred_x_mm, spot_pred_y_mm, detector_distance_mm, wavelength, crystal_init_orientation = inputs
    else:
      txt_exception += txt_organize_input + '\n'
      return None, txt_exception
    #2. Select data for post-refinement (only select indices that are common with the reference set
    observations_non_polar, index_basis_name = self.get_observations_non_polar(observations_original, pickle_filename, iparams)
    matches = miller.match_multi_indices(
                  miller_indices_unique=miller_array_ref.indices(),
                  miller_indices=observations_non_polar.indices())
    pair_0 = flex.size_t([pair[0] for pair in matches.pairs()])
    pair_1 = flex.size_t([pair[1] for pair in matches.pairs()])
    references_sel = miller_array_ref.select(pair_0)
    observations_original_sel = observations_original.select(pair_1)
    observations_non_polar_sel = observations_non_polar.select(pair_1)
    alpha_angle_set = alpha_angle.select(pair_1)
    spot_pred_x_mm_set = spot_pred_x_mm.select(pair_1)
    spot_pred_y_mm_set = spot_pred_y_mm.select(pair_1)
    #4. Do least-squares refinement
    lsqrh = leastsqr_handler()
    try:
      refined_params, stats, n_refl_postrefined = lsqrh.optimize(references_sel.data(),
                                                                   observations_original_sel, wavelength,
                                                                   crystal_init_orientation, alpha_angle_set,
                                                                   spot_pred_x_mm_set, spot_pred_y_mm_set,
                                                                   iparams,
                                                                   pres_in,
                                                                   observations_non_polar_sel,
                                                                   detector_distance_mm)
    except Exception:
      txt_exception += 'optimization failed.\n'
      return None, txt_exception
    #caculate partiality for output (with target_anomalous check)
    G_fin, B_fin, rotx_fin, roty_fin, ry_fin, rz_fin, r0_fin, re_fin, voigt_nu_fin, \
        a_fin, b_fin, c_fin, alpha_fin, beta_fin, gamma_fin = refined_params
    inputs, txt_organize_input = self.organize_input(observations_pickle, iparams, avg_mode, pickle_filename=pickle_filename)
    observations_original, alpha_angle, spot_pred_x_mm, spot_pred_y_mm, detector_distance_mm, wavelength, crystal_init_orientation = inputs
    observations_non_polar, index_basis_name = self.get_observations_non_polar(observations_original, pickle_filename, iparams)
    from cctbx.uctbx import unit_cell
    uc_fin = unit_cell((a_fin, b_fin, c_fin, alpha_fin, beta_fin, gamma_fin))
    if pres_in is not None:
      crystal_init_orientation = pres_in.crystal_orientation
    two_theta = observations_original.two_theta(wavelength=wavelength).data()
    ph = partiality_handler()
    partiality_fin, dummy, rs_fin, rh_fin = ph.calc_partiality_anisotropy_set(uc_fin, rotx_fin, roty_fin,
                                                           observations_original.indices(),
                                                           ry_fin, rz_fin, r0_fin, re_fin, voigt_nu_fin,
                                                           two_theta, alpha_angle, wavelength,
                                                           crystal_init_orientation,
                                                           spot_pred_x_mm, spot_pred_y_mm,
                                                           detector_distance_mm,
                                                           iparams.partiality_model,
                                                           iparams.flag_beam_divergence)
    #calculate the new crystal orientation
    O = sqr(uc_fin.orthogonalization_matrix()).transpose()
    R = sqr(crystal_init_orientation.crystal_rotation_matrix()).transpose()
    from cctbx.crystal_orientation import crystal_orientation, basis_type
    CO = crystal_orientation(O*R, basis_type.direct)
    crystal_fin_orientation = CO.rotate_thru((1,0,0), rotx_fin
                               ).rotate_thru((0,1,0), roty_fin)
    #remove reflections with partiality below threshold
    i_sel = partiality_fin > iparams.merge.partiality_min
    partiality_fin_sel = partiality_fin.select(i_sel)
    rs_fin_sel = rs_fin.select(i_sel)
    rh_fin_sel = rh_fin.select(i_sel)
    observations_non_polar_sel = observations_non_polar.customized_copy(\
        indices=observations_non_polar.indices().select(i_sel),
        data=observations_non_polar.data().select(i_sel),
        sigmas=observations_non_polar.sigmas().select(i_sel))
    observations_original_sel = observations_original.customized_copy(\
        indices=observations_original.indices().select(i_sel),
        data=observations_original.data().select(i_sel),
        sigmas=observations_original.sigmas().select(i_sel))
    pres = postref_results()
    pres.set_params(observations = observations_non_polar_sel,
            observations_original = observations_original_sel,
            refined_params=refined_params,
            stats=stats,
            partiality=partiality_fin_sel,
            rs_set=rs_fin_sel,
            rh_set=rh_fin_sel,
            frame_no=frame_no,
            pickle_filename=pickle_filename,
            wavelength=wavelength,
            crystal_orientation=crystal_fin_orientation,
            detector_distance_mm=detector_distance_mm)
    r_change = ((pres.R_final - pres.R_init)/pres.R_init)*100
    r_xy_change = ((pres.R_xy_final - pres.R_xy_init)/pres.R_xy_init)*100
    cc_change = ((pres.CC_final - pres.CC_init)/pres.CC_init)*100
    txt_postref= '{0:40} => RES:{1:5.2f} NREFL:{2:5d} R:{3:6.1f}% RXY:{4:5.1f}% CC:{5:5.1f}% G:{6:6.4f} B:{7:5.1f} CELL:{8:6.1f}{9:6.1f} {10:6.1f} {11:5.1f} {12:5.1f} {13:5.1f}'.format(img_filename_only+' ('+index_basis_name+')', observations_original_sel.d_min(), len(observations_original_sel.data()), r_change, r_xy_change, cc_change, pres.G, pres.B, a_fin, b_fin, c_fin, alpha_fin, beta_fin, gamma_fin)
    print(txt_postref)
    txt_postref += '\n'
    return pres, txt_postref

  def calc_mean_intensity(self, pickle_filename, iparams, avg_mode):
    observations_pickle = read_frame(pickle_filename)
    pickle_filepaths = pickle_filename.split('/')
    txt_exception = ' {0:40} ==> '.format(pickle_filepaths[len(pickle_filepaths)-1])
    if observations_pickle is None:
      txt_exception += 'empty or bad input file\n'
      return None, txt_exception
    inputs, txt_organize_input = self.organize_input(observations_pickle, iparams, avg_mode, pickle_filename=pickle_filename)
    if inputs is not None:
      observations_original, alpha_angle_obs, spot_pred_x_mm, spot_pred_y_mm, detector_distance_mm, wavelength, crystal_init_orientation = inputs
    else:
      txt_exception += txt_organize_input + '\n'
      return None, txt_exception
    #filter resolution
    observations_sel = observations_original.resolution_filter(d_min=iparams.scale.d_min, d_max=iparams.scale.d_max)
    #filer sigma
    i_sel = (observations_sel.data()/observations_sel.sigmas()) > iparams.scale.sigma_min
    if len(observations_sel.data().select(i_sel)) == 0:
      return None, txt_exception
    mean_I = flex.median(observations_sel.data().select(i_sel))
    return mean_I, txt_exception+'ok'

  def scale_frame_by_mean_I(self, frame_no, pickle_filename, iparams, mean_of_mean_I, avg_mode):
    observations_pickle = read_frame(pickle_filename)
    pickle_filepaths = pickle_filename.split('/')
    img_filename_only = pickle_filepaths[len(pickle_filepaths)-1]
    txt_exception = ' {0:40} ==> '.format(img_filename_only)
    if observations_pickle is None:
      txt_exception += 'empty or bad input file\n'
      return None, txt_exception
    inputs, txt_organize_input = self.organize_input(observations_pickle, iparams, avg_mode, pickle_filename=pickle_filename)
    if inputs is not None:
      observations_original, alpha_angle, spot_pred_x_mm, spot_pred_y_mm, detector_distance_mm, wavelength, crystal_init_orientation = inputs
    else:
      txt_exception += txt_organize_input + '\n'
      return None, txt_exception
    #select only reflections matched with scale input params.
    #filter by resolution
    i_sel_res = observations_original.resolution_filter_selection(d_min=iparams.scale.d_min,
                                                                  d_max=iparams.scale.d_max)
    observations_original_sel = observations_original.select(i_sel_res)
    alpha_angle_sel = alpha_angle.select(i_sel_res)
    spot_pred_x_mm_sel = spot_pred_x_mm.select(i_sel_res)
    spot_pred_y_mm_sel = spot_pred_y_mm.select(i_sel_res)
    #filter by sigma
    i_sel_sigmas = (observations_original_sel.data()/observations_original_sel.sigmas()) > iparams.scale.sigma_min
    observations_original_sel = observations_original_sel.select(i_sel_sigmas)
    alpha_angle_sel = alpha_angle_sel.select(i_sel_sigmas)
    spot_pred_x_mm_sel = spot_pred_x_mm_sel.select(i_sel_sigmas)
    spot_pred_y_mm_sel = spot_pred_y_mm_sel.select(i_sel_sigmas)
    observations_non_polar_sel, index_basis_name = self.get_observations_non_polar(observations_original_sel, pickle_filename, iparams)
    observations_non_polar, index_basis_name = self.get_observations_non_polar(observations_original, pickle_filename, iparams)
    uc_params = observations_original.unit_cell().parameters()
    ph = partiality_handler()
    r0 = ph.calc_spot_radius(sqr(crystal_init_orientation.reciprocal_matrix()),
                                 observations_original_sel.indices(), wavelength)
    #calculate first G
    (G, B) = (1,0)
    stats = (0,0,0,0,0,0,0,0,0,0)
    if mean_of_mean_I > 0:
      G = flex.median(observations_original_sel.data())/mean_of_mean_I
    if iparams.flag_apply_b_by_frame:
      try:
        mxh = mx_handler()
        asu_contents = mxh.get_asu_contents(iparams.n_residues)
        observations_as_f = observations_non_polar_sel.as_amplitude_array()
        binner_template_asu = observations_as_f.setup_binner(auto_binning=True)
        wp = statistics.wilson_plot(observations_as_f, asu_contents, e_statistics=True)
        G = wp.wilson_intensity_scale_factor * 1e2
        B = wp.wilson_b
      except Exception:
        txt_exception += 'warning B-factor calculation failed.\n'
        return None, txt_exception
    two_theta = observations_original.two_theta(wavelength=wavelength).data()
    sin_theta_over_lambda_sq = observations_original.two_theta(wavelength=wavelength).sin_theta_over_lambda_sq().data()
    ry, rz, re, voigt_nu, rotx, roty = (0, 0, iparams.gamma_e, iparams.voigt_nu, 0, 0)
    partiality_init, delta_xy_init, rs_init, rh_init = ph.calc_partiality_anisotropy_set(\
                                                          crystal_init_orientation.unit_cell(),
                                                          rotx, roty, observations_original.indices(),
                                                          ry, rz, r0, re, voigt_nu,
                                                          two_theta, alpha_angle, wavelength,
                                                          crystal_init_orientation, spot_pred_x_mm, spot_pred_y_mm,
                                                          detector_distance_mm, iparams.partiality_model,
                                                          iparams.flag_beam_divergence)
    if iparams.flag_plot_expert:
      n_bins = 20
      binner = observations_original.setup_binner(n_bins=n_bins)
      binner_indices = binner.bin_indices()
      avg_partiality_init = flex.double()
      avg_rs_init = flex.double()
      avg_rh_init = flex.double()
      one_dsqr_bin = flex.double()
      for i in range(1,n_bins+1):
        i_binner = (binner_indices == i)
        if len(observations_original.data().select(i_binner)) > 0:
          print(binner.bin_d_range(i)[1], flex.mean(partiality_init.select(i_binner)), flex.mean(rs_init.select(i_binner)), flex.mean(rh_init.select(i_binner)), len(partiality_init.select(i_binner)))
    #monte-carlo merge
    if iparams.flag_monte_carlo:
      G = 1
      B = 0
      partiality_init=flex.double([1]*len(partiality_init))
    #save results
    refined_params = flex.double([G,B,rotx,roty,ry,rz,r0,re,voigt_nu,uc_params[0],uc_params[1],uc_params[2],uc_params[3],uc_params[4],uc_params[5]])
    pres = postref_results()
    pres.set_params(observations = observations_non_polar,
            observations_original = observations_original,
            refined_params=refined_params,
            stats=stats,
            partiality=partiality_init,
            rs_set=rs_init,
            rh_set=rh_init,
            frame_no=frame_no,
            pickle_filename=pickle_filename,
            wavelength=wavelength,
            crystal_orientation=crystal_init_orientation,
            detector_distance_mm=detector_distance_mm)
    txt_scale_frame_by_mean_I = ' {0:40} ==> RES:{1:5.2f} NREFL:{2:5d} G:{3:6.4f} B:{4:6.1f} CELL:{5:6.2f} {6:6.2f} {7:6.2f} {8:6.2f} {9:6.2f} {10:6.2f}'.format(img_filename_only+' ('+index_basis_name+')', observations_original.d_min(), len(observations_original_sel.data()), G, B, uc_params[0],uc_params[1],uc_params[2],uc_params[3],uc_params[4],uc_params[5])
    print(txt_scale_frame_by_mean_I)
    txt_scale_frame_by_mean_I += '\n'
    return pres, txt_scale_frame_by_mean_I
