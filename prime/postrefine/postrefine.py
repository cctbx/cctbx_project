from __future__ import division
from cctbx.array_family import flex
from cctbx import miller
import numpy as np
import cPickle as pickle
from mod_util import intensities_scaler
from mod_leastsqr import leastsqr_handler
from mod_results import postref_results
from cctbx.crystal import symmetry
import math
from scitbx.matrix import sqr

class postref_handler(object):
  '''
  handle post-refinement
  - read-in and store input in input_handler object
  - generate a mean-intensity-scaled mtz file as a reference set
  - perform post-refinement
  '''
  def __init__(self):
    '''
    Constructor
    '''

  def organize_input(self, observations_pickle, iparams, pickle_filename=None):

    """Given the pickle file, extract and prepare observations object and
    the alpha angle (meridional to equatorial).
    """
    observations = observations_pickle["observations"][0]

    detector_distance_mm = observations_pickle['distance']
    mm_predictions = iparams.pixel_size_mm*(observations_pickle['mapped_predictions'][0])
    xbeam = observations_pickle["xbeam"]
    ybeam = observations_pickle["ybeam"]
    alpha_angle_obs = flex.double([math.atan(abs(pred[0]-xbeam)/abs(pred[1]-ybeam)) \
                                   for pred in mm_predictions])
    spot_pred_x_mm = flex.double([pred[0]-xbeam for pred in mm_predictions])
    spot_pred_y_mm = flex.double([pred[1]-ybeam for pred in mm_predictions])
    assert len(alpha_angle_obs)==len(observations.indices()), 'Size of alpha angles and observations are not equal %6.0f, %6.0f'%(len(alpha_angle_obs),len(observations.indices()))

    #Lorentz-polarization correction
    wavelength = observations_pickle["wavelength"]

    if iparams.flag_LP_correction:
      fx = iparams.polarization_horizontal_fraction
      fy = 1 - fx
      if fx > 1.0 or fx < 0:
        print 'Horizontal polarization fraction is not correct. The value must be >= 0 and <= 1'
        print 'No polarization correction. Continue with post-refinement'
      else:
        phi_angle_obs = flex.double([math.atan2(pred[1]-ybeam, pred[0]-xbeam) \
                                         for pred in mm_predictions])
        bragg_angle_obs = observations.two_theta(wavelength).data()
        P = ((fx*((np.sin(phi_angle_obs)**2)+((np.cos(phi_angle_obs)**2)*np.cos(bragg_angle_obs)**2)))+\
          (fy*((np.cos(phi_angle_obs)**2)+((np.sin(phi_angle_obs)**2)*np.cos(bragg_angle_obs)**2))))
        I_prime = observations.data()/P
        sigI_prime =observations.sigmas()/P
        observations_LP_corrected = observations.customized_copy(data=flex.double(I_prime),
                                                    sigmas=flex.double(sigI_prime))

    #set observations with target space group - !!! required for correct
    #merging due to map_to_asu command.
    miller_set = symmetry(
        unit_cell=observations.unit_cell().parameters(),
        space_group_symbol=iparams.target_space_group
      ).build_miller_set(
        anomalous_flag=iparams.target_anomalous_flag,
        d_min=iparams.merge.d_min)

    observations = observations.customized_copy(anomalous_flag=iparams.target_anomalous_flag,
                    crystal_symmetry=miller_set.crystal_symmetry())

    import os.path
    if os.path.isfile(iparams.run_no+'/rejections.txt'):
      txt_out = pickle_filename + ' \nN_before_rejection: ' + str(len(observations.data())) + '\n'
      #remove observations from rejection list
      file_reject = open(iparams.run_no+'/rejections.txt', 'r')
      data_reject=file_reject.read().split("\n")
      miller_indices_ori_rejected = flex.miller_index()
      for row_reject in data_reject:
        col_reject = row_reject.split()
        if len(col_reject) > 0:
          if col_reject[0].strip() == pickle_filename:
            miller_indices_ori_rejected.append((int(col_reject[1].strip()), int(col_reject[2].strip()), int(col_reject[3].strip())))

      if len(miller_indices_ori_rejected) > 0:
        i_sel_flag = flex.bool([True]*len(observations.data()))
        for miller_index_ori_rejected in miller_indices_ori_rejected:
          i_index_ori = 0
          for miller_index_ori in observations.indices():
            if miller_index_ori_rejected == miller_index_ori:
              i_sel_flag[i_index_ori] = False
              txt_out += ' -Discard:' + str(miller_index_ori[0]) + \
          ','+str(miller_index_ori[1])+','+str(miller_index_ori[2]) + '\n'
            i_index_ori += 1


        observations = observations.customized_copy(indices=observations.indices().select(i_sel_flag),
            data=observations.data().select(i_sel_flag),
            sigmas=observations.sigmas().select(i_sel_flag))
        alpha_angle_obs = alpha_angle_obs.select(i_sel_flag)
        spot_pred_x_mm = spot_pred_x_mm.select(i_sel_flag)
        spot_pred_y_mm = spot_pred_y_mm.select(i_sel_flag)
        txt_out += 'N_after_rejection: ' + str(len(observations.data())) + '\n'

      if iparams.flag_output_verbose:
        print txt_out
    #filter resolution
    i_sel_res = observations.resolution_filter_selection(d_max=iparams.merge.d_max, d_min=iparams.merge.d_min)
    observations = observations.customized_copy(indices=observations.indices().select(i_sel_res),
        data=observations.data().select(i_sel_res),
        sigmas=observations.sigmas().select(i_sel_res))
    alpha_angle_obs = alpha_angle_obs.select(i_sel_res)
    spot_pred_x_mm = spot_pred_x_mm.select(i_sel_res)
    spot_pred_y_mm = spot_pred_y_mm.select(i_sel_res)


    #Filter weak
    i_sel = (observations.data()/observations.sigmas()) > iparams.merge.sigma_min
    observations = observations.customized_copy(indices=observations.indices().select(i_sel),
        data=observations.data().select(i_sel),
        sigmas=observations.sigmas().select(i_sel)
        )
    alpha_angle_obs = alpha_angle_obs.select(i_sel)
    spot_pred_x_mm = spot_pred_x_mm.select(i_sel)
    spot_pred_y_mm = spot_pred_y_mm.select(i_sel)

    #calculate Wilson scale and B-factor only when data extend to high reslution
    wilson_b = 0

    observations_as_f = observations.as_amplitude_array()
    observations_as_f.setup_binner(auto_binning=True)
    from cctbx import statistics

    pdb_asu_contents = {'C':iparams.n_residues*5,'H':0,'N':iparams.n_residues*2,'O':int(round(iparams.n_residues*1.5)),'S':0}
    wilson_plot = statistics.wilson_plot(observations_as_f, pdb_asu_contents)
    wilson_k_b = (wilson_plot.wilson_intensity_scale_factor, wilson_plot.wilson_b)
    wilson_b = wilson_k_b[1]
    if iparams.flag_plot_expert:
      binner = observations.setup_binner(n_bins=iparams.n_bins)
      binner_indices = binner.bin_indices()
      avg_I_by_bin = flex.double()
      one_dsqr_by_bin = flex.double()
      for i in range(1,iparams.n_bins+1):
        i_binner = (binner_indices == i)
        if len(observations.data().select(i_binner)) > 0:
          avg_I_by_bin.append(np.mean(observations.data().select(i_binner)))
          one_dsqr_by_bin.append(1/binner.bin_d_range(i)[1]**2)

      import matplotlib.pyplot as plt
      x_axis = one_dsqr_by_bin
      plt.plot(x_axis, avg_I_by_bin, linestyle='-', linewidth=2.0, c='b', label='<I>')
      plt.title('Intensity plot by resolutions (B-factor=%6.2f)'%(wilson_k_b[1]))
      plt.xlabel('1/d^2')
      plt.ylabel('<I>')
      plt.show()

    if wilson_k_b[1] > iparams.wilson_b_max:
      print pickle_filename, ' - discarded (high B-factor: %6.2f)'%wilson_k_b[1]
      return None, 0, 0


    if iparams.flag_apply_b_by_frame:
      #apply Wilson B-factor to the intensity
      sin_theta_over_lambda_sq = observations.two_theta(wavelength=wavelength).sin_theta_over_lambda_sq().data()
      I_b = flex.exp(-2*wilson_b*sin_theta_over_lambda_sq) * observations.data()
      sigI_b = flex.exp(-2*wilson_b*sin_theta_over_lambda_sq) * observations.sigmas()
      observations = observations.customized_copy(data=I_b, sigmas=sigI_b)

    return observations, alpha_angle_obs, spot_pred_x_mm, spot_pred_y_mm, wilson_b, detector_distance_mm


  def determine_polar(self, observations_original, iparams, pickle_filename):

    """
    Determine polarity based on input data.
    The function still needs isomorphous reference so, if flag_polar is True,
    miller_array_iso must be supplied in input file.
    """
    if iparams.indexing_ambiguity.flag_on == False:
      return 'h,k,l', 0 , 0

    pickle_filename_arr = pickle_filename.split('/')
    if len(pickle_filename_arr) == 1:
      pickle_filename_only = pickle_filename_arr[0]
    else:
      pickle_filename_only = pickle_filename_arr[len(pickle_filename_arr)-1]

    #use basis in the given input file
    polar_hkl = 'h,k,l'
    basis_pickle = pickle.load(open(iparams.indexing_ambiguity.index_basis_in,"rb"))
    keys = basis_pickle.viewkeys()
    for key in keys:
      if key.find(pickle_filename_only) > 0:
        polar_hkl = basis_pickle[key]
        break

    return polar_hkl, 0, 0


  def get_observations_non_polar(self, observations_original, polar_hkl):
    #return observations with correct polarity
    observations_asu = observations_original.map_to_asu()
    assert len(observations_original.indices())==len(observations_asu.indices()), 'No. of original and asymmetric-unit indices are not equal %6.0f, %6.0f'%(len(observations_original.indices()), len(observations_asu.indices()))

    if polar_hkl == 'h,k,l':
      return observations_asu
    else:
      from cctbx import sgtbx
      cb_op = sgtbx.change_of_basis_op(polar_hkl)
      observations_rev = observations_asu.change_basis(cb_op).map_to_asu()
      assert len(observations_original.indices())==len(observations_rev.indices()), 'No. of original and inversed asymmetric-unit indices are not equal %6.0f, %6.0f'%(len(observations_original.indices()), len(observations_rev.indices()))
      return observations_rev

  def postrefine_by_frame(self, frame_no, pickle_filename, iparams, miller_array_ref, pres_in):

    #1. Prepare data
    observations_pickle = pickle.load(open(pickle_filename,"rb"))
    crystal_init_orientation = observations_pickle["current_orientation"][0]
    wavelength = observations_pickle["wavelength"]

    #grab img. name
    imgname = pickle_filename

    try:
      observations_original, alpha_angle, spot_pred_x_mm, spot_pred_y_mm,  wilson_b, detector_distance_mm = self.organize_input(observations_pickle, iparams, pickle_filename=pickle_filename)
    except Exception:
      observations_original = None

    if observations_original is None:
      print frame_no, '-fail obs is none'
      return None

    #2. Determine polarity - always do this even if flag_polar = False
    #the function will take care of it.
    polar_hkl, cc_iso_raw_asu, cc_iso_raw_rev = self.determine_polar(observations_original, iparams, pickle_filename)

    #3. Select data for post-refinement (only select indices that are common with the reference set
    observations_non_polar = self.get_observations_non_polar(observations_original, polar_hkl)
    matches = miller.match_multi_indices(
                  miller_indices_unique=miller_array_ref.indices(),
                  miller_indices=observations_non_polar.indices())

    I_ref_match = flex.double([miller_array_ref.data()[pair[0]] for pair in matches.pairs()])
    miller_indices_ref_match = flex.miller_index((miller_array_ref.indices()[pair[0]] for pair in matches.pairs()))
    I_obs_match = flex.double([observations_non_polar.data()[pair[1]] for pair in matches.pairs()])
    sigI_obs_match = flex.double([observations_non_polar.sigmas()[pair[1]] for pair in matches.pairs()])
    miller_indices_original_obs_match = flex.miller_index((observations_original.indices()[pair[1]] \
                                                           for pair in matches.pairs()))
    miller_indices_non_polar_obs_match = flex.miller_index((observations_non_polar.indices()[pair[1]] \
                                                           for pair in matches.pairs()))
    alpha_angle_set = flex.double([alpha_angle[pair[1]] for pair in matches.pairs()])
    spot_pred_x_mm_set = flex.double([spot_pred_x_mm[pair[1]] for pair in matches.pairs()])
    spot_pred_y_mm_set = flex.double([spot_pred_y_mm[pair[1]] for pair in matches.pairs()])
    references_sel = miller_array_ref.customized_copy(data=I_ref_match, indices=miller_indices_ref_match)
    observations_original_sel = observations_original.customized_copy(data=I_obs_match,
                                                                      sigmas=sigI_obs_match,
                                                                      indices=miller_indices_original_obs_match)

    observations_non_polar_sel = observations_non_polar.customized_copy(data=I_obs_match,
                                                                       sigmas=sigI_obs_match,
                                                                       indices=miller_indices_non_polar_obs_match)

    #4. Do least-squares refinement
    lsqrh = leastsqr_handler()
    try:
      refined_params, stats, n_refl_postrefined, spot_radius = lsqrh.optimize(I_ref_match,
                                                                 observations_original_sel, wavelength,
                                                                 crystal_init_orientation, alpha_angle_set,
                                                                 spot_pred_x_mm_set, spot_pred_y_mm_set,
                                                                 iparams,
                                                                 pres_in,
                                                                 observations_non_polar_sel,
                                                                 detector_distance_mm)
    except Exception:
      return None

    #caculate partiality for output (with target_anomalous check)
    G_fin, B_fin, rotx_fin, roty_fin, ry_fin, rz_fin, re_fin, \
        a_fin, b_fin, c_fin, alpha_fin, beta_fin, gamma_fin = refined_params
    observations_original, alpha_angle, spot_pred_x_mm, spot_pred_y_mm,  wilson_b, detector_distance_mm = \
        self.organize_input(observations_pickle, iparams, pickle_filename=pickle_filename)
    observations_non_polar = self.get_observations_non_polar(observations_original, polar_hkl)

    from cctbx.uctbx import unit_cell
    uc_fin = unit_cell((a_fin, b_fin, c_fin, alpha_fin, beta_fin, gamma_fin))
    if pres_in is not None:
      crystal_init_orientation = pres_in.crystal_orientation

    two_theta = observations_original.two_theta(wavelength=wavelength).data()
    from mod_leastsqr import calc_partiality_anisotropy_set
    partiality_fin, dummy, rs_fin, rh_fin = calc_partiality_anisotropy_set(uc_fin, rotx_fin, roty_fin,
                                                           observations_original.indices(),
                                                           ry_fin, rz_fin, spot_radius, re_fin,
                                                           two_theta, alpha_angle, wavelength, crystal_init_orientation,
                                                           spot_pred_x_mm, spot_pred_y_mm,
                                                           detector_distance_mm,
                                                           iparams.partiality_model)

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
            spot_radius=spot_radius)

    _tmp_pickle_filename = pickle_filename.split('/')
    print '%6.0f %5.2f %8.0f %8.0f %8.2f %8.2f %8.2f %8.2f %9.2f %7.2f %10.2f %10.2f   '%( \
      pres.frame_no, observations_non_polar.d_min(), len(observations_non_polar.indices()), \
      n_refl_postrefined, pres.R_init, pres.R_final, pres.R_xy_init, pres.R_xy_final, \
      pres.CC_init*100, pres.CC_final*100,pres.CC_iso_init*100, pres.CC_iso_final*100), \
      _tmp_pickle_filename[len(_tmp_pickle_filename)-1]

    return pres

  def calc_mean_intensity(self, pickle_filename, iparams):
    observations_pickle = pickle.load(open(pickle_filename,"rb"))
    wavelength = observations_pickle["wavelength"]

    try:
      observations_original, alpha_angle_obs, spot_pred_x_mm, spot_pred_y_mm, wilson_b, detector_distance_mm = self.organize_input(observations_pickle, iparams, pickle_filename=pickle_filename)
    except Exception:
      observations_original = None

    if observations_original is None:
      return None

    #filter resolution
    observations_sel = observations_original.resolution_filter(d_min=iparams.scale.d_min, d_max=iparams.scale.d_max)
    #filer sigma
    i_sel = (observations_sel.data()/observations_sel.sigmas()) > iparams.scale.sigma_min

    if len(observations_sel.data().select(i_sel)) == 0:
      return None
    mean_I = np.median(observations_sel.data().select(i_sel))

    return mean_I


  def scale_frame_by_mean_I(self, frame_no, pickle_filename, iparams, mean_of_mean_I):

    observations_pickle = pickle.load(open(pickle_filename,"rb"))

    try:
      observations_original, alpha_angle, spot_pred_x_mm, spot_pred_y_mm, wilson_b, detector_distance_mm = self.organize_input(observations_pickle, iparams, pickle_filename=pickle_filename)
    except Exception:
      observations_original = None

    if observations_original is None:
      return None
    wavelength = observations_pickle["wavelength"]
    crystal_init_orientation = observations_pickle["current_orientation"][0]
    crystal_pointgroup = observations_pickle["pointgroup"]

    if iparams.target_pointgroup is not None and crystal_pointgroup != iparams.target_pointgroup:
      print 'frame %6.0f'%frame_no, ' - wrong pointgroup', crystal_pointgroup
      return None

    #select only reflections matched with scale input params.
    #filter by resolution
    i_sel_res = observations_original.resolution_filter_selection(d_min=iparams.scale.d_min,
                                                                  d_max=iparams.scale.d_max)
    observations_original_sel = observations_original.customized_copy( \
      indices=observations_original.indices().select(i_sel_res),
      data=observations_original.data().select(i_sel_res),
      sigmas=observations_original.sigmas().select(i_sel_res))
    alpha_angle_sel = alpha_angle.select(i_sel_res)
    spot_pred_x_mm_sel = spot_pred_x_mm.select(i_sel_res)
    spot_pred_y_mm_sel = spot_pred_y_mm.select(i_sel_res)

    #filter by sigma
    i_sel_sigmas = (observations_original_sel.data()/observations_original_sel.sigmas()) > iparams.scale.sigma_min
    observations_original_sel = observations_original_sel.customized_copy(\
      indices=observations_original_sel.indices().select(i_sel_sigmas),
      data=observations_original_sel.data().select(i_sel_sigmas),
      sigmas=observations_original_sel.sigmas().select(i_sel_sigmas))
    alpha_angle_sel = alpha_angle_sel.select(i_sel_sigmas)
    spot_pred_x_mm_sel = spot_pred_x_mm_sel.select(i_sel_sigmas)
    spot_pred_y_mm_sel = spot_pred_y_mm_sel.select(i_sel_sigmas)

    polar_hkl, cc_iso_raw_asu, cc_iso_raw_rev = self.determine_polar(observations_original, iparams, pickle_filename)
    observations_non_polar_sel = self.get_observations_non_polar(observations_original_sel,
                                                                 polar_hkl)
    observations_non_polar = self.get_observations_non_polar(observations_original, polar_hkl)
    uc_params = observations_original.unit_cell().parameters()
    from mod_leastsqr import calc_spot_radius
    spot_radius = calc_spot_radius(sqr(crystal_init_orientation.reciprocal_matrix()),
                                          observations_original_sel.indices(), wavelength)

    G = mean_of_mean_I/np.median(observations_original_sel.data())
    B = 0
    stats = (0,0,0,0,0,0,0,0,0,0)


    from mod_leastsqr import calc_partiality_anisotropy_set
    two_theta = observations_original.two_theta(wavelength=wavelength).data()
    sin_theta_over_lambda_sq = observations_original.two_theta(wavelength=wavelength).sin_theta_over_lambda_sq().data()
    ry = 0
    rz = 0
    re = iparams.gamma_e
    rotx = 0
    roty = 0
    partiality_init, delta_xy_init, rs_init, rh_init = calc_partiality_anisotropy_set(crystal_init_orientation.unit_cell(), rotx, roty, observations_original.indices(), ry, rz, spot_radius, re, two_theta, alpha_angle, wavelength, crystal_init_orientation, spot_pred_x_mm, spot_pred_y_mm, detector_distance_mm,iparams.partiality_model)

    refined_params = np.array([G,B,rotx,roty,ry,rz,re,uc_params[0],uc_params[1],uc_params[2],uc_params[3],uc_params[4],uc_params[5]])

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
            spot_radius=spot_radius)

    _tmp_pickle_filename = pickle_filename.split('/')
    print '%6.0f %6.2f %6.0f %6.0f %14.2f %14.2f %14.2f %6.2f %6.2f'%(frame_no,
      observations_original.d_min(), len(observations_original.data()), len(observations_original_sel.data()),
      np.sum(observations_original_sel.data()), np.mean(observations_original_sel.data()),
      np.median(observations_original_sel.data()), G, B), _tmp_pickle_filename[len(_tmp_pickle_filename)-1]

    return pres


def prepare_output(results,
        iparams,
        avg_mode):

  #results is a list of postref_results objects
  #lenght of this list equals to number of input frames

  inten_scaler = intensities_scaler()
  prep_output = inten_scaler.prepare_output(results, iparams, avg_mode)
  return prep_output

def calc_avg_I(group_no, miller_index, miller_indices_ori, I, sigI, G, B, p_set, rs_set, wavelength_set,
               sin_theta_over_lambda_sq, SE, avg_mode, iparams, pickle_filename):

  #results is a list of postref_results objects
  #lenght of this list equals to number of input frames
  inten_scaler = intensities_scaler()
  avg_I_result = inten_scaler.calc_avg_I(group_no, miller_index, miller_indices_ori, I, sigI, G, B,
                     p_set, rs_set, wavelength_set, sin_theta_over_lambda_sq, SE, avg_mode,
                     iparams, pickle_filename)
  return avg_I_result

def write_output(miller_indices_merge, I_merge, sigI_merge, stat_all,
                   I_even, I_odd, iparams, uc_mean, wavelength_mean,
                   output_mtz_file_prefix, avg_mode):

  #results is a list of postref_results objects
  #lenght of this list equals to number of input frames
  inten_scaler = intensities_scaler()
  miller_array_merge, txt_out, csv_out = inten_scaler.write_output(miller_indices_merge, \
                                            I_merge, sigI_merge, stat_all, \
                                            I_even, I_odd, iparams, uc_mean, \
                                            wavelength_mean, output_mtz_file_prefix, avg_mode)
  return miller_array_merge, txt_out, csv_out

def read_input(args):
  from mod_input import process_input
  iparams, txt_out_input = process_input(args)
  return iparams, txt_out_input
