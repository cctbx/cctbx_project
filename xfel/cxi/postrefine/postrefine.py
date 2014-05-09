"""
Read pickle files output from integration process to determine the polarity and
refine the rotation matrix
Usage:
cxi.postrefine input=input.inp
"""
from __future__ import division
from cctbx.array_family import flex
from cctbx import miller
from scitbx.matrix import col
import numpy as np
import cPickle as pickle
from mod_util import intensities_scaler, file_handler, input_handler
from mod_leastsqr import leastsqr_handler
from mod_results import postref_results
from cctbx.crystal import symmetry
import math

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

  def read_input_parameters(self, file_name_input):
    iph = input_handler()
    iph.read_input(file_name_input)
    return iph

  def organize_input(self, observations_pickle, iph):

    """Given the pickle file, extract and prepare observations object and
    the alpha angle (meridional to equatorial).
    """
    observations = observations_pickle["observations"][0]
    mm_predictions = iph.pixel_size_mm*(observations_pickle['mapped_predictions'][0])
    xbeam = observations_pickle["xbeam"]
    ybeam = observations_pickle["ybeam"]
    alpha_angle_obs = flex.double([math.atan(abs(pred[0]-xbeam)/abs(pred[1]-ybeam)) for pred in mm_predictions])
    assert len(alpha_angle_obs)==len(observations.indices()), 'Size of alpha angles and observations are not equal %6.0f, %6.0f'%(len(alpha_angle_obs),len(observations.indices()))

    #Lorentz-polarization correction
    wavelength = observations_pickle["wavelength"]
    two_theta = observations.two_theta(wavelength=wavelength).data()
    one_over_LP = (2 * flex.sin(two_theta))/(1 + (flex.cos(two_theta)**2))
    one_over_P = 2/(1 + (flex.cos(two_theta)**2))
    observations = observations.customized_copy(data=observations.data()*one_over_P)

    #set observations with target space group - !!! required for correct
    #merging due to map_to_asu command.
    miller_set = symmetry(
      unit_cell=iph.target_unit_cell,
      space_group_symbol=iph.target_space_group
    ).build_miller_set(
      anomalous_flag=iph.target_anomalous_flag,
      d_min=iph.d_min)

    #Filter negative intensities
    i_I_positive = (observations.data() > 0)
    miller_indices_positive = observations.indices().select(i_I_positive)
    I_positive = observations.data().select(i_I_positive)
    sigI_positive = observations.sigmas().select(i_I_positive)
    alpha_angle_obs = alpha_angle_obs.select(i_I_positive)

    observations = observations.customized_copy(indices=miller_indices_positive,
        data=I_positive,
        sigmas=sigI_positive,
        anomalous_flag=iph.target_anomalous_flag,
        crystal_symmetry=miller_set.crystal_symmetry())

    if observations.crystal_symmetry().is_compatible_unit_cell() == False:
      return None

    #filter out weak data
    I_over_sigi = observations.data()/ observations.sigmas()
    i_I_obs_sel = (I_over_sigi > iph.sigma_max)
    observations = observations.customized_copy(indices=observations.indices().select(i_I_obs_sel),
        data=observations.data().select(i_I_obs_sel),
        sigmas=observations.sigmas().select(i_I_obs_sel),
        )
    alpha_angle_obs = alpha_angle_obs.select(i_I_obs_sel)

    #filter resolution
    i_sel_res = observations.resolution_filter_selection(d_max=iph.d_max, d_min=iph.d_min)
    observations = observations.customized_copy(indices=observations.indices().select(i_sel_res),
        data=observations.data().select(i_sel_res),
        sigmas=observations.sigmas().select(i_sel_res),
        )
    alpha_angle_obs = alpha_angle_obs.select(i_sel_res)

    assert len(alpha_angle_obs)==len(observations.indices()), 'Size of alpha angles and observations are not equal %6.0f, %6.0f'%(len(alpha_angle_obs),len(observations.indices()))

    return observations, alpha_angle_obs


  def determine_polar(self, observations_original, iph, pickle_filename):

    """
    Determine polarity based on input data.
    The function still needs isomorphous reference so, if flag_polar is True,
    miller_array_iso must be supplied in input file.
    """
    if iph.flag_polar == False:
      return 'h,k,l', 0 , 0

    #use basis in the given input file
    basis_pickle = pickle.load(open(iph.index_basis_in,"rb"))
    polar_hkl = basis_pickle[pickle_filename]
    return polar_hkl, 0, 0


  def calc_spot_radius(self, a_star_matrix, miller_indices, wavelength):
    #calculate spot_radius based on rms delta_S for all spots
    delta_S_all = flex.double()
    for miller_index in miller_indices:
      S0 = -1*col((0,0,1./wavelength))
      h = col(miller_index)
      x = a_star_matrix * h
      S = x + S0
      delta_S = S.length() - (1./wavelength)
      delta_S_all.append(delta_S)

    spot_radius = math.sqrt(flex.mean(delta_S_all*delta_S_all))

    return spot_radius

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

  def postrefine_by_frame(self, frame_no, pickle_filename, iph, miller_array_ref):

    #1. Prepare data
    observations_pickle = pickle.load(open(pickle_filename,"rb"))
    crystal_init_orientation = observations_pickle["current_orientation"][0]
    wavelength = observations_pickle["wavelength"]

    #grab img. name
    imgname = pickle_filename
    if iph.file_name_in_img != '':
      fh = file_handler()
      imgname = fh.get_imgname_from_pickle_filename(iph.file_name_in_img, pickle_filename)

    observations_original, alpha_angle_obs = self.organize_input(observations_pickle, iph)
    if observations_original is None:
      print frame_no, '-fail obs is none'
      return None

    #2. Determine polarity - always do this even if flag_polar = False
    #the function will take care of it.
    polar_hkl, cc_iso_raw_asu, cc_iso_raw_rev = self.determine_polar(observations_original, iph, pickle_filename)

    #3. Select data for post-refinement (only select indices that are common with the reference set
    observations_non_polar = self.get_observations_non_polar(observations_original, polar_hkl)
    matches = miller.match_multi_indices(
                  miller_indices_unique=miller_array_ref.indices(),
                  miller_indices=observations_non_polar.indices())

    I_ref_match = flex.double([miller_array_ref.data()[pair[0]] for pair in matches.pairs()])
    miller_indices_ref_match = flex.miller_index((miller_array_ref.indices()[pair[0]] for pair in matches.pairs()))
    I_obs_match = flex.double([observations_non_polar.data()[pair[1]] for pair in matches.pairs()])
    sigI_obs_match = flex.double([observations_non_polar.sigmas()[pair[1]] for pair in matches.pairs()])
    miller_indices_original_obs_match = flex.miller_index((observations_original.indices()[pair[1]] for pair in matches.pairs()))
    alpha_angle_set = flex.double([alpha_angle_obs[pair[1]] for pair in matches.pairs()])

    references_sel = miller_array_ref.customized_copy(data=I_ref_match, indices=miller_indices_ref_match)
    observations_original_sel = observations_original.customized_copy(data=I_obs_match,
          sigmas=sigI_obs_match,
          indices=miller_indices_original_obs_match)

    #4. Do least-squares refinement
    lsqrh = leastsqr_handler()
    refined_params, se_params, stats, partiality_sel, SE_I, var_I_p, var_k, var_p = lsqrh.optimize(I_ref_match, observations_original_sel,
              wavelength, crystal_init_orientation, alpha_angle_set, iph)

    if SE_I is None:
      print 'frame', frame_no, ' - failed'
      return None
    else:
      pres = postref_results()
      observations_non_polar_sel = self.get_observations_non_polar(observations_original_sel, polar_hkl)
      pres.set_params(observations = observations_non_polar_sel,
            refined_params=refined_params,
            se_params=se_params,
            stats=stats,
            partiality=partiality_sel,
            frame_no=frame_no,
            pickle_filename=pickle_filename,
            wavelength=wavelength,
            SE_I=SE_I,
            var_I_p=var_I_p,
            var_k=var_k,
            var_p=var_p)
      print 'frame %6.0f'%pres.frame_no, polar_hkl, ' SE=%7.2f R-sq=%7.2f CC=%7.2f'%pres.stats

    return pres

  def calc_mean_intensity(self, pickle_filename, iph):
    observations_pickle = pickle.load(open(pickle_filename,"rb"))
    observations_original, alpha_angle_obs = self.organize_input(observations_pickle, iph)
    if observations_original is None:
      return None

    mean_I = flex.mean(observations_original.data())
    return mean_I


  def scale_frame_by_mean_I(self, frame_no, pickle_filename, iph, mean_of_mean_I):

    observations_pickle = pickle.load(open(pickle_filename,"rb"))
    observations_original, alpha_angle_obs = self.organize_input(observations_pickle, iph)
    wavelength = observations_pickle["wavelength"]
    crystal_init_orientation = observations_pickle["current_orientation"][0]
    crystal_pointgroup = observations_pickle["pointgroup"]
    if observations_original is None:
      return None

    if iph.target_pointgroup != '' and crystal_pointgroup != iph.target_pointgroup:
      print 'frame %6.0f'%frame_no, ' - wrong pointgroup', crystal_pointgroup
      return None

    polar_hkl, cc_iso_raw_asu, cc_iso_raw_rev = self.determine_polar(observations_original, iph, pickle_filename)
    observations_non_polar = self.get_observations_non_polar(observations_original, polar_hkl)
    uc_params = crystal_init_orientation.unit_cell().parameters()

    G = np.mean(observations_non_polar.data())/mean_of_mean_I
    refined_params = np.array([G,0,0,0,0,0,0,uc_params[0],uc_params[1],uc_params[2],uc_params[3],uc_params[4],uc_params[5]])
    se_params = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0])
    stats = (0,0,0)
    partiality_sel = flex.double([1]*len(observations_non_polar.data()))
    SE_I = flex.double([1]*len(observations_non_polar.data()))
    var_I_p = flex.double([0]*len(observations_non_polar.data()))
    var_k = flex.double([0]*len(observations_non_polar.data()))
    var_p = flex.double([0]*len(observations_non_polar.data()))

    pres = postref_results()
    pres.set_params(observations = observations_non_polar,
            refined_params=refined_params,
            se_params=se_params,
            stats=stats,
            partiality=partiality_sel,
            frame_no=frame_no,
            pickle_filename=pickle_filename,
            wavelength=wavelength,
            SE_I=SE_I,
            var_I_p=var_I_p,
            var_k=var_k,
            var_p=var_p)

    print 'frame %6.0f'%frame_no, polar_hkl, '<I>=%9.2f <G>=%9.2f G=%9.2f'%(np.mean(observations_non_polar.data()), mean_of_mean_I, G)
    return pres


def merge_observations(results,
        iph,
        output_mtz_file_prefix,
        avg_mode):

  #results is a list of postref_results objects
  #lenght of this list equals to number of input frames

  inten_scaler = intensities_scaler()
  miller_array_merge, txt_out = inten_scaler.output_mtz_files(results, iph, output_mtz_file_prefix, avg_mode)
  return miller_array_merge, txt_out
