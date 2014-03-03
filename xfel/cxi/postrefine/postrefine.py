"""
Read pickle files output from integration process to determine the polarity and
refine the rotation matrix
Usage:
cxi.postrefine input=input.inp
"""
from __future__ import division

from cctbx.array_family import flex
from cctbx import miller
from cctbx.crystal import symmetry
from scitbx.matrix import sqr, col

from cctbx.crystal_orientation import crystal_orientation
from scitbx.lstbx import normal_eqns_solving

import matplotlib.pyplot as plt

import numpy as np
import os,cPickle as pickle,math

from mod_partiality import partiality_handler
from mod_util import intensities_scaler, file_handler, input_handler
from mod_normal_eqns import normal_eqns_handler

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

  def organize_input(self, observations,
      d_min,
      d_max,
      sigma_max,
      target_unit_cell,
      target_space_group,
      target_anomalous_flag,
      flag_polar,
      wavelength,
      miller_array_iso=None):

    """Given the supplied observations and Miller array, return a tuple of
    the original observations, the ASU-mapped observations, and the
    observations reindexed using (k, h, -l).  Merge Bijvoet mates unless
    anomalous_flag is True.
    """
    observations = observations.resolution_filter(d_min=d_min, d_max=d_max)

    #Lorentz-polarization correction
    two_theta = observations.two_theta(wavelength=wavelength).data()
    one_over_LP = (2 * flex.sin(two_theta))/(1 + (flex.cos(two_theta)**2))
    one_over_P = 2/(1 + (flex.cos(two_theta)**2))
    observations = observations.customized_copy(data=observations.data()*one_over_P)

    miller_set = symmetry(
        unit_cell=target_unit_cell,
        space_group_symbol=target_space_group
      ).build_miller_set(
        anomalous_flag=target_anomalous_flag,
        d_min=d_min)

    #Filter negative intensities
    i_I_positive = (observations.data() > 0)
    miller_indices_positive = observations.indices().select(i_I_positive)
    I_positive = observations.data().select(i_I_positive)
    sigI_positive = observations.sigmas().select(i_I_positive)

    observations = observations.customized_copy(indices=miller_indices_positive,
        data=I_positive,
        sigmas=sigI_positive,
        anomalous_flag=target_anomalous_flag,
        crystal_symmetry=miller_set.crystal_symmetry()
        )

    #Filter possible outliers (use keyword sigma_max=value to control the outliers)
    I_over_sigi = observations.data()/ observations.sigmas()
    i_I_obs_sel = (I_over_sigi > sigma_max)

    observations = observations.customized_copy(indices=observations.indices().select(i_I_obs_sel),
        data=observations.data().select(i_I_obs_sel),
        sigmas=observations.sigmas().select(i_I_obs_sel),
        anomalous_flag=target_anomalous_flag,
        crystal_symmetry=miller_set.crystal_symmetry()
        )

    #Prepare original, asu, and rev (for polar case)
    observations_asu = observations.customized_copy(
        anomalous_flag=target_anomalous_flag,
        crystal_symmetry=miller_set.crystal_symmetry()
        ).map_to_asu()

    observations_original = observations.customized_copy(
        anomalous_flag=target_anomalous_flag,
        crystal_symmetry=miller_set.crystal_symmetry()
        )

    if flag_polar:
      from cctbx import sgtbx
      cb_op = sgtbx.change_of_basis_op('k,h,-l')
      observations_rev = observations_asu.change_basis(cb_op).map_to_asu()
    else:
      observations_rev = observations_asu

    observations_original = observations_original.resolution_filter(d_min=d_min, d_max=d_max)
    observations_asu = observations_asu.resolution_filter(d_min=d_min, d_max=d_max)
    observations_rev = observations_rev.resolution_filter(d_min=d_min, d_max=d_max)

    return observations_original, observations_asu, observations_rev


  def determine_polar(self, observations_original, observations_asu, observations_rev, miller_array_iso, flag_polar):

    """
    Determine polarity based on input data.
    The function still needs isomorphous reference so, if flag_polar is True,
    miller_array_iso must be supplied in input file.
    """
    matches_asu = miller.match_multi_indices(
                  miller_indices_unique=miller_array_iso.indices(),
                  miller_indices=observations_asu.indices())

    I_ref_asu = flex.double([miller_array_iso.data()[pair[0]] for pair in matches_asu.pairs()])
    miller_indices_ref_asu = flex.miller_index((miller_array_iso.indices()[pair[0]] for pair in matches_asu.pairs()))
    I_obs_asu = flex.double([observations_asu.data()[pair[1]] for pair in matches_asu.pairs()])
    sigI_obs_asu = flex.double([observations_asu.sigmas()[pair[1]] for pair in matches_asu.pairs()])
    miller_indices_ori_asu = flex.miller_index((observations_original.indices()[pair[1]] for pair in matches_asu.pairs()))

    corr_raw_asu, slope_raw_asu = get_overall_correlation(I_obs_asu, I_ref_asu)

    matches_rev = miller.match_multi_indices(
                  miller_indices_unique=miller_array_iso.indices(),
                  miller_indices=observations_rev.indices())

    I_ref_rev = flex.double([miller_array_iso.data()[pair[0]] for pair in matches_rev.pairs()])
    miller_indices_ref_rev = flex.miller_index((miller_array_iso.indices()[pair[0]] for pair in matches_rev.pairs()))
    I_obs_rev = flex.double([observations_rev.data()[pair[1]] for pair in matches_rev.pairs()])
    sigI_obs_rev = flex.double([observations_rev.sigmas()[pair[1]] for pair in matches_rev.pairs()])
    miller_indices_ori_rev = flex.miller_index((observations_original.indices()[pair[1]] for pair in matches_rev.pairs()))

    corr_raw_rev, slope_raw_rev = get_overall_correlation(I_obs_rev, I_ref_rev)


    if corr_raw_asu >= corr_raw_rev:
      polar_hkl = 'h,k,l'
    else:
      polar_hkl = 'k,h,-l'

    if flag_polar == False:
      polar_hkl = 'h,k,l'



    return polar_hkl, corr_raw_asu, corr_raw_rev


  def get_pickle_info(self, pickle_filename):
    """
    #in case there is a reference matrix, calculate the angular error.
    """

    trial_results = pickle.load(open(pickle_filename,"rb"))
    observations = trial_results["observations"][0]
    current_hexagonal_ori = trial_results["current_orientation"][0]
    wavelength = trial_results["wavelength"]

    current_cb_op_to_primitive = trial_results["current_cb_op_to_primitive"][0]
    current_triclinic_ori = current_hexagonal_ori.change_basis(current_cb_op_to_primitive)

    crystal_pointgroup = trial_results["pointgroup"]

    CONTAINER_SZ=1000
    mosflm = sqr((1,0,0,0,1,0,0,0,1)) # mosflm lab vectors in their own frame
    labelit = sqr((0,0,-1,-1,0,0,0,1,0)) # mosflm basis vectors in labelit frame
    LM = mosflm * labelit.inverse() # converts labelit frame coords to mosflm frame
    SWAPXY = sqr((0,1,0,1,0,0,0,0,1)) # in labelit frame
    SWAPZ  = sqr((1,0,0,0,1,0,0,0,-1)) # in labelit frame
    R90    = sqr((0,-1,0,1,0,0,0,0,1)) # in labelit frame, rotation 90 on beam axis

    i_found_item = pickle_filename.find('/int-data_')
    item_no = int(pickle_filename[i_found_item+10:i_found_item+15].strip())
    container_no = item_no//CONTAINER_SZ


    filename = "/net/viper/raid1/sauter/fake/holton/mosflm_matrix"
    filename = os.path.join( filename, "%02d"%container_no, "%05d.mat"%item_no)


    lines = open(filename).readlines()

    A0 = lines[0].strip().split()
    A1 = lines[1].strip().split()
    A2 = lines[2].strip().split()
    A = sqr((float(A0[0]), float(A0[1]), float(A0[2]),
               float(A1[0]), float(A1[1]), float(A1[2]),
               float(A2[0]), float(A2[1]), float(A2[2])))


    A = A/wavelength
    Holton_hexagonal_ori = crystal_orientation(SWAPZ*SWAPXY*R90*LM*A,True)
    Holton_triclinic_ori = Holton_hexagonal_ori.change_basis(current_cb_op_to_primitive)

    c_inv_r_best = Holton_triclinic_ori.best_similarity_transformation(
                   other=current_triclinic_ori,
                   fractional_length_tolerance=50.,
                   unimodular_generator_range=1)
    c_inv_r_int = tuple([int(round(ij,0)) for ij in c_inv_r_best])
    from cctbx import sgtbx
    c_inv = sgtbx.rt_mx(sgtbx.rot_mx(c_inv_r_int))
    cb_op = sgtbx.change_of_basis_op(c_inv)
    comparison_triclinic = Holton_triclinic_ori.change_basis(cb_op)
    comparison_hexagonal = comparison_triclinic.change_basis(current_cb_op_to_primitive.inverse())

    return current_hexagonal_ori, comparison_hexagonal, wavelength, observations, crystal_pointgroup

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


  def get_refined_observations(self, refine_mode,
        refined_parameters,
        observations,
        iph,
        wavelength,
        polar_hkl,
        miller_array_ref,
        crystal_init_orientation=None):


    observations_original, observations_asu, observations_rev = self.organize_input(observations,
        iph.d_min_merge,
        iph.d_max,
        iph.sigma_max,
        iph.target_unit_cell,
        iph.target_space_group,
        iph.target_anomalous_flag,
        iph.flag_polar,
        wavelength)

    observations_sel = observations_asu.deep_copy()
    if iph.flag_polar:
      if polar_hkl == 'k,h,-l':
        observations_sel = observations_rev.deep_copy()
    observations_sel_two_theta = observations_sel.two_theta(wavelength=wavelength)
    two_theta = observations_sel_two_theta.data()
    sin_theta_over_lambda_sq = observations_sel_two_theta.sin_theta_over_lambda_sq().data()

    if refine_mode == 'scale_factor':
      G_best = refined_parameters[0]
      B_factor_best = refined_parameters[1]
      I_obs_refined = (G_best * (flex.exp(-2*B_factor_best*sin_theta_over_lambda_sq)) * observations_sel.data())/ observations_sel.sigmas()
      sigI_obs_refined = (G_best * (flex.exp(-2*B_factor_best*sin_theta_over_lambda_sq)) * observations_sel.sigmas())

      observations_weight = 1/observations_sel.sigmas()
      observations_refined = observations_sel.customized_copy(data=I_obs_refined,
                  sigmas=sigI_obs_refined)

    elif (refine_mode == 'crystal_rotation' or refine_mode == 'reflecting_range'):
      G_best = refined_parameters[0]
      B_factor_best = refined_parameters[1]
      rotx_best = refined_parameters[2]
      roty_best = refined_parameters[3]
      ry_best = refined_parameters[4]
      rz_best = refined_parameters[5]

      effective_orientation = crystal_init_orientation.rotate_thru((1,0,0), rotx_best
               ).rotate_thru((0,1,0), roty_best)
      effective_a_star = sqr(effective_orientation.reciprocal_matrix())
      ph = partiality_handler(wavelength, 0)
      partiality_final = ph.calc_partiality_anisotropy_set(effective_a_star, observations_original.indices(), ry_best, rz_best, two_theta)

      I_obs_refined = (G_best * (flex.exp(-2*B_factor_best*sin_theta_over_lambda_sq)) * observations_sel.data())/ (partiality_final * observations_sel.sigmas())
      sigI_obs_refined = (G_best * (flex.exp(-2*B_factor_best*sin_theta_over_lambda_sq)) * observations_sel.sigmas())/partiality_final
      observations_weight = 1/observations_sel.sigmas()
      observations_refined = observations_sel.customized_copy(data=I_obs_refined,
                sigmas=sigI_obs_refined)

    #calculate cc_ref, r_ref
    matches = miller.match_multi_indices(
                  miller_indices_unique=miller_array_ref.indices(),
                  miller_indices=observations_refined.indices())

    I_ref = flex.double([miller_array_ref.data()[pair[0]] for pair in matches.pairs()])
    I_obs = flex.double([observations_sel.data()[pair[1]] for pair in matches.pairs()])
    I_obs_refined_match = flex.double([I_obs_refined[pair[1]] for pair in matches.pairs()])
    I_weight_match = flex.double([observations_weight[pair[1]] for pair in matches.pairs()])
    I_obs_weight_match = I_obs_refined_match/ I_weight_match

    cc_ref_init, slope_ref_init = get_overall_correlation(I_obs, I_ref)
    cc_ref_final, slope_ref_final = get_overall_correlation(I_obs_weight_match, I_ref)
    r_ref_init = sum(flex.abs(I_ref - (I_obs*slope_ref_init)))/sum(I_obs*slope_ref_init)
    r_ref_final = sum(flex.abs(I_ref - (I_obs_weight_match*slope_ref_final)))/sum(I_obs_weight_match*slope_ref_final)

    """
    plt.subplot(211)
    plt.scatter(I_ref, I_obs,s=10, marker='x', c='r')
    plt.title('CC=%6.5f Slope=%6.5f R=%6.5f'%(cc_ref_init, slope_ref_init, r_ref_init))
    plt.xlabel('Reference intensity')
    plt.ylabel('Observed intensity')
    plt.subplot(212)
    plt.scatter(I_ref, I_obs_weight_match,s=10, marker='x', c='r')
    plt.title('CC=%6.5f Slope=%6.5f R=%6.5f Refined G=%6.5f'%(cc_ref_final, slope_ref_final, r_ref_final, G_best))
    plt.xlabel('Reference intensity')
    plt.ylabel('Observed intensity')
    plt.show()
    """

    #calculate cc_iso, r_iso
    cc_iso_init = 0
    cc_iso_final = 0
    r_iso_init = 0
    r_iso_final = 0
    if iph.miller_array_iso is not None:
      matches = miller.match_multi_indices(
                    miller_indices_unique=iph.miller_array_iso.indices(),
                    miller_indices=observations_refined.indices())
      I_ref = flex.double([iph.miller_array_iso.data()[pair[0]] for pair in matches.pairs()])
      I_obs = flex.double([observations_sel.data()[pair[1]] for pair in matches.pairs()])
      I_obs_refined_match = flex.double([I_obs_refined[pair[1]] for pair in matches.pairs()])
      I_weight_match = flex.double([observations_weight[pair[1]] for pair in matches.pairs()])
      I_obs_weight_match = I_obs_refined_match/ I_weight_match

      cc_iso_init, slope_iso_init = get_overall_correlation(I_obs, I_ref)
      cc_iso_final, slope_iso_final = get_overall_correlation(I_obs_weight_match, I_ref)
      r_iso_init = sum(flex.abs(I_ref - (I_obs*slope_iso_init)))/sum(I_obs*slope_iso_init)
      r_iso_final = sum(flex.abs(I_ref - (I_obs_weight_match*slope_iso_final)))/sum(I_obs_weight_match*slope_iso_final)

    cc_r_results = (cc_ref_init, cc_ref_final, cc_iso_init, cc_iso_final,
        r_ref_init, r_ref_final, r_iso_init, r_iso_final)


    return observations_refined, observations_weight, cc_r_results


  def postrefine_by_frame(self, refine_mode, pickle_filename, iph, miller_array_ref,
        scale_factors=None, rotations=None, reflecting_ranges=None):
    """
    Main module to support multi-processor
    """
    #1. Prepare data
    if iph.flag_reference_a_matrix:
      crystal_init_orientation, crystal_compare_orientation, wavelength, observations, crystal_pointgroup = self.get_pickle_info(pickle_filename)
      target_unit_cell = iph.target_unit_cell
    else:
      trial_results = pickle.load(open(pickle_filename,"rb"))
      crystal_init_orientation = trial_results["current_orientation"][0]
      crystal_compare_orientation = trial_results["current_orientation"][0]
      wavelength = trial_results["wavelength"]
      crystal_pointgroup = trial_results["pointgroup"]
      unit_cell = trial_results["current_orientation"][0].unit_cell()
      target_unit_cell = unit_cell.parameters()
      observations = trial_results["observations"][0]


    #grab img. name
    imgname = pickle_filename
    if iph.file_name_in_img != '':
      fh = file_handler()
      imgname = fh.get_imgname_from_pickle_filename(iph.file_name_in_img, pickle_filename)

    observations_original, observations_asu, observations_rev = self.organize_input(observations,
        iph.d_min,
        iph.d_max,
        iph.sigma_max,
        target_unit_cell,
        iph.target_space_group,
        iph.target_anomalous_flag,
        iph.flag_polar,
        wavelength,
        miller_array_iso=iph.miller_array_iso)

    #2. Determine polarity - always do this even if flag_polar = False
    #the function will take care of it.
    polar_hkl, cc_iso_raw_asu, cc_iso_raw_rev = self.determine_polar(observations_original,
          observations_asu, observations_rev, iph.miller_array_iso, iph.flag_polar)


    if polar_hkl == 'h,k,l':
      observations_non_polar = observations_asu
    elif polar_hkl == 'k,h,-l':
      observations_non_polar = observations_rev

    #3. Select data for post-refinement (only select indices that are common with the reference set
    matches = miller.match_multi_indices(
                  miller_indices_unique=miller_array_ref.indices(),
                  miller_indices=observations_non_polar.indices())

    I_ref_match = flex.double([miller_array_ref.data()[pair[0]] for pair in matches.pairs()])
    miller_indices_ref_match = flex.miller_index((miller_array_ref.indices()[pair[0]] for pair in matches.pairs()))
    I_obs_match = flex.double([observations_non_polar.data()[pair[1]] for pair in matches.pairs()])
    sigI_obs_match = flex.double([observations_non_polar.sigmas()[pair[1]] for pair in matches.pairs()])
    miller_indices_original_obs_match = flex.miller_index((observations_original.indices()[pair[1]] for pair in matches.pairs()))

    references_sel = miller_array_ref.customized_copy(data=I_ref_match, indices=miller_indices_ref_match)
    observations_original_sel = observations_original.customized_copy(data=I_obs_match,
          sigmas=sigI_obs_match,
          indices=miller_indices_original_obs_match,
          anomalous_flag=iph.target_anomalous_flag)

    I_ref = references_sel.data()
    I_obs = observations_original_sel.data()
    sigI_obs = observations_original_sel.sigmas()

    #for B-factor weighting, precalculate (sin(2theta)/lambda)^2
    observations_original_sel_two_theta = observations_original_sel.two_theta(wavelength=wavelength)
    two_theta = observations_original_sel_two_theta.data()
    sin_theta_over_lambda_sq = observations_original_sel_two_theta.sin_theta_over_lambda_sq().data()

    #4. Do least-squares refinement using different refine_mode:
    #scale_factor    -> G0 and B-factor
    #crystal_setting -> rotx, roty, ry, rz, uc (a,b,c,alpha,beta,gamma)

    cc_ref_init, G_init = get_overall_correlation(I_obs, I_ref)
    corr_threshold_postref = 0.25
    if (cc_ref_init < corr_threshold_postref):
      return None
    else:
      if refine_mode == 'scale_factor':
        B_factor_init = 0
        grad_thres = 1e-5
        parameters = (G_init, B_factor_init)

        neh = normal_eqns_handler()
        helper = neh.get_helper_refine_scale(I_ref, observations_original_sel,
                    wavelength, parameters)
        helper.restart()
        iterations = normal_eqns_solving.naive_iterations(
                     non_linear_ls = helper,
                     gradient_threshold = grad_thres)
        lstsqr_results =  helper.x

        observations_refined, observations_weight, cc_r_results = self.get_refined_observations(refine_mode,
            lstsqr_results,
            observations,
            iph,
            wavelength,
            polar_hkl,
            miller_array_ref)

      elif refine_mode == 'crystal_rotation':
        rotx_init = 0
        roty_init = 0
        grad_thres = 1e-5
        a_star_init = sqr(crystal_init_orientation.reciprocal_matrix())
        spot_radius = self.calc_spot_radius(a_star_init, observations_original_sel.indices(), wavelength)
        parameters = (rotx_init, roty_init)
        neh = normal_eqns_handler()
        helper = neh.get_helper_refine_crystal_rotation(I_ref, observations_original_sel,
                  wavelength, parameters, scale_factors, spot_radius, crystal_init_orientation)
        iterations = normal_eqns_solving.naive_iterations(
                     non_linear_ls = helper,
                     gradient_threshold = grad_thres)
        lstsqr_results =  helper.x

        observations_refined, observations_weight, cc_r_results = self.get_refined_observations(refine_mode,
            (scale_factors[0], scale_factors[1], lstsqr_results[0], lstsqr_results[1], spot_radius, spot_radius),
            observations,
            iph,
            wavelength,
            polar_hkl,
            miller_array_ref,
            crystal_init_orientation=crystal_init_orientation)

      elif refine_mode == 'reflecting_range':
        grad_thres = 1e-5
        a_star_init = sqr(crystal_init_orientation.reciprocal_matrix())
        spot_radius = self.calc_spot_radius(a_star_init, observations_original_sel.indices(), wavelength)
        parameters = (spot_radius, spot_radius)
        neh = normal_eqns_handler()
        helper = neh.get_helper_refine_reflecting_range(I_ref, observations_original_sel,
                  wavelength, parameters, scale_factors, rotations, crystal_init_orientation)
        iterations = normal_eqns_solving.naive_iterations(
                     non_linear_ls = helper,
                     gradient_threshold = grad_thres)
        lstsqr_results =  helper.x

        observations_refined, observations_weight, cc_r_results = self.get_refined_observations(refine_mode,
            (scale_factors[0], scale_factors[1], rotations[0], rotations[1], lstsqr_results[0], lstsqr_results[1]),
            observations,
            iph,
            wavelength,
            polar_hkl,
            miller_array_ref,
            crystal_init_orientation=crystal_init_orientation)


    cc_ref_init, cc_ref_final, cc_iso_init, cc_iso_final, r_ref_init, r_ref_final, r_iso_init, r_iso_final = cc_r_results

    _tmp_imgname = imgname.split('/')
    print _tmp_imgname[len(_tmp_imgname)-1], '%4.0f %4.0f %2.4f %2.4f %2.4f %2.4f'%(len(observations.indices()), len(observations_refined.indices()), cc_ref_init, cc_ref_final, r_ref_init, r_ref_final)

    return pickle_filename, observations_refined, observations_weight, cc_r_results, lstsqr_results


  def calc_mean_intensity(self, pickle_filename, iph):
    trial_results = pickle.load(open(pickle_filename,"rb"))
    observations = trial_results["observations"][0]
    wavelength = trial_results["wavelength"]

    observations_original, observations_asu, observations_rev = self.organize_input(observations,
        iph.d_min_merge,
        iph.d_max,
        iph.sigma_max,
        iph.target_unit_cell,
        iph.target_space_group,
        iph.target_anomalous_flag,
        iph.flag_polar,
        wavelength)

    mean_I = flex.mean(observations_original.data())

    return mean_I


  def scale_frame_by_mean_I(self, pickle_filename, iph, mean_of_mean_I):

    trial_results = pickle.load(open(pickle_filename,"rb"))
    observations = trial_results["observations"][0]
    wavelength = trial_results["wavelength"]
    observations_original, observations_asu, observations_rev = self.organize_input(observations,
        iph.d_min_merge,
        iph.d_max,
        iph.sigma_max,
        iph.target_unit_cell,
        iph.target_space_group,
        iph.target_anomalous_flag,
        iph.flag_polar,
        wavelength)
    polar_hkl, cc_iso_raw_asu, cc_iso_raw_rev = self.determine_polar(observations_original,
          observations_asu, observations_rev, iph.miller_array_iso, iph.flag_polar)

    if polar_hkl == 'h,k,l':
      observations_sel = observations_asu
    elif polar_hkl == 'k,h,-l':
      observations_sel = observations_rev


    this_G = mean_of_mean_I/ np.mean(observations_sel.data())
    observations_scaled = observations_sel.customized_copy(
        data = this_G * (observations_sel.data()/ observations_sel.sigmas()),
        sigmas = this_G * observations_sel.sigmas())

    observations_weight = this_G**2/ observations_sel.sigmas()
    return pickle_filename, observations_sel, observations_weight




class ScoringContainer:
  """
  Compare two rotation matrix and report angular error (Nick's legacy code)
  """
  pass

  def angular_rotation(self):
    """Compare the two rotation operators self.model to self.reference and
    return the angle (in radians) between them.
    """

    #print "model"
    model_unit = self.direct_matrix_as_unit_vectors(self.model)
    #print "reference"
    reference_unit = self.direct_matrix_as_unit_vectors(self.reference)
    rotation = model_unit*reference_unit.inverse()
    UQ = rotation.r3_rotation_matrix_as_unit_quaternion()
    UQ = UQ.normalize() # bugfix; without this many evaluations come to 0.00000 degrees
    angle, axis = UQ.unit_quaternion_as_axis_and_angle()
    #print "axis length",axis.length()
    #print "axis %7.4f %7.4f %7.4f"%(axis[0],axis[1],axis[2])
    return angle

  def direct_matrix_as_unit_vectors(self, ori):
    """Return ori's direct_matrix(), with normalized rows.
    """
    direct = sqr(ori.direct_matrix())
    A = col((direct[0],direct[1],direct[2])).normalize()
    B = col((direct[3],direct[4],direct[5])).normalize()
    C = col((direct[6],direct[7],direct[8])).normalize()
    #print "gamma deg",math.acos(A.dot(B))*180./math.pi
    #print "alpha deg",math.acos(B.dot(C))*180./math.pi
    #print " beta deg",math.acos(C.dot(A))*180./math.pi
    direct_as_unit = sqr((A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2]))
    return direct_as_unit

def get_overall_correlation (data_a, data_b) :
  """
  Correlate the averaged intensities to the intensities from the
  reference data set.
  """

  assert len(data_a) == len(data_b)
  corr = 0
  slope = 0
  try:
    sum_xx = 0
    sum_xy = 0
    sum_yy = 0
    sum_x  = 0
    sum_y  = 0
    N      = 0
    for i in xrange(len(data_a)):
      I_r       = data_a[i]
      I_o       = data_b[i]
      N      += 1
      sum_xx += I_r**2
      sum_yy += I_o**2
      sum_xy += I_r * I_o
      sum_x  += I_r
      sum_y  += I_o
    slope = (N * sum_xy - sum_x * sum_y) / (N * sum_xx - sum_x**2)
    corr  = (N * sum_xy - sum_x * sum_y) / (math.sqrt(N * sum_xx - sum_x**2) *
               math.sqrt(N * sum_yy - sum_y**2))
  except:
    pass

  return corr, slope

def get_observations (dir_name,data_subset):
  file_names = []
  for file_name in os.listdir(dir_name):
    if (file_name.endswith("_00000.pickle")):
      if data_subset==0 or \
        (data_subset==1 and (int(os.path.basename(file_name).split("_00000.pickle")[0][-1])%2==1)) or \
        (data_subset==2 and (int(os.path.basename(file_name).split("_00000.pickle")[0][-1])%2==0)):
        file_names.append(os.path.join(dir_name, file_name))
    elif (file_name.endswith(".pickle")):
      if data_subset==0 or \
        (data_subset==1 and (int(os.path.basename(file_name).split(".pickle")[0][-1])%2==1)) or \
        (data_subset==2 and (int(os.path.basename(file_name).split(".pickle")[0][-1])%2==0)):
        file_names.append(os.path.join(dir_name, file_name))
  print "Number of pickle files found:", len(file_names)
  return file_names

def merge_observations(observations_set,
        observations_weight_set,
        iph,
        output_mtz_file_prefix):

  inten_scaler = intensities_scaler()
  miller_array_merge, cc_merge, slope_merge, txt_output_mtz = inten_scaler.output_mtz_files(observations_set,
    observations_weight_set, iph, output_mtz_file_prefix)

  return miller_array_merge, cc_merge, slope_merge, txt_output_mtz
