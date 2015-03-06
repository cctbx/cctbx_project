'''
Author      : Uervirojnangkoorn, M.
Created     : 7/13/2014
Description :

The leastsqr_handler class refines post-refinement parameters for each frame m in:

J = sum(I_hi - (G0 * exp(-2B*sin_theta_over_lambda_sq) * p * I_h))^2

where I_hi is observed intensitiy, G0 and B are scale factors, and p
is the partiality function of these parameters:
rot_x, rot_y: crystal orientation angles around x- and y- axes
gamma_y, gamma_z: parameters to model anisotropy of crystal mosaicity
gamma_e: spectral dispersion
unit-cell parameters: a,b,c,alpha,beta,gamma.

The class implements Leveberg-Marquardt algorithms to scale the refined parameters
using the lamda updates. The unit-cell parameters are refined with restraints based
on the 7 crystal systems (6 conditions).
'''
from __future__ import division
from scipy import optimize
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from cctbx.array_family import flex
from scitbx.matrix import sqr, col
from cctbx.uctbx import unit_cell
from cctbx.crystal_orientation import crystal_orientation, basis_type

def calc_full_refl(I_o_p_set, sin_theta_over_lambda_sq_set,
                   G, B, p_set, rs_set, flag_volume_correction):
  I_o_full_set = ((G * np.exp(-2*B*sin_theta_over_lambda_sq_set) * I_o_p_set)/p_set)
  if flag_volume_correction:
    I_o_full_set = I_o_full_set * (4/3) * (rs_set)

  I_o_full_set = flex.double(I_o_full_set)

  return I_o_full_set

def calc_spot_radius(a_star_matrix, miller_indices, wavelength):
  #calculate spot_radius based on rms delta_S for all spots
  delta_S_all = flex.double()
  for miller_index in miller_indices:
    S0 = -1*col((0,0,1./wavelength))
    h = col(miller_index)
    x = a_star_matrix * h
    S = x + S0
    delta_S = S.length() - (1./wavelength)
    delta_S_all.append(delta_S)

  #spot_radius = math.sqrt(flex.mean(delta_S_all*delta_S_all))
  spot_radius = np.std(delta_S_all)*1

  return spot_radius

def coefficient_of_determination(y, y_model):
  mean_y = np.mean(y)
  r_sqr = np.sum((y_model - mean_y)**2)/np.sum((y - mean_y)**2)
  return r_sqr

def standard_error_of_the_estimate(y, y_model, n_params):
  s = np.sqrt(np.sum((y - y_model)**2)/(len(y) - n_params))
  return s

def calc_partiality_anisotropy_set(my_uc, rotx, roty, miller_indices, ry, rz, r0, re,
    bragg_angle_set, alpha_angle_set, wavelength, crystal_init_orientation,
    spot_pred_x_mm_set, spot_pred_y_mm_set, detector_distance_mm,
    partiality_model, flag_beam_divergence):
  #use III.4 in Winkler et al 1979 (A35; P901) for set of miller indices
  O = sqr(my_uc.orthogonalization_matrix()).transpose()
  R = sqr(crystal_init_orientation.crystal_rotation_matrix()).transpose()
  CO = crystal_orientation(O*R, basis_type.direct)
  CO_rotate = CO.rotate_thru((1,0,0), rotx
               ).rotate_thru((0,1,0), roty)
  A_star = sqr(CO_rotate.reciprocal_matrix())
  S0 = -1*col((0,0,1./wavelength))
  partiality_set = flex.double()
  delta_xy_set = flex.double()
  rs_set = flex.double()
  rh_set = flex.double()
  for miller_index, bragg_angle, alpha_angle, spot_pred_x_mm, spot_pred_y_mm in \
      zip(miller_indices, bragg_angle_set, alpha_angle_set,
          spot_pred_x_mm_set, spot_pred_y_mm_set):
    if flag_beam_divergence:
      rs = math.sqrt((ry * math.cos(alpha_angle))**2 + (rz * math.sin(alpha_angle))**2) + \
        (r0 + (abs(re)*math.tan(bragg_angle)))
    else:
      rs = r0 + (abs(re)*math.tan(bragg_angle))
    h = col(miller_index)
    x = A_star * h
    S = x + S0
    rh = S.length() - (1/wavelength)


    if partiality_model == 'Lorentzian':
      #Lorenzian
      spot_partiality = ((rs**2)/((2*(rh**2))+(rs**2)))
    elif partiality_model == 'Disc':
      #Disc
      if abs(rs) - abs(rh) > 0:
        spot_partiality = 1-(rh**2/rs**2)
      else:
        spot_partiality = 0.02
    elif partiality_model == 'Kabsch':
      #Kabsch
      t = rh/(math.sqrt(2)*rs)
      spot_partiality = math.exp(-(t**2))
    elif partiality_model == 'Rossmann':
      if abs(rs) - abs(rh) > 0:
        q=(rs-rh)/(2*rs);
        spot_partiality = (3*(q**2)) - (2*(q**3));
      else:
        spot_partiality = 0.02

    partiality_set.append(spot_partiality)
    rs_set.append(rs)
    rh_set.append(rh)

    #finding coordinate x,y on the detector
    d_ratio = -detector_distance_mm/S[2]
    dx_mm = S[0]*d_ratio
    dy_mm = S[1]*d_ratio
    pred_xy = col((spot_pred_x_mm, spot_pred_y_mm))
    calc_xy = col((dx_mm, dy_mm))
    diff_xy = pred_xy - calc_xy
    delta_xy_set.append(diff_xy.length())

  return partiality_set, delta_xy_set, rs_set, rh_set

def prep_input(params, cs):
  #From crystal system cs, determine refined parameters
  a, b, c, alpha, beta, gamma = params
  if cs == 'Triclinic':
    x0 = params
  elif cs == 'Monoclinic':
    x0 = np.array([a,b,c,beta])
  elif cs == 'Orthorhombic':
    x0 = np.array([a,b,c])
  elif cs == 'Tetragonal':
    x0 = np.array([a,c])
  elif cs == 'Trigonal' or cs == 'Hexagonal':
    x0 = np.array([a,c])
  elif cs == 'Cubic':
    x0 = np.array([a])

  return x0

def prep_output(params, cs):
  if cs == 'Triclinic':
    xopt = params
  elif cs == 'Monoclinic':
    xopt = np.array([params[0],params[1],params[2],90,params[3],90])
  elif cs == 'Orthorhombic':
    xopt = np.array([params[0],params[1],params[2],90,90,90])
  elif cs == 'Tetragonal':
    xopt = np.array([params[0],params[0],params[1],90,90,90])
  elif cs == 'Trigonal' or cs == 'Hexagonal':
    xopt = np.array([params[0],params[0],params[1],90,90,120])
  elif cs == 'Cubic':
    xopt = np.array([params[0],params[0],params[0],90,90,90])

  return xopt

def func(params, *args):
  I_r = args[0]
  miller_array_o = args[1]
  wavelength = args[2]
  alpha_angle_set = args[3]
  crystal_init_orientation = args[4]
  spot_pred_x_mm_set = args[5]
  spot_pred_y_mm_set = args[6]
  detector_distance_mm = args[7]
  refine_mode = args[8]
  const_params = args[9]
  b0 = args[10]
  miller_array_iso = args[11]
  iparams = args[12]

  partiality_model = iparams.partiality_model
  flag_volume_correction = iparams.flag_volume_correction
  flag_beam_divergence = iparams.flag_beam_divergence
  b_refine_d_min = iparams.b_refine_d_min
  I_o = miller_array_o.data().as_numpy_array()
  sigI_o = miller_array_o.sigmas().as_numpy_array()
  miller_indices_original = miller_array_o.indices()
  sin_theta_over_lambda_sq = miller_array_o.two_theta(wavelength=wavelength).sin_theta_over_lambda_sq().data().as_numpy_array()
  two_theta_flex = miller_array_o.two_theta(wavelength=wavelength).data()
  cs = miller_array_o.crystal_symmetry().space_group().crystal_system()

  if refine_mode == 'scale_factor':
    G, B = params
    rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma = const_params
  elif refine_mode == 'crystal_orientation':
    rotx, roty = params
    G, B, ry, rz, r0, re, a, b, c, alpha, beta, gamma = const_params
  elif refine_mode == 'reflecting_range':
    ry, rz, r0, re = params
    G, B, rotx, roty, a, b, c, alpha, beta, gamma = const_params
  elif refine_mode == 'unit_cell':
    a, b, c, alpha, beta, gamma = prep_output(params, cs)
    G, B, rotx, roty, ry, rz, r0, re = const_params
  elif refine_mode == 'allparams':
    a, b, c, alpha, beta, gamma = prep_output(params[8:], cs)
    G, B, rotx, roty, ry, rz, r0, re = params[0:8]

  try:
    uc = unit_cell((a,b,c,alpha,beta,gamma))
  except Exception:
    return None

  p_calc_flex, delta_xy_flex, rs_set, dummy = calc_partiality_anisotropy_set(\
    uc, rotx, roty, miller_indices_original, ry, rz, r0, re, two_theta_flex,
    alpha_angle_set, wavelength, crystal_init_orientation, spot_pred_x_mm_set,
    spot_pred_y_mm_set, detector_distance_mm, partiality_model, flag_beam_divergence)
  p_calc_set = p_calc_flex.as_numpy_array()

  if miller_array_o.d_min() < b_refine_d_min:
    I_o_full = calc_full_refl(I_o, sin_theta_over_lambda_sq,
                              G, B, p_calc_set, rs_set, flag_volume_correction)

  else:
    I_o_full = calc_full_refl(I_o, sin_theta_over_lambda_sq,
                              G, b0, p_calc_set, rs_set, flag_volume_correction)

  if refine_mode == 'unit_cell':
    error = delta_xy_flex
  else:
    error = ((I_r - I_o_full)/sigI_o)

  CC_ref = np.corrcoef(I_r, I_o_full)[0,1]
  CC_iso = 0
  I_o_match = flex.double()
  if miller_array_iso is not None:
    miller_array_o_asu = miller_array_o.map_to_asu()

    from cctbx import miller
    matches = miller.match_multi_indices(
                      miller_indices_unique=miller_array_iso.indices(),
                      miller_indices=miller_array_o_asu.indices())
    I_iso_match = flex.double([miller_array_iso.data()[pair[0]] for pair in matches.pairs()])
    I_o_match = flex.double([I_o_full[pair[1]] for pair in matches.pairs()])

    CC_iso = np.corrcoef(I_iso_match, I_o_match)[0,1]

  #print refine_mode, 'G=%.4g B=%.4g rotx=%.4g roty=%.4g ry=%.4g rz=%.4g r0=%.4g re=%.4g a=%.4g b=%.4g c=%.4g alp=%.4g beta=%.4g gam=%.4g fpr=%.4g fxy=%.4g CCref=%5.2f CCiso=%5.2f dd_mm=%6.2f n_refl=%5.0f n_iso=%5.0f'%(G, B, rotx*180/math.pi, roty*180/math.pi, ry, rz, r0, re, a, b, c, alpha, beta, gamma, np.sum(((I_r - I_o_full)/sigI_o)**2), np.sum(delta_xy_flex**2), CC_ref*100, CC_iso*100, detector_distance_mm, len(I_o_full), len(I_o_match))
  #print refine_mode+' %10.2f %10.2f'%(np.sum(((I_r - I_o_full)/sigI_o)**2), np.sum(delta_xy_flex**2))
  return error

def good_unit_cell(uc_params, iparams, uc_tol):

  flag_good_uc = False
  if (abs(uc_params[0]-iparams.target_unit_cell.parameters()[0]) \
      <= (uc_tol*iparams.target_unit_cell.parameters()[0]/100) \
                and abs(uc_params[1]-iparams.target_unit_cell.parameters()[1]) \
                <= (uc_tol*iparams.target_unit_cell.parameters()[1]/100) \
                and abs(uc_params[2]-iparams.target_unit_cell.parameters()[2]) \
                <= (uc_tol*iparams.target_unit_cell.parameters()[2]/100) \
                and abs(uc_params[3]-iparams.target_unit_cell.parameters()[3]) \
                <= (uc_tol*iparams.target_unit_cell.parameters()[3]/100) \
                and abs(uc_params[4]-iparams.target_unit_cell.parameters()[4]) \
                <= (uc_tol*iparams.target_unit_cell.parameters()[4]/100) \
                and abs(uc_params[5]-iparams.target_unit_cell.parameters()[5]) \
                <= (uc_tol*iparams.target_unit_cell.parameters()[5]/100)):
    flag_good_uc = True
  return flag_good_uc

class leastsqr_handler(object):
  '''
  A wrapper class for least-squares refinement
  '''

  def __init__(self):
    '''
    Intialitze parameters
    '''

  def get_filtered_data(self, filter_mode, filter_params,
                        observations_in, alpha_angle_in,
                        spot_pred_x_mm_in, spot_pred_y_mm_in,
                        I_ref_in, partiality_in=False):
    if filter_mode == 'resolution':
      i_sel = observations_in.resolution_filter_selection(d_min=filter_params[0],
                                                                    d_max=filter_params[1])
    elif filter_mode == 'sigma':
      i_sel = (observations_in.data()/observations_in.sigmas()) > filter_params[0]
    elif filter_mode == 'partiality':
      i_sel = partiality_in > filter_params[0]

    observations_out = observations_in.customized_copy(\
      indices=observations_in.indices().select(i_sel),
      data=observations_in.data().select(i_sel),
      sigmas=observations_in.sigmas().select(i_sel))
    alpha_angle_out = alpha_angle_in.select(i_sel)
    spot_pred_x_mm_out = spot_pred_x_mm_in.select(i_sel)
    spot_pred_y_mm_out = spot_pred_y_mm_in.select(i_sel)
    I_ref_out = I_ref_in.select(i_sel)

    return observations_out, alpha_angle_out, spot_pred_x_mm_out, spot_pred_y_mm_out, I_ref_out



  def optimize_scalefactors(self, I_r_flex, observations_original,
              wavelength, crystal_init_orientation,
              alpha_angle, spot_pred_x_mm, spot_pred_y_mm, iparams,
              pres_in, observations_non_polar, detector_distance_mm, const_params, G=1, B=0):

    self.gamma_e = iparams.gamma_e
    pr_d_min = iparams.postref.scale.d_min
    pr_d_max = iparams.postref.scale.d_max
    pr_sigma_min = iparams.postref.scale.sigma_min
    pr_partiality_min = iparams.postref.scale.partiality_min

    #filter by resolution
    observations_original_sel, alpha_angle_sel, spot_pred_x_mm_sel, \
      spot_pred_y_mm_sel, I_ref_sel = self.get_filtered_data(\
        'resolution', [pr_d_min, pr_d_max], observations_original, alpha_angle,
        spot_pred_x_mm, spot_pred_y_mm, I_r_flex)

    #filter by sigma
    observations_original_sel, alpha_angle_sel, spot_pred_x_mm_sel, \
          spot_pred_y_mm_sel, I_ref_sel = self.get_filtered_data(\
            'sigma', [pr_sigma_min], observations_original_sel, alpha_angle_sel,
            spot_pred_x_mm_sel, spot_pred_y_mm_sel, I_ref_sel)

    I_r_true = I_ref_sel.as_numpy_array()
    I_o_true = observations_original_sel.data().as_numpy_array()

    if pres_in is not None:
      b0 = pres_in.B
      r0 = pres_in.r0
    else:
      b0 = 0
      r0 = calc_spot_radius(sqr(crystal_init_orientation.reciprocal_matrix()),
                                                     observations_original_sel.indices(), wavelength)

    refine_mode = 'scale_factor'
    xinp = np.array([G,B])

    xopt, cov_x, infodict, mesg, ier = optimize.leastsq(func, xinp,
                                                          args=(I_r_true, observations_original_sel,
                                                                wavelength, alpha_angle_sel,
                                                                crystal_init_orientation,
                                                                spot_pred_x_mm_sel, spot_pred_y_mm_sel,
                                                                detector_distance_mm,
                                                                refine_mode, const_params,
                                                                b0, None, iparams),
                                                          full_output=True, maxfev=100)
    G_fin, B_fin = xopt
    rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma = const_params
    two_theta = observations_original.two_theta(wavelength=wavelength)
    sin_theta_over_lambda_sq = two_theta.sin_theta_over_lambda_sq().data()
    uc = unit_cell((a,b,c,alpha,beta,gamma))

    partiality_init, delta_xy_init, rs_init, dummy = calc_partiality_anisotropy_set(uc, rotx, roty,
                                                                             observations_original.indices(),
                                                                             ry, rz, r0, re,
                                                                             two_theta.data(),
                                                                             alpha_angle,
                                                                             wavelength,
                                                                             crystal_init_orientation, spot_pred_x_mm,
                                                                             spot_pred_y_mm,
                                                                             detector_distance_mm,
                                                                             iparams.partiality_model,
                                                                             iparams.flag_beam_divergence)

    I_o_init = calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                              G, B, partiality_init, rs_init, iparams.flag_volume_correction)
    I_o_fin = calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                              G_fin, B_fin, partiality_init, rs_init, iparams.flag_volume_correction)

    SE_of_the_estimate = standard_error_of_the_estimate(I_r_flex, I_o_fin, 2)
    R_sq = coefficient_of_determination(I_r_flex,I_o_fin)*100

    CC_init = np.corrcoef(I_r_flex, I_o_init)[0,1]
    CC_final = np.corrcoef(I_r_flex, I_o_fin)[0,1]
    err_init = (I_r_flex - I_o_init)/observations_original.sigmas()
    R_init = np.sqrt(np.sum(err_init**2))
    err_final = (I_r_flex - I_o_fin)/observations_original.sigmas()
    R_final = np.sqrt(np.sum(err_final**2))
    R_xy_init = 0
    R_xy_final = 0
    CC_iso_init = 0
    CC_iso_final = 0

    return xopt, (SE_of_the_estimate, R_sq, CC_init, CC_final, R_init, R_final, R_xy_init, R_xy_final, CC_iso_init, CC_iso_final)

  def prepare_data_microcycle(self, refine_mode, iparams,
                        observations_original, alpha_angle,
                        spot_pred_x_mm, spot_pred_y_mm,
                        I_r_flex, init_params, crystal_init_orientation,
                        wavelength, detector_distance_mm):
    #prepare data
    if refine_mode == 'scale_factor':
      pr_d_min = iparams.postref.scale.d_min
      pr_d_max = iparams.postref.scale.d_max
      pr_sigma_min = iparams.postref.scale.sigma_min
      pr_partiality_min = iparams.postref.scale.partiality_min
      pr_uc_tol = iparams.postref.unit_cell.uc_tolerance
    elif refine_mode == 'crystal_orientation':
      pr_d_min = iparams.postref.crystal_orientation.d_min
      pr_d_max = iparams.postref.crystal_orientation.d_max
      pr_sigma_min = iparams.postref.crystal_orientation.sigma_min
      pr_partiality_min = iparams.postref.crystal_orientation.partiality_min
      pr_uc_tol = iparams.postref.unit_cell.uc_tolerance
    elif refine_mode == 'reflecting_range':
      pr_d_min = iparams.postref.reflecting_range.d_min
      pr_d_max = iparams.postref.reflecting_range.d_max
      pr_sigma_min = iparams.postref.reflecting_range.sigma_min
      pr_partiality_min = iparams.postref.reflecting_range.partiality_min
      pr_uc_tol = iparams.postref.unit_cell.uc_tolerance
    elif refine_mode == 'unit_cell':
      pr_d_min = iparams.postref.unit_cell.d_min
      pr_d_max = iparams.postref.unit_cell.d_max
      pr_sigma_min = iparams.postref.unit_cell.sigma_min
      pr_partiality_min = iparams.postref.unit_cell.partiality_min
      pr_uc_tol = iparams.postref.unit_cell.uc_tolerance
    elif refine_mode == 'allparams':
      pr_d_min = iparams.postref.allparams.d_min
      pr_d_max = iparams.postref.allparams.d_max
      pr_sigma_min = iparams.postref.allparams.sigma_min
      pr_partiality_min = iparams.postref.allparams.partiality_min
      pr_uc_tol = iparams.postref.unit_cell.uc_tolerance


    #filter by resolution
    observations_original_sel, alpha_angle_sel, spot_pred_x_mm_sel, \
          spot_pred_y_mm_sel, I_ref_sel = self.get_filtered_data(\
            'resolution', [pr_d_min, pr_d_max], observations_original, alpha_angle,
            spot_pred_x_mm, spot_pred_y_mm, I_r_flex)

    #filter by sigma
    observations_original_sel, alpha_angle_sel, spot_pred_x_mm_sel, \
              spot_pred_y_mm_sel, I_ref_sel = self.get_filtered_data(\
                'sigma', [pr_sigma_min], observations_original_sel, alpha_angle_sel,
                spot_pred_x_mm_sel, spot_pred_y_mm_sel, I_ref_sel)

    #extract refined parameters
    G, B, rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma = init_params

    #filter by partiality
    two_theta = observations_original_sel.two_theta(wavelength=wavelength).data()
    uc = unit_cell((a,b,c,alpha,beta,gamma))
    partiality_init, delta_xy_init, rs_init, dummy = calc_partiality_anisotropy_set(uc, rotx, roty,
                                                                    observations_original_sel.indices(),
                                                                    ry, rz, r0, re, two_theta,
                                                                    alpha_angle_sel, wavelength,
                                                                    crystal_init_orientation,
                                                                    spot_pred_x_mm_sel,
                                                                    spot_pred_y_mm_sel,
                                                                    detector_distance_mm,
                                                                    iparams.partiality_model,
                                                                    iparams.flag_beam_divergence)

    observations_original_sel, alpha_angle_sel, spot_pred_x_mm_sel, \
              spot_pred_y_mm_sel, I_ref_sel = self.get_filtered_data(\
                'partiality', [pr_partiality_min], observations_original_sel,
                alpha_angle_sel, spot_pred_x_mm_sel, spot_pred_y_mm_sel, I_ref_sel,
                partiality_in=partiality_init)

    return observations_original_sel, alpha_angle_sel, spot_pred_x_mm_sel, \
              spot_pred_y_mm_sel, I_ref_sel

  def optimize(self, I_r_flex, observations_original,
              wavelength, crystal_init_orientation,
              alpha_angle, spot_pred_x_mm, spot_pred_y_mm, iparams,
              pres_in, observations_non_polar, detector_distance_mm):

    self.gamma_e = iparams.gamma_e

    if iparams.postref.allparams.flag_on:
      refine_steps = ['allparams']
    else:
      refine_steps = ['scale_factor']
      if iparams.postref.crystal_orientation.flag_on:
        refine_steps.append('crystal_orientation')
      if iparams.postref.reflecting_range.flag_on:
        refine_steps.append('reflecting_range')
      if iparams.postref.unit_cell.flag_on:
        refine_steps.append('unit_cell')


    #get miller array iso, if given.
    miller_array_iso = None
    if iparams.hklisoin is not None:
      from iotbx import reflection_file_reader
      reflection_file_iso = reflection_file_reader.any_reflection_file(iparams.hklisoin)
      miller_arrays_iso=reflection_file_iso.as_miller_arrays()
      is_found_iso_as_intensity_array = False
      is_found_iso_as_amplitude_array = False
      for miller_array in miller_arrays_iso:
        if miller_array.is_xray_intensity_array():
          miller_array_iso = miller_array.deep_copy()
          is_found_iso_as_intensity_array = True
          break
        elif miller_array.is_xray_amplitude_array():
          is_found_iso_as_amplitude_array = True
          miller_array_converted_to_intensity = miller_array.as_intensity_array()
      if is_found_iso_as_intensity_array == False:
        if is_found_iso_as_amplitude_array:
          miller_array_iso = miller_array_converted_to_intensity.deep_copy()


    #prepare data
    pr_d_min = iparams.postref.allparams.d_min
    pr_d_max = iparams.postref.allparams.d_max
    pr_sigma_min = iparams.postref.allparams.sigma_min
    pr_partiality_min = iparams.postref.allparams.partiality_min
    pr_uc_tol = iparams.postref.allparams.uc_tolerance
    cs = observations_original.crystal_symmetry().space_group().crystal_system()

    #filter by resolution
    observations_original_sel, alpha_angle_sel, spot_pred_x_mm_sel, \
          spot_pred_y_mm_sel, I_ref_sel = self.get_filtered_data(\
            'resolution', [pr_d_min, pr_d_max], observations_original, alpha_angle,
            spot_pred_x_mm, spot_pred_y_mm, I_r_flex)

    #filter by sigma
    observations_original_sel, alpha_angle_sel, spot_pred_x_mm_sel, \
              spot_pred_y_mm_sel, I_ref_sel = self.get_filtered_data(\
                'sigma', [pr_sigma_min], observations_original_sel, alpha_angle_sel,
                spot_pred_x_mm_sel, spot_pred_y_mm_sel, I_ref_sel)

    #initialize values only in the first sub cycle and the first refine step.
    spot_radius = calc_spot_radius(sqr(crystal_init_orientation.reciprocal_matrix()),
                                                     observations_original_sel.indices(), wavelength)
    if pres_in is None:
      ry, rz, r0, re, rotx, roty = 0, 0, spot_radius, self.gamma_e, 0.0, 0.0

      #apply constrain on the unit cell using crystal system
      uc_scale_inp = prep_input(observations_original.unit_cell().parameters(), cs)
      uc_scale_constrained = prep_output(uc_scale_inp, cs)
      a,b,c,alpha,beta,gamma = uc_scale_constrained
      const_params_scale = (rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma)
      xopt_scalefactors, stats = self.optimize_scalefactors(I_r_flex,
                                                                  observations_original,
                                                                  wavelength, crystal_init_orientation,
                                                                  alpha_angle,
                                                                  spot_pred_x_mm,
                                                                  spot_pred_y_mm,
                                                                  iparams,
                                                                  pres_in,
                                                                  observations_non_polar,
                                                                  detector_distance_mm,
                                                                  const_params_scale)
      G, B = xopt_scalefactors

    else:
      G, B, ry, rz, r0, re, rotx, roty = pres_in.G, pres_in.B, pres_in.ry, pres_in.rz, pres_in.r0, pres_in.re, 0.0 , 0.0
      a,b,c,alpha,beta,gamma = pres_in.unit_cell.parameters()
      crystal_init_orientation = pres_in.crystal_orientation


    #filter by partiality
    two_theta = observations_original_sel.two_theta(wavelength=wavelength).data()
    uc = unit_cell((a,b,c,alpha,beta,gamma))
    partiality_init, delta_xy_init, rs_init, dummy = calc_partiality_anisotropy_set(uc, rotx, roty,
                                                                    observations_original_sel.indices(),
                                                                    ry, rz, r0, re, two_theta,
                                                                    alpha_angle_sel, wavelength,
                                                                    crystal_init_orientation,
                                                                    spot_pred_x_mm_sel,
                                                                    spot_pred_y_mm_sel,
                                                                    detector_distance_mm,
                                                                    iparams.partiality_model,
                                                                    iparams.flag_beam_divergence)

    observations_original_sel, alpha_angle_sel, spot_pred_x_mm_sel, \
              spot_pred_y_mm_sel, I_ref_sel = self.get_filtered_data(\
                'partiality', [pr_partiality_min], observations_original_sel,
                alpha_angle_sel, spot_pred_x_mm_sel, spot_pred_y_mm_sel, I_ref_sel,
                partiality_in=partiality_init)

    I_r_true = I_ref_sel.as_numpy_array()
    I_o_true = observations_original_sel.data().as_numpy_array()


    #calculate initial residual and residual_xy error
    const_params_uc = (G, B, rotx, roty, ry, rz, r0, re)
    xinp_uc = prep_input((a,b,c,alpha,beta,gamma), cs)
    uc_params_err = func(xinp_uc,
                I_r_true,
                observations_original_sel, wavelength, alpha_angle_sel,
                crystal_init_orientation,
                spot_pred_x_mm_sel, spot_pred_y_mm_sel,
                detector_distance_mm,
                'unit_cell', const_params_uc,
                B, miller_array_iso, iparams)
    init_residual_xy_err = np.sum(uc_params_err**2)

    const_params_all = None
    xinp_all = np.array([G, B, rotx, roty, ry, rz, r0, re])
    xinp_all_uc = prep_input((a,b,c,alpha,beta,gamma), cs)
    xinp_all = np.append(xinp_all, xinp_all_uc)
    all_params_err = func(xinp_all,
                I_r_true,
                observations_original_sel, wavelength, alpha_angle_sel,
                crystal_init_orientation,
                spot_pred_x_mm_sel, spot_pred_y_mm_sel,
                detector_distance_mm,
                'allparams', const_params_all,
                B, miller_array_iso, iparams)
    init_residual_err = np.sum(all_params_err**2)

    t_pr_list = [init_residual_err]
    t_xy_list = [init_residual_xy_err]
    refined_params_hist = [(G, B, rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma)]
    txt_out = ''

    for i_sub_cycle in range(iparams.n_postref_sub_cycle):

      for j_refine_step in range(len(refine_steps)):
        refine_mode = refine_steps[j_refine_step]

        #prepare data
        init_params = (G, B, rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma)
        observations_original_sel, alpha_angle_sel, spot_pred_x_mm_sel, \
              spot_pred_y_mm_sel, I_ref_sel = self.prepare_data_microcycle(refine_mode, iparams,
                        observations_original, alpha_angle,
                        spot_pred_x_mm, spot_pred_y_mm,
                        I_r_flex, init_params, crystal_init_orientation,
                        wavelength, detector_distance_mm)

        I_r_true = I_ref_sel.as_numpy_array()
        I_o_true = observations_original_sel.data().as_numpy_array()

        if refine_mode == 'scale_factor':
          xinp = np.array([G, B])
          const_params = (rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma)
        elif refine_mode == 'crystal_orientation':
          xinp = np.array([rotx, roty])
          const_params = (G, B, ry, rz, r0, re, a, b, c, alpha, beta, gamma)
        elif refine_mode == 'reflecting_range':
          xinp = np.array([ry, rz, r0, re])
          const_params = (G, B, rotx, roty, a, b, c, alpha, beta, gamma)
        elif refine_mode == 'unit_cell':
          xinp = prep_input((a,b,c,alpha,beta,gamma), cs)
          const_params = (G, B, rotx, roty, ry, rz, r0, re)
        elif refine_mode == 'allparams':
          xinp = np.array([G, B, rotx, roty, ry, rz, r0, re])
          xinp_uc = prep_input((a,b,c,alpha,beta,gamma), cs)
          xinp = np.append(xinp, xinp_uc)
          const_params = None

        xopt, cov_x, infodict, mesg, ier = optimize.leastsq(func, xinp,
                                                            args=(I_r_true,
                                                                  observations_original_sel,
                                                                  wavelength, alpha_angle_sel,
                                                                  crystal_init_orientation,
                                                                  spot_pred_x_mm_sel,
                                                                  spot_pred_y_mm_sel,
                                                                  detector_distance_mm,
                                                                  refine_mode, const_params,
                                                                  B, miller_array_iso, iparams),
                                                            full_output=True, maxfev=100)


        if refine_mode == 'scale_factor' or refine_mode == 'crystal_orientation' or \
            refine_mode == 'reflecting_range' or refine_mode == 'allparams':

          current_residual_err = np.sum(infodict['fvec']**2)

          #calculate residual_xy_error (for refine_mode = SF, CO, RR, and all params)
          xinp_uc = prep_input((a,b,c,alpha,beta,gamma), cs)
          if refine_mode == 'scale_factor':
            G, B = xopt
          elif refine_mode == 'crystal_orientation':
            rotx, roty = xopt
          elif refine_mode == 'reflecting_range':
            ry, rz, r0, re = xopt
          elif refine_mode == 'allparams':
            G, B, rotx, roty, ry, rz, r0, re = xopt[:8]
            xinp_uc = xopt[8:]
            a, b, c, alpha, beta, gamma = prep_output(xinp_uc, cs)

          const_params_uc = (G, B, rotx, roty, ry, rz, r0, re)
          uc_params_err = func(xinp_uc,
                  I_r_true,
                  observations_original_sel, wavelength, alpha_angle_sel,
                  crystal_init_orientation,
                  spot_pred_x_mm_sel, spot_pred_y_mm_sel,
                  detector_distance_mm,
                  'unit_cell', const_params_uc,
                  B, miller_array_iso, iparams)
          current_residual_xy_err = np.sum(uc_params_err**2)
        elif refine_mode == 'unit_cell':
          current_residual_xy_err = np.sum(infodict['fvec']**2)
          xopt_uc = prep_output(xopt, cs)
          a, b, c, alpha, beta, gamma = xopt_uc

          #check the unit-cell with the reference intensity
          xinp = np.array([G, B, rotx, roty, ry, rz, r0, re])
          xinp_uc = prep_input((a, b, c, alpha, beta, gamma), cs)
          xinp = np.append(xinp, xinp_uc)
          const_params = None
          all_params_err = func(xinp,
                I_r_true,
                observations_original_sel, wavelength, alpha_angle_sel,
                crystal_init_orientation,
                spot_pred_x_mm_sel, spot_pred_y_mm_sel,
                detector_distance_mm,
                'allparams', const_params,
                B, miller_array_iso, iparams)
          current_residual_err = np.sum(all_params_err**2)


        flag_success = False
        if refine_mode == 'allparams':
          #if allparams refinement, only check the post-refine target function
          if current_residual_err < (t_pr_list[len(t_pr_list)-1] + \
              (t_pr_list[len(t_pr_list)-1]*iparams.postref.residual_threshold/100)):
            t_pr_list.append(current_residual_err)
            t_xy_list.append(current_residual_xy_err)
            refined_params_hist.append((G, B, rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma))
            flag_success = True
        else:
          #print '---'
          #print '%.4g %.4g %.4g %.4g %.4g %.4g'%(current_residual_err, t_pr_list[len(t_pr_list)-1], t_pr_list[len(t_pr_list)-1]*iparams.postref.residual_threshold/100, current_residual_xy_err, t_xy_list[len(t_xy_list)-1], t_xy_list[len(t_xy_list)-1]*iparams.postref.residual_threshold_xy/100)
          #print '---'
          if current_residual_err < (t_pr_list[len(t_pr_list)-1] + \
                (t_pr_list[len(t_pr_list)-1]*iparams.postref.residual_threshold/100)):
            if current_residual_xy_err < (t_xy_list[len(t_xy_list)-1] + \
                (t_xy_list[len(t_xy_list)-1]*iparams.postref.residual_threshold_xy/100)):
              t_pr_list.append(current_residual_err)
              t_xy_list.append(current_residual_xy_err)
              refined_params_hist.append((G, B, rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma))
              flag_success = True
              #print 'accecpt', G, B, rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma

        if flag_success is False:
          G,B,rotx,roty,ry,rz,r0,re,a,b,c,alpha,beta,gamma = refined_params_hist[len(refined_params_hist)-1]
          #print 'reverted to', G, B, rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma

        tmp_txt_out = refine_mode + ' %3.0f %6.4f %6.4f %6.4f %6.4f %10.8f %10.8f %10.8f %10.8f %6.3f %6.3f %.4g %6.3f\n'%(i_sub_cycle,G,B,rotx*180/math.pi,roty*180/math.pi,ry,rz,r0,re,a,c,t_pr_list[len(t_pr_list)-1],t_xy_list[len(t_pr_list)-1])
        txt_out += tmp_txt_out



    #apply the refined parameters on the full (original) reflection set
    two_theta = observations_original.two_theta(wavelength=wavelength).data()
    sin_theta_over_lambda_sq = observations_original.two_theta(wavelength=wavelength).sin_theta_over_lambda_sq().data()

    if pres_in is None:
      partiality_init, delta_xy_init, rs_init, rh_init = calc_partiality_anisotropy_set(observations_original.unit_cell(),
                                                                      0.0, 0.0,
                                                                      observations_original.indices(),
                                                                      0, 0, spot_radius, self.gamma_e,
                                                                      two_theta, alpha_angle, wavelength,
                                                                      crystal_init_orientation,
                                                                      spot_pred_x_mm, spot_pred_y_mm,
                                                                      detector_distance_mm,
                                                                      iparams.partiality_model,
                                                                      iparams.flag_beam_divergence)

      I_o_init = calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                                1, 0, partiality_init, rs_init, iparams.flag_volume_correction)

    else:
      partiality_init, delta_xy_init, rs_init, rh_init = calc_partiality_anisotropy_set(pres_in.unit_cell,
                                                                      0.0, 0.0,
                                                                      observations_original.indices(),
                                                                      pres_in.ry, pres_in.rz,
                                                                      pres_in.r0, pres_in.re,
                                                                      two_theta, alpha_angle, wavelength,
                                                                      crystal_init_orientation,
                                                                      spot_pred_x_mm, spot_pred_y_mm,
                                                                      detector_distance_mm,
                                                                      iparams.partiality_model,
                                                                      iparams.flag_beam_divergence)
      I_o_init = calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                              pres_in.G, pres_in.B, partiality_init, rs_init, iparams.flag_volume_correction)



    partiality_fin, delta_xy_fin, rs_fin, rh_fin = calc_partiality_anisotropy_set(unit_cell((a,b,c,alpha,beta,gamma)),
                                                                              rotx, roty,
                                                                              observations_original.indices(),
                                                                              ry, rz, r0, re,
                                                                              two_theta, alpha_angle, wavelength,
                                                                              crystal_init_orientation,
                                                                              spot_pred_x_mm, spot_pred_y_mm,
                                                                              detector_distance_mm,
                                                                              iparams.partiality_model,
                                                                              iparams.flag_beam_divergence)

    I_o_fin = calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                              G, B, partiality_fin, rs_fin, iparams.flag_volume_correction)
    SE_of_the_estimate = standard_error_of_the_estimate(I_r_flex,I_o_fin, 13)
    R_sq = coefficient_of_determination(I_r_flex,I_o_fin)*100

    CC_init = np.corrcoef(I_r_flex, I_o_init)[0,1]
    CC_final = np.corrcoef(I_r_flex, I_o_fin)[0,1]
    err_init = (I_r_flex - I_o_init)/observations_original.sigmas()
    R_init = np.sqrt(np.sum(err_init**2))
    err_final = (I_r_flex - I_o_fin)/observations_original.sigmas()
    R_final = np.sqrt(np.sum(err_final**2))
    R_xy_init = math.sqrt(np.sum(delta_xy_init**2))
    R_xy_final = math.sqrt(np.sum(delta_xy_fin**2))

    if R_init < R_final:
      CC_final = CC_init
      R_final = R_init
      R_xy_final = R_xy_init
      if pres_in is None:
        G = 1.0
        B = 0.0
        r0 = spot_radius
        ry = 0.0
        rz = 0.0
        re = self.gamma_e
        rotx = 0.0
        roty = 0.0
        a,b,c,alpha,beta,gamma = observations_original.unit_cell().parameters()
      else:
        G = pres_in.G
        B = pres_in.B
        r0 = pres_in.r0
        ry = pres_in.ry
        rz = pres_in.rz
        re = pres_in.re
        rotx = 0.0
        roty = 0.0
        a,b,c,alpha,beta,gamma = pres_in.unit_cell.parameters()
        crystal_init_orientation = pres_in.crystal_orientation

    #calculate CCiso if hklisoin is given
    CC_iso_init = 0
    CC_iso_final = 0
    if iparams.hklisoin is not None:
      if miller_array_iso is not None:
        from cctbx import miller
        matches = miller.match_multi_indices(
                          miller_indices_unique=miller_array_iso.indices(),
                          miller_indices=observations_non_polar.indices())
        I_iso_match = flex.double([miller_array_iso.data()[pair[0]] for pair in matches.pairs()])
        I_o_init_match = flex.double([I_o_init[pair[1]] for pair in matches.pairs()])
        I_o_fin_match = flex.double([I_o_fin[pair[1]] for pair in matches.pairs()])

        CC_iso_init = np.corrcoef(I_iso_match, I_o_init_match)[0,1]
        CC_iso_final = np.corrcoef(I_iso_match, I_o_fin_match)[0,1]

    if iparams.flag_plot_expert:
      print 'Predictor'
      print 'Resolution: %6.2f Angstroms'%(observations_original_sel.d_min())
      print 'G %.4g'%(G)
      print 'B-factor %.4g'%(B)
      print 'rotx %.4g (%.4g degrees)'%(rotx, rotx*180/math.pi)
      print 'roty %.4g (%.4g degrees)'%(roty, roty*180/math.pi)
      print 'ry %.4g'%(ry)
      print 'rz %.4g'%(rz)
      print 'r0 %.4g'%(r0)
      print 're %.4g'%(re)
      print 'uc', a,b,c,alpha,beta,gamma
      print 'S = %.4g'%SE_of_the_estimate
      print 'Target R = %.4g'%(R_final)
      print 'Target R (x,y) = %.4g mm.'%(R_xy_final)
      print 'rh_init mean = %8.6f max = %8.6f min = %8.6f'%(np.mean(flex.abs(rh_init)), np.max(flex.abs(rh_init)), np.min(flex.abs(rh_init)))
      print 'rh_final mean = %8.6f max = %8.6f min = %8.6f'%(np.mean(flex.abs(rh_fin)), np.max(flex.abs(rh_fin)), np.min(flex.abs(rh_fin)))
      print 'CC = %.4g'%(CC_final)

      plt.subplot(221)
      plt.scatter(I_r_flex, I_o_init,s=10, marker='x', c='r')
      plt.title('CCinit=%5.2f Rinit=%5.2f'%(CC_init, R_init))
      plt.xlabel('I_ref')
      plt.ylabel('I_obs')
      plt.subplot(222)
      plt.scatter(I_r_flex, I_o_fin,s=10, marker='o', c='b')
      plt.title('CCfinal=%5.2f Rfinal=%5.2f'%(CC_final, R_final))
      plt.xlabel('I_ref')
      plt.ylabel('I_obs')

      n_bins = 16
      binner = observations_original.setup_binner(n_bins=n_bins)
      binner_indices = binner.bin_indices()
      avg_delta_xy_init = flex.double()
      avg_delta_xy_fin = flex.double()
      avg_partiality_init = flex.double()
      avg_partiality_fin = flex.double()
      avg_rh_init = flex.double()
      avg_rh_fin = flex.double()
      one_dsqr_bin = flex.double()
      for i in range(1,n_bins+1):
        i_binner = (binner_indices == i)
        if len(observations_original.data().select(i_binner)) > 0:
          avg_delta_xy_init.append(np.mean(delta_xy_init.select(i_binner)))
          avg_delta_xy_fin.append(np.mean(delta_xy_fin.select(i_binner)))

          avg_partiality_init.append(np.mean(partiality_init.select(i_binner)))
          avg_partiality_fin.append(np.mean(partiality_fin.select(i_binner)))

          avg_rh_init.append(np.mean(flex.abs(rh_init.select(i_binner))))
          avg_rh_fin.append(np.mean(flex.abs(rh_fin.select(i_binner))))

          one_dsqr_bin.append(1/binner.bin_d_range(i)[1]**2)

      plt.subplot(223)
      plt.plot(one_dsqr_bin, avg_delta_xy_init, linestyle='-', linewidth=2.0, c='r', label='Initial <delta_xy>=%4.2f'%np.mean(delta_xy_init))
      plt.plot(one_dsqr_bin, avg_delta_xy_fin, linestyle='-', linewidth=2.0, c='b', label='Final <delta_xy>=%4.2f'%np.mean(delta_xy_fin))
      legend = plt.legend(loc='upper left', shadow=False)
      for label in legend.get_texts():
        label.set_fontsize('medium')
      plt.title('Distance (delta_xy) from spot centroid')
      plt.xlabel('1/(d^2)')
      plt.ylabel('(mm)')

      plt.subplot(224)
      one_dsqr = 1/(observations_original.d_spacings().data()**2)
      perm = flex.sort_permutation(one_dsqr, reverse=True)
      one_dsqr_sort = one_dsqr.select(perm)
      rs_fin_sort = rs_fin.select(perm)
      rs_init_sort = rs_init.select(perm)
      plt.plot(one_dsqr_sort, rs_init_sort, linestyle='-', linewidth=2.0, c='r',
                           label='Initial ry=%.4g ry=%.4g r0=%.4g re=%.4g'%(0,0,
                                                                            spot_radius,self.gamma_e))
      plt.plot(one_dsqr_sort, rs_fin_sort, linestyle='-', linewidth=2.0, c='b',
               label='Final ry=%.4g ry=%.4g r0=%.4g re=%.4g'%(ry,rz,r0,re))
      legend = plt.legend(loc='lower right', shadow=False)
      for label in legend.get_texts():
        label.set_fontsize('medium')
      plt.title('Reciprocal lattice radius (r_s)')
      plt.xlabel('1/(d^2)')
      plt.ylabel('1/Angstroms')
      plt.show()

      #plot partiality (histogram and function of resolutions).
      plt.subplot(221)
      x = partiality_init.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('Partiality before\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(222)
      x = partiality_fin.as_numpy_array()
      mu = np.mean(x)
      med = np.median(x)
      sigma = np.std(x)
      num_bins = 10
      n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
      y = mlab.normpdf(bins, mu, sigma)
      plt.plot(bins, y, 'r--')
      plt.ylabel('Frequencies')
      plt.title('Partiality after\nmean %5.3f median %5.3f sigma %5.3f' %(mu, med, sigma))

      plt.subplot(223)
      plt.plot(one_dsqr_bin, avg_partiality_init, linestyle='-', linewidth=2.0, c='r', label='Initial ')
      plt.plot(one_dsqr_bin, avg_partiality_fin, linestyle='-', linewidth=2.0, c='b', label='Final ')
      legend = plt.legend(loc='lower left', shadow=False)
      for label in legend.get_texts():
        label.set_fontsize('medium')
      plt.title('Partiality')
      plt.xlabel('1/(d^2)')
      plt.ylabel('(mm)')

      plt.subplot(224)
      plt.plot(one_dsqr_bin, avg_rh_init, linestyle='-', linewidth=2.0, c='r', label='Initial ')
      plt.plot(one_dsqr_bin, avg_rh_fin, linestyle='-', linewidth=2.0, c='b', label='Final ')
      legend = plt.legend(loc='lower left', shadow=False)
      for label in legend.get_texts():
        label.set_fontsize('medium')
      plt.title('Ewald-sphere offset (r_h)')
      plt.xlabel('1/(d^2)')
      plt.ylabel('1/Angstrom')
      plt.show()

    xopt = (G, B, rotx, roty, ry, rz, r0, re,a,b,c,alpha,beta,gamma)

    return xopt, (SE_of_the_estimate, R_sq, CC_init, CC_final, R_init, R_final, R_xy_init, R_xy_final, CC_iso_init, CC_iso_final), len(I_ref_sel)
