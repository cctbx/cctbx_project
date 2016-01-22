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
import math
from cctbx.array_family import flex
from scitbx.matrix import sqr, col
from cctbx.uctbx import unit_cell
from cctbx.crystal_orientation import crystal_orientation, basis_type

def calc_full_refl(I_o_p_set, sin_theta_over_lambda_sq_set,
                   G, B, p_set, rs_set, flag_volume_correction):
  I_o_full_set = I_o_p_set/(G * flex.exp(-2*B*sin_theta_over_lambda_sq_set) * p_set)
  return I_o_full_set

def calc_spot_radius(a_star_matrix, miller_indices, wavelength):
  #calculate spot_radius based on rms delta_S for all spots
  S0 = -1*col((0,0,1./wavelength))
  sd_array = a_star_matrix.elems * miller_indices.as_vec3_double() + S0.elems
  rh_set = sd_array.norms() - (1/wavelength)
  return rh_set.standard_deviation_of_the_sample()

def coefficient_of_determination(y, y_model):
  mean_y = flex.mean(y)
  r_sqr = flex.sum((y_model - mean_y)**2)/flex.sum((y - mean_y)**2)
  return r_sqr

def standard_error_of_the_estimate(y, y_model, n_params):
  s = math.sqrt(flex.sum((y - y_model)**2)/(len(y) - n_params))
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
  #caculate rs
  rs_set = r0 + (re * flex.tan(bragg_angle_set))
  if flag_beam_divergence:
    rs_set += ((ry * flex.cos(alpha_angle_set))**2 + (rz * flex.sin(alpha_angle_set))**2)**(1/2)
  #calculate rh
  sd_array = A_star.elems * miller_indices.as_vec3_double() + S0.elems
  rh_set = sd_array.norms() - (1/wavelength)
  #calculate partiality
  partiality_set = ((rs_set**2)/((2*(rh_set**2))+(rs_set**2)))
  #calculate delta_xy
  d_ratio = -detector_distance_mm/sd_array.parts()[2]
  calc_xy_array = flex.vec3_double(sd_array.parts()[0]*d_ratio, \
      sd_array.parts()[1]*d_ratio, flex.double([0]*len(d_ratio)))
  pred_xy_array = flex.vec3_double(spot_pred_x_mm_set, spot_pred_y_mm_set, flex.double([0]*len(d_ratio)))
  delta_xy_set = (pred_xy_array - calc_xy_array).norms()
  return partiality_set, delta_xy_set, rs_set, rh_set

def prep_input(params, cs):
  #From crystal system cs, determine refined parameters
  a, b, c, alpha, beta, gamma = params
  if cs == 'Triclinic':
    x0 = params
  elif cs == 'Monoclinic':
    x0 = flex.double([a,b,c,beta])
  elif cs == 'Orthorhombic':
    x0 = flex.double([a,b,c])
  elif cs == 'Tetragonal':
    x0 = flex.double([a,c])
  elif cs == 'Trigonal' or cs == 'Hexagonal':
    x0 = flex.double([a,c])
  elif cs == 'Cubic':
    x0 = flex.double([a])
  return x0

def prep_output(params, cs):
  if cs == 'Triclinic':
    xopt = params
  elif cs == 'Monoclinic':
    xopt = flex.double([params[0],params[1],params[2],90,params[3],90])
  elif cs == 'Orthorhombic':
    xopt = flex.double([params[0],params[1],params[2],90,90,90])
  elif cs == 'Tetragonal':
    xopt = flex.double([params[0],params[0],params[1],90,90,90])
  elif cs == 'Trigonal' or cs == 'Hexagonal':
    xopt = flex.double([params[0],params[0],params[1],90,90,120])
  elif cs == 'Cubic':
    xopt = flex.double([params[0],params[0],params[0],90,90,90])
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
  I_o = miller_array_o.data()
  sigI_o = miller_array_o.sigmas()
  miller_indices_original = miller_array_o.indices()
  sin_theta_over_lambda_sq = miller_array_o.two_theta(wavelength=wavelength).sin_theta_over_lambda_sq().data()
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
    a, b, c, alpha, beta, gamma = prep_output(params[6:], cs)
    rotx, roty, ry, rz, r0, re = params[0:6]
    G,B = const_params
  try:
    uc = unit_cell((a,b,c,alpha,beta,gamma))
  except Exception:
    return None
  G,B,rotx,roty,ry,rz,r0,re,a,b,c,alpha,beta,gamma = flex.double([G,B,rotx,roty,ry,rz,r0,re,a,b,c,alpha,beta,gamma])
  p_calc_flex, delta_xy_flex, rs_set, dummy = calc_partiality_anisotropy_set(\
    uc, rotx, roty, miller_indices_original, ry, rz, r0, re, two_theta_flex,
    alpha_angle_set, wavelength, crystal_init_orientation, spot_pred_x_mm_set,
    spot_pred_y_mm_set, detector_distance_mm, partiality_model, flag_beam_divergence)
  if miller_array_o.d_min() < b_refine_d_min:
    I_o_full = calc_full_refl(I_o, sin_theta_over_lambda_sq,
                              G, B, p_calc_flex, rs_set, flag_volume_correction)
  else:
    I_o_full = calc_full_refl(I_o, sin_theta_over_lambda_sq,
                              G, b0, p_calc_flex, rs_set, flag_volume_correction)
  if refine_mode == 'unit_cell':
    error = delta_xy_flex
  else:
    error = ((I_r - I_o_full)/sigI_o)
  #print refine_mode, 'G=%.4g B=%.4g rotx=%.4g roty=%.4g ry=%.4g rz=%.4g r0=%.4g re=%.4g a=%.4g b=%.4g c=%.4g alp=%.4g beta=%.4g gam=%.4g fpr=%.4g fxy=%.4g dd_mm=%6.2f n_refl=%5.0f'%(G, B, rotx*180/math.pi, roty*180/math.pi, ry, rz, r0, re, a, b, c, alpha, beta, gamma, flex.sum(((I_r - I_o_full)/sigI_o)**2), flex.sum(delta_xy_flex**2), detector_distance_mm, len(I_o_full))
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
  """
  A wrapper class for least-squares refinement
  """
  def __init__(self):
    """
    Intialitze parameters
    """

  def get_filtered_data(self, filter_mode, filter_params,
                        observations_in, alpha_angle_in,
                        spot_pred_x_mm_in, spot_pred_y_mm_in,
                        I_ref_in, partiality_in=False):
    if filter_mode == 'resolution':
      i_sel = observations_in.resolution_filter_selection(d_min=filter_params[0],\
          d_max=filter_params[1])
    elif filter_mode == 'sigma':
      i_sel = (observations_in.data()/observations_in.sigmas()) > filter_params[0]
    elif filter_mode == 'partiality':
      i_sel = partiality_in > filter_params[0]
    observations_out = observations_in.select(i_sel)
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
        'resolution', [pr_d_min, pr_d_max], observations_original, alpha_angle,\
        spot_pred_x_mm, spot_pred_y_mm, I_r_flex)
    #filter by sigma
    observations_original_sel, alpha_angle_sel, spot_pred_x_mm_sel, \
        spot_pred_y_mm_sel, I_ref_sel = self.get_filtered_data(\
        'sigma', [pr_sigma_min], observations_original_sel, alpha_angle_sel,\
        spot_pred_x_mm_sel, spot_pred_y_mm_sel, I_ref_sel)
    I_r_true = I_ref_sel[:]
    I_o_true = observations_original_sel.data()[:]
    if pres_in is not None:
      b0 = pres_in.B
      r0 = pres_in.r0
    else:
      b0 = 0
      r0 = calc_spot_radius(sqr(crystal_init_orientation.reciprocal_matrix()),\
               observations_original_sel.indices(), wavelength)
    refine_mode = 'scale_factor'
    xinp = flex.double([G,B])
    xopt, cov_x, infodict, mesg, ier = optimize.leastsq(func, xinp,
                                                          args=(I_r_true, observations_original_sel,
                                                                wavelength, alpha_angle_sel,
                                                                crystal_init_orientation,
                                                                spot_pred_x_mm_sel, spot_pred_y_mm_sel,
                                                                detector_distance_mm,
                                                                refine_mode, const_params,
                                                                b0, None, iparams),
                                                          full_output=True, maxfev=100)
    G_fin, B_fin = flex.double(xopt)
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
                                                                             crystal_init_orientation,
                                                                             spot_pred_x_mm,
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
    CC_init = flex.linear_correlation(I_r_flex, I_o_init).coefficient()
    CC_final = flex.linear_correlation(I_r_flex, I_o_fin).coefficient()
    err_init = (I_r_flex - I_o_init)/observations_original.sigmas()
    R_init = math.sqrt(flex.sum(err_init**2))
    err_final = (I_r_flex - I_o_fin)/observations_original.sigmas()
    R_final = math.sqrt(flex.sum(err_final**2))
    R_xy_init = 0
    R_xy_final = 0
    CC_iso_init = 0
    CC_iso_final = 0
    return flex.double(xopt), (SE_of_the_estimate, R_sq, CC_init, CC_final, R_init, R_final, R_xy_init, R_xy_final, CC_iso_init, CC_iso_final)

  def prepare_data_microcycle(self, refine_mode, iparams,
                        observations_original, alpha_angle,
                        spot_pred_x_mm, spot_pred_y_mm,
                        I_r_flex, init_params, crystal_init_orientation,
                        wavelength, detector_distance_mm):
    #prepare data
    if refine_mode == 'crystal_orientation':
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
      refine_steps = ['crystal_orientation']
      if iparams.postref.reflecting_range.flag_on:
        refine_steps.append('reflecting_range')
      if iparams.postref.unit_cell.flag_on:
        refine_steps.append('unit_cell')
    #get miller array iso, if given.
    miller_array_iso = None

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
    I_r_true = I_ref_sel[:]
    I_o_true = observations_original_sel.data()[:]
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
    init_residual_xy_err = flex.sum(uc_params_err**2)
    const_params_all = None
    xinp_all = flex.double([rotx, roty, ry, rz, r0, re])
    xinp_all.extend(prep_input((a,b,c,alpha,beta,gamma), cs))
    const_params_all= (G,B)
    all_params_err = func(xinp_all,
                I_r_true,
                observations_original_sel, wavelength, alpha_angle_sel,
                crystal_init_orientation,
                spot_pred_x_mm_sel, spot_pred_y_mm_sel,
                detector_distance_mm,
                'allparams', const_params_all,
                B, miller_array_iso, iparams)
    init_residual_err = flex.sum(all_params_err**2)
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
        I_r_true = I_ref_sel[:]
        I_o_true = observations_original_sel.data()
        if refine_mode == 'crystal_orientation':
          xinp = flex.double([rotx, roty])
          const_params = (G, B, ry, rz, r0, re, a, b, c, alpha, beta, gamma)
        elif refine_mode == 'reflecting_range':
          xinp = flex.double([ry, rz, r0, re])
          const_params = (G, B, rotx, roty, a, b, c, alpha, beta, gamma)
        elif refine_mode == 'unit_cell':
          xinp = prep_input((a,b,c,alpha,beta,gamma), cs)
          const_params = (G, B, rotx, roty, ry, rz, r0, re)
        elif refine_mode == 'allparams':
          xinp = flex.double([rotx, roty, ry, rz, r0, re])
          xinp.extend(prep_input((a,b,c,alpha,beta,gamma), cs))
          const_params = (G,B)
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
        xopt = flex.double(xopt)
        if refine_mode == 'crystal_orientation' or \
            refine_mode == 'reflecting_range' or refine_mode == 'allparams':
          current_residual_err = flex.sum(flex.double(infodict['fvec'])**2)
          #calculate residual_xy_error (for refine_mode = SF, CO, RR, and all params)
          xinp_uc = prep_input((a,b,c,alpha,beta,gamma), cs)
          if refine_mode == 'crystal_orientation':
            rotx, roty = xopt
          elif refine_mode == 'reflecting_range':
            ry, rz, r0, re = xopt
          elif refine_mode == 'allparams':
            rotx, roty, ry, rz, r0, re = xopt[:6]
            xinp_uc = xopt[6:]
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
          current_residual_xy_err = flex.sum(flex.double(uc_params_err)**2)
        elif refine_mode == 'unit_cell':
          current_residual_xy_err = flex.sum(flex.double(infodict['fvec'])**2)
          xopt_uc = prep_output(xopt, cs)
          a, b, c, alpha, beta, gamma = xopt_uc
          #check the unit-cell with the reference intensity
          xinp = flex.double([rotx, roty, ry, rz, r0, re])
          xinp.extend(prep_input((a, b, c, alpha, beta, gamma), cs))
          const_params = (G,B)
          all_params_err = func(xinp,
                I_r_true,
                observations_original_sel, wavelength, alpha_angle_sel,
                crystal_init_orientation,
                spot_pred_x_mm_sel, spot_pred_y_mm_sel,
                detector_distance_mm,
                'allparams', const_params,
                B, miller_array_iso, iparams)
          current_residual_err = flex.sum(flex.double(all_params_err)**2)
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
              #print 'Refine ', refine_mode, '- accecpted', G, B, rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma
        if flag_success is False:
          G,B,rotx,roty,ry,rz,r0,re,a,b,c,alpha,beta,gamma = refined_params_hist[len(refined_params_hist)-1]
          #print 'Refine ', refine_mode, '- reverted to', G, B, rotx, roty, ry, rz, r0, re, a, b, c, alpha, beta, gamma

        tmp_txt_out = refine_mode + ' %3.0f %6.4f %6.4f %6.4f %6.4f %10.8f %10.8f %10.8f %10.8f %6.3f %6.3f %.4g %6.3f\n'%(i_sub_cycle,G,B,rotx*180/math.pi,roty*180/math.pi,ry,rz,r0,re,a,c,t_pr_list[len(t_pr_list)-1],t_xy_list[len(t_pr_list)-1])
        txt_out += tmp_txt_out
    #apply the refined parameters on the full (original) reflection set
    two_theta = observations_original.two_theta(wavelength=wavelength).data()
    sin_theta_over_lambda_sq = observations_original.two_theta(wavelength=wavelength).sin_theta_over_lambda_sq().data()
    if pres_in is None:
      partiality_init, delta_xy_init, rs_init, rh_init = calc_partiality_anisotropy_set(\
          observations_original.unit_cell(),0.0, 0.0,observations_original.indices(),
          0, 0, spot_radius, self.gamma_e, two_theta, alpha_angle, wavelength,
          crystal_init_orientation,spot_pred_x_mm, spot_pred_y_mm,detector_distance_mm,
          iparams.partiality_model,iparams.flag_beam_divergence)
      I_o_init = calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                                1, 0, partiality_init, rs_init, iparams.flag_volume_correction)
    else:
      partiality_init, delta_xy_init, rs_init, rh_init = calc_partiality_anisotropy_set(\
          pres_in.unit_cell,0.0, 0.0,observations_original.indices(),
          pres_in.ry, pres_in.rz,pres_in.r0, pres_in.re,two_theta, alpha_angle, wavelength,
          crystal_init_orientation,spot_pred_x_mm, spot_pred_y_mm,detector_distance_mm,
          iparams.partiality_model,iparams.flag_beam_divergence)
      I_o_init = calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                              pres_in.G, pres_in.B, partiality_init, rs_init, iparams.flag_volume_correction)
    partiality_fin, delta_xy_fin, rs_fin, rh_fin = calc_partiality_anisotropy_set(\
        unit_cell((a,b,c,alpha,beta,gamma)),rotx, roty,observations_original.indices(),
        ry, rz, r0, re,two_theta, alpha_angle, wavelength,crystal_init_orientation,
        spot_pred_x_mm, spot_pred_y_mm,detector_distance_mm,
        iparams.partiality_model,iparams.flag_beam_divergence)
    I_o_fin = calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                              G, B, partiality_fin, rs_fin, iparams.flag_volume_correction)
    SE_of_the_estimate = standard_error_of_the_estimate(I_r_flex,I_o_fin, 13)
    R_sq = coefficient_of_determination(I_r_flex,I_o_fin)*100
    CC_init = flex.linear_correlation(I_r_flex, I_o_init).coefficient()
    CC_final = flex.linear_correlation(I_r_flex, I_o_fin).coefficient()
    err_init = (I_r_flex - I_o_init)/observations_original.sigmas()
    R_init = math.sqrt(flex.sum(err_init**2))
    err_final = (I_r_flex - I_o_fin)/observations_original.sigmas()
    R_final = math.sqrt(flex.sum(err_final**2))
    R_xy_init = math.sqrt(flex.sum(delta_xy_init**2))
    R_xy_final = math.sqrt(flex.sum(delta_xy_fin**2))
    if R_init < R_final:
      CC_final = CC_init
      R_final = R_init
      R_xy_final = R_xy_init
      if pres_in is None:
        G,B,r0,ry,rz,re,rotx,roty = (1.0,0.0,spot_radius,0.0,0.0,self.gamma_e,0.0,0.0)
        a,b,c,alpha,beta,gamma = observations_original.unit_cell().parameters()
      else:
        G,B,r0,ry,rz,re,rotx,roty = (pres_in.G,pres_in.B,pres_in.r0,pres_in.ry,pres_in.rz,pres_in.re,0.0,0.0)
        a,b,c,alpha,beta,gamma = pres_in.unit_cell.parameters()
        crystal_init_orientation = pres_in.crystal_orientation
    #calculate CCiso if hklisoin is given
    CC_iso_init,CC_iso_final = (0,0)
    if iparams.hklisoin is not None:
      if miller_array_iso is not None:
        from cctbx import miller
        matches = miller.match_multi_indices(
                          miller_indices_unique=miller_array_iso.indices(),
                          miller_indices=observations_non_polar.indices())
        I_iso_match = flex.double([miller_array_iso.data()[pair[0]] for pair in matches.pairs()])
        I_o_init_match = flex.double([I_o_init[pair[1]] for pair in matches.pairs()])
        I_o_fin_match = flex.double([I_o_fin[pair[1]] for pair in matches.pairs()])
        CC_iso_init = flex.linear_correlation(I_iso_match, I_o_init_match).coefficient()
        CC_iso_final = flex.linear_correlation(I_iso_match, I_o_fin_match).coefficient()
    xopt = (G, B, rotx, roty, ry, rz, r0, re,a,b,c,alpha,beta,gamma)
    return xopt, (SE_of_the_estimate, R_sq, CC_init, CC_final, R_init, R_final, R_xy_init, R_xy_final, CC_iso_init, CC_iso_final), len(I_ref_sel)
