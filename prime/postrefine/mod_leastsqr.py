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
from cctbx.array_family import flex
from scitbx.matrix import sqr, col
from cctbx.uctbx import unit_cell
from cctbx.crystal_orientation import crystal_orientation, basis_type
import logging

def calc_full_refl(I_o_p_set, sin_theta_over_lambda_sq_set, G, B, p_set, rs_set, flag_volume_correction):
  if flag_volume_correction:
    I_o_full_set = flex.double(((G * np.exp(-2*B*sin_theta_over_lambda_sq_set) * I_o_p_set)/p_set)*((4/3)*rs_set))
  else:
    I_o_full_set = flex.double(((G * np.exp(-2*B*sin_theta_over_lambda_sq_set) * I_o_p_set)/p_set))
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
    if delta_S > 0.1:
      logging.warning("Rh is woryingly large: {}".format(delta_S))
    logging.debug("Delta S: {}".format(delta_S))

  spot_radius = math.sqrt(flex.mean(delta_S_all*delta_S_all))

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
    partiality_model):
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
  for miller_index, bragg_angle, alpha_angle, spot_pred_x_mm, spot_pred_y_mm in \
      zip(miller_indices, bragg_angle_set, alpha_angle_set,
          spot_pred_x_mm_set, spot_pred_y_mm_set):
    rs = math.sqrt((ry * math.cos(alpha_angle))**2 + (rz * math.sin(alpha_angle))**2) + (re*math.tan(bragg_angle))
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

    #finding coordinate x,y on the detector
    d_ratio = -detector_distance_mm/S[2]
    dx_mm = S[0]*d_ratio
    dy_mm = S[1]*d_ratio
    pred_xy = col((spot_pred_x_mm, spot_pred_y_mm))
    calc_xy = col((dx_mm, dy_mm))
    diff_xy = pred_xy - calc_xy
    delta_xy_set.append(diff_xy.length())

  return partiality_set, delta_xy_set, rs_set

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
  b_refine_d_min = args[4]
  crystal_init_orientation = args[5]
  spot_pred_x_mm_set = args[6]
  spot_pred_y_mm_set = args[7]
  detector_distance_mm = args[8]
  refine_mode = args[9]
  const_params = args[10]
  partiality_model = args[11]
  flag_volume_correction = args[12]
  r0 = args[13]
  I_o = miller_array_o.data().as_numpy_array()
  sigI_o = miller_array_o.sigmas().as_numpy_array()
  miller_indices_original = miller_array_o.indices()
  sin_theta_over_lambda_sq = miller_array_o.two_theta(wavelength=wavelength).sin_theta_over_lambda_sq().data().as_numpy_array()
  two_theta_flex = miller_array_o.two_theta(wavelength=wavelength).data()

  if refine_mode == 'scale_factor':
    G, B = params
    rotx, roty, ry, rz, re, a, b, c, alpha, beta, gamma = const_params
  elif refine_mode == 'crystal_orientation':
    rotx, roty = params
    G, B, ry, rz, re, a, b, c, alpha, beta, gamma = const_params
  elif refine_mode == 'reflecting_range':
    ry, rz, re = params
    G, B, rotx, roty, a, b, c, alpha, beta, gamma = const_params
  elif refine_mode == 'unit_cell':
    cs = miller_array_o.crystal_symmetry().space_group().crystal_system()
    a, b, c, alpha, beta, gamma = prep_output(params, cs)
    G, B, rotx, roty, ry, rz, re = const_params

  uc = unit_cell((a,b,c,alpha,beta,gamma))

  p_calc_flex, delta_xy_flex, rs_set = calc_partiality_anisotropy_set(uc, rotx, roty,
                                                                                miller_indices_original,
                                                                                ry, rz, r0, re,
                                                                                two_theta_flex,
                                                                                alpha_angle_set,
                                                                                wavelength, crystal_init_orientation,
                                                                                spot_pred_x_mm_set,
                                                                                spot_pred_y_mm_set,
                                                                                detector_distance_mm,
                                                                                partiality_model)
  p_calc_set = p_calc_flex.as_numpy_array()

  if miller_array_o.d_min() < b_refine_d_min:
    I_o_full = calc_full_refl(I_o, sin_theta_over_lambda_sq,
                              G, B, p_calc_set, rs_set, flag_volume_correction)

  else:
    I_o_full = calc_full_refl(I_o, sin_theta_over_lambda_sq,
                              G, 0, p_calc_set, rs_set, flag_volume_correction)


  if refine_mode == 'unit_cell':
    error = delta_xy_flex
  else:
    error = ((I_r - I_o_full)/sigI_o)
  #print refine_mode, 'G=%.4g B=%.4g rotx=%.4g roty=%.4g ry=%.4g rz=%.4g re=%.4g a=%.4g b=%.4g c=%.4g alp=%.4g beta=%.4g gam=%.4g f=%.4g'%(G, B, rotx*180/math.pi, roty*180/math.pi, ry, rz, re, a, b, c, alpha, beta, gamma, np.sum(error**2))
  return error



class leastsqr_handler(object):
  '''
  A wrapper class for least-squares refinement
  '''

  def __init__(self):
    '''
    Intialitze parameters
    '''
    self.gamma_e = 0.003

  def optimize_scalefactors(self, I_r_flex, observations_original,
              wavelength, crystal_init_orientation,
              alpha_angle, spot_pred_x_mm, spot_pred_y_mm, iparams,
              pres_in, observations_non_polar, detector_distance_mm):

    pr_d_min = iparams.postref.scale.d_min
    pr_d_max = iparams.postref.scale.d_max
    pr_sigma_min = iparams.postref.scale.sigma_min
    pr_partiality_min = iparams.postref.scale.partiality_min
    #filter by resolution
    i_sel_res = observations_original.resolution_filter_selection(d_min=pr_d_min, d_max=pr_d_max)
    observations_original_sel = observations_original.customized_copy( \
      indices=observations_original.indices().select(i_sel_res),
      data=observations_original.data().select(i_sel_res),
      sigmas=observations_original.sigmas().select(i_sel_res))
    alpha_angle_sel = alpha_angle.select(i_sel_res)
    spot_pred_x_mm_sel = spot_pred_x_mm.select(i_sel_res)
    spot_pred_y_mm_sel = spot_pred_y_mm.select(i_sel_res)
    I_ref_sel = I_r_flex.select(i_sel_res)

    #filter by sigma
    i_sel_sigmas = (observations_original_sel.data()/observations_original_sel.sigmas()) > pr_sigma_min
    observations_original_sel = observations_original_sel.customized_copy(indices=observations_original_sel.indices().select(i_sel_sigmas),
        data=observations_original_sel.data().select(i_sel_sigmas),
        sigmas=observations_original_sel.sigmas().select(i_sel_sigmas))
    alpha_angle_sel = alpha_angle_sel.select(i_sel_sigmas)
    spot_pred_x_mm_sel = spot_pred_x_mm_sel.select(i_sel_sigmas)
    spot_pred_y_mm_sel = spot_pred_y_mm_sel.select(i_sel_sigmas)
    I_ref_sel = I_ref_sel.select(i_sel_sigmas)
    I_r_true = I_ref_sel.as_numpy_array()
    I_o_true = observations_original_sel.data().as_numpy_array()

    spot_radius = calc_spot_radius(sqr(crystal_init_orientation.reciprocal_matrix()),
                                  observations_original_sel.indices(), wavelength)
    G = np.median(I_ref_sel)/np.median(observations_original_sel.data())
    B = 0
    ry = 0
    rz = 0
    re = self.gamma_e
    rotx = 0.0
    roty = 0.0
    a,b,c,alpha,beta,gamma = crystal_init_orientation.unit_cell().parameters()

    refine_mode = 'scale_factor'
    xinp = np.array([G,B])

    const_params = (rotx, roty, ry, rz, re, a, b, c, alpha, beta, gamma)
    xopt, cov_x, infodict, mesg, ier = optimize.leastsq(func, xinp,
                                                          args=(I_r_true, observations_original_sel,
                                                                wavelength, alpha_angle_sel,
                                                                iparams.b_refine_d_min,
                                                                crystal_init_orientation,
                                                                spot_pred_x_mm_sel, spot_pred_y_mm_sel,
                                                                detector_distance_mm,
                                                                refine_mode, const_params,
                                                                iparams.partiality_model,
                                                                iparams.flag_volume_correction,
                                                                spot_radius),
                                                          full_output=True, maxfev=100)
    G_fin, B_fin = xopt


    two_theta = observations_original.two_theta(wavelength=wavelength)
    sin_theta_over_lambda_sq = two_theta.sin_theta_over_lambda_sq().data()
    uc = unit_cell((a,b,c,alpha,beta,gamma))
    partiality_init, delta_xy_init, rs_init = calc_partiality_anisotropy_set(uc, rotx, roty,
                                                                             observations_original.indices(), ry, rz,
                                                                             spot_radius, re,
                                                                             two_theta.data(),
                                                                             alpha_angle,
                                                                             wavelength,
                                                                             crystal_init_orientation, spot_pred_x_mm,
                                                                             spot_pred_y_mm,
                                                                             detector_distance_mm,
                                                                             iparams.partiality_model)

    I_o_init = calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                              G, B, partiality_init, rs_init, iparams.flag_volume_correction)
    I_o_fin = calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                              G_fin, B_fin, partiality_init, rs_init, iparams.flag_volume_correction)

    SE_of_the_estimate = standard_error_of_the_estimate(I_r_flex/observations_original.sigmas(),
                                                        I_o_fin/observations_original.sigmas(), 2)
    R_sq = coefficient_of_determination(I_r_flex/observations_original.sigmas(),
                                        I_o_fin/observations_original.sigmas())*100

    CC_init = np.corrcoef(I_r_flex, I_o_init)[0,1]
    CC_final = np.corrcoef(I_r_flex, I_o_fin)[0,1]
    R_init = np.sum(((I_r_flex - I_o_init)/observations_original.sigmas())**2)
    R_final = np.sum(((I_r_flex - I_o_fin)/observations_original.sigmas())**2)
    R_xy_init = 0
    R_xy_final = 0
    CC_iso_init = 0
    CC_iso_final = 0

    return xopt, (SE_of_the_estimate, R_sq, CC_init, CC_final, R_init, R_final, R_xy_init, R_xy_final, CC_iso_init, CC_iso_final)

  def optimize(self, I_r_flex, observations_original,
              wavelength, crystal_init_orientation,
              alpha_angle, spot_pred_x_mm, spot_pred_y_mm, iparams,
              pres_in, observations_non_polar, detector_distance_mm):

    refine_steps = ['scale_factor']

    if iparams.postref.crystal_orientation.flag_on:
      refine_steps.append('crystal_orientation')
    if iparams.postref.reflecting_range.flag_on:
      refine_steps.append('reflecting_range')
    if iparams.postref.unit_cell.flag_on:
      refine_steps.append('unit_cell')

    residuals = []
    for i_sub_cycle in range(iparams.n_postref_sub_cycle):
      residual = [0,0,0,0]
      for j_refine_step in range(len(refine_steps)):
        refine_mode = refine_steps[j_refine_step]
        if refine_mode == 'scale_factor':
          pr_d_min = iparams.postref.scale.d_min
          pr_d_max = iparams.postref.scale.d_max
          pr_sigma_min = iparams.postref.scale.sigma_min
          pr_partiality_min = iparams.postref.scale.partiality_min
        elif refine_mode == 'crystal_orientation':
          pr_d_min = iparams.postref.crystal_orientation.d_min
          pr_d_max = iparams.postref.crystal_orientation.d_max
          pr_sigma_min = iparams.postref.crystal_orientation.sigma_min
          pr_partiality_min = iparams.postref.crystal_orientation.partiality_min
        elif refine_mode == 'reflecting_range':
          pr_d_min = iparams.postref.reflecting_range.d_min
          pr_d_max = iparams.postref.reflecting_range.d_max
          pr_sigma_min = iparams.postref.reflecting_range.sigma_min
          pr_partiality_min = iparams.postref.reflecting_range.partiality_min
        elif refine_mode == 'unit_cell':
          pr_d_min = iparams.postref.unit_cell.d_min
          pr_d_max = iparams.postref.unit_cell.d_max
          pr_sigma_min = iparams.postref.unit_cell.sigma_min
          pr_partiality_min = iparams.postref.unit_cell.partiality_min
          pr_uc_tol = iparams.postref.unit_cell.uc_tolerance


        #filter by resolution
        i_sel_res = observations_original.resolution_filter_selection(d_min=pr_d_min, d_max=pr_d_max)
        observations_original_sel = observations_original.customized_copy(indices=observations_original.indices().select(i_sel_res),
              data=observations_original.data().select(i_sel_res),
              sigmas=observations_original.sigmas().select(i_sel_res))
        alpha_angle_sel = alpha_angle.select(i_sel_res)
        spot_pred_x_mm_sel = spot_pred_x_mm.select(i_sel_res)
        spot_pred_y_mm_sel = spot_pred_y_mm.select(i_sel_res)
        I_ref_sel = I_r_flex.select(i_sel_res)

        #filter by sigma
        i_sel_sigmas = (observations_original_sel.data()/observations_original_sel.sigmas()) > pr_sigma_min
        observations_original_sel = observations_original_sel.customized_copy(indices=observations_original_sel.indices().select(i_sel_sigmas),
            data=observations_original_sel.data().select(i_sel_sigmas),
            sigmas=observations_original_sel.sigmas().select(i_sel_sigmas))
        alpha_angle_sel = alpha_angle_sel.select(i_sel_sigmas)
        spot_pred_x_mm_sel = spot_pred_x_mm_sel.select(i_sel_sigmas)
        spot_pred_y_mm_sel = spot_pred_y_mm_sel.select(i_sel_sigmas)
        I_ref_sel = I_ref_sel.select(i_sel_sigmas)

        #initialize values only in the first sub cycle and the first refine step.
        if j_refine_step == 0 and i_sub_cycle == 0:
          spot_radius = calc_spot_radius(sqr(crystal_init_orientation.reciprocal_matrix()),
                                                     observations_original_sel.indices(), wavelength)
          if pres_in is None:
            xopt_scalefactors, stats = self.optimize_scalefactors(I_r_flex,
                                                                  observations_original,
                                                                  wavelength, crystal_init_orientation,
                                                                  alpha_angle,
                                                                  spot_pred_x_mm,
                                                                  spot_pred_y_mm,
                                                                  iparams,
                                                                  pres_in,
                                                                  observations_non_polar,
                                                                  detector_distance_mm)
            G, B = xopt_scalefactors
            ry = 0
            rz = 0
            re = self.gamma_e
            rotx = 0.0
            roty = 0.0
            a,b,c,alpha,beta,gamma = crystal_init_orientation.unit_cell().parameters()
          else:
            G = pres_in.G
            B = pres_in.B
            ry = pres_in.ry
            rz = pres_in.rz
            re = pres_in.re
            rotx = 0.0
            roty = 0.0
            a,b,c,alpha,beta,gamma = pres_in.unit_cell.parameters()
            crystal_init_orientation = pres_in.crystal_orientation


        #filter by partiality
        two_theta = observations_original_sel.two_theta(wavelength=wavelength).data()
        uc = unit_cell((a,b,c,alpha,beta,gamma))
        partiality_init, delta_xy_init, rs_init = calc_partiality_anisotropy_set(uc, rotx, roty,
                                                                    observations_original_sel.indices(),
                                                                    ry, rz, spot_radius, re, two_theta,
                                                                    alpha_angle_sel, wavelength,
                                                                    crystal_init_orientation,
                                                                    spot_pred_x_mm_sel,
                                                                    spot_pred_y_mm_sel,
                                                                    detector_distance_mm,
                                                                    iparams.partiality_model)
        i_sel_p = partiality_init > pr_partiality_min
        observations_original_sel = observations_original_sel.customized_copy(
            indices=observations_original_sel.indices().select(i_sel_p),
            data=observations_original_sel.data().select(i_sel_p),
            sigmas=observations_original_sel.sigmas().select(i_sel_p))
        alpha_angle_sel = alpha_angle_sel.select(i_sel_p)
        spot_pred_x_mm_sel = spot_pred_x_mm_sel.select(i_sel_p)
        spot_pred_y_mm_sel = spot_pred_y_mm_sel.select(i_sel_p)
        I_ref_sel = I_ref_sel.select(i_sel_p)

        I_r_true = I_ref_sel.as_numpy_array()
        I_o_true = observations_original_sel.data().as_numpy_array()
        sigI_o_true = observations_original_sel.sigmas().as_numpy_array()
        cs = observations_original_sel.crystal_symmetry().space_group().crystal_system()


        if refine_mode == 'scale_factor':
          xinp = np.array([G,B])
          const_params = (rotx, roty, ry, rz, re, a, b, c, alpha, beta, gamma)
        elif refine_mode == 'crystal_orientation':
          xinp = np.array([rotx, roty])
          const_params = (G, B, ry, rz, re, a, b, c, alpha, beta, gamma)
        elif refine_mode == 'reflecting_range':
          xinp = np.array([ry, rz, re])
          const_params = (G, B, rotx, roty, a, b, c, alpha, beta, gamma)
        elif refine_mode == 'unit_cell':
          xinp = prep_input((a,b,c,alpha,beta,gamma), cs)
          const_params = (G, B, rotx, roty, ry, rz, re)

        xopt, cov_x, infodict, mesg, ier = optimize.leastsq(func, xinp,
                                                            args=(I_r_true,
                                                                  observations_original_sel, wavelength, alpha_angle_sel, iparams.b_refine_d_min, crystal_init_orientation, spot_pred_x_mm_sel, spot_pred_y_mm_sel, detector_distance_mm,
                                                                  refine_mode, const_params,
                                                                  iparams.partiality_model,
                                                                  iparams.flag_volume_correction,
                                                                  spot_radius),
                                                            full_output=True, maxfev=100)

        if refine_mode == 'scale_factor':
          G, B = xopt
          residual[0] = np.sum(infodict['fvec']**2)
        elif refine_mode == 'crystal_orientation':
          if (abs(xopt[0]*180/math.pi) <= 2) and (abs(xopt[1]*180/math.pi) <= 2) :
            rotx, roty = xopt
            residual[1] = np.sum(infodict['fvec']**2)
        elif refine_mode == 'reflecting_range':
          ry, rz, re = xopt
          residual[2] = np.sum(infodict['fvec']**2)
        elif refine_mode == 'unit_cell':
          xopt_uc = prep_output(xopt, cs)
          if (abs(xopt_uc[0]-iparams.target_unit_cell.parameters()[0]) <= (pr_uc_tol*iparams.target_unit_cell.parameters()[0]/100) \
              and abs(xopt_uc[1]-iparams.target_unit_cell.parameters()[1]) <= (pr_uc_tol*iparams.target_unit_cell.parameters()[1]/100) \
              and abs(xopt_uc[2]-iparams.target_unit_cell.parameters()[2]) <= (pr_uc_tol*iparams.target_unit_cell.parameters()[2]/100) \
              and abs(xopt_uc[3]-iparams.target_unit_cell.parameters()[3]) <= (pr_uc_tol*iparams.target_unit_cell.parameters()[3]/100) \
              and abs(xopt_uc[4]-iparams.target_unit_cell.parameters()[4]) <= (pr_uc_tol*iparams.target_unit_cell.parameters()[4]/100) \
              and abs(xopt_uc[5]-iparams.target_unit_cell.parameters()[5]) <= (pr_uc_tol*iparams.target_unit_cell.parameters()[5]/100)):
            a,b,c,alpha,beta,gamma = xopt_uc
            residual[3] = np.sum(infodict['fvec']**2)
          else:
            a,b,c,alpha,beta,gamma = crystal_init_orientation.unit_cell().parameters()

      residuals.append(residual)
      if i_sub_cycle > 0:
        if ((residuals[i_sub_cycle-1][2] - residuals[i_sub_cycle][2]) < 0.1) or \
           ((residuals[i_sub_cycle-1][3] - residuals[i_sub_cycle][3]) < 0.0001):
          break


    #do final round of the scale factors
    refine_mode = 'scale_factor'
    xinp = np.array([G,B])
    const_params = (rotx, roty, ry, rz, re, a, b, c, alpha, beta, gamma)
    xopt, cov_x, infodict, mesg, ier = optimize.leastsq(func, xinp,
                                                          args=(I_r_true, observations_original_sel,
                                                                wavelength, alpha_angle_sel,
                                                                iparams.b_refine_d_min,
                                                                crystal_init_orientation,
                                                                spot_pred_x_mm_sel, spot_pred_y_mm_sel,
                                                                detector_distance_mm,
                                                                refine_mode, const_params,
                                                                iparams.partiality_model,
                                                                iparams.flag_volume_correction,
                                                                spot_radius),
                                                          full_output=True, maxfev=100)



    #apply the refined parameters on the full (original) reflection set
    two_theta = observations_original.two_theta(wavelength=wavelength).data()
    sin_theta_over_lambda_sq = observations_original.two_theta(wavelength=wavelength).sin_theta_over_lambda_sq().data()

    if pres_in is None:
      partiality_init, delta_xy_init, rs_init = calc_partiality_anisotropy_set(crystal_init_orientation.unit_cell(),
                                                                      0.0, 0.0,
                                                                      observations_original.indices(),
                                                                      spot_radius, spot_radius,
                                                                      spot_radius, self.gamma_e,
                                                                      two_theta, alpha_angle, wavelength,
                                                                      crystal_init_orientation,
                                                                      spot_pred_x_mm, spot_pred_y_mm,
                                                                      detector_distance_mm,
                                                                      iparams.partiality_model)

      I_o_init = calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                                1, 0, partiality_init, rs_init, iparams.flag_volume_correction)
    else:
      partiality_init, delta_xy_init, rs_init = calc_partiality_anisotropy_set(pres_in.unit_cell,
                                                                      0.0, 0.0,
                                                                      observations_original.indices(),
                                                                      pres_in.ry, pres_in.rz,
                                                                      spot_radius, pres_in.re,
                                                                      two_theta, alpha_angle, wavelength,
                                                                      crystal_init_orientation,
                                                                      spot_pred_x_mm, spot_pred_y_mm,
                                                                      detector_distance_mm,
                                                                      iparams.partiality_model)
      I_o_init = calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                              pres_in.G, pres_in.B, partiality_init, rs_init, iparams.flag_volume_correction)



    partiality_fin, delta_xy_fin, rs_fin = calc_partiality_anisotropy_set(unit_cell((a,b,c,alpha,beta,gamma)),
                                                                              rotx, roty,
                                                                              observations_original.indices(),
                                                                              ry, rz, spot_radius, re,
                                                                              two_theta, alpha_angle, wavelength,
                                                                              crystal_init_orientation,
                                                                              spot_pred_x_mm, spot_pred_y_mm,
                                                                              detector_distance_mm,
                                                                              iparams.partiality_model)
    I_o_fin = calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                              G, B, partiality_fin, rs_fin, iparams.flag_volume_correction)
    SE_of_the_estimate = standard_error_of_the_estimate(I_r_flex/observations_original.sigmas(),
                                                        I_o_fin/observations_original.sigmas(), 13)
    R_sq = coefficient_of_determination(I_r_flex/observations_original.sigmas(),
                                        I_o_fin/observations_original.sigmas())*100

    CC_init = np.corrcoef(I_r_flex, I_o_init)[0,1]
    CC_final = np.corrcoef(I_r_flex, I_o_fin)[0,1]
    R_init = np.sum(((I_r_flex - I_o_init)/observations_original.sigmas())**2)/np.sum(I_r_flex)
    R_final = np.sum(((I_r_flex - I_o_fin)/observations_original.sigmas())**2)/np.sum(I_r_flex)
    R_xy_init = np.sum(delta_xy_init**2)
    R_xy_final = np.sum(delta_xy_fin**2)

    if R_init < R_final:
      CC_final = CC_init
      R_final = R_init
      R_xy_final = R_xy_init
      if pres_in is None:
        spot_radius = calc_spot_radius(sqr(crystal_init_orientation.reciprocal_matrix()),
                                       observations_original_sel.indices(), wavelength)
        G = 1.0
        B = 0.0
        ry = 0
        rz = 0
        re = self.gamma_e
        rotx = 0.0
        roty = 0.0
        a,b,c,alpha,beta,gamma = crystal_init_orientation.unit_cell().parameters()
      else:
        G = pres_in.G
        B = pres_in.B
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
      flag_hklisoin_found = True
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
        else:
          flag_hklisoin_found = False

      if flag_hklisoin_found:
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
      print 'r0 %.4g'%(spot_radius)
      print 're %.4g'%(re)
      print 'uc', a,b,c,alpha,beta,gamma
      print 'S = %.4g'%SE_of_the_estimate
      print 'Target R = %.4g%%'%(R_final)
      print 'Target R (x,y) = %.4g mm.'%(R_xy_final)
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
      one_dsqr = flex.double()
      for i in range(1,n_bins+1):
        i_binner = (binner_indices == i)
        if len(observations_original.data().select(i_binner)) > 0:
          avg_delta_xy_init.append(np.mean(delta_xy_init.select(i_binner)))
          avg_delta_xy_fin.append(np.mean(delta_xy_fin.select(i_binner)))
          one_dsqr.append(1/binner.bin_d_range(i)[1]**2)

      plt.subplot(223)
      plt.plot(one_dsqr, avg_delta_xy_init, linestyle='-', linewidth=2.0, c='r', label='Initial <delta_xy>=%4.2f'%np.mean(delta_xy_init))
      plt.plot(one_dsqr, avg_delta_xy_fin, linestyle='-', linewidth=2.0, c='b', label='Final <delta_xy>=%4.2f'%np.mean(delta_xy_fin))
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
               label='Final ry=%.4g ry=%.4g r0=%.4g re=%.4g'%(ry,rz,spot_radius,re))
      legend = plt.legend(loc='lower right', shadow=False)
      for label in legend.get_texts():
        label.set_fontsize('medium')
      plt.title('Reflecting range')
      plt.xlabel('1/(d^2)')
      plt.ylabel('1/Angstroms')
      plt.show()

    xopt = (G, B, rotx, roty, ry, rz, re,a,b,c,alpha,beta,gamma)

    return xopt, (SE_of_the_estimate, R_sq, CC_init, CC_final, R_init, R_final, R_xy_init, R_xy_final, CC_iso_init, CC_iso_final), len(I_ref_sel)
