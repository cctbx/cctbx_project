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
from __future__ import absolute_import, division, print_function
import math
from cctbx.array_family import flex
from scitbx.matrix import sqr
from cctbx.uctbx import unit_cell
from cctbx.crystal_orientation import crystal_orientation
from .mod_lbfgs import lbfgs_handler
from .mod_lbfgs_partiality import lbfgs_partiality_handler
from .mod_partiality import partiality_handler

def coefficient_of_determination(y, y_model):
  mean_y = flex.mean(y)
  r_sqr = flex.sum((y_model - mean_y)**2)/flex.sum((y - mean_y)**2)
  return r_sqr

def standard_error_of_the_estimate(y, y_model, n_params):
  s = math.sqrt(flex.sum((y - y_model)**2)/(len(y) - n_params))
  return s

def good_unit_cell(uc_params, iparams, uc_tol, target_unit_cell=None):
  if iparams is None:
    expected_unit_cell = target_unit_cell
  else:
    expected_unit_cell = iparams.target_unit_cell
  flag_good_uc = False
  if (abs(uc_params[0]-expected_unit_cell.parameters()[0]) \
      <= (uc_tol*expected_unit_cell.parameters()[0]/100) \
                and abs(uc_params[1]-expected_unit_cell.parameters()[1]) \
                <= (uc_tol*expected_unit_cell.parameters()[1]/100) \
                and abs(uc_params[2]-expected_unit_cell.parameters()[2]) \
                <= (uc_tol*expected_unit_cell.parameters()[2]/100) \
                and abs(uc_params[3]-expected_unit_cell.parameters()[3]) \
                <= (uc_tol*expected_unit_cell.parameters()[3]/100) \
                and abs(uc_params[4]-expected_unit_cell.parameters()[4]) \
                <= (uc_tol*expected_unit_cell.parameters()[4]/100) \
                and abs(uc_params[5]-expected_unit_cell.parameters()[5]) \
                <= (uc_tol*expected_unit_cell.parameters()[5]/100)):
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
              pres_in, observations_non_polar, detector_distance_mm, const_params):
    ph = partiality_handler()
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
      G, B, b0 = pres_in.G, pres_in.B, pres_in.B
    else:
      G,B,b0 = (1,0,0)
    refine_mode = 'scale_factor'
    xinp = flex.double([G,B])
    args = (I_r_true, observations_original_sel, wavelength, alpha_angle_sel,
            crystal_init_orientation, spot_pred_x_mm_sel, spot_pred_y_mm_sel,
            detector_distance_mm, refine_mode, const_params, b0, None, iparams)
    lh = lbfgs_handler(current_x=xinp, args=args)
    G_fin, B_fin = (lh.x[0], lh.x[1])
    rotx, roty, ry, rz, r0, re, voigt_nu, a, b, c, alpha, beta, gamma = const_params
    two_theta = observations_original.two_theta(wavelength=wavelength)
    sin_theta_over_lambda_sq = two_theta.sin_theta_over_lambda_sq().data()
    uc = unit_cell((a,b,c,alpha,beta,gamma))
    ph = partiality_handler()
    partiality_init, delta_xy_init, rs_init, dummy = ph.calc_partiality_anisotropy_set(uc, rotx, roty,
                                                                             observations_original.indices(),
                                                                             ry, rz, r0, re, voigt_nu,
                                                                             two_theta.data(),
                                                                             alpha_angle,
                                                                             wavelength,
                                                                             crystal_init_orientation,
                                                                             spot_pred_x_mm,
                                                                             spot_pred_y_mm,
                                                                             detector_distance_mm,
                                                                             iparams.partiality_model,
                                                                             iparams.flag_beam_divergence)
    I_o_init = ph.calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                              G, B, partiality_init, rs_init)
    I_o_fin = ph.calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                              G_fin, B_fin, partiality_init, rs_init)
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
    return flex.double(list(lh.x)), (SE_of_the_estimate, R_sq, CC_init, CC_final, R_init, R_final, R_xy_init, R_xy_final, CC_iso_init, CC_iso_final)

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
    G, B, rotx, roty, ry, rz, r0, re, voigt_nu, a, b, c, alpha, beta, gamma = init_params
    #filter by partiality
    two_theta = observations_original_sel.two_theta(wavelength=wavelength).data()
    uc = unit_cell((a,b,c,alpha,beta,gamma))
    ph = partiality_handler()
    partiality_init, delta_xy_init, rs_init, dummy = ph.calc_partiality_anisotropy_set(uc, rotx, roty,
                                                                    observations_original_sel.indices(),
                                                                    ry, rz, r0, re, voigt_nu, two_theta,
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
    ph = partiality_handler()
    lph = lbfgs_partiality_handler()
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
    spot_radius = ph.calc_spot_radius(sqr(crystal_init_orientation.reciprocal_matrix()),
                                                     observations_original_sel.indices(), wavelength)
    if pres_in is None:
      ry, rz, r0, re, voigt_nu, rotx, roty = 0, 0, spot_radius, iparams.gamma_e, iparams.voigt_nu, 0.0, 0.0
      #apply constrain on the unit cell using crystal system
      uc_scale_inp = lph.prep_input(observations_original.unit_cell().parameters(), cs)
      uc_scale_constrained = lph.prep_output(uc_scale_inp, cs)
      a,b,c,alpha,beta,gamma = uc_scale_constrained
      const_params_scale = (rotx, roty, ry, rz, r0, re, voigt_nu, a, b, c, alpha, beta, gamma)
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
      G, B, ry, rz, r0, re, voigt_nu, rotx, roty = pres_in.G, pres_in.B, pres_in.ry, pres_in.rz, pres_in.r0, pres_in.re, pres_in.voigt_nu, 0.0 , 0.0
      a,b,c,alpha,beta,gamma = pres_in.unit_cell.parameters()
      crystal_init_orientation = pres_in.crystal_orientation
    #filter by partiality
    two_theta = observations_original_sel.two_theta(wavelength=wavelength).data()
    uc = unit_cell((a,b,c,alpha,beta,gamma))
    partiality_init, delta_xy_init, rs_init, dummy = ph.calc_partiality_anisotropy_set(uc, rotx, roty,
                                                                    observations_original_sel.indices(),
                                                                    ry, rz, r0, re, voigt_nu, two_theta,
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
    #calculate initial residual_xy error
    const_params_uc = (G, B, rotx, roty, ry, rz, r0, re, voigt_nu)
    xinp_uc = lph.prep_input((a,b,c,alpha,beta,gamma), cs)
    args_uc = (I_r_true, observations_original_sel, wavelength, alpha_angle_sel,
            crystal_init_orientation, spot_pred_x_mm_sel, spot_pred_y_mm_sel,
            detector_distance_mm, 'unit_cell', const_params_uc, B, miller_array_iso, iparams)
    uc_params_err = lph.func(xinp_uc, args_uc)
    init_residual_xy_err = flex.sum(uc_params_err**2)
    #calculate initial residual_pr error
    const_params_all= (G,B)
    xinp_all = flex.double([rotx, roty, ry, rz, r0, re, voigt_nu])
    xinp_all.extend(lph.prep_input((a,b,c,alpha,beta,gamma), cs))
    args_all = (I_r_true, observations_original_sel, wavelength, alpha_angle_sel,
            crystal_init_orientation, spot_pred_x_mm_sel, spot_pred_y_mm_sel,
            detector_distance_mm, 'allparams', const_params_all, B, miller_array_iso, iparams)
    all_params_err = lph.func(xinp_all, args_all)
    init_residual_err = flex.sum(all_params_err**2)
    #keep in list
    t_pr_list = [init_residual_err]
    t_xy_list = [init_residual_xy_err]
    refined_params_hist = [(G, B, rotx, roty, ry, rz, r0, re, voigt_nu, a, b, c, alpha, beta, gamma)]
    txt_out = ''
    for i_sub_cycle in range(iparams.n_postref_sub_cycle):
      for j_refine_step in range(len(refine_steps)):
        refine_mode = refine_steps[j_refine_step]
        #prepare data
        init_params = (G, B, rotx, roty, ry, rz, r0, re, voigt_nu, a, b, c, alpha, beta, gamma)
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
          const_params = (G, B, ry, rz, r0, re, voigt_nu, a, b, c, alpha, beta, gamma)
        elif refine_mode == 'reflecting_range':
          xinp = flex.double([ry, rz, r0, re, voigt_nu])
          const_params = (G, B, rotx, roty, a, b, c, alpha, beta, gamma)
        elif refine_mode == 'unit_cell':
          xinp = lph.prep_input((a,b,c,alpha,beta,gamma), cs)
          const_params = (G, B, rotx, roty, ry, rz, r0, re, voigt_nu)
        elif refine_mode == 'allparams':
          xinp = flex.double([rotx, roty, ry, rz, r0, re, voigt_nu])
          xinp.extend(lph.prep_input((a,b,c,alpha,beta,gamma), cs))
          const_params = (G,B)
        args=(I_r_true, observations_original_sel, wavelength, alpha_angle_sel,
              crystal_init_orientation, spot_pred_x_mm_sel, spot_pred_y_mm_sel,
              detector_distance_mm, refine_mode, const_params, B, miller_array_iso, iparams)
        lh = lbfgs_handler(current_x=xinp, args=args)
        xopt = flex.double(list(lh.x))
        if refine_mode == 'crystal_orientation' or \
            refine_mode == 'reflecting_range' or refine_mode == 'allparams':
          current_residual_err = lh.f
          #calculate residual_xy_error (for refine_mode = SF, CO, RR, and all params)
          xinp_uc = lph.prep_input((a,b,c,alpha,beta,gamma), cs)
          if refine_mode == 'crystal_orientation':
            rotx, roty = xopt
          elif refine_mode == 'reflecting_range':
            ry, rz, r0, re, voigt_nu = xopt
          elif refine_mode == 'allparams':
            rotx, roty, ry, rz, r0, re, voigt_nu = xopt[:7]
            xinp_uc = xopt[7:]
            a, b, c, alpha, beta, gamma = lph.prep_output(xinp_uc, cs)
          const_params_uc = (G, B, rotx, roty, ry, rz, r0, re, voigt_nu)
          xinp_uc = lph.prep_input((a,b,c,alpha,beta,gamma), cs)
          args_uc = (I_r_true, observations_original_sel, wavelength, alpha_angle_sel,
                  crystal_init_orientation, spot_pred_x_mm_sel, spot_pred_y_mm_sel,
                  detector_distance_mm, 'unit_cell', const_params_uc, B, miller_array_iso, iparams)
          uc_params_err = lph.func(xinp_uc, args_uc)
          current_residual_xy_err = flex.sum(uc_params_err**2)
        elif refine_mode == 'unit_cell':
          current_residual_xy_err = lh.f
          xopt_uc = lph.prep_output(xopt, cs)
          a, b, c, alpha, beta, gamma = xopt_uc
          #check the unit-cell with the reference intensity
          xinp = flex.double([rotx, roty, ry, rz, r0, re, voigt_nu])
          xinp.extend(lph.prep_input((a, b, c, alpha, beta, gamma), cs))
          const_params_all = (G,B)
          args_all = (I_r_true, observations_original_sel, wavelength, alpha_angle_sel,
            crystal_init_orientation, spot_pred_x_mm_sel, spot_pred_y_mm_sel,
            detector_distance_mm, 'allparams', const_params_all, B, miller_array_iso, iparams)
          all_params_err = lph.func(xinp_all, args_all)
          current_residual_err = flex.sum(all_params_err**2)
        flag_success = False
        if refine_mode == 'allparams':
          #if allparams refinement, only check the post-refine target function
          if current_residual_err < (t_pr_list[len(t_pr_list)-1] + \
              (t_pr_list[len(t_pr_list)-1]*iparams.postref.residual_threshold/100)):
            t_pr_list.append(current_residual_err)
            t_xy_list.append(current_residual_xy_err)
            refined_params_hist.append((G, B, rotx, roty, ry, rz, r0, re, voigt_nu, a, b, c, alpha, beta, gamma))
            flag_success = True
        else:
          if current_residual_err < (t_pr_list[len(t_pr_list)-1] + \
                (t_pr_list[len(t_pr_list)-1]*iparams.postref.residual_threshold/100)):
            if current_residual_xy_err < (t_xy_list[len(t_xy_list)-1] + \
                (t_xy_list[len(t_xy_list)-1]*iparams.postref.residual_threshold_xy/100)):
              t_pr_list.append(current_residual_err)
              t_xy_list.append(current_residual_xy_err)
              refined_params_hist.append((G, B, rotx, roty, ry, rz, r0, re, voigt_nu, a, b, c, alpha, beta, gamma))
              flag_success = True
        if flag_success is False:
          G,B,rotx,roty,ry,rz,r0,re,voigt_nu,a,b,c,alpha,beta,gamma = refined_params_hist[len(refined_params_hist)-1]
        tmp_txt_out = refine_mode + ' %3.0f %6.4f %6.4f %6.4f %6.4f %10.8f %10.8f %10.8f %10.8f %10.8f %6.3f %6.3f %.4g %6.3f\n'%(i_sub_cycle,G,B,rotx*180/math.pi,roty*180/math.pi,ry,rz,r0,re,voigt_nu,a,c,t_pr_list[len(t_pr_list)-1],t_xy_list[len(t_pr_list)-1])
        txt_out += tmp_txt_out
    #apply the refined parameters on the full (original) reflection set
    two_theta = observations_original.two_theta(wavelength=wavelength).data()
    sin_theta_over_lambda_sq = observations_original.two_theta(wavelength=wavelength).sin_theta_over_lambda_sq().data()
    if pres_in is None:
      partiality_init, delta_xy_init, rs_init, rh_init = ph.calc_partiality_anisotropy_set(\
          observations_original.unit_cell(),0.0, 0.0,observations_original.indices(),
          0, 0, spot_radius, iparams.gamma_e, iparams.voigt_nu,
          two_theta, alpha_angle, wavelength,
          crystal_init_orientation,spot_pred_x_mm, spot_pred_y_mm,detector_distance_mm,
          iparams.partiality_model,iparams.flag_beam_divergence)
      I_o_init = ph.calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                                1, 0, partiality_init, rs_init)
    else:
      partiality_init, delta_xy_init, rs_init, rh_init = ph.calc_partiality_anisotropy_set(\
          pres_in.unit_cell,0.0, 0.0,observations_original.indices(),
          pres_in.ry, pres_in.rz,pres_in.r0, pres_in.re, pres_in.voigt_nu,
          two_theta, alpha_angle, wavelength,
          crystal_init_orientation,spot_pred_x_mm, spot_pred_y_mm,detector_distance_mm,
          iparams.partiality_model,iparams.flag_beam_divergence)
      I_o_init = ph.calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                              pres_in.G, pres_in.B, partiality_init, rs_init)
    partiality_fin, delta_xy_fin, rs_fin, rh_fin = ph.calc_partiality_anisotropy_set(\
        unit_cell((a,b,c,alpha,beta,gamma)),rotx, roty,observations_original.indices(),
        ry, rz, r0, re, voigt_nu, two_theta, alpha_angle, wavelength,crystal_init_orientation,
        spot_pred_x_mm, spot_pred_y_mm,detector_distance_mm,
        iparams.partiality_model,iparams.flag_beam_divergence)
    I_o_fin = ph.calc_full_refl(observations_original.data(), sin_theta_over_lambda_sq,
                              G, B, partiality_fin, rs_fin)
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
    if R_init < R_final or re > (iparams.gamma_e * 3):
      CC_final = CC_init
      R_final = R_init
      R_xy_final = R_xy_init
      if pres_in is None:
        G,B,r0,ry,rz,re,rotx,roty = (1.0,0.0,spot_radius,0.0,0.0,iparams.gamma_e,0.0,0.0)
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
    xopt = (G, B, rotx, roty, ry, rz, r0, re, voigt_nu, a,b,c,alpha,beta,gamma)
    return xopt, (SE_of_the_estimate, R_sq, CC_init, CC_final, R_init, R_final, R_xy_init, R_xy_final, CC_iso_init, CC_iso_final), len(I_ref_sel)
