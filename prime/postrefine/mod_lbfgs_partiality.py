from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.uctbx import unit_cell
from .mod_partiality import partiality_handler
"""
lbfgs_partiality_handler
calculate the minimization function according to given arguments and return sum(error^2)
"""
class lbfgs_partiality_handler(object):
  """
  lbfgs_handler optimizes set of parameters (params) by fitting
  data[0] to data [1] using given function (func).
  """
  def __init__(self):
    """
    Intialitze parameters
    """

  def prep_input(self, params, cs):
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

  def prep_output(self, params, cs):
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

  def func(self, params, args):
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
      rotx, roty, ry, rz, r0, re, nu, a, b, c, alpha, beta, gamma = const_params
    elif refine_mode == 'crystal_orientation':
      rotx, roty = params
      G, B, ry, rz, r0, re, nu, a, b, c, alpha, beta, gamma = const_params
    elif refine_mode == 'reflecting_range':
      ry, rz, r0, re, nu = params
      G, B, rotx, roty, a, b, c, alpha, beta, gamma = const_params
    elif refine_mode == 'unit_cell':
      a, b, c, alpha, beta, gamma = self.prep_output(params, cs)
      G, B, rotx, roty, ry, rz, r0, re, nu = const_params
    elif refine_mode == 'allparams':
      a, b, c, alpha, beta, gamma = self.prep_output(params[7:], cs)
      rotx, roty, ry, rz, r0, re, nu = params[0:7]
      G,B = const_params
    try:
      uc = unit_cell((a,b,c,alpha,beta,gamma))
    except Exception:
      return None
    G,B,rotx,roty,ry,rz,r0,re,nu,a,b,c,alpha,beta,gamma = flex.double([G,B,rotx,roty,ry,rz,r0,re,nu,a,b,c,alpha,beta,gamma])
    ph = partiality_handler()
    p_calc_flex, delta_xy_flex, rs_set, dummy = ph.calc_partiality_anisotropy_set(\
      uc, rotx, roty, miller_indices_original, ry, rz, r0, re, nu, two_theta_flex,
      alpha_angle_set, wavelength, crystal_init_orientation, spot_pred_x_mm_set,
      spot_pred_y_mm_set, detector_distance_mm, partiality_model, flag_beam_divergence)
    if miller_array_o.d_min() < b_refine_d_min:
      I_o_full = ph.calc_full_refl(I_o, sin_theta_over_lambda_sq, G, B, p_calc_flex, rs_set, flag_volume_correction)
    else:
      I_o_full = ph.calc_full_refl(I_o, sin_theta_over_lambda_sq, G, b0, p_calc_flex, rs_set, flag_volume_correction)
    if refine_mode == 'unit_cell':
      error = delta_xy_flex
    else:
      error = ((I_r - I_o_full)/sigI_o)
    """
    import math
    print refine_mode, 'G:%.4g B:%.4g rotx:%.4g roty:%.4g r0:%.4g re:%.4g nu:%6.4f a:%.4g b:%.4g c:%.4g fpr:%.4g fxy:%.4g n_refl=%5.0f'%(G, B, rotx*180/math.pi, roty*180/math.pi, r0, re, nu, a, b, c, flex.sum(((I_r - I_o_full)/sigI_o)**2), flex.sum(delta_xy_flex**2), len(I_o_full))
    """
    return error
