"""
leastsqr_handler class refines post-refinement parameters for each frame m in:

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
"""
from __future__ import division
from scipy import optimize, stats, misc
import numpy as np
import math
import matplotlib.pyplot as plt
from cctbx.array_family import flex
from scitbx.matrix import sqr, col
from cctbx.uctbx import unit_cell
from cctbx.crystal_orientation import crystal_orientation, basis_type
from mod_partiality import partiality_handler

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

  spot_radius = math.sqrt(flex.mean(delta_S_all*delta_S_all))

  return spot_radius

def get_crystal_orientation(ortho_matrix, rot_matrix):
  #From orthogonalization matrix and rotation matrix,
  #generate and return crystal orientation
  O = sqr(ortho_matrix).transpose()
  R = sqr(rot_matrix).transpose()
  X = O*R
  co = crystal_orientation(X, basis_type.direct)
  return co

def coefficient_of_determination(y, y_model):
  mean_y = np.mean(y)
  r_sqr = np.sum((y_model - mean_y)**2)/np.sum((y - mean_y)**2)
  return r_sqr

def standard_error_of_the_estimate(y, y_model, n_params):
  s = np.sqrt(np.sum((y - y_model)**2)/(len(y) - n_params))
  return s

def func_scale(params, *args):
  I_r = args[0]
  miller_array_o = args[1]
  wavelength = args[2]
  I_o = miller_array_o.data().as_numpy_array()
  sigI_o = miller_array_o.sigmas().as_numpy_array()
  sin_theta_over_lambda_sq = miller_array_o.two_theta(wavelength=wavelength).sin_theta_over_lambda_sq().data().as_numpy_array()

  G, B = params

  I_o_model = G * np.exp(-2*B*sin_theta_over_lambda_sq) * I_r
  error = (I_o - I_o_model)/sigI_o

  #print 'G=%.4g B=%.4g f=%.4g'%(G, B, np.sum(error**2))
  return error

def func(params, *args):
  I_r = args[0]
  miller_array_o = args[1]
  ph = args[2]
  crystal_init_orientation = args[3]
  alpha_angle_set = args[4]
  crystal_rotation_matrix = crystal_init_orientation.crystal_rotation_matrix()
  I_o = miller_array_o.data().as_numpy_array()
  sigI_o = miller_array_o.sigmas().as_numpy_array()
  miller_indices_original = miller_array_o.indices()
  sin_theta_over_lambda_sq = miller_array_o.two_theta(wavelength=ph.wavelength).sin_theta_over_lambda_sq().data().as_numpy_array()
  two_theta_flex = miller_array_o.two_theta(wavelength=ph.wavelength).data()

  #determine which of 6 uc parameters will be refined, based on the crystal system
  cs = miller_array_o.crystal_symmetry().space_group().crystal_system()
  params_all = prep_output(params, cs)
  G, B, rotx, roty, ry, rz, re, a, b, c, alpha, beta, gamma = params_all

  uc = unit_cell((a,b,c,alpha,beta,gamma))
  crystal_init_orientation = get_crystal_orientation(uc.orthogonalization_matrix(), crystal_rotation_matrix)
  crystal_orientation_model = crystal_init_orientation.rotate_thru((1,0,0), rotx
               ).rotate_thru((0,1,0), roty)
  a_star_model = sqr(crystal_orientation_model.reciprocal_matrix())
  partiality_model_flex = ph.calc_partiality_anisotropy_set(a_star_model, miller_indices_original, ry, rz, re, two_theta_flex, alpha_angle_set)
  partiality_model = partiality_model_flex.as_numpy_array()

  I_o_model = G * np.exp(-2*B*sin_theta_over_lambda_sq) * partiality_model * I_r
  error = (I_o - I_o_model)/sigI_o

  #print 'G=%.4g B=%.4g rotx=%.4g roty=%.4g ry=%.4g rz=%.4g re=%.4g a=%.4g b=%.4g c=%.4g alp=%.4g beta=%.4g gam=%.4g f=%.4g'%(G, B, rotx*180/math.pi, roty*180/math.pi, ry, rz, re, a, b, c, alpha, beta, gamma, np.sum(error**2))
  return error

def func_partiality(x, *args):
  miller_index = args[0]
  crystal_rotation_matrix = args[1]
  wavelength = args[2]
  bragg_angle, alpha_angle = args[3]
  const_params = args[4]
  fmode = args[5]

  G, B, rotx, roty, ry, rz, re, a, b, c, alpha, beta, gamma = const_params
  if fmode == 'G':
    G = x
  elif fmode== 'B':
    B = x
  elif fmode== 'rotx':
    rotx = x
  elif fmode== 'roty':
    roty = x
  elif fmode== 'ry':
    ry = x
  elif fmode== 'rz':
    rz = x
  elif fmode== 're':
    re = x
  elif fmode== 'a':
    a = x
  elif fmode== 'b':
    b = x
  elif fmode== 'c':
    c = x
  elif fmode== 'alpha':
    alpha = x
  elif fmode== 'beta':
    beta = x
  elif fmode== 'gamma':
    gamma = x

  uc = unit_cell((a,b,c,alpha,beta,gamma))
  crystal_init_orientation = get_crystal_orientation(uc.orthogonalization_matrix(), crystal_rotation_matrix)
  crystal_orientation_model = crystal_init_orientation.rotate_thru((1,0,0), rotx
                 ).rotate_thru((0,1,0), roty)
  a_star_model = sqr(crystal_orientation_model.reciprocal_matrix())
  ph = partiality_handler(wavelength, 0)
  partiality, dummy, dummy = ph.calc_partiality_anisotropy(a_star_model, miller_index, ry, rz, re, bragg_angle, alpha_angle)

  return partiality

def prep_input(params, cs):
  #From crystal system cs, determine refined parameters
  #0)G, 1)B, 2)rotx, 3)roty, 4)ry, 5)rz, 6)re, 7)a, 8)b, 9)c, 10)alpha, 11)beta, 12)gamma
  G, B, rotx, roty, ry, rz, re, a, b, c, alpha, beta, gamma = params
  if cs == 'Triclinic':
    x0 = params
  elif cs == 'Monoclinic':
    x0 = np.array([G, B, rotx, roty, ry, rz, re,a,b,c,beta])
  elif cs == 'Orthorhombic':
    x0 = np.array([G, B, rotx, roty, ry, rz, re,a,b,c])
  elif cs == 'Tetragonal':
    x0 = np.array([G, B, rotx, roty, ry, rz, re,a,c])
  elif cs == 'Trigonal' or cs == 'Hexagonal':
    x0 = np.array([G, B, rotx, roty, ry, rz, re,a,c])
  elif cs == 'Cubic':
    x0 = np.array([G, B, rotx, roty, ry, rz, re,a])

  return x0

def prep_output(params, cs):
  if cs == 'Triclinic':
    xopt = params
  elif cs == 'Monoclinic':
    xopt = np.array([params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8],params[9],90,params[10],90])
  elif cs == 'Orthorhombic':
    xopt = np.array([params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8],params[9],90,90,90])
  elif cs == 'Tetragonal':
    xopt = np.array([params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[7],params[8],90,90,90])
  elif cs == 'Trigonal' or cs == 'Hexagonal':
    xopt = np.array([params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[7],params[8],90,90,120])
  elif cs == 'Cubic':
    xopt = np.array([params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[7],params[7],90,90,90])

  return xopt

def prep_variance(params, cs):
  if cs == 'Triclinic':
    se_xopt = params
  elif cs == 'Monoclinic':
    se_xopt = np.array([params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8],params[9],0,params[10],0])
  elif cs == 'Orthorhombic':
    se_xopt = np.array([params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8],params[9],0,0,0])
  elif cs == 'Tetragonal':
    se_xopt = np.array([params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],0,params[8],0,0,0])
  elif cs == 'Trigonal' or cs == 'Hexagonal':
    se_xopt = np.array([params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],0,params[8],0,0,0])
  elif cs == 'Cubic':
    se_xopt = np.array([params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],0,0,0,0,0])

  return se_xopt

class leastsqr_handler(object):
  '''
  A wrapper class for least-squares refinement
  '''

  def __init__(self):
    '''
    Do nothing
    '''

  def optimize(self, I_r_flex, observations_original,
              wavelength, crystal_init_orientation, alpha_angle_set, iph):

    uc_len_tol = 3.5
    uc_angle_tol = 2.0

    assert len(alpha_angle_set)==len(observations_original.indices()), 'Size of alpha angles and observations are not equal %6.0f, %6.0f'%(len(alpha_angle_set),len(observations_original.indices()))

    uc_init = crystal_init_orientation.unit_cell()
    uc_init_params = uc_init.parameters()
    I_r_true = I_r_flex.as_numpy_array()
    I_o_true = observations_original.data().as_numpy_array()
    sigI_o_true = observations_original.sigmas().as_numpy_array()
    sin_theta_over_lambda_sq = observations_original.two_theta(wavelength=wavelength).sin_theta_over_lambda_sq().data().as_numpy_array()
    two_theta_flex = observations_original.two_theta(wavelength=wavelength).data()
    cs = observations_original.crystal_symmetry().space_group().crystal_system()

    #calculate spot_radius
    a_star_true = sqr(crystal_init_orientation.reciprocal_matrix())
    spot_radius = calc_spot_radius(a_star_true, observations_original.indices(), wavelength)
    ph = partiality_handler(wavelength, spot_radius)

    #1. first optain best G and k from linregress
    x0 = np.array([1, 0])
    xopt_scale, success = optimize.leastsq(func_scale, x0, args=(I_r_true, observations_original, wavelength))

    #2. optimize rot, rs, uc
    x0_all = np.array([xopt_scale[0], xopt_scale[1], 0*math.pi/180, 0*math.pi/180, spot_radius, spot_radius, 0.0026,
        uc_init_params[0], uc_init_params[1],uc_init_params[2], uc_init_params[3],uc_init_params[4], uc_init_params[5]])
    x0 = prep_input(x0_all, cs)
    xopt_limit, cov_x, infodict, errmsg, success = optimize.leastsq(func, x0, args=(I_r_true, observations_original, ph, crystal_init_orientation, alpha_angle_set), full_output=True)
    xopt = prep_output(xopt_limit, cs)
    G, B, rotx, roty, ry, rz, re, a, b, c, alpha, beta, gamma = xopt

    #3. decide wheter to take the refined parameters
    if (abs(a-uc_init_params[0]) > uc_len_tol or abs(b-uc_init_params[1]) > uc_len_tol or abs(c-uc_init_params[2]) > uc_len_tol \
        or abs(alpha-uc_init_params[3]) > uc_angle_tol or abs(beta-uc_init_params[4]) > uc_angle_tol or abs(gamma-uc_init_params[5]) > uc_angle_tol):
      print 'Refinement failed - unit-cell parameters exceed the limits (%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f)'%(a,b,c,alpha,beta,gamma)
      if iph.flag_force_accept_all_frames:
        a, b, c, alpha, beta, gamma = uc_init_params
        print ' flag_force_accept_all_frames is on, reset unit cell to (%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f)'%(a,b,c,alpha,beta,gamma)
        rotx = 0
        roty = 0
        ry = spot_radius
        rz = spot_radius
        re = 0.0026
        G, B = xopt_scale
        xopt = np.array([G, B, rotx, roty, ry, rz, re, a, b, c, alpha, beta, gamma])

    #caclculate stats
    uc_opt = unit_cell((a,b,c,alpha,beta,gamma))
    crystal_orientation_opt = get_crystal_orientation(uc_opt.orthogonalization_matrix(), crystal_init_orientation.crystal_rotation_matrix())
    crystal_orientation_model = crystal_orientation_opt.rotate_thru((1,0,0), rotx
                   ).rotate_thru((0,1,0), roty)
    a_star_model = sqr(crystal_orientation_model.reciprocal_matrix())
    partiality_model_flex = ph.calc_partiality_anisotropy_set(a_star_model, observations_original.indices(), ry, rz, re, two_theta_flex, alpha_angle_set)
    partiality_model = partiality_model_flex.as_numpy_array()

    I_o_model = G * np.exp(-2*B*sin_theta_over_lambda_sq) * partiality_model * I_r_true

    SE_of_the_estimate = standard_error_of_the_estimate(I_o_true/sigI_o_true, I_o_model/sigI_o_true, len(x0))
    R_sq = coefficient_of_determination(I_o_true/sigI_o_true, I_o_model/sigI_o_true)*100
    CC = np.corrcoef(I_o_true/sigI_o_true, I_o_model/sigI_o_true)[0,1]
    var_I_p = ((observations_original.sigmas()/observations_original.data())**2).as_numpy_array()

    #calculate standard error for the parameters
    if cov_x is None:
      se_xopt = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0])
      var_k = flex.double([0]*len(I_o_model))
      var_p = flex.double([0]*len(I_o_model))
      SE_I = np.sqrt(var_I_p)*(I_o_model)
    else:
      pcov = cov_x * (SE_of_the_estimate**2)
      var_xopt_limit = np.array([pcov[i,i] for i in range(len(x0))])
      var_xopt = prep_variance(var_xopt_limit, cs)
      se_xopt = np.sqrt(var_xopt)

      #calculate error propagation in I_full
      se_G, se_B, se_rotx, se_roty, se_ry, se_rz, se_re, se_a, se_b, se_c, se_alpha, se_beta, se_gamma = se_xopt

      ###k = G_0 * exp(-2*B*sin_theta_over_lambda_sq)
      se_k_sq = (np.exp(-2*B*sin_theta_over_lambda_sq)*se_G)**2 + (-2 * G * sin_theta_over_lambda_sq * np.exp(-2*B*sin_theta_over_lambda_sq) * se_B)**2
      var_k = se_k_sq/((G*np.exp(-2*B*sin_theta_over_lambda_sq))**2)

      ###use finite differences to propagate errors in partility function p.
      var_p_sq = flex.double()
      fmode_arr = ('G','B','rotx','roty','ry','rz','re','a','b','c','alpha','beta','gamma')
      for miller_index, bragg_angle, alpha_angle in zip(observations_original.indices(), two_theta_flex, alpha_angle_set):
        Dp = flex.double()
        cn_fmode = 0
        for fmode in fmode_arr:
          dp = misc.derivative(func_partiality, xopt[cn_fmode], args=(miller_index,
              crystal_init_orientation.crystal_rotation_matrix(), wavelength, (bragg_angle, alpha_angle), xopt, fmode))
          cn_fmode += 1
          Dp.append(dp)

        se_p_sq = 0
        for dp, se in zip(Dp, se_xopt):
          se_p_sq += (dp*se)**2

        var_p_sq.append(se_p_sq)


      var_p = (var_p_sq/(partiality_model_flex**2)).as_numpy_array()
      SE_I = np.sqrt(var_I_p + var_k + var_p)*(I_o_model)

    """
    print 'Predictor'
    print 'G %.4g'%(G)
    print 'B-factor %.4g'%(B)
    print 'rotx %.4g'%(rotx*180/math.pi)
    print 'roty %.4g'%(roty*180/math.pi)
    print 'ry %.4g'%(ry)
    print 'rz %.4g'%(rz)
    print 're %.4g'%(re)
    print 'uc', uc_opt
    print 'S = %.4g'%SE_of_the_estimate
    print 'R-Sq = %.4g%%'%(R_sq)
    print 'CC = %.4g'%(CC)

    plt.scatter(I_r_true, I_o_true,s=10, marker='x', c='r')
    plt.scatter(I_r_true, I_o_model,s=10, marker='o', c='b')
    plt.title('G=%.4g B=%.4g rotx=%.4g roty=%.4g R-sq=%.4g%%'%(xopt[0], xopt[1], xopt[2]*180/math.pi, xopt[3]*180/math.pi, R_sq))
    plt.xlabel('I_ref')
    plt.ylabel('I_obs')
    plt.show()
    """

    return xopt, se_xopt, (SE_of_the_estimate, R_sq, CC), partiality_model, SE_I, var_I_p, var_k, var_p
