from __future__ import division
from cctbx.array_family import flex
from scitbx.matrix import sqr, col
from cctbx.crystal_orientation import crystal_orientation, basis_type

class partiality_handler(object):
  """
  mod_partiality:
  1. Calculate partiality for given
  miller indices, crystal orientation, unit cell, wavelength.
  2. Cacluate spot centroid delta distance
  """
  def __init__(self):
    """
    Intialitze parameters
    """

  def calc_full_refl(self, I_o_p_set, sin_theta_over_lambda_sq_set,
                   G, B, p_set, rs_set, flag_volume_correction=True):
    I_o_full_set = I_o_p_set/(G * flex.exp(-2*B*sin_theta_over_lambda_sq_set) * p_set)
    return I_o_full_set

  def calc_spot_radius(self, a_star_matrix, miller_indices, wavelength):
    #calculate spot_radius based on rms delta_S for all spots
    S0 = -1*col((0,0,1./wavelength))
    sd_array = a_star_matrix.elems * miller_indices.as_vec3_double() + S0.elems
    rh_set = sd_array.norms() - (1/wavelength)
    return rh_set.standard_deviation_of_the_sample()

  def calc_partiality_anisotropy_set(self, my_uc, rotx, roty, miller_indices, ry, rz, r0, re,
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
    x = A_star.elems * miller_indices.as_vec3_double()
    sd_array = x + S0.elems
    rh_set = sd_array.norms() - (1/wavelength)
    #calculate partiality
    partiality_set = ((rs_set**2)/((2*(rh_set**2))+(rs_set**2)))
    #partiality_set = 1 - ((rh_set**2)/(rs_set**2))
    #calculate delta_xy
    d_ratio = -detector_distance_mm/sd_array.parts()[2]
    calc_xy_array = flex.vec3_double(sd_array.parts()[0]*d_ratio, \
        sd_array.parts()[1]*d_ratio, flex.double([0]*len(d_ratio)))
    pred_xy_array = flex.vec3_double(spot_pred_x_mm_set, spot_pred_y_mm_set, flex.double([0]*len(d_ratio)))
    delta_xy_set = (pred_xy_array - calc_xy_array).norms()
    return partiality_set, delta_xy_set, rs_set, rh_set
