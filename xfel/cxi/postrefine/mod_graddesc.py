from __future__ import division
import math
from cctbx.array_family import flex
from scitbx.matrix import col, sqr
from mod_partiality import partiality_handler

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


class graddesc_handler(object):

  def __init__(self):
    '''
    Constructor
    '''

  def compute_cost(self, I_ref, I_obs, sigI_obs, crystal_init_orientation,
            miller_indices_ori, wavelength, current_values):
    G = current_values[0]
    rotx = current_values[1]
    roty = current_values[2]
    rotz = current_values[3]
    rs = current_values[4]
    m = len(I_ref)

    effective_orientation = crystal_init_orientation.rotate_thru((1,0,0),rotx
             ).rotate_thru((0,1,0),roty
             ).rotate_thru((0,0,1),rotz)

    effective_a_star = sqr(effective_orientation.reciprocal_matrix())
    ph = partiality_handler(wavelength, rs)
    partiality = flex.double([ph.calc_partiality(effective_a_star, miller_index)[0] for miller_index in miller_indices_ori])

    J = flex.sum(((((G * I_obs)/partiality)- I_ref)/ sigI_obs)**2)/(2*m)

    return J


  def gradient_descend(self,I_ref, I_obs, sigI_obs, crystal_init_orientation,
        miller_indices_ori, wavelength, parameters, alpha, n_iters):

    m = len(I_ref)

    G = parameters[0]
    rotx = parameters[1]
    roty = parameters[2]
    rotz = parameters[3]
    rs = parameters[4]

    delta_rot = 0.00001*math.pi/180
    J_history = flex.double([0]*n_iters)

    converge_iter = n_iters - 1
    for i in range(n_iters):
      #0. Calculate current partiality based on current rotation
      crystal_current_orientation = crystal_init_orientation.rotate_thru((1,0,0),rotx
             ).rotate_thru((0,1,0),roty
             ).rotate_thru((0,0,1),rotz)
      a_star = sqr(crystal_current_orientation.reciprocal_matrix())
      ph = partiality_handler(wavelength, rs)
      partiality = flex.double([ph.calc_partiality(a_star, miller_index)[0] for miller_index in miller_indices_ori])

      #1. Calculate partial derivatives of J function
      Ai = sqr(crystal_init_orientation.reciprocal_matrix())
      Rx = col((1,0,0)).axis_and_angle_as_r3_rotation_matrix(rotx)
      Ry = col((0,1,0)).axis_and_angle_as_r3_rotation_matrix(roty)
      Rz = col((0,0,1)).axis_and_angle_as_r3_rotation_matrix(rotz)
      dRx_drotx = col((1,0,0)).axis_and_angle_as_r3_derivative_wrt_angle(rotx)
      dRy_droty = col((0,1,0)).axis_and_angle_as_r3_derivative_wrt_angle(roty)
      dRz_drotz = col((0,0,1)).axis_and_angle_as_r3_derivative_wrt_angle(rotz)
      dA_drotx = Rz * Ry * dRx_drotx * Ai
      dA_droty = Rz * dRy_droty * Rx * Ai
      dA_drotz = dRz_drotz * Ry * Rx * Ai

      s0 = -1*col((0,0,1./wavelength))
      pDJ_pDrotx = 0
      pDJ_pDroty = 0
      pDJ_pDrotz = 0
      pDJ_pDG = 0
      pDJ_pDrs = 0
      for j in range(len(miller_indices_ori)):
        miller_index = miller_indices_ori[j]
        hvec = col(miller_index)
        xvec = a_star * hvec
        svec = s0 + xvec
        p, rh = ph.calc_partiality(a_star, miller_index)


        #for rotx
        pDxvec_pDrotx = dA_drotx * hvec
        pDrh_pDrotx = (svec.dot(pDxvec_pDrotx))/svec.length()

        """
        #for finite differences
        #Rx_deltax = col((1,0,0)).axis_and_angle_as_r3_rotation_matrix(rotx+delta_rot)
        #a_star_deltax = Ry * Rx_deltax * Ai
        crystal_deltax_orientation = crystal_init_orientation.rotate_thru((1,0,0),rotx+delta_rot
             ).rotate_thru((0,1,0),roty
             ).rotate_thru((0,0,1),rotz)
        a_star_deltax = sqr(crystal_deltax_orientation.reciprocal_matrix())
        p_deltax, rh_deltax = ph.calc_partiality(a_star_deltax, miller_index)
        xvec_deltax = a_star_deltax * hvec
        svec_deltax = s0 + xvec_deltax
        pDxvec_pDrotx_fd = (xvec_deltax - xvec)/delta_rot
        pDrh_pDrotx_fd = (rh_deltax - rh)/delta_rot
        if j==1:
          print pDrh_pDrotx, pDrh_pDrotx_fd, pDrh_pDrotx-pDrh_pDrotx_fd
        """

        #for roty
        pDxvec_pDroty = dA_droty * hvec
        pDrh_pDroty = (svec.dot(pDxvec_pDroty))/svec.length()


        """
        #for finite differences
        crystal_deltay_orientation = crystal_init_orientation.rotate_thru((1,0,0),rotx
             ).rotate_thru((0,1,0),roty+delta_rot
             ).rotate_thru((0,0,1),rotz)
        a_star_deltay = sqr(crystal_deltay_orientation.reciprocal_matrix())
        p_deltay, rh_deltay = ph.calc_partiality(a_star_deltay, miller_index)
        xvec_deltay = a_star_deltay * hvec
        svec_deltay = s0 + xvec_deltay
        pDxvec_pDroty_fd = (xvec_deltay - xvec)/delta_rot
        pDrh_pDroty_fd = (rh_deltay - rh)/delta_rot
        if j==1:
          print pDrh_pDroty, pDrh_pDroty_fd, pDrh_pDroty-pDrh_pDroty_fd
        """

        #for rotz
        pDxvec_pDrotz = dA_drotz * hvec
        pDrh_pDrotz = (svec.dot(pDxvec_pDrotz))/svec.length()



        #for finite differences
        crystal_deltaz_orientation = crystal_init_orientation.rotate_thru((1,0,0),rotx
             ).rotate_thru((0,1,0),roty
             ).rotate_thru((0,0,1),rotz+delta_rot)
        a_star_deltaz = sqr(crystal_deltaz_orientation.reciprocal_matrix())
        p_deltaz, rh_deltaz = ph.calc_partiality(a_star_deltaz, miller_index)
        xvec_deltaz = a_star_deltaz * hvec
        svec_deltaz = s0 + xvec_deltaz
        pDxvec_pDrotz = (xvec_deltaz - xvec)/delta_rot
        pDrh_pDrotz = (rh_deltaz - rh)/delta_rot


        #shared derivatives on rotx and roty
        pDP_pDrh = (-4*(rs**2)*rh)/((2*(rh**2))+(rs**2))**2
        pDg_pDP = (-1 * G  * I_obs[j])/(sigI_obs[j]*(p**2))

        #for G
        pDg_pDG = I_obs[j] / (sigI_obs[j] * p)

        #for rs
        pDP_pDrs = (4*rs*(rh**2))/(((2*(rh**2))+(rs**2))**2)

        pDJ = ((((G * I_obs[j])/p) - I_ref[j])/sigI_obs[j])/m

        pDJ_pDrotx += pDJ * pDg_pDP * pDP_pDrh * pDrh_pDrotx
        pDJ_pDroty += pDJ * pDg_pDP * pDP_pDrh * pDrh_pDroty
        pDJ_pDrotz += pDJ * pDg_pDP * pDP_pDrh * pDrh_pDrotz
        pDJ_pDG += pDJ * pDg_pDG
        pDJ_pDrs += pDJ * pDg_pDP * pDP_pDrs



      """
      #finite differences for PDJ
      #Calculate finite differences of J function over rotation angle x and y
      J_at_rotx_and_delta = compute_cost(I_ref, I_obs, crystal_init_orientation, miller_indices, ph, (k, G, rotx+delta_rot, roty))
      J_at_roty_and_delta = compute_cost(I_ref, I_obs, crystal_init_orientation, miller_indices, ph, (k, G, rotx, roty+delta_rot))
      J_at_rot = compute_cost(I_ref, I_obs, crystal_init_orientation, miller_indices, ph, (k, G, rotx, roty))
      pDJ_pDrotx_fd = (J_at_rotx_and_delta - J_at_rot)/delta_rot
      pDJ_pDroty_fd = (J_at_roty_and_delta - J_at_rot)/delta_rot

      print i, pDJ_pDrotx, pDJ_pDrotx_fd, pDJ_pDroty, pDJ_pDroty_fd
      """

      #update parameters
      rotx = rotx - (alpha*pDJ_pDrotx)
      roty = roty - (alpha*pDJ_pDroty)
      rotz = rotz - (alpha*pDJ_pDrotz)
      G = G - (alpha*pDJ_pDG)
      rs = rs - (alpha*pDJ_pDrs)

      J_history[i] = self.compute_cost(I_ref, I_obs, sigI_obs, crystal_init_orientation, miller_indices_ori, wavelength, (G, rotx, roty, rotz, rs))

      print i, G, rotx*180/math.pi, roty*180/math.pi, rotz*180/math.pi, rs

      #check convergence
      gradient_threshold = 0.05
      if i > 0:
        if abs((J_history[i]-J_history[i-1])/delta_rot) < gradient_threshold:
          converge_iter = i
          break

    return G, rotx, roty, rotz, rs
