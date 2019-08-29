from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.uctbx import unit_cell

class postref_results(object):
  '''
  Author      : Uervirojnangkoorn, M.
  Created     : 7/13/2014
  A wrapper class to store refinement result
  '''

  def __init__(self):
    '''
    Do nothing
    '''

  def set_params(self, observations=None,
      observations_original=None,
      refined_params=None,
      stats=None,
      partiality=None,
      rs_set=None,
      rh_set=None,
      frame_no=None,
      pickle_filename=None,
      wavelength=None,
      crystal_orientation=None,
      detector_distance_mm=None):

    self.observations = observations
    self.observations_original = observations_original
    self.refined_params = refined_params
    self.partiality = partiality
    self.rs_set = rs_set
    self.rh_set = rh_set
    self.frame_no = frame_no
    self.pickle_filename = pickle_filename
    self.wavelength = wavelength
    self.detector_distance_mm = detector_distance_mm

    #refined_params
    #note params = G,B,rotx,roty,ry,a,b,c,alpha,beta,gamma
    self.G = refined_params[0]
    self.B = refined_params[1]
    self.rotx = refined_params[2]
    self.roty = refined_params[3]
    self.ry = refined_params[4]
    self.rz = refined_params[5]
    self.r0 = refined_params[6]
    self.re = refined_params[7]
    self.voigt_nu = refined_params[8]
    self.uc_params = flex.double([refined_params[9], refined_params[10], refined_params[11],
          refined_params[12], refined_params[13], refined_params[14]])
    self.unit_cell = unit_cell((refined_params[9], refined_params[10], refined_params[11],
          refined_params[12], refined_params[13], refined_params[14]))
    self.crystal_orientation = crystal_orientation

    #SE, R_sq, CC_init, CC_final, R_init, R_final, R_xy_init, R_xy_final
    self.SE = stats[0]
    self.R_sq = stats[1]
    self.CC_init = stats[2]
    self.CC_final = stats[3]
    self.R_init = stats[4]
    self.R_final = stats[5]
    self.R_xy_init = stats[6]
    self.R_xy_final = stats[7]
    self.CC_iso_init = stats[8]
    self.CC_iso_final = stats[9]

  def get_full_observations(self):
    sin_theta_over_lambda_sq = self.observations.two_theta(wavelength=self.wavelength).sin_theta_over_lambda_sq().data()
    return self.observations.customized_copy( \
        data=self.observations.data()/(self.G * flex.exp(-2*self.B*sin_theta_over_lambda_sq) * self.partiality), \
        sigmas=self.observations.sigmas()/(self.G * flex.exp(-2*self.B*sin_theta_over_lambda_sq) * self.partiality))
