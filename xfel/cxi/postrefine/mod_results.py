from __future__ import division
from cctbx.array_family import flex

class postref_results(object):
  '''
  A wrapper class to store refinement result
  '''

  def __init__(self):
    '''
    Do nothing
    '''

  def set_params(self, observations=None,
      refined_params=None,
      se_params=None,
      stats=None,
      partiality=None,
      frame_no=None,
      pickle_filename=None,
      wavelength=None,
      SE_I=None,
      var_I_p=None,
      var_k=None,
      var_p=None):
    
    self.observations = observations
    self.refined_params = refined_params
    self.se_params = se_params
    self.stats = stats
    self.partiality = partiality
    self.frame_no = frame_no
    self.pickle_filename = pickle_filename
    self.wavelength = wavelength
    self.SE_I = SE_I 
    self.var_I_p = var_I_p
    self.var_k = var_k
    self.var_p = var_p
    
    #refined_params
    #note params = G,B,rotx,roty,ry,a,b,c,alpha,beta,gamma
    self.G = refined_params[0]
    self.B = refined_params[1]
    self.rotx = refined_params[2]
    self.roty = refined_params[3]
    self.ry = refined_params[4]
    self.rz = refined_params[5]
    self.re = refined_params[6]
    self.uc_params = flex.double([refined_params[7], refined_params[8], refined_params[9], 
          refined_params[10], refined_params[11], refined_params[12]])
          
    self.se_G = se_params[0]
    self.se_B = se_params[1]
    self.se_rotx = se_params[2]
    self.se_roty = se_params[3]
    self.se_ry = se_params[4]
    self.se_rz = se_params[5]
    self.se_re = se_params[6]
    self.se_uc_params = flex.double([se_params[7], se_params[8], se_params[9], 
          se_params[10], se_params[11], se_params[12]])
