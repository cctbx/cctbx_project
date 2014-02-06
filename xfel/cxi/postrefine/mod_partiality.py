import math
from scitbx.matrix import col
from cctbx.array_family import flex

class partiality_handler(object):
  '''
  calculate partiality based on off-the-Ewald sphere error
  '''

  def __init__(self, wavelength, spot_radius):
    '''
    Constructor
    '''
    self.reciprocal_wavelength = 1/wavelength
    self.S0 = -1*col((0,0,self.reciprocal_wavelength)) #NOTE: S0 has negative sign!!!
    self.spot_radius = spot_radius
   
   
  def calc_partiality(self, a_star_matrix, miller_index):
    
    h = col(miller_index)
      
    #get the pickle S and calculate delta offset from Ewald sphere
    x = a_star_matrix * h
    S = x + self.S0
    rh = S.length() - self.reciprocal_wavelength
      
    #calculate partiality
    spot_partiality = (pow(self.spot_radius,2)/((2*pow(rh,2))+pow(self.spot_radius,2)))
    
    
    return spot_partiality, rh
  
  def calc_partiality_anisotropy(self, a_star_matrix, miller_index, ry, rz, alpha_angle):
    
    #use III.4 in Winkler et al 1979 (A35; P901)
    rs = math.sqrt((ry * math.cos(alpha_angle))**2 + (rz * math.sin(alpha_angle))**2)
    h = col(miller_index)
      
    #get the pickle S and calculate delta offset from Ewald sphere
    x = a_star_matrix * h
    S = x + self.S0
    rh = S.length() - self.reciprocal_wavelength
      
    #calculate partiality
    spot_partiality = pow(rs,2)/((2*pow(rh,2))+pow(rs,2))
    
    return spot_partiality, rh, rs
  
  def calc_partiality_anisotropy_set(self, a_star_matrix, miller_indices, ry, rz, alpha_angle_set):
    
    #use III.4 in Winkler et al 1979 (A35; P901)
    partiality_set = flex.double()
    for miller_index, alpha_angle in zip(miller_indices, alpha_angle_set):
      rs = math.sqrt((ry * math.cos(alpha_angle))**2 + (rz * math.sin(alpha_angle))**2)
      h = col(miller_index)
        
      #get the pickle S and calculate delta offset from Ewald sphere
      x = a_star_matrix * h
      S = x + self.S0
      rh = S.length() - self.reciprocal_wavelength
        
      #calculate partiality
      spot_partiality = pow(rs,2)/((2*pow(rh,2))+pow(rs,2))
      
      partiality_set.append(spot_partiality)
    
    return partiality_set
