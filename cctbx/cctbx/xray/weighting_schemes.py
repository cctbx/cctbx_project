from cctbx.array_family import flex

class unit_weighting(object):
  """ Mere unit weights for F and quasi-unit  weights 1/(4 Fo^2) for F^2. 
      
      For the last case, the weights is replaced by 1/(4 sigma(Fo^2)^2) 
      for weak Fo^2, which are defined as Fo^2 < n_sigma * sigma(Fo^2)
  """
  
  depends_only_on_obs = True
  
  def __init__(self, n_sigma=1.0):
    self.n_sigma = n_sigma
  
  def amplitude_weights(self, f_obs, f_calc):
    return None
    
  def intensity_weights(self, f_obs_square, f_calc):
    f_sqr = f_obs_square.data()
    sig_f_sqr = f_obs_square.sigmas()
    result = flex.double(f_sqr.size())
    if self.n_sigma == 1:
      strongs = f_sqr > sig_f_sqr
    else:
      strongs = f_sqr > self.n_sigma * sig_f_sqr
    weaks = ~strongs
    result.set_selected(strongs, 0.25/f_sqr.select(strongs))
    result.set_selected(weaks,   0.25/flex.pow2( sig_f_sqr.select(weaks) ))
    return result
  

class pure_statistical_weighting(object):
  """ 1/sigma^2 weights """
  
  depends_only_on_obs = True
  
  def amplitude_weights(self, f_obs, f_calc):
    return self.common_weights(y)
    
  def intensity_weights(self, f_obs_square, f_calc):
    return self.common_weights(y)
  
  def common_weights(self, y):
    sigmas_squared = flex.pow2(y.sigmas())
    assert sigmas_squared.all_gt(0)
    return 1 / sigmas_squared
  

