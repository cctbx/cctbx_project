from scitbx import lbfgs
from scitbx.array_family import flex
from scitbx.math import chebyshev_lsq
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_base
import math


class chebyshev_lsq_fit:
  def __init__(self,
               n_terms,
               x_obs,
               y_obs,
               w_obs=None,
               free_flags=None,
               low_limit=None,
               high_limit=None,
               randomise=False):
    self.n = n_terms
    self.x_obs = x_obs
    self.y_obs = y_obs
    self.free_flags = free_flags
    if self.free_flags is None:
      self.free_flags = flex.bool(x_obs.size(), True)
    self.w_obs = flex.double(x_obs.size(), 1.0)
    if w_obs is not None:
      self.w_obs = w_obs
    self.x = flex.double(self.n, 0)
    if randomise:
      self.x = (flex.random_double(self.n)-0.5)*10.0
    self.low_limit = flex.min(self.x_obs)
    self.high_limit = flex.max(self.x_obs)
    self.f = None
    if low_limit is not None:
      self.low_limit = low_limit
    if high_limit is not None:
      self.high_limit = high_limit

    ## Set the first term equal to twice mean of the data points.
    ## Although not really needed, seems like a good idea anyway.
    ## It should speed up convergence.
    self.x[0] = flex.mean(self.y_obs)*2.0
    self.lsq_object = chebyshev_lsq(self.n,
                                    self.low_limit,
                                    self.high_limit,
                                    self.x_obs,
                                    self.y_obs,
                                    self.w_obs,
                                    self.free_flags)
    self.lsq_object.replace(self.x)
    self.minimizer = lbfgs.run(target_evaluator=self)
    self.coefs = self.lsq_object.coefs()
    self.f = self.lsq_object.residual()
    self.free_f = self.lsq_object.free_residual()
    del self.x

  def __call__(self):
    self.lsq_object.replace(self.x)
    f = self.lsq_object.residual()
    g = self.lsq_object.gradient()
    self.f = f
    return self.x, f ,g


def cross_validate_to_determine_number_of_terms(x_obs,y_obs,w_obs=None,
                                                min_terms=10,max_terms=25,
                                                n_free=100, n_goes = 5):
  if (n_goes==None):
    if (min_terms<2):
      min_terms=2

    free_residuals = []

    free_flags = flex.bool(x_obs.size(),True)
    free_permut = flex.random_permutation(x_obs.size())
    for ii in range(n_free):
      free_flags[free_permut[ii]]=False

    for count in range(min_terms,max_terms):
      fit = chebyshev_lsq_fit(count,x_obs,y_obs,w_obs,free_flags)
      free_residuals.append(fit.free_f)
    return(flex.double(free_residuals))

  else :
    if w_obs is None:
      w_obs = flex.double(x_obs.size(),1)
    free_resid = flex.double(max_terms-min_terms,0)
    for jj in range(n_goes):
      free_resid += cross_validate_to_determine_number_of_terms(
          x_obs,y_obs,w_obs,
          min_terms=min_terms,max_terms=max_terms,
          n_free=n_free, n_goes = None)
    return( min_terms + flex.min_index( free_resid ) )
