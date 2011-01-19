from cctbx.array_family import flex
from libtbx import adopt_init_args
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit

class manager(object):
  def __init__(self,
               fmodel):
    adopt_init_args(self, locals())

  def tst_call(self):
    print "TST CALL FROM TA_UTILS"

  #Takes array of x,y where x is d*^3 and fits smoothing chebyshev function
  #over all resolution range, returning a miller array (based on sigmaa est)

  def chebyshev_smoothing(self,
                          x_data,
                          y_data,
                          x_miller,
                          n_chebyshev_terms = 10,
                          ):
#    reparam_y = -flex.log( 1.0/y_data -1.0 )

    max_y = max(y_data)
    reparam_y = -flex.log(y_data)
    fit_lsq = chebyshev_lsq_fit.chebyshev_lsq_fit(
      n_terms = n_chebyshev_terms,
      x_obs   = x_data,
      y_obs   = reparam_y,
      w_obs   = None)

    min_h = flex.min(x_miller)
    max_h = flex.max(x_miller)

    cheb_pol = chebyshev_polynome(
      n_chebyshev_terms,
      min_h,
      max_h,
      fit_lsq.coefs)

    def reverse_reparam(values):
      return flex.exp(-values)
#      return 1.0/(1.0 + flex.exp(-values))
    y_fitted = reverse_reparam(cheb_pol.f(x_data))
    fitted_miller_array = reverse_reparam(cheb_pol.f(x_miller))
    fitted_miller_array = self.fmodel.f_obs.array(data=fitted_miller_array)
    return fitted_miller_array
