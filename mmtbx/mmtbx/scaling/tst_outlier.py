from mmtbx import scaling
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import math


def tst_outliers_find_posterior_mode():
  # first check if we can find the posterior mode for the acentric
  # when alpha is close to 1, this should be very close to fcalc
  fobs  =  flex.double( [3]*10 )
  fcalc =  flex.double( range(10) )*10 + 10
  sigmas = flex.double( [0]*10 )
  epsilon = flex.double( [1]*10 )
  centric = flex.bool( [False]*10 )
  beta = flex.double( [1]*10 ) 
  alpha = flex.double( [0.99]*10 )
  tmp_object = scaling.likelihood_ratio_outlier_test(
     fobs,
     sigmas,
     fcalc,
     epsilon,
     centric,
     alpha,
     beta)
  posterior_mode = tmp_object.posterior_mode()
  for f, m in zip(fcalc,posterior_mode):
    assert approx_equal(f/m, 1, eps=0.05)
  
  # have a look at centrics
  fobs  =  flex.double( [3]*10 )
  fcalc =  flex.double( range(10) )*100 + 100
  sigmas = flex.double( [0]*10 )
  epsilon = flex.double( [1]*10 )
  centric = flex.bool( [True]*10 )
  beta = flex.double( [1]*10 )
  alpha = flex.double( [0.099]*10 )
  tmp_object = scaling.likelihood_ratio_outlier_test(
     fobs,
     sigmas,
     fcalc,
     epsilon,
     centric,
     alpha,
     beta)
  posterior_mode = tmp_object.posterior_mode()
  for f, m in zip(fcalc,posterior_mode):
    assert approx_equal(m/f, 0.099, eps=0.001)

def tst_loglikelihoods():
  fobs  =  flex.double( range(1000) )/50
  fcalc =  flex.double( [10]*1000 )
  sigmas = flex.double( [0]*1000 )
  epsilon = flex.double( [1]*1000 )
  centric = flex.bool( [False]*1000 )
  beta = flex.double( [1]*1000 )
  alpha = flex.double( [0.99]*1000 )
  tmp_object = scaling.likelihood_ratio_outlier_test(
     fobs,
     sigmas,
     fcalc,
     epsilon,
     centric,
     alpha,
     beta)
  cur_lik = tmp_object.log_likelihood()
  pm_lik  = tmp_object.posterior_mode_log_likelihood()

  level=4.5
  flags   = tmp_object.flag_potential_outliers( level )
  
  for fl,pl,l in zip(flags,pm_lik,cur_lik):
    if (pl-l)<=level:
      assert fl
    else:
      assert not fl

def run():
  tst_outliers_find_posterior_mode()
  tst_loglikelihoods() 
  print "OK"

if (__name__ == "__main__"):
  run()


