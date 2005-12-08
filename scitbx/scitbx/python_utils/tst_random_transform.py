import random
import os
import time
import math
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from scitbx.python_utils import random_transform as rt

def exersize_gauss():
  data = rt.normal_variate(mu=0,sigma=1,N=1000000)
  mu1 = flex.mean(data)
  mu2 = flex.mean(data*data)
  mu3 = flex.mean(data*data*data)
  assert approx_equal(mu1,0,eps=0.02)
  assert approx_equal(mu2,1,eps=0.02)
  assert approx_equal(mu3,0,eps=0.02)

def exersize_t_variate():
  data = rt.t_variate(a=6, mu=0,sigma=1,N=1000000)
  mu1 = flex.mean(data)
  mu2 = flex.mean(data*data)
  assert approx_equal(mu1,0,eps=0.02)
  assert approx_equal(mu2,1.5,eps=0.04)

def exersize_wilson_amplitude_variate():
  data = rt.wilson_amplitude_variate(N=1000000)
  mu1 = flex.mean(data)
  mu2 = flex.mean(data*data)
  assert approx_equal(mu1,0.886,eps=0.02)
  assert approx_equal(mu2,1.000,eps=0.02)

def exersize_wilson_intensity_variate():
  data = rt.wilson_intensity_variate(N=1000000)
  mu1 = flex.mean(data)
  mu2 = flex.mean(data*data)
  assert approx_equal(mu1,1.000,eps=0.02)
  assert approx_equal(mu2,2.000,eps=0.04)

def exersize_pseudo_normalized_abs_delta_i():
  data = rt.pseudo_normalized_abs_delta_i(N=1000000)
  mu1 = flex.mean( data )
  assert approx_equal(mu1,0.5,eps=0.02)


def run():
  exersize_gauss()
  exersize_t_variate()
  exersize_wilson_amplitude_variate()
  exersize_wilson_intensity_variate()
  exersize_pseudo_normalized_abs_delta_i()

  print 'OK'

if (__name__ == "__main__"):
  run()
