from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from scitbx.python_utils import random_transform as rt

def exercise_gauss():
  data = rt.normal_variate(mu=0,sigma=1,N=1000000)
  mu1 = flex.mean(data)
  mu2 = flex.mean(data*data)
  mu3 = flex.mean(data*data*data)
  assert approx_equal(mu1,0,eps=0.02)
  assert approx_equal(mu2,1,eps=0.02)
  assert approx_equal(mu3,0,eps=0.02)

def exercise_t_variate():
  data = rt.t_variate(a=6, mu=0,sigma=1,N=1000000)
  mu1 = flex.mean(data)
  mu2 = flex.mean(data*data)
  assert approx_equal(mu1,0,eps=0.02)
  assert approx_equal(mu2,1.5,eps=0.04)

def exercise_wilson_amplitude_variate():
  data = rt.wilson_amplitude_variate(N=1000000)
  mu1 = flex.mean(data)
  mu2 = flex.mean(data*data)
  assert approx_equal(mu1,0.886,eps=0.02)
  assert approx_equal(mu2,1.000,eps=0.02)

def exercise_wilson_intensity_variate():
  data = rt.wilson_intensity_variate(N=1000000)
  mu1 = flex.mean(data)
  mu2 = flex.mean(data*data)
  assert approx_equal(mu1,1.000,eps=0.02)
  assert approx_equal(mu2,2.000,eps=0.04)

def exercise_pseudo_normalized_abs_delta_i():
  data = rt.pseudo_normalized_abs_delta_i(N=1000000)
  mu1 = flex.mean( data )
  assert approx_equal(mu1,0.5,eps=0.02)

def exercise_poisson():
  import random
  random.seed(0) # failure rate is fairly high without fixed seed
  a = rt.poisson_variate(100000,1).as_double()
  m = flex.mean(a)
  v = flex.mean(a*a)-m*m
  assert approx_equal(m,1.0,eps=0.05)
  assert approx_equal(v,1.0,eps=0.05)

  a = rt.poisson_variate(100000,10).as_double()
  m = flex.mean(a)
  v = flex.mean(a*a)-m*m
  assert approx_equal(m,10.0,eps=0.05)
  assert approx_equal(v,10.0,eps=0.05)

def run():
  exercise_poisson()
  exercise_gauss()
  exercise_t_variate()
  exercise_wilson_amplitude_variate()
  exercise_wilson_intensity_variate()
  exercise_pseudo_normalized_abs_delta_i()
  print('OK')

if (__name__ == "__main__"):
  run()
