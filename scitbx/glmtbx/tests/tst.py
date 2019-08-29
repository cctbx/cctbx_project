
from __future__ import absolute_import, division, print_function
from six.moves import range

def tst_glm():
  from scitbx.glmtbx import glm
  from scitbx.array_family import flex
  from numpy.random import poisson, seed
  from math import exp

  seed(0)

  n_obs = 100

  for c in range(1, 100):

    # Test for a constant value
    X = flex.double([1 for i in range(n_obs)])
    X.reshape(flex.grid(n_obs, 1))
    Y = flex.double(list(poisson(c, n_obs)))
    B = flex.double([0])
    P = flex.double([1 for i in range(n_obs)])
    result = glm(X, Y, B, P, family="poisson", max_iter=100)
    assert(abs(c - exp(result.parameters()[0])) < 0.1*c)


  print('OK')

def tst_robust_glm():
  from scitbx.glmtbx import robust_glm
  from scitbx.array_family import flex
  from numpy.random import poisson, seed
  from math import exp

  seed(0)

  n_obs = 100

  for c in range(1, 100):

    # Test for a constant value
    X = flex.double([1 for i in range(n_obs)])
    X.reshape(flex.grid(n_obs, 1))
    Y = flex.double(list(poisson(c, n_obs)))
    B = flex.double([0])
    result = robust_glm(X, Y, B, family="poisson", max_iter=100)
    assert(abs(c - exp(result.parameters()[0])) < 0.1*c)

  # Now test with a massive outlier
  for c in range(1, 100):

    # Test for a constant value
    X = flex.double([1 for i in range(n_obs)])
    X.reshape(flex.grid(n_obs, 1))
    Y = flex.double(list(poisson(c, n_obs)))
    Y[n_obs//2] = c * 100
    B = flex.double([0])
    result = robust_glm(X, Y, B, family="poisson", max_iter=100)
    assert(abs(c - exp(result.parameters()[0])) < 0.1*c)

  print('OK')

def run():
  tst_glm()
  tst_robust_glm()

if __name__ == '__main__':

  run()
