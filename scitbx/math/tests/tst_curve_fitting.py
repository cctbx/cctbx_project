from stdlib import random

from libtbx.test_utils import approx_equal
from libtbx.utils import frange
from scitbx.array_family import flex
from scitbx.math import curve_fitting

def run():

  def do_polynomial_fit(x, params):
    n_terms = len(params)
    y = flex.double(x.size())
    for i in range(len(params)):
      y += params[i] * flex.pow(x, i)
    fit = curve_fitting.univariate_polynomial_fit(x, y, degree=n_terms-1)
    assert approx_equal(params, fit.params, eps=1e-4)

  x = flex.double(range(-50,50))
  do_polynomial_fit(x, (2,3,5)) # y = 2 + 3x + 5x^2
  do_polynomial_fit(x, (-0.0002, -1000)) # y = -0.0002 -1000x

  for n_terms in range(1, 6):
    params = [100*random.random() for i in range(n_terms)]
    x = flex.double(frange(-random.randint(1,10), random.randint(1,10), 0.1))
    do_polynomial_fit(x, params)

  # test fitting of a gaussian
  def do_gaussian_fit(scale, mu, sigma):
    start = mu - 6 * sigma
    stop = mu + 6 * sigma
    step = (stop - start)/1000
    x = flex.double(frange(start, stop, step))
    y = scale * flex.exp(-flex.pow2(x - mu) / (2 * sigma**2))
    fit = curve_fitting.single_gaussian_fit(x, y)
    assert approx_equal(fit.scale, scale, 1e-4)
    assert approx_equal(fit.mu, mu, eps=1e-4)
    assert approx_equal(fit.sigma, sigma, eps=1e-4)

  for i in range(10):
    scale = random.random() * 1000
    sigma = (random.random() + 0.0001) * 10
    mu = (-1)**random.randint(0,1) * random.random() * 1000
    do_gaussian_fit(scale, mu, sigma)

  scale = 123
  mu = 3.2
  sigma = 0.1
  # if we take the log of a gaussian we can fit a parabola
  y = scale * flex.exp(-flex.pow2(x - mu)/(2 * sigma**2))
  # need to be careful to only use values of y > 0
  eps = 1e-15
  x = flex.double([x[i] for i in range(x.size()) if y[i] > eps])
  y = flex.double([y[i] for i in range(y.size()) if y[i] > eps])
  fit = curve_fitting.univariate_polynomial_fit(x, flex.log(y), degree=2)
  c, b, a = fit.params
  assert approx_equal(mu, -b/(2*a))
  assert approx_equal(sigma*sigma, -1/(2*a))

if (__name__ == "__main__"):
  run()
  print "OK"
