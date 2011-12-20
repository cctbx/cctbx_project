from stdlib import random

from libtbx.test_utils import approx_equal
from libtbx.utils import frange
from scitbx.array_family import flex
from scitbx.math import curve_fitting
from scitbx import matrix
import scitbx.lbfgs

if (1): # fixed random seed to avoid rare failures
  random.seed(0)
  flex.set_random_seed(0)

def finite_differences(functor, x, eps=1e-4):
  grads = []
  for i in range(len(functor.params)):
    params = flex.double(functor.params)
    params[i] += eps
    f = functor.__class__(*params)
    qm = matrix.col(f(x))
    params[i] -= 2 * eps
    f = functor.__class__(*params)
    qp = matrix.col(f(x))
    dq = (qm-qp)/(2*eps)
    grads.append(list(dq))
  return grads

def exercise_polynomial_fit():

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
    functor = curve_fitting.univariate_polynomial(*params)
    fd_grads = finite_differences(functor, x)
    assert approx_equal(functor.partial_derivatives(x), fd_grads, 1e-4)
    do_polynomial_fit(x, params)

def exercise_gaussian_fit():

  # test fitting of a gaussian
  def do_gaussian_fit(scale, mu, sigma):
    start = mu - 6 * sigma
    stop = mu + 6 * sigma
    step = (stop - start)/1000
    x = flex.double(frange(start, stop, step))
    y = scale * flex.exp(-flex.pow2(x - mu) / (2 * sigma**2))
    fit = curve_fitting.single_gaussian_fit(x, y)
    assert approx_equal(fit.a, scale, 1e-4)
    assert approx_equal(fit.b, mu, eps=1e-4)
    assert approx_equal(fit.c, sigma, eps=1e-4)

  for i in range(10):
    scale = random.random() * 1000
    sigma = (random.random() + 0.0001) * 10
    mu = (-1)**random.randint(0,1) * random.random() * 1000
    functor = curve_fitting.gaussian(scale, mu, sigma)
    start = mu - 6 * sigma
    stop = mu + 6 * sigma
    step = (stop - start)/1000
    x = flex.double(frange(start, stop, step))
    fd_grads = finite_differences(functor, x)
    assert approx_equal(functor.partial_derivatives(x), fd_grads, 1e-4)
    do_gaussian_fit(scale, mu, sigma)

  # if we take the log of a gaussian we can fit a parabola
  scale = 123
  mu = 3.2
  sigma = 0.1
  x = flex.double(frange(2, 4, 0.01))
  y = scale * flex.exp(-flex.pow2(x - mu) / (2 * sigma**2))
  # need to be careful to only use values of y > 0
  eps = 1e-15
  x = flex.double([x[i] for i in range(x.size()) if y[i] > eps])
  y = flex.double([y[i] for i in range(y.size()) if y[i] > eps])
  fit = curve_fitting.univariate_polynomial_fit(x, flex.log(y), degree=2)
  c, b, a = fit.params
  assert approx_equal(mu, -b/(2*a))
  assert approx_equal(sigma*sigma, -1/(2*a))

  # test multiple gaussian fits
  gaussians = [curve_fitting.gaussian(0.3989538, 3.7499764, 0.7500268),
               curve_fitting.gaussian(0.7978957, 6.0000004, 0.5000078)]
  x = flex.double(frange(0, 10, 0.1))
  y = flex.double(x.size())
  for i in range(len(gaussians)):
    g = gaussians[i]
    scale, mu, sigma = g.a, g.b, g.c
    y += g(x)

  starting_gaussians = [
    curve_fitting.gaussian(1, 4, 1),
    curve_fitting.gaussian(1, 5, 1)]
  fit = curve_fitting.gaussian_fit(x, y, starting_gaussians)
  for g1, g2 in zip(gaussians, fit.gaussians):
    assert approx_equal(g1.a, g2.a, eps=1e-4)
    assert approx_equal(g1.b, g2.b, eps=1e-4)
    assert approx_equal(g1.c, g2.c, eps=1e-4)

  # use example of 5-gaussian fit from here:
  # http://research.stowers-institute.org/efg/R/Statistics/MixturesOfDistributions/index.htm
  gaussians = [curve_fitting.gaussian(0.10516252, 23.32727, 2.436638),
               curve_fitting.gaussian(0.46462715, 33.09053, 2.997594),
               curve_fitting.gaussian(0.29827916, 41.27244, 4.274585),
               curve_fitting.gaussian(0.08986616, 51.24468, 5.077521),
               curve_fitting.gaussian(0.04206501, 61.31818, 7.070303)]

  x = flex.double(frange(0, 80, 0.1))
  y = flex.double(x.size())
  for i in range(len(gaussians)):
    g = gaussians[i]
    scale, mu, sigma = g.a, g.b, g.c
    y += g(x)

  termination_params = scitbx.lbfgs.termination_parameters(
    min_iterations=500)
  starting_gaussians = [curve_fitting.gaussian(1, 21, 2.1),
                        curve_fitting.gaussian(1, 30, 2.8),
                        curve_fitting.gaussian(1, 40, 2.2),
                        curve_fitting.gaussian(1, 51, 1.2),
                        curve_fitting.gaussian(1, 60, 2.3)]
  fit = curve_fitting.gaussian_fit(
    x, y, starting_gaussians, termination_params=termination_params)
  y_calc = fit.compute_y_calc()
  assert approx_equal(y, y_calc, eps=1e-2)


def exercise_skew_normal_fit():
  shape, location, scale = 8.0, 4.0, 2.0
  x_obs = flex.double(frange(0, 10, 0.1))
  f = curve_fitting.skew_normal(shape, location, scale)
  y_obs = f(x_obs)

  if 0:
    from matplotlib import pyplot
    pyplot.plot(x_obs, y_obs)
    pyplot.show()

  fd_grads = finite_differences(f, x_obs)
  ana_grads = f.partial_derivatives(x_obs)
  assert approx_equal(ana_grads, fd_grads, 1e-4)
  shape, location, scale = 1, 0, 8
  starting_f = curve_fitting.skew_normal(shape, location, scale)
  termination_params = scitbx.lbfgs.termination_parameters(
    min_iterations=100)
  fit = curve_fitting.generic_minimiser(
    [starting_f], x_obs, y_obs, termination_params=termination_params)
  assert approx_equal(fit.functions[0].params, f.params)



def run():
  exercise_skew_normal_fit()
  exercise_polynomial_fit()
  exercise_gaussian_fit()

if (__name__ == "__main__"):
  run()
  print "OK"
