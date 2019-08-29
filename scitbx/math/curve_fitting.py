from __future__ import absolute_import, division, print_function
import math

import libtbx
import libtbx.load_env
from libtbx import adopt_init_args
from scitbx.array_family import flex
import scitbx.lbfgs
import scitbx.math
from scitbx import matrix
from six.moves import range


class function_base(object):

  def __call__(self, x_obs):
    raise NotImplementedError

  def partial_derivatives(self, x_obs):
    """ This default implementation returns the finite difference partial
    derivatives. Override this function to calculate the derivatives
    analytically if required.
    """
    return [flex.double(g) for g in self.finite_differences(x_obs)]

  def finite_differences(self, x_obs, eps=1e-4):
    grads = []
    for i in range(len(self.params)):
      params = flex.double(self.params)
      params[i] += eps
      f = self.__class__(*params)
      qm = matrix.col(f(x_obs))
      params[i] -= 2 * eps
      f = self.__class__(*params)
      qp = matrix.col(f(x_obs))
      dq = (qm-qp)/(2*eps)
      grads.append(list(dq))
    return grads


class univariate_polynomial(function_base):

  def __init__(self, *params):
    """A polynomial of degree n:
         f(x) = a[0] + a[1] x**1 + ... * a[n] x**n.
    """
    self.params = params
    self.n_terms = len(params)
    self.degree = self.n_terms - 1

  def __call__(self, x_obs):
    y_calc = flex.double(x_obs.size())
    for n in range(self.n_terms):
      y_calc += self.params[n] * flex.pow(x_obs, n)
    return y_calc

  def partial_derivatives(self, x_obs):
    g = []
    for n in range(self.n_terms):
      g.append(flex.pow(x_obs, n))
    return g


class gaussian(function_base):

  def __init__(self, a, b, c):
    """Simple wrapper for the parameters associated with a gaussian
         f(x) = a * exp(-(x - b)^2 / (2 * c^2))
    """
    adopt_init_args(self, locals())
    self.params = (a, b, c)

  def __call__(self, x_obs):
    a, b, c = self.params
    y_calc = a * flex.exp(-flex.pow2(x_obs - b) / (2 * c**2))
    return y_calc

  def partial_derivatives(self, x_obs):
    a, b, c = self.params
    exponential_part = flex.exp(-flex.pow2(x_obs - b) / (2 * c**2))
    return(exponential_part,
           a * (x_obs - b) / c**2 * exponential_part,
           a * flex.pow2(x_obs - b) / c**3 * exponential_part)

  @property
  def sigma(self):
    return abs(self.params[2])

class skew_normal(function_base):

  def __init__(self, shape, location, scale):
    adopt_init_args(self, locals())
    self.params = (shape, location, scale)

  def __call__(self, x_obs):
    shape, location, scale = self.params
    normal_part = (2 / (scale * math.sqrt(2 * math.pi))
                   * flex.exp(- flex.pow2(x_obs - location)/(2 * scale**2)))
    cdf_part = 0.5 * (
      1 + scitbx.math.erf(shape * (x_obs - location)/ (math.sqrt(2) * scale)))
    y_calc = normal_part * cdf_part
    return y_calc

  def partial_derivatives(self, x_obs):
    shape, location, scale = self.params
    exponential_part = (1/(math.sqrt(2 * math.pi))
                        * flex.exp(- flex.pow2(x_obs - location)/(2 * scale**2)))
    normal_part = 2 / scale * exponential_part
    cdf_part = 0.5 * (
      1 + scitbx.math.erf(shape * (x_obs - location)/ (math.sqrt(2) * scale)))
    d_normal_part_d_location = 2 / scale**3 * (x_obs - location) * exponential_part
    d_normal_part_d_scale = \
      2 / scale**4 * (flex.pow2(x_obs - location) - scale**2) * exponential_part
    exponential_part_with_shape = (
      1 / (math.sqrt(math.pi)) *
      flex.exp(-shape**2 * flex.pow2(x_obs - location)/(2 * scale**2)))
    d_cdf_d_shape = \
      (x_obs - location) / (math.sqrt(2) * scale) * exponential_part_with_shape
    d_cdf_d_location = \
      -shape / (math.sqrt(2) * scale) * exponential_part_with_shape
    d_cdf_d_scale = (-shape * (x_obs - location) * exponential_part_with_shape /
                     (math.sqrt(2) * scale**2))
    # product rule
    return (d_cdf_d_shape * normal_part,
            d_normal_part_d_location * cdf_part + d_cdf_d_location * normal_part,
            d_normal_part_d_scale * cdf_part + d_cdf_d_scale * normal_part)

class tanh(function_base):

  def __init__(self, *params):
    """
  Curve fitting as suggested by Ed Pozharski to a tanh function
  of the form (1/2)(1 - tanh(z)) where z = (s - s0)/r,
  s0 is the value of s at the half-falloff value, and r controls the
  steepness of falloff
    """
    self.params = params

  def __call__(self, x_obs):
    s = x_obs
    r, s0 = self.params
    z = (s - s0)/r
    return 0.5 * (1 - flex.tanh(z))

class tanh_fit(object):

  def __init__(self, x_obs, y_obs, r=1, s0=1,
        min_iterations=0,
        max_iterations=None):
    """Curve fitting as suggested by Ed Pozharski to a tanh function
       of the form (1/2)(1 - tanh(z)) where z = (s - s0)/r,
       s0 is the value of s at the half-falloff value, and r controls the
       steepness of falloff

       :param x_obs: x-coordinates of the data
       :type x_obs: flex.double
       :param y_obs: y-coordinates of the data
       :type y_obs: flex.double
       :param s0: s0 is the value of s at the half-falloff value
       :type s0: float
       :param r: r controls the steepness of falloff
       :type r: float
    """
    self.x_obs = x_obs
    self.y_obs = y_obs
    assert r >= 0
    f = tanh(r, s0)
    fit = lbfgs_minimiser(
      functions=[f],
      x_obs=x_obs,
      y_obs=self.y_obs,
      termination_params=scitbx.lbfgs.termination_parameters(
        min_iterations=min_iterations,
        max_iterations=max_iterations))
    self.params = fit.functions[0].params

class univariate_polynomial_fit(object):

  def __init__(self, x_obs, y_obs, degree,
        min_iterations=0,
        max_iterations=None,
        number_of_cycles=1):
    """Fit a polynomial of degree n to points (x_obs, y_obs)
         f(x) = a[0] + a[1] x**1 + ... * a[n] x**n.

       :param x_obs: x-coordinates of the data
       :type x_obs: flex.double
       :param y_obs: y-coordinates of the data
       :type y_obs: flex.double
       :param degree: the degree of the polynomial - the largest power of x
       :type degree: int
    """
    self.x_obs = x_obs
    self.y_obs = y_obs
    assert isinstance(degree, int)
    assert degree >= 0
    self.degree = degree
    self.n_terms = degree + 1
    params = flex.double([1] * self.n_terms)
    for cycle in range(number_of_cycles):
      polynomial = univariate_polynomial(*params)
      fit = lbfgs_minimiser(
        functions=[polynomial],
        x_obs=x_obs,
        y_obs=self.y_obs,
        termination_params=scitbx.lbfgs.termination_parameters(
          min_iterations=min_iterations,
          max_iterations=max_iterations))
      self.params = fit.functions[0].params
      params = self.params


class single_gaussian_fit(object):

  def __init__(self, x_obs, y_obs):
    """Fit a gaussian to points (x_obs, y_obs):
         f(x) = A exp(-(x - mu)**2 / (2 * sigma**2))

       :param x_obs: x-coordinates of the data
       :type x_obs: flex.double
       :param y_obs: y-coordinates of the data
       :type y_obs: flex.double
    """
    self.x_obs = x_obs
    self.y_obs = y_obs
    max_i = flex.max_index(y_obs)
    # quick estimate of scale and mean to give the optimiser a helping hand
    scale = y_obs[max_i]
    mu = x_obs[max_i]
    sigma = 1 # can we make a simple estimate of sigma too?
    fit = gaussian_fit(x_obs, y_obs, [gaussian(scale, mu, sigma)])
    self.a = fit.gaussians[0].a
    self.b = fit.gaussians[0].b
    self.c = fit.gaussians[0].c


class gaussian_fit(object):

  def __init__(self, x_obs, y_obs, starting_gaussians, termination_params=None):
    """Fit one or more gaussians to points (x_obs, y_obs):
         f(x) = sum_i(A_i exp(-(x - mu_i)**2 / (2 * sigma_i**2)))

       :param x_obs: x-coordinates of the data
       :type x_obs: flex.double
       :param y_obs: y-coordinates of the data
       :type y_obs: flex.double
       :param gaussian: a list or tuple of gaussian objects
       :type gaussian: list
    """
    self.n_cycles = 0
    self.x_obs = x_obs
    self.y_obs = y_obs
    self.n_gaussians = len(starting_gaussians)
    assert self.n_gaussians > 0
    fit = lbfgs_minimiser(
      functions=starting_gaussians, x_obs=x_obs, y_obs=self.y_obs)
    self.gaussians = fit.functions

  def compute_y_calc(self):
    y_calc = flex.double(self.x_obs.size())
    for i in range(self.n_gaussians):
      y_calc += self.gaussians[i](self.x_obs)
    return y_calc

  def pyplot(self):
    from matplotlib import pyplot
    pyplot.plot(self.x_obs, self.y_obs)
    pyplot.plot(self.x_obs, self.compute_y_calc())
    for i in range(self.n_gaussians):
      scale, mu, S = tuple(self.x[i*3:i*3+3])
      y_calc = scale * flex.exp(-flex.pow2(self.x_obs-mu) * S**2)
      pyplot.plot(self.x_obs, y_calc)
    pyplot.show()

class minimiser_base(object):

  def __init__(self, functions, x_obs, y_obs):
    self.n_cycles = 0
    self.x_obs = x_obs
    self.y_obs = y_obs
    self.n_functions = len(functions)
    self.functions = functions

  def compute_functional(self, params):
    self.x = params
    y_calc = self.compute_y_calc()
    delta_y = self.y_obs - y_calc
    f = flex.sum(flex.pow2(delta_y))
    return f

  def compute_y_calc(self):
    y_calc = flex.double(self.x_obs.size())
    for f in self.functions:
      y_calc += f(self.x_obs)
    return y_calc

  def compute_functional_and_gradients(self):
    y_calc = self.compute_y_calc()
    delta_y = self.y_obs - y_calc
    f = flex.sum(flex.pow2(delta_y))
    g = flex.double()
    for funct in self.functions:
      partial_ders = funct.partial_derivatives(self.x_obs)
      for i, partial in enumerate(partial_ders):
        g.append(-2 * flex.sum(delta_y * partial))
    return f, g

  def callback_after_step(self, minimizer):
    self.n_cycles += 1

  def pyplot(self):
    from matplotlib import pyplot
    pyplot.plot(self.x_obs, self.y_obs)
    pyplot.plot(self.x_obs, self.compute_y_calc())
    for f in self.functions:
      y_calc = f(self.x_obs)
      pyplot.plot(self.x_obs, y_calc)
    pyplot.show()

  @property
  def functions(self):
    x = self.x.deep_copy()
    for i, f in enumerate(self._functions):
      f = self._functions[i]
      self._functions[i] = f.__class__(*x[:len(f.params)])
      x = x[len(f.params):]
    return self._functions

  @functions.setter
  def functions(self, functions):
    self._functions = functions
    x = []
    for f in self._functions:
      x.extend(f.params)
    self.x = flex.double(x)


class lbfgs_minimiser(minimiser_base):

  def __init__(self, functions, x_obs, y_obs, termination_params=None):
    super(lbfgs_minimiser, self).__init__(functions, x_obs, y_obs)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self, termination_params=termination_params)


have_cma_es = libtbx.env.has_module("cma_es")

if have_cma_es:
  from cma_es import cma_es_interface

  class cma_es_minimiser(minimiser_base):

    def __init__(self, functions, x_obs, y_obs):
      super(cma_es_minimiser, self).__init__(functions, x_obs, y_obs)
      sigma = flex.double(self.x.size(), 1)
      self.minimizer = cma_es_interface.cma_es_driver(len(self.x), self.x.deep_copy(), sigma, self.compute_functional)


generic_minimiser = lbfgs_minimiser # XXX backward compatibility 2012-02-07
