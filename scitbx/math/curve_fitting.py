import itertools
import math

from libtbx import adopt_init_args
import scitbx.lbfgs
from scitbx.array_family import flex

class univariate_polynomial_fit(object):

  def __init__(self, x_obs, y_obs, degree):
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
    self.x = flex.double([1] * self.n_terms)
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self)
    self.params = self.x
    del self.x

  def compute_functional_and_gradients(self):
    y_calc = flex.double(self.x_obs.size())
    for n in range(self.n_terms):
      y_calc += self.x[n] * flex.pow(self.x_obs, n)
    delta_y = self.y_obs - y_calc
    f = flex.sum(flex.pow2(delta_y))
    g = flex.double()
    for n in range(self.n_terms):
      g.append(-2 * flex.sum(delta_y * flex.pow(self.x_obs, n)))
    return f, g

  #def callback_after_step(self, minimizer):
    #print list(self.x)


class gaussian(object):

  def __init__(self, scale, mu, sigma):
    """Simple wrapper for the parameters associated with a gaussian
         f(x) = scale * exp(-(x - mu)**2 / (2 * sigma**2))
    """
    adopt_init_args(self, locals())

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
    self.scale = fit.gaussians[0].scale
    self.mu = fit.gaussians[0].mu
    self.sigma = fit.gaussians[0].sigma

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
    max_i = flex.max_index(y_obs)
    self.x = flex.double(list(itertools.chain(
      *[(g.scale, g.mu, math.sqrt(1/(2*g.sigma**2))) for g in starting_gaussians])))
    # run the minimizer
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self,
                                      termination_params=termination_params)
    # prepare the results
    self.gaussians = []
    for i in range(self.n_gaussians):
      scale, mu, S = self.x[i*3:i*3+3]
      sigma = math.sqrt(1/(2 * S**2))
      self.gaussians.append(gaussian(scale, mu, sigma))

  def compute_y_calc(self):
    y_calc = flex.double(self.x_obs.size())
    for i in range(self.n_gaussians):
      scale, mu, S = tuple(self.x[i*3:i*3+3])
      # reparametrisation:
      #   S**2 = 1/(2 * sigma**2)
      y_calc += scale * flex.exp(-flex.pow2(self.x_obs-mu) * S**2)
    return y_calc

  def compute_functional_and_gradients(self):
    y_calc = self.compute_y_calc()
    delta_y = self.y_obs - y_calc
    f = flex.sum(flex.pow2(delta_y))
    g = flex.double()
    for i in range(self.n_gaussians):
      scale, mu, S = tuple(self.x[i*3:i*3+3])
      exponential_part = flex.exp(-flex.pow2(self.x_obs-mu) * S**2)
      g.extend(flex.double((
        -2 * flex.sum(delta_y * exponential_part),
        -2 * flex.sum(delta_y * 2 * S**2 * scale * exponential_part * (self.x_obs - mu)),
        -2 * flex.sum(delta_y * - 2 * scale * S * exponential_part * flex.pow2(mu - self.x_obs))
      )))
    return f, g

  def callback_after_step(self, minimizer):
    self.n_cycles += 1
    #print self.n_cycles
    #for i in range(self.n_gaussians):
      #scale, mu, S = tuple(self.x[i*3:i*3+3])
      #sigma = math.sqrt(1/(2 * S**2))
      #print scale, mu, sigma

  def pyplot(self):
    from matplotlib import pyplot
    pyplot.plot(self.x_obs, self.y_obs)
    pyplot.plot(self.x_obs, self.compute_y_calc())
    for i in range(self.n_gaussians):
      scale, mu, S = tuple(self.x[i*3:i*3+3])
      y_calc = scale * flex.exp(-flex.pow2(self.x_obs-mu) * S**2)
      pyplot.plot(self.x_obs, y_calc)
    pyplot.show()
