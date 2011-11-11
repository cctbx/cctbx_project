import math

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
    self.x = flex.double([scale, mu, 1]) # can we make a simple estimate of sigma too?
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self)
    self.scale, self.mu, S = tuple(self.x)
    self.sigma = math.sqrt(1/(2 * S**2))
    del self.x

  def compute_functional_and_gradients(self):
    scale, mu, S = tuple(self.x)
    # reparametrisation:
    #   S**2 = 1/(2 * sigma**2)
    exponential_part = flex.exp(-flex.pow2(self.x_obs-mu) * S**2)
    y_calc = scale * exponential_part
    delta_y = self.y_obs - y_calc
    f = flex.sum(flex.pow2(delta_y))
    g = flex.double((
      -2 * flex.sum(delta_y * exponential_part),
      -2 * flex.sum(delta_y * 2 * S**2 * scale * exponential_part * (self.x_obs - mu)),
      -2 * flex.sum(delta_y * - 2 * scale * S * exponential_part * flex.pow2(mu - self.x_obs))
    ))
    return f, g

  #def callback_after_step(self, minimizer):
    #print list(self.x)

