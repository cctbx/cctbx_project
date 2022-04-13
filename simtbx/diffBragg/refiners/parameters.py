
from __future__ import absolute_import, division, print_function
import numpy as np
from numpy import sin, cos, arcsin
from collections import OrderedDict


class Parameters(OrderedDict):
  def add(self, p):
    p.xpos = len(self)
    self[p.name] = p


class RangedParameter:
  # TODO, make setting attributes named 'max' and 'min' attributes illegal
  """
  simple tool for managing parameters during refinement

  We re-parameterize to create a psuedo-free parameter
  See https://lmfit.github.io/lmfit-py/bounds.html
  """

  def __init__(self, init=0, minval=-1, maxval=1, sigma=1, fix=False, center=None, beta=None,
               name="param", is_global=False):
    """

    :param init: initial value for parameter
    :param minval: min value
    :param maxval: max value
    :param sigma: refinement sensitivity factor
    :param fix: whether to fix the parameter
    :param center: restraint center
    :param beta: restraint variance (smaller values give rise to tighter restraints)
    :param name: str, an optional name for this parameter, for bookkeeping
    :param is_global: bool, useful flag, if True, the parameter is common to all shots in a multi-shot refinement
    """
    self.misc_data = None  #custom field for user data
    self.minval = minval
    self.maxval = maxval
    self.sigma = sigma
    self.init = init
    self.fix = fix
    self.center = center
    self.beta = beta
    self.name = name
    self.is_global = is_global  # useful flag for potential APIs, specifies if the parameter is shared amongst all shots ?
    self._current_val = None  # place holder for the last value (output of get_val)
    # TODO use _rescaled_val in get_restraint_term and get_deriv in order to limit the valls to get_val
    self.xpos = 0  # position of parameter in list of params
    if fix:
      self.minval = init - 1
      self.maxval = init + 1
    self._arcsin_term = None

  def get_restraint_deriv(self, reparam_val):
    val = self.get_val(reparam_val)
    delta = self.center - val
    deriv = self.get_deriv(reparam_val, -delta/self.beta)
    return deriv

  def get_restraint_val(self, reparam_val):
    if not self.refine:
      return 0
    val = self.get_val(reparam_val)
    dist = self.center - val
    restraint_term = .5*(np.log(2*np.pi*self.beta) + dist**2/self.beta)
    return restraint_term

  @property
  def refine(self):
    return not self.fix

  @property
  def maxval(self):
    return self._maxval

  @maxval.setter
  def maxval(self, val):
    self._maxval = val

  @property
  def minval(self):
    return self._minval

  @minval.setter
  def minval(self, val):
    self._minval = val

  @property
  def rng(self):
    if self.minval >= self.maxval:
      raise ValueError("minval (%f) for RangedParameter must be less than the maxval (%f)" % (self.minval, self.maxval))
    return self.maxval - self.minval

  @property
  def arcsin_term(self):
    if self._arcsin_term is None:
      self._arcsin_term = arcsin(2 * (self.init - self.minval) / self.rng - 1)
    return self._arcsin_term

  def get_val(self, x_current):
    sin_arg = self.sigma * (x_current - 1) + self.arcsin_term
    val = (sin(sin_arg) + 1) * self.rng / 2 + self.minval
    self.value = val
    self.xvalue = x_current
    return val

  def get_deriv(self, x_current, deriv):
    cos_arg = self.sigma * (x_current - 1) + self.arcsin_term #arcsin(2 * (self.init - self.minval) / self.rng - 1)
    dtheta_dx = self.rng / 2 * cos(cos_arg) * self.sigma
    return deriv*dtheta_dx

  def get_second_deriv(self, x_current, deriv, second_deriv):
    sin_arg = self.sigma * (x_current - 1) + arcsin(2 * (self.init - self.minval) / self.rng - 1)
    cos_arg = self.sigma * (x_current - 1) + arcsin(2 * (self.init - self.minval) / self.rng - 1)
    dtheta_dx = self.rng / 2 * cos(cos_arg) * self.sigma
    d2theta_dx2 = -sin(sin_arg)*self.sigma*self.sigma * self.rng / 2.
    return dtheta_dx*dtheta_dx*second_deriv + d2theta_dx2*deriv

  @property
  def misc_data(self):
    return self._misc_data

  @misc_data.setter
  def misc_data(self, val):
    self._misc_data = val


class PositiveParameter(RangedParameter):

  def get_val(self, x_current):
    return self.init*np.exp(self.sigma*(x_current-1))

  def get_deriv(self, x_current, deriv, x_is_theta=False):
    """
    :param x_current: unscaled or rescaled parameter (see x_is_theta description)
    :param deriv: gradient array
    :param x_is_theta: optional, can skip recomputing the exponential if x_current is the unscaled parameter
    :return:
    """
    if x_is_theta:
      dtheta_dx = x_current *self.sigma
    else:
      dtheta_dx = self.get_val(x_current) * self.sigma
    return deriv*dtheta_dx


class NormalParameter(RangedParameter):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)

  def get_val(self, x):
    return x

  def get_deriv(self, x, deriv):
    return deriv
