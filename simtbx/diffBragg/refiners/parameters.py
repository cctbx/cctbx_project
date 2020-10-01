
from numpy import sin, cos, arcsin


class RangedParameter:

  def __init__(self):
    self.minval = 0
    self.maxval = 1
    self.sigma = None
    self.init = None

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

  def get_val(self, x_current):
    sin_arg = self.sigma * (x_current - 1) + arcsin(2 * (self.init - self.minval) / self.rng - 1)
    val = (sin(sin_arg) + 1) * self.rng / 2 + self.minval
    return val

  def get_deriv(self, x_current, deriv):
    cos_arg = self.sigma * (x_current - 1) + arcsin(2 * (self.init - self.minval) / self.rng - 1)
    dtheta_dx = self.rng / 2 * cos(cos_arg) * self.sigma
    return deriv*dtheta_dx
