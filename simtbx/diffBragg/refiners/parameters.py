
from __future__ import absolute_import, division, print_function
from numpy import sin, cos, arcsin


class RangedParameter:
    # TODO, make setting attributes named 'max' and 'min' attributes illegal

  def __init__(self):
    self.minval = 0
    self.maxval = 1
    self.sigma = None
    self.init = None
    self.fix = False
    self._arcsin_term = None

  #@property
  #def init(self):
  #  return self.__init

  #@init.setter
  #def init(self, val):
  #  if val is not None:
  #    if val < self.minval:
  #      raise ValueError("Parameter cannot be initialized to less than the minimum")
  #    if val >self.maxval:
  #      raise ValueError("Parameter cannot be initialized to more than the maximum")
  #  self._init = val
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
    #if self.fix:
    #  return self.init
    #else:

    #sin_arg = self.sigma * (x_current - 1) + arcsin(2 * (self.init - self.minval) / self.rng - 1)
    sin_arg = self.sigma * (x_current - 1) + self.arcsin_term
    val = (sin(sin_arg) + 1) * self.rng / 2 + self.minval
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


class NormalParameter(RangedParameter):

  def __init__(self):
    super().__init__()

  def get_val(self, x):
    return x

  def get_deriv(self, x, deriv):
    return deriv


class Parameters:

  def __init__(self):
    self.Ncells_abc ={}
    self.Ncells_def ={}
    self.rotXYZ ={}
    self.Bmatrix ={}
    self.spot_scale ={}
    self.eta ={}
    self.wavelen_offset ={}
    self.wavelen_scale ={}
    self.panelX = []
    self.panelY = []
    self.panelZ = []
    self.panelO = []
    self.panelF = []
    self.panelS = []
    self.panelOrig = []
    self.panelFast =[]
    self.panelSlow = []
    self.keys = []

  def safe_append(self, some_dict,name,val):
    if name not in some_dict:
      some_dict[name] = [val]
    else:
      some_dict[name].append(val)
    if name not in self.keys:
      self.keys.append(name)

  def add_panelOrig(self, vals):
    self.panelOrig.append(vals)
  def add_panelFast(self, vals):
    self.panelFast.append(vals)
  def add_panelSlow(self, vals):
    self.panelSlow.append(vals)

  def add_panelX(self, vals):
    self.panelX.append(vals)
  def add_panelY(self, vals):
    self.panelY.append(vals)
  def add_panelZ(self, vals):
    self.panelZ.append(vals)

  def add_panelO(self, vals):
    self.panelO.append(vals)
  def add_panelF(self, vals):
    self.panelF.append(vals)
  def add_panelS(self, vals):
    self.panelS.append(vals)

  def add_Ncells_abc(self,name, val):
    #if len(val) != 3:
    #  raise ValueError("Ncells abc must be a 3-tuple")
    self.safe_append(self.Ncells_abc, name, val)

  def add_Ncells_def(self,name, val):
    self.safe_append(self.Ncells_def, name, val)

  def add_spot_scale(self,name, val):
    self.safe_append(self.spot_scale,name,val)

  def add_rotXYZ(self,name, val):
    self.safe_append(self.rotXYZ, name, val)

  def add_Bmatrix(self,name, val):
    self.safe_append(self.Bmatrix,name, val.elems)

  def add_eta(self,name,val):
    self.safe_append(self.eta, name, val)

  def add_wavelen_offset(self,name, val):
    self.safe_append(self.wavelen_offset, name, val)

  def add_wavelen_scale(self,name, val):
    self.safe_append(self.wavelen_scale, name, val)
