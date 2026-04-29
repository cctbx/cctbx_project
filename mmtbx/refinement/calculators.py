from __future__ import division
from libtbx import adopt_init_args
from cctbx import adptbx

class individual(object):
  def __init__(self,
               data,
               x,
               data_weight=1,
               restraints=None,
               restraints_weight=None):
    adopt_init_args(self, locals())
    self.t = None
    self.g = None
    self.d = None
    self.use_curvatures=False
    self.n = x.size()
    self.lower_bound = restraints.lower_bound()
    self.upper_bound = restraints.upper_bound()
    self.bound_flags = restraints.bound_flags()
    #
    self.update(x = self.x)
    self.f_start = self.t

  def __call__(self):
    f, g = self.target_and_gradients()
    return self.x, f, g

  def target(self): return self.t

  def gradients(self): return self.g

  def update(self, x):
    self.data.update(x = x)
    if self.restraints is not None:
      self.restraints.update(x = x)
    self.x = x
    self.t, self.g = 0, 0
    for p in [[self.data,       self.data_weight],
              [self.restraints, self.restraints_weight]]:
      source, weight = p
      if not None in [source, weight]:
        t,g = source.target(), source.gradients()
        if t is not None:
          self.t += t*weight
          self.g += g*weight

  def target_and_gradients(self):
    self.update(x = self.x)
    return self.t, self.g

  def compute_functional_and_gradients(self):
    return self.target_and_gradients()

class xyz(object):
  def __init__(self,
               data              = None,
               restraints        = None,
               selection         = None,
               data_weight       = 1.,
               restraints_weight = 1.,
               max_shift         = None):
    adopt_init_args(self, locals())
    assert isinstance(max_shift, float)
    assert [data, restraints].count(None) != 2
    if data is not None:
      data.set_refine_sites(selection = selection)
      x = data.get_x()
    if restraints is not None:
      restraints.set_use_xyz(selection = selection, max_shift = max_shift)
      x = restraints.get_x()
    self._calculator = individual(
      data              = data,
      restraints        = restraints,
      data_weight       = data_weight,
      restraints_weight = restraints_weight,
      x                 = x)

  def calculator(self):
    return self._calculator

class adp(object):
  def __init__(self,
               data              = None,
               restraints        = None,
               selection         = None,
               data_weight       = 1,
               restraints_weight = None,
               u_min             = None,
               u_max             = None):
    adopt_init_args(self, locals())
    assert [data, restraints].count(None) != 2
    if data is not None:
      data.set_refine_u_iso(selection = selection)
      x = data.get_x()
    if restraints is not None:
      restraints.set_use_adp(
        selection = selection,
        b_min     = adptbx.u_as_b(u_min),
        b_max     = adptbx.u_as_b(u_max))
      x = restraints.get_x()
    self._calculator = individual(
      data              = data,
      restraints        = restraints,
      data_weight       = data_weight,
      restraints_weight = restraints_weight,
      x                 = x)

  def calculator(self):
    return self._calculator

class occ(object):
  def __init__(self,
               data              = None,
               restraints        = None,
               selection         = None,
               data_weight       = 1,
               restraints_weight = None,
               q_min             = None,
               q_max             = None):
    adopt_init_args(self, locals())
    assert [data, restraints].count(None) != 2
    if data is not None:
      data.set_refine_occupancy(selection = selection)
      x = data.get_x()
    if restraints is not None:
      restraints.set_use_occ(
        selection = selection,
        q_min     = q_min,
        q_max     = q_max)
      x = restraints.get_x()
    self._calculator = individual(
      data              = data,
      restraints        = restraints,
      data_weight       = data_weight,
      restraints_weight = restraints_weight,
      x                 = x)

  def calculator(self):
    return self._calculator
