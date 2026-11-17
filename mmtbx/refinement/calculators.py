from __future__ import division
from libtbx import adopt_init_args
from cctbx import adptbx
from cctbx.array_family import flex

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

class scaled(object):
  """
  Present a calculator to the minimizer in scaled variables x' = scale * x.

  The minimizer sees x', the bounds scaled the same way, the target of
  the wrapped calculator at x = x' / scale and its gradient divided by
  the scale. A diagonal scale keeps box bounds as box bounds, so this
  works with L-BFGS-B. Scaling by the square root of a curvature
  diagonal makes the target's curvature close to one in every variable.

  Parameters
  ----------
  calculator : individual
      The calculator to wrap; its x, bounds and model are updated as the
      minimizer moves.
  scale : flex.double
      One positive factor per variable.
  """
  def __init__(self, calculator, scale):
    adopt_init_args(self, locals())
    assert scale.size() == calculator.x.size()
    assert flex.min(scale) > 0
    self.x           = calculator.x * scale
    self.lower_bound = calculator.lower_bound * scale
    self.upper_bound = calculator.upper_bound * scale
    self.bound_flags = calculator.bound_flags
    self.n           = self.x.size()

  def update(self, x):
    self.x = x
    self.calculator.update(x = x / self.scale)

  def target(self):
    return self.calculator.target()

  def gradients(self):
    return self.calculator.gradients() / self.scale

  def target_and_gradients(self):
    self.update(x = self.x)
    return self.target(), self.gradients()

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

class xyz_adp(object):
  """
  Calculator for joint refinement of sites and u_iso.

  The variables are in block order, the sites_cart of all atoms followed
  by their u_iso; see mmtbx.refinement.data.fs.set_refine_sites_and_u_iso.
  Only the data term is available (restraints provide the box bounds).

  Parameters
  ----------
  data : mmtbx.refinement.data.fs
  restraints : mmtbx.refinement.restraints.manager
  selection : flex.bool, optional
  data_weight : float
  max_shift : float
      Half-width of the coordinate box, in Angstrom.
  u_min, u_max : float
      Bounds on u_iso, in Angstrom squared.
  """
  def __init__(self,
               data,
               restraints,
               selection   = None,
               data_weight = 1.,
               max_shift   = None,
               u_min       = None,
               u_max       = None):
    adopt_init_args(self, locals())
    assert isinstance(max_shift, float)
    data.set_refine_sites_and_u_iso(selection = selection)
    restraints.set_use_xyz_adp(
      selection = selection,
      max_shift = max_shift,
      b_min     = adptbx.u_as_b(u_min),
      b_max     = adptbx.u_as_b(u_max))
    x = restraints.get_x()
    assert x.size() == data.get_x().size()
    self._calculator = individual(
      data              = data,
      restraints        = restraints,
      data_weight       = data_weight,
      restraints_weight = None,
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
