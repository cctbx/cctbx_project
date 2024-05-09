from __future__ import division
from libtbx import adopt_init_args
from scitbx.array_family import flex

class individual(object):
  def __init__(self,
               data,
               x,
               data_weight=1,
               restraints=None,
               restraints_weight=None,
               lower_bound=None,
               upper_bound=None,
               bound_flags=None):
    adopt_init_args(self, locals())
    self.t = None
    self.g = None
    self.d = None
    self.use_curvatures=False
    self.n = x.size()
    #
    self.update(x = self.x)
    self.f_start = self.t

  def __call__(self):
    f, g = self.target_and_gradients()
    return self.x, f, g

  def update(self, x):
    self.data.update(x = x)
    self.x = x
    self.t, self.g = 0, 0
    for p in [[self.data,       self.data_weight],
              [self.restraints, self.restraints_weight]]:
      source, weight = p
      if source is not None:
        self.t += source.target()*weight
        self.g += source.gradients()*weight

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
               data_weight       = 1,
               restraints_weight = None,
               max_shift         = None):
    adopt_init_args(self, locals())
    assert isinstance(max_shift, float)
    data.set_refine_sites(selection = selection)
    x = flex.vec3_double(data.get_x())
    shift = flex.vec3_double(x.size(), [max_shift, max_shift, max_shift])
    lower_bound = x.deep_copy()
    lower_bound.set_selected(selection, x-shift)
    upper_bound = x.deep_copy()
    upper_bound.set_selected(selection, x+shift)
    self._calculator = individual(
      data        = data,
      x           = x.as_double(),
      lower_bound = lower_bound.as_double(),
      upper_bound = upper_bound.as_double(),
      bound_flags = flex.int(x.size()*3, 2))

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
    data.set_refine_u_iso(selection = selection)
    x = data.get_x()
    lower_bound = x.deep_copy()
    lower_bound.set_selected(selection, u_min)
    upper_bound = x.deep_copy()
    upper_bound.set_selected(selection, u_max)
    self._calculator = individual(
      data        = data,
      x           = x,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      bound_flags = flex.int(x.size(), 2))

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
    data.set_refine_occupancy(selection = selection)
    x = data.get_x()
    lower_bound = x.deep_copy()
    lower_bound.set_selected(selection, q_min)
    upper_bound = x.deep_copy()
    upper_bound.set_selected(selection, q_max)
    self._calculator = individual(
      data        = data,
      x           = x,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      bound_flags = flex.int(x.size(), 2))

  def calculator(self):
    return self._calculator
