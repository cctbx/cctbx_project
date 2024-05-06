from __future__ import division
from libtbx import adopt_init_args

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
