from cctbx.uctbx.reduction_base import iteration_limit_exceeded
from cctbx.uctbx.reduction_base import reduction_base
from cctbx.uctbx.reduction_base import minimal_reduction_mixin
from scitbx.python_utils.math_utils import ifloor
import math

def entier(x):
  result = ifloor(x)
  if (x-result >= 1): return result + 1 # work around rounding errors
  assert 0 <= x-result < 1, "entier error"
  return result

class reduction(reduction_base):

  def __init__(self, unit_cell, relative_epsilon=None, iteration_limit=None):
    reduction_base.__init__(self, unit_cell, relative_epsilon, iteration_limit)
    while (self.step()): pass

  def _name(self):
    return "Gruber"

  def step(s):
    eq = s.eps_eq
    gt = s.eps_gt
    # N1
    if (gt(s.a, s.b) or (eq(s.a, s.b) and gt(abs(s.d), abs(s.e)))):
      s.n1_action()
    # N2
    if (gt(s.b, s.c) or (eq(s.b, s.c) and gt(abs(s.e), abs(s.f)))):
      s.n2_action()
      return 0001
    # N3
    if (s.def_gt_0()):
      s.n3_true_action()
    else:
      s.n3_false_action()
    if (s.b2_action()): return 0001
    if (s.b3_action()): return 0001
    if (s.b4_action()): return 0001
    if (s.b5_action()): return 0001
    return 00000

  def b2_action(s):
    if (not s.eps_gt(abs(s.d), s.b)): return 00000
    j = entier((s.d+s.b)/(2*s.b))
    if (j == 0): return 00000
    s.c += j*j*s.b - j*s.d
    s.d -= 2*j*s.b
    s.e -= j*s.f
    s.cb_update((1,0,0,0,1,-j,0,0,1))
    assert s.c > 0
    return 0001

  def b3_action(s):
    if (not s.eps_gt(abs(s.e), s.a)): return 00000
    j = entier((s.e+s.a)/(2*s.a))
    if (j == 0): return 00000
    s.c += j*j*s.a - j*s.e
    s.d -= j*s.f
    s.e -= 2*j*s.a
    s.cb_update((1,0,-j,0,1,0,0,0,1))
    assert s.c > 0
    return 0001

  def b4_action(s):
    if (not s.eps_gt(abs(s.f), s.a)): return 00000
    j = entier((s.f+s.a)/(2*s.a))
    if (j == 0): return 00000
    s.b += j*j*s.a - j*s.f
    s.d -= j*s.e
    s.f -= 2*j*s.a
    s.cb_update((1,-j,0,0,1,0,0,0,1))
    assert s.b > 0
    return 0001

  def b5_action(s):
    de = s.d + s.e
    fab = s.f + s.a + s.b
    if (not s.eps_lt(de+fab, 0)): return 00000
    j = entier((de+fab)/(2*fab))
    if (j == 0): return 00000
    s.c += j*j*fab-j*de
    s.d -= j*(2*s.b+s.f)
    s.e -= j*(2*s.a+s.f)
    s.cb_update((1,0,-j,0,1,-j,0,0,1))
    assert s.c > 0
    return 0001

class minimal_reduction(minimal_reduction_mixin, reduction):

  def __init__(self, unit_cell, expected_cycle_limit=None,
                                iteration_limit=None):
    minimal_reduction_mixin.__init__(self,
      unit_cell, expected_cycle_limit, iteration_limit)

  def n3_false_action(self):
    self.current_cycle_id = 1
    return reduction.n3_false_action(self)

  def b5_action(self):
    self.current_cycle_id = 2
    return reduction.b5_action(self)
