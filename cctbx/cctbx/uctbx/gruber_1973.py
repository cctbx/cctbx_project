from cctbx.uctbx.reduction_base import iteration_limit_exceeded
from cctbx.uctbx.reduction_base import reduction_base
from cctbx.uctbx.reduction_base import minimal_reduction_mixin
from cctbx import uctbx
from cctbx import matrix
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
    s.cb_update((1,0,-j,0,1,-j,0,0,1))
    s.c += j*j*fab-j*de
    s.d -= j*(2*s.b+s.f)
    s.e -= j*(2*s.a+s.f)
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

class fast_minimal_reduction:

  def __init__(self, unit_cell, expected_cycle_limit=None,
                                iteration_limit=None):
    if (expected_cycle_limit is None):
      self.expected_cycle_limit = 2
    else:
      self.expected_cycle_limit = expected_cycle_limit
    if (iteration_limit is None):
      self._iteration_limit = 100
    else:
      self._iteration_limit = iteration_limit
    sym_mat3 = unit_cell.metrical_matrix()
    self.a = sym_mat3[0]
    self.b = sym_mat3[1]
    self.c = sym_mat3[2]
    self.d = 2*sym_mat3[5]
    self.e = 2*sym_mat3[4]
    self.f = 2*sym_mat3[3]
    self._r_inv = matrix.sqr((1,0,0,0,1,0,0,0,1))
    self._n_iterations = 0
    self.current_cycle_id = 0
    self.last_cycle_id = 0
    self.n_expected_cycles = 0
    self.had_expected_cycle = 00000
    while (self.step()): pass

  def as_gruber_matrix(self):
    return (self.a, self.b, self.c, self.d, self.e, self.f)

  def as_niggli_matrix(self):
    return (self.a, self.b, self.c, self.d/2, self.e/2, self.f/2)

  def as_sym_mat3(self):
    return (self.a, self.b, self.c, self.f/2, self.e/2, self.d/2)

  def as_unit_cell(self):
    return uctbx.unit_cell(self.as_sym_mat3(), is_metrical_matrix=0001)

  def iteration_limit(self):
    return self._iteration_limit

  def r_inv(self):
    return self._r_inv

  def n_iterations(self):
    return self._n_iterations

  def def_test(s):
    n_zero = 0
    n_positive = 0
    if (0 < s.d): n_positive += 1
    elif (not s.d < 0): n_zero += 1
    if (0 < s.e): n_positive += 1
    elif (not s.e < 0): n_zero += 1
    if (0 < s.f): n_positive += 1
    elif (not s.f < 0): n_zero += 1
    return n_zero, n_positive

  def def_gt_0(self):
    n_zero, n_positive = self.def_test()
    return n_positive == 3 or (n_zero == 0 and n_positive == 1)

  def type(self):
    n_zero, n_positive = self.def_test()
    if (n_positive == 3): return 1
    if (n_positive == 0): return 2
    return 0

  def step(s):
    # N1
    if (s.b < s.a):
      s.n1_action()
    # N2
    if (s.c < s.b):
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

  def cb_update(self, m_elems):
    if (self.current_cycle_id == 1):
      if (self.last_cycle_id != 2):
        self.n_expected_cycles = 0
    elif (self.current_cycle_id == 2):
      if (self.last_cycle_id != 1):
        self.n_expected_cycles = 0
      else:
        self.n_expected_cycles += 1
        if (self.n_expected_cycles == self.expected_cycle_limit):
          self.had_expected_cycle = 0001
          return
    else:
      self.n_expected_cycles = 0
    self.last_cycle_id = self.current_cycle_id
    self.current_cycle_id = 0
    if (self._n_iterations == self._iteration_limit):
      raise iteration_limit_exceeded(
        "Gruber iteration limit exceeded (limit=%d)."
        % self._iteration_limit)
    self._r_inv *= matrix.sqr(m_elems)
    self._n_iterations += 1

  def n1_action(self):
    self.a, self.b = self.b, self.a
    self.d, self.e = self.e, self.d
    self.cb_update((0,-1,0, -1,0,0, 0,0,-1))

  def n2_action(self):
    self.b, self.c = self.c, self.b
    self.e, self.f = self.f, self.e
    self.cb_update((-1,0,0, 0,0,-1, 0,-1,0))

  def n3_true_action(s):
    i,j,k = 1,1,1
    if (s.d < 0): i = -1
    if (s.e < 0): j = -1
    if (s.f < 0): k = -1
    s.d = abs(s.d)
    s.e = abs(s.e)
    s.f = abs(s.f)
    s.cb_update((i,0,0, 0,j,0, 0,0,k))

  def n3_false_action(s):
    s.current_cycle_id = 1
    f = [1,1,1]
    z = -1
    if (0 < s.d): f[0] = -1
    elif (not s.d < 0): z = 0
    if (0 < s.e): f[1] = -1
    elif (not s.e < 0): z = 1
    if (0 < s.f): f[2] = -1
    elif (not s.f < 0): z = 2
    if (f[0]*f[1]*f[2] < 0):
      assert z != -1
      f[z] = -1
    s.d = -abs(s.d)
    s.e = -abs(s.e)
    s.f = -abs(s.f)
    s.cb_update((f[0],0,0, 0,f[1],0, 0,0,f[2]))

  def b2_action(s):
    if (not s.b < abs(s.d)): return 00000
    j = entier((s.d+s.b)/(2*s.b))
    if (j == 0): return 00000
    s.c += j*j*s.b - j*s.d
    s.d -= 2*j*s.b
    s.e -= j*s.f
    s.cb_update((1,0,0,0,1,-j,0,0,1))
    assert 0 < s.c
    return 0001

  def b3_action(s):
    if (not s.a < abs(s.e)): return 00000
    j = entier((s.e+s.a)/(2*s.a))
    if (j == 0): return 00000
    s.c += j*j*s.a - j*s.e
    s.d -= j*s.f
    s.e -= 2*j*s.a
    s.cb_update((1,0,-j,0,1,0,0,0,1))
    assert 0 < s.c
    return 0001

  def b4_action(s):
    if (not s.a < abs(s.f)): return 00000
    j = entier((s.f+s.a)/(2*s.a))
    if (j == 0): return 00000
    s.b += j*j*s.a - j*s.f
    s.d -= j*s.e
    s.f -= 2*j*s.a
    s.cb_update((1,-j,0,0,1,0,0,0,1))
    assert 0 < s.b
    return 0001

  def b5_action(s):
    s.current_cycle_id = 2
    de = s.d + s.e
    fab = s.f + s.a + s.b
    if (not de+fab < 0): return 00000
    j = entier((de+fab)/(2*fab))
    if (j == 0): return 00000
    s.cb_update((1,0,-j,0,1,-j,0,0,1))
    if (s.had_expected_cycle): return 00000
    s.c += j*j*fab-j*de
    s.d -= j*(2*s.b+s.f)
    s.e -= j*(2*s.a+s.f)
    assert 0 < s.c
    return 0001
