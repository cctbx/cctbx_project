from cctbx.uctbx.reduction_base import iteration_limit_exceeded
from cctbx.uctbx.reduction_base import reduction_base
from cctbx.uctbx.reduction_base import minimum_reduction_mixin
from cctbx import uctbx
from scitbx import matrix

def entier(x):
  "greatest integer which is not greater than x"
  result = int(x)
  if (x-result < 0): result -= 1
  if (not (x-result < 1)): result += 1 # work around rounding errors
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
      return True
    # N3
    if (s.def_gt_0()):
      s.n3_true_action()
    else:
      s.n3_false_action()
    if (s.b2_action()): return True
    if (s.b3_action()): return True
    if (s.b4_action()): return True
    if (s.b5_action()): return True
    return False

  def b2_action(s):
    if (not s.eps_gt(abs(s.d), s.b)): return False
    j = entier((s.d+s.b)/(2*s.b))
    if (j == 0): return False
    s.cb_update((1,0,0,0,1,-j,0,0,1))
    s.c += j*j*s.b - j*s.d
    s.d -= 2*j*s.b
    s.e -= j*s.f
    assert s.c > 0
    return True

  def b3_action(s):
    if (not s.eps_gt(abs(s.e), s.a)): return False
    j = entier((s.e+s.a)/(2*s.a))
    if (j == 0): return False
    s.cb_update((1,0,-j,0,1,0,0,0,1))
    s.c += j*j*s.a - j*s.e
    s.d -= j*s.f
    s.e -= 2*j*s.a
    assert s.c > 0
    return True

  def b4_action(s):
    if (not s.eps_gt(abs(s.f), s.a)): return False
    j = entier((s.f+s.a)/(2*s.a))
    if (j == 0): return False
    s.cb_update((1,-j,0,0,1,0,0,0,1))
    s.b += j*j*s.a - j*s.f
    s.d -= j*s.e
    s.f -= 2*j*s.a
    assert s.b > 0
    return True

  def b5_action(s):
    de = s.d + s.e
    fab = s.f + s.a + s.b
    if (not s.eps_lt(de+fab, 0)): return False
    j = entier((de+fab)/(2*fab))
    if (j == 0): return False
    s.cb_update((1,0,-j,0,1,-j,0,0,1))
    s.c += j*j*fab-j*de
    s.d -= j*(2*s.b+s.f)
    s.e -= j*(2*s.a+s.f)
    assert s.c > 0
    return True

class minimum_reduction(minimum_reduction_mixin, reduction):
  """Development and regression test code. Do not use for applications.
     Use uctbx.fast_minimum_reduction instead.
  """

  def __init__(self, unit_cell, iteration_limit=None):
    minimum_reduction_mixin.__init__(self, unit_cell, iteration_limit)

  def n3_false_action(self):
    reduction.n3_false_action(self)
    if (not self.significant_change_test()):
      raise StopIteration

class fast_minimum_reduction(object):
  """Development and regression test code. Do not use for applications.
     Use uctbx.fast_minimum_reduction instead.
  """

  def __init__(self, unit_cell, iteration_limit=None,
                                multiplier_significant_change_test=16,
                                min_n_no_significant_change=2):
    if (iteration_limit is None):
      self._iteration_limit = 100
    else:
      self._iteration_limit = iteration_limit
    self.multiplier_significant_change_test=multiplier_significant_change_test
    self.min_n_no_significant_change=min_n_no_significant_change
    sym_mat3 = unit_cell.metrical_matrix()
    self.a = sym_mat3[0]
    self.b = sym_mat3[1]
    self.c = sym_mat3[2]
    self.d = 2*sym_mat3[5]
    self.e = 2*sym_mat3[4]
    self.f = 2*sym_mat3[3]
    self._r_inv = matrix.sqr((1,0,0,0,1,0,0,0,1))
    self._n_iterations = 0
    self._n_no_significant_change = 0
    self._last_abc_significant_change_test = (-self.a,-self.b,-self.c)
    while (self.step()): pass

  def as_gruber_matrix(self):
    return (self.a, self.b, self.c, self.d, self.e, self.f)

  def as_niggli_matrix(self):
    return (self.a, self.b, self.c, self.d/2, self.e/2, self.f/2)

  def as_sym_mat3(self):
    return (self.a, self.b, self.c, self.f/2, self.e/2, self.d/2)

  def as_unit_cell(self):
    return uctbx.unit_cell(metrical_matrix=self.as_sym_mat3())

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
      return True
    # N3
    if (s.def_gt_0()):
      s.n3_true_action()
    else:
      s.n3_false_action()
      if (not s.significant_change_test()):
        return False
    if (s.b2_action()): return True
    if (s.b3_action()): return True
    if (s.b4_action()): return True
    if (s.b5_action()): return True
    return False

  def significant_change_test(self):
    abc = (self.a,self.b,self.c)
    m = self.multiplier_significant_change_test
    change = tuple([(new*m+(new-last))-new*m
      for new,last in zip(abc, self._last_abc_significant_change_test)])
    if (change == (0,0,0)):
      self._n_no_significant_change += 1
      if (self._n_no_significant_change >= self.min_n_no_significant_change):
        return False
    else:
      self._n_no_significant_change = 0
    self._last_abc_significant_change_test = abc
    return True

  def cb_update(self, m_elems):
    if (self._n_iterations == self._iteration_limit):
      raise iteration_limit_exceeded(
        "Gruber iteration limit exceeded (limit=%d)."
        % self._iteration_limit)
    self._r_inv *= matrix.sqr(m_elems)
    self._n_iterations += 1

  def n1_action(self):
    self.cb_update((0,-1,0, -1,0,0, 0,0,-1))
    self.a, self.b = self.b, self.a
    self.d, self.e = self.e, self.d

  def n2_action(self):
    self.cb_update((-1,0,0, 0,0,-1, 0,-1,0))
    self.b, self.c = self.c, self.b
    self.e, self.f = self.f, self.e

  def n3_true_action(s):
    i,j,k = 1,1,1
    if (s.d < 0): i = -1
    if (s.e < 0): j = -1
    if (s.f < 0): k = -1
    s.cb_update((i,0,0, 0,j,0, 0,0,k))
    s.d = abs(s.d)
    s.e = abs(s.e)
    s.f = abs(s.f)

  def n3_false_action(s):
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
    s.cb_update((f[0],0,0, 0,f[1],0, 0,0,f[2]))
    s.d = -abs(s.d)
    s.e = -abs(s.e)
    s.f = -abs(s.f)

  def b2_action(s):
    if (not s.b < abs(s.d)): return False
    j = entier((s.d+s.b)/(2*s.b))
    if (j == 0): return False
    s.cb_update((1,0,0,0,1,-j,0,0,1))
    s.c += j*j*s.b - j*s.d
    s.d -= 2*j*s.b
    s.e -= j*s.f
    assert 0 < s.c
    return True

  def b3_action(s):
    if (not s.a < abs(s.e)): return False
    j = entier((s.e+s.a)/(2*s.a))
    if (j == 0): return False
    s.cb_update((1,0,-j,0,1,0,0,0,1))
    s.c += j*j*s.a - j*s.e
    s.d -= j*s.f
    s.e -= 2*j*s.a
    assert 0 < s.c
    return True

  def b4_action(s):
    if (not s.a < abs(s.f)): return False
    j = entier((s.f+s.a)/(2*s.a))
    if (j == 0): return False
    s.cb_update((1,-j,0,0,1,0,0,0,1))
    s.b += j*j*s.a - j*s.f
    s.d -= j*s.e
    s.f -= 2*j*s.a
    assert 0 < s.b
    return True

  def b5_action(s):
    de = s.d + s.e
    fab = s.f + s.a + s.b
    if (not de+fab < 0): return False
    j = entier((de+fab)/(2*fab))
    if (j == 0): return False
    s.cb_update((1,0,-j,0,1,-j,0,0,1))
    s.c += j*j*fab-j*de
    s.d -= j*(2*s.b+s.f)
    s.e -= j*(2*s.a+s.f)
    assert 0 < s.c
    return True
