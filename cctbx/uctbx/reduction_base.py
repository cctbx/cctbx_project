from __future__ import absolute_import, division, print_function
from cctbx import uctbx
from scitbx import matrix
from six.moves import zip

class gruber_parameterization(object):

  def __init__(self, unit_cell, relative_epsilon=None):
    if (relative_epsilon is None): relative_epsilon = 1.e-5
    sym_mat3 = unit_cell.metrical_matrix()
    self.a = sym_mat3[0]
    self.b = sym_mat3[1]
    self.c = sym_mat3[2]
    self.d = 2*sym_mat3[5]
    self.e = 2*sym_mat3[4]
    self.f = 2*sym_mat3[3]
    self.epsilon = unit_cell.volume()**(1/3.) * relative_epsilon

  def as_gruber_matrix(self):
    return (self.a, self.b, self.c, self.d, self.e, self.f)

  def as_niggli_matrix(self):
    return (self.a, self.b, self.c, self.d/2, self.e/2, self.f/2)

  def as_sym_mat3(self):
    return (self.a, self.b, self.c, self.f/2, self.e/2, self.d/2)

  def as_unit_cell(self):
    return uctbx.unit_cell(metrical_matrix=self.as_sym_mat3())

  def eps_lt(self, x, y):
    return x < y - self.epsilon

  def eps_gt(self, x, y):
    return self.eps_lt(y, x)

  def eps_eq(self, x, y):
    return not (self.eps_lt(x, y) or self.eps_lt(y, x))

  def def_test(self):
    lt = self.eps_lt
    d,e,f = (self.d,self.e,self.f)
    n_zero = 0
    n_positive = 0
    if (lt(0, d)): n_positive += 1
    elif (not lt(d, 0)): n_zero += 1
    if (lt(0, e)): n_positive += 1
    elif (not lt(e, 0)): n_zero += 1
    if (lt(0, f)): n_positive += 1
    elif (not lt(f, 0)): n_zero += 1
    return n_zero, n_positive

  def def_gt_0(self):
    n_zero, n_positive = self.def_test()
    return n_positive == 3 or (n_zero == 0 and n_positive == 1)

  def type(self):
    n_zero, n_positive = self.def_test()
    if (n_positive == 3): return 1
    if (n_positive == 0): return 2
    return 0

  def meets_primary_conditions(self):
    gt = self.eps_gt
    s = self
    if (gt(s.a, s.b)): return False
    if (gt(s.b, s.c)): return False
    if (gt(abs(s.d), s.b)): return False
    if (gt(abs(s.e), s.a)): return False
    if (gt(abs(s.f), s.a)): return False
    return True

  def meets_main_conditions(self):
    if (not self.meets_primary_conditions()): return False
    type = self.type()
    if (type == 0): return False
    if (type == 2):
      lt = self.eps_lt
      s = self
      if (lt(s.d+s.e+s.f+s.a+s.b, 0)): return False
    return True

  def is_buerger_cell(self):
    if (not self.meets_main_conditions()): return False
    eq = self.eps_eq
    gt = self.eps_gt
    s = self
    if (eq(s.a, s.b)):
      if (gt(abs(s.d), abs(s.e))): return False
    if (eq(s.b, s.c)):
      if (gt(abs(s.e), abs(s.f))): return False
    return True

  def is_niggli_cell(self):
    if (not self.is_buerger_cell()): return False
    eq = self.eps_eq
    gt = self.eps_gt
    s = self
    if (eq(s.d, s.b)):
      if (gt(s.f, s.e+s.e)): return False
    if (eq(s.e, s.a)):
      if (gt(s.f, s.d+s.d)): return False
    if (eq(s.f, s.a)):
      if (gt(s.e, s.d+s.d)): return False
    if (eq(s.d, -s.b)):
      if (not eq(s.f, 0)): return False
    if (eq(s.e, -s.a)):
      if (not eq(s.f, 0)): return False
    if (eq(s.f, -s.a)):
      if (not eq(s.e, 0)): return False
    if (eq(s.d+s.e+s.f+s.a+s.b, 0)):
      if (gt(s.a+s.a+s.e+s.e+s.f, 0)): return False
    return True

class iteration_limit_exceeded(RuntimeError): pass

class reduction_base(gruber_parameterization):

  def __init__(self, unit_cell, relative_epsilon, iteration_limit):
    if (iteration_limit is None):
      self._iteration_limit = 1000
    else:
      self._iteration_limit = iteration_limit
    gruber_parameterization.__init__(self, unit_cell, relative_epsilon)
    self._r_inv = matrix.sqr((1,0,0,0,1,0,0,0,1))
    self._n_iterations = 0
    self._last_after_all_obtuse_action = (-1,-1,-1)

  def iteration_limit(self):
    return self._iteration_limit

  def r_inv(self):
    return self._r_inv

  def change_of_basis_op(self):
    from cctbx import sgtbx
    return sgtbx.change_of_basis_op(
      sgtbx.rt_mx(sgtbx.rot_mx(self._r_inv.elems, 1))).inverse()

  def n_iterations(self):
    return self._n_iterations

  def cb_update(self, m_elems):
    if (self._n_iterations == self._iteration_limit):
      raise iteration_limit_exceeded(
        "%s iteration limit exceeded (limit=%d)."
        % (self._name(), self._iteration_limit))
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

  def n3_true_action(self):
    lt = self.eps_lt
    i,j,k = 1,1,1
    if (lt(self.d, 0)): i = -1
    if (lt(self.e, 0)): j = -1
    if (lt(self.f, 0)): k = -1
    self.cb_update((i,0,0, 0,j,0, 0,0,k))
    self.d = abs(self.d)
    self.e = abs(self.e)
    self.f = abs(self.f)

  def n3_false_action(self):
    lt = self.eps_lt
    gt = self.eps_gt
    f = [1,1,1]
    z = -1
    if (gt(self.d, 0)): f[0] = -1
    elif (not lt(self.d, 0)): z = 0
    if (gt(self.e, 0)): f[1] = -1
    elif (not lt(self.e, 0)): z = 1
    if (gt(self.f, 0)): f[2] = -1
    elif (not lt(self.f, 0)): z = 2
    if (f[0]*f[1]*f[2] < 0):
      assert z != -1
      f[z] = -1
    self.cb_update((f[0],0,0, 0,f[1],0, 0,0,f[2]))
    self.d = -abs(self.d)
    self.e = -abs(self.e)
    self.f = -abs(self.f)

class minimum_reduction_mixin(object):
  """Development and regression test code. Do not use for applications.
     Use uctbx.fast_minimum_reduction instead.
  """

  def __init__(self, unit_cell, iteration_limit=None,
                                multiplier_significant_change_test=16,
                                min_n_no_significant_change=2):
    self.multiplier_significant_change_test=multiplier_significant_change_test
    self.min_n_no_significant_change=min_n_no_significant_change
    self._n_no_significant_change = 0
    a,b,c = unit_cell.metrical_matrix()[:3]
    self._last_abc_significant_change_test = (-a,-b,-c)
    self.reduction = self.__class__.__bases__[1]
    try:
      self.reduction.__init__(self,
        unit_cell=unit_cell,
        relative_epsilon=0,
        iteration_limit=iteration_limit)
      self.termination_due_to_significant_change_test = False
    except StopIteration:
      self.termination_due_to_significant_change_test = True

  def eps_eq(self, x, y):
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
