from cctbx import uctbx
from cctbx import matrix

class gruber_parameterization:

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
    return uctbx.unit_cell(self.as_sym_mat3(), is_metrical_matrix=0001)

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

  def meets_main_conditions(self):
    lt = self.eps_lt
    gt = self.eps_gt
    a,b,c,d,e,f = (self.a,self.b,self.c,self.d,self.e,self.f)
    if (gt(a, b)): return 00000
    if (gt(b, c)): return 00000
    if (gt(abs(d), b)): return 00000
    if (gt(abs(e), a)): return 00000
    if (gt(abs(f), a)): return 00000
    type = self.type()
    if (type == 0): return 00000
    if (type == 2):
      if (lt(d+e+f+a+b, 0)): return 00000
    return 0001

  def is_buerger_cell(self):
    if (not self.meets_main_conditions()): return 00000
    eq = self.eps_eq
    gt = self.eps_gt
    a,b,c,d,e,f = (self.a,self.b,self.c,self.d,self.e,self.f)
    if (eq(a, b)):
      if (gt(abs(d), abs(e))): return 00000
    if (eq(b, c)):
      if (gt(abs(e), abs(f))): return 00000
    return 0001

  def is_niggli_cell(self):
    if (not self.is_buerger_cell()): return 00000
    eq = self.eps_eq
    gt = self.eps_gt
    a,b,c,d,e,f = (self.a,self.b,self.c,self.d,self.e,self.f)
    if (eq(d, b)):
      if (gt(f, e+e)): return 00000
    if (eq(e, a)):
      if (gt(f, d+d)): return 00000
    if (eq(f, a)):
      if (gt(e, d+d)): return 00000
    if (eq(d, -b)):
      if (not eq(f, 0)): return 00000
    if (eq(e, -a)):
      if (not eq(f, 0)): return 00000
    if (eq(f, -a)):
      if (not eq(e, 0)): return 00000
    if (eq(d+e+f+a+b, 0)):
      if (gt(a+a+e+e+f, 0)): return 00000
    return 0001
