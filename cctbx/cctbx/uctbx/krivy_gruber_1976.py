from cctbx.uctbx import reduction_base
from cctbx import uctbx
from cctbx import matrix

class reduction(reduction_base.gruber_parameterization):

  def __init__(self, unit_cell, relative_epsilon=None, iteration_limit=1000):
    if (iteration_limit is None): iteration_limit = 1000
    reduction_base.gruber_parameterization.__init__(
      self, unit_cell, relative_epsilon)
    self._r_inv = matrix.sqr((1,0,0,0,1,0,0,0,1))
    self._n_iterations = 0
    while 1:
      if (self._n_iterations == iteration_limit):
        raise RuntimeError("Krivy-Gruber iteration limit exceeded (limit=%d)."
          % iteration_limit)
      stop = self.step()
      self._n_iterations += 1
      if (stop): break

  def r_inv(self):
    return self._r_inv

  def n_iterations(self):
    return self._n_iterations

  def step(self):
    s = self
    eq = s.eps_eq
    lt = s.eps_lt
    gt = s.eps_gt
    # A1
    if (gt(s.a, s.b) or (eq(s.a, s.b) and gt(abs(s.d), abs(s.e)))):
      s.a1_action()
    # A2
    if (gt(s.b, s.c) or (eq(s.b, s.c) and gt(abs(s.e), abs(s.f)))):
      s.a2_action()
      return 00000
    # A3
    if (s.def_gt_0()):
      s.a3_action()
    # A4
    else:
      s.a4_action()
    # A5
    if (gt(abs(s.d), s.b)
        or (eq(s.d, s.b) and lt(2*s.e, s.f))
        or (eq(s.d, -s.b) and lt(s.f, 0))):
      s.a5_action()
      return 00000
    # A6
    if (gt(abs(s.e), s.a)
        or (eq(s.e, s.a) and lt(2*s.d, s.f))
        or (eq(s.e, -s.a) and lt(s.f, 0))):
      s.a6_action()
      return 00000
    # A7
    if (gt(abs(s.f), s.a)
        or (eq(s.f, s.a) and lt(2*s.d, s.e))
        or (eq(s.f, -s.a) and lt(s.e, 0))):
      s.a7_action()
      return 00000
    # A8
    if (lt(s.d+s.e+s.f+s.a+s.b, 0)
        or (eq(s.d+s.e+s.f+s.a+s.b, 0) and gt(2*(s.a+s.e)+s.f, 0))):
      s.a8_action()
      return 00000
    return 0001

  def cb_update(self, m_elems):
    self._r_inv *= matrix.sqr(m_elems)

  def sign(self, x):
    assert abs(x) > self.epsilon
    if (x > 0): return 1
    return -1

  def a1_action(self):
    self.a, self.b = self.b, self.a
    self.d, self.e = self.e, self.d
    self.cb_update((0,-1,0, -1,0,0, 0,0,-1))

  def a2_action(self):
    self.b, self.c = self.c, self.b
    self.e, self.f = self.f, self.e
    self.cb_update((-1,0,0, 0,0,-1, 0,-1,0))

  def a3_action(self):
    lt = self.eps_lt
    i,j,k = 1,1,1
    if (lt(self.d, 0)): i = -1
    if (lt(self.e, 0)): j = -1
    if (lt(self.f, 0)): k = -1
    self.d = abs(self.d)
    self.e = abs(self.e)
    self.f = abs(self.f)
    self.cb_update((i,0,0, 0,j,0, 0,0,k))

  def a4_action(self):
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
    self.d = -abs(self.d)
    self.e = -abs(self.e)
    self.f = -abs(self.f)
    self.cb_update((f[0],0,0, 0,f[1],0, 0,0,f[2]))

  def a5_action(self):
    j = self.sign(self.d)
    self.c += self.b - j*self.d
    self.d -= 2*j*self.b
    self.e -= j*self.f
    self.cb_update(((1,0,0,0,1,-j,0,0,1)))

  def a6_action(self):
    j = self.sign(self.e)
    self.c += self.a - j*self.e
    self.d -= j*self.f
    self.e -= 2*j*self.a
    self.cb_update(((1,0,-j,0,1,0,0,0,1)))

  def a7_action(self):
    j = self.sign(self.f)
    self.b += self.a - j*self.f
    self.d -= j*self.e
    self.f -= 2*j*self.a
    self.cb_update(((1,-j,0,0,1,0,0,0,1)))

  def a8_action(self):
    self.c += self.a+self.b+self.d+self.e+self.f
    self.d += 2*self.b+self.f
    self.e += 2*self.a+self.f
    self.cb_update((1,0,1,0,1,1,0,0,1))
