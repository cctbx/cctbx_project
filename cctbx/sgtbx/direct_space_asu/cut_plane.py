from libtbx import slots_getstate_setstate
from scitbx import matrix
from boost import rational
import string

class cut_expr_ops(slots_getstate_setstate):

  __slots__ = []

  def __and__(self, other): return cut_expression("&", self, other)
  def __or__(self, other): return cut_expression("|", self, other)

class cut_expression(cut_expr_ops):

  __slots__ = ["op", "lhs", "rhs"]

  def __init__(self, op, lhs, rhs):
    self.op = op
    self.lhs = lhs
    self.rhs = rhs

  def __repr__(self):
    if (self.op == "&"):
      if (isinstance(self.lhs, cut) or self.lhs.op == "&"):
        lhs = str(self.lhs)
      else:
        lhs = "(" + str(self.lhs) + ")"
      if (isinstance(self.rhs, cut) or self.rhs.op == "&"):
        return lhs + " & " + str(self.rhs)
      return lhs + " & (" + str(self.rhs) + ")"
    if (self.op == "|"):
      return str(self.lhs) + " | " + str(self.rhs)
    raise RuntimeError

  def as_xyz(self):
    if (self.op == "&"):
      if (isinstance(self.lhs, cut) or self.lhs.op == "&"):
        lhs = self.lhs.as_xyz()
      else:
        lhs = "(" + self.lhs.as_xyz() + ")"
      if (isinstance(self.rhs, cut) or self.rhs.op == "&"):
        return lhs + " & " + self.rhs.as_xyz()
      return lhs + " & (" + self.rhs.as_xyz() + ")"
    if (self.op == "|"):
      return self.lhs.as_xyz() + " | " + self.rhs.as_xyz()
    raise RuntimeError

  def is_inside(self, point):
    if (self.op == "&"):
      return self.lhs.is_inside(point) and self.rhs.is_inside(point)
    if (self.op == "|"):
      return self.lhs.is_inside(point) or self.rhs.is_inside(point)
    raise RuntimeError

  def extract_all_cuts(self, result):
    self.lhs.extract_all_cuts(result)
    self.rhs.extract_all_cuts(result)

  def change_basis(self, cb_op):
    return cut_expression(
      self.op,
      self.lhs.change_basis(cb_op),
      self.rhs.change_basis(cb_op))

class cut(cut_expr_ops):

  __slots__ = ["n", "c", "inclusive", "cut_expr"]

  def __init__(self, n, c, inclusive=True, cut_expr=None):
    assert inclusive in (True, False)
    assert cut_expr is None or isinstance(cut_expr, cut_expr_ops)
    self.n = tuple(n)
    self.c = c
    self.inclusive = inclusive
    self.cut_expr = cut_expr

  def __repr__(self):
    s = self.base_symbol()
    if (not self.inclusive): s = "+" + s
    if (self.has_cuts()):
      if (("*" in s or "/" in s) and s[-1] in string.digits):
        s = "(" + s + ")"
      s += "(" + str(self.cut_expr) + ")"
    return s

  def as_xyz(self):
    n_negative = 0
    n_non_zero = 0
    for n in self.n:
      if (n < 0): n_negative += 1
      if (n != 0): n_non_zero += 1
    if (n_non_zero == 1):
      if (n_negative == 0):
        f = 1
      else:
        f = -1
    else:
      if (self.c > 0): n_negative += 1
      if (n_negative <= n_non_zero//2):
        f = 1
      else:
        f = -1
    s = ""
    for n,x in zip(self.n, "xyz"):
      nf = n * f
      if (nf == 0): continue
      if (nf > 0): s += "+"
      if (abs(nf) != 1):
        s += "%s*%s" % (str(nf), x)
      else:
        if (nf < 0): s += "-"
        s += x
    if (s[0] == "+"):
      s = s[1:]
    if (f > 0):
      s += ">"
    else:
      s += "<"
    if (self.inclusive):
      s += "="
    s += str(-self.c*f)
    if (self.has_cuts()):
      s += " [" + self.cut_expr.as_xyz() + "]"
    return s

  def base_symbol(self):
    from cctbx.sgtbx.direct_space_asu import short_cuts
    n = tuple(self.n)
    minus_n = tuple([-e for e in n])
    matching_n = None
    for key,value in short_cuts.__dict__.items():
      if (isinstance(value, cut)):
        if (value.n == n):
          if (value.c == self.c):
            return key
          elif (value.c == -self.c):
            return "-~"+key
          elif (value.c == 1):
            assert matching_n is None
            if (self.c < 0):
              matching_n = "-~"+key
            else:
              matching_n = key
        elif (value.n == minus_n):
          if (value.c == -self.c):
            return "-"+key
          elif (value.c == self.c):
            return "~"+key
          elif (value.c == 1):
            assert matching_n is None
            if (self.c < 0):
              matching_n = "-"+key
            else:
              matching_n = "~"+key
    c = short_cuts.r1 * self.c
    num = c.numerator()
    abs_num = abs(num)
    den = c.denominator()
    if (matching_n is None):
      if (num == 0):
        s = "0"
      else:
        s = "r1"
        if (num < 0): s = "-"+s
        if (abs_num != 1): s += "*"+str(abs_num)
        if (den != 1): s += "/"+str(den)
      return "cut(" + str(self.n).replace(" ", "") + "," + s + ")"
    s = matching_n
    assert num != 0
    if (abs_num != 1): s += "*"+str(abs_num)
    if (den != 1): s += "/"+str(den)
    return s

  def __pos__(self):
    "unsets inclusive flag"
    assert self.inclusive == True
    assert self.cut_expr is None
    return cut(self.n, self.c, inclusive=False)

  def __neg__(self):
    "-n, -c: flips inside/outside"
    return cut(n=[-e for e in self.n], c=-self.c,
               inclusive=self.inclusive, cut_expr=self.cut_expr)

  def __invert__(self):
    "-n, c"
    return cut(n=[-e for e in self.n], c=self.c,
               inclusive=self.inclusive, cut_expr=self.cut_expr)

  def __mul__(self, other):
    assert isinstance(other, int)
    return cut(n=self.n, c=other*self.c,
               inclusive=self.inclusive, cut_expr=self.cut_expr)

  def __truediv__(self, other):
    assert isinstance(other, int)
    assert other != 0
    assert self.c != 0
    return cut(n=self.n, c=rational.int(1)*self.c/other,
               inclusive=self.inclusive, cut_expr=self.cut_expr)

  def __div__(self, other):
    return self.__truediv__(other)

  def one(self):
    return cut(n=self.n, c=1,
               inclusive=self.inclusive, cut_expr=self.cut_expr)

  def __call__(self, expr):
    assert self.inclusive == True
    assert self.cut_expr is None
    return cut(self.n, self.c, cut_expr=expr)

  def has_cuts(self):
    return self.cut_expr is not None

  def evaluate(self, point):
    result = self.c
    for i in xrange(3):
      result += self.n[i] * point[i]
    return result

  def is_inside(self, point):
    i = self.evaluate(point)
    if (i < 0): return False
    if (i > 0): return True
    if (not self.has_cuts()): return self.inclusive
    return self.cut_expr.is_inside(point)

  def extract_all_cuts(self, result):
    result.append(self)
    if (self.has_cuts()):
      self.cut_expr.extract_all_cuts(result)

  def strip(self):
    return cut(self.n, self.c)

  def change_basis(self, cb_op):
    r_inv_tp = cb_op.c_inv().r().as_rational().transpose()
    t = cb_op.c().t().as_rational()
    np = r_inv_tp * matrix.col(self.n)
    cp = self.c - np.dot(t)
    if (self.has_cuts()):
      return cut(np.elems, cp, inclusive=self.inclusive,
        cut_expr=self.cut_expr.change_basis(cb_op))
    return cut(np.elems, cp, inclusive=self.inclusive)

  def lcm_of_denominators(self, start_lcm=1):
    result = start_lcm
    for e in self.n:
      result = rational.lcm(result, rational.int(e).denominator())
    result = rational.lcm(result, rational.int(self.c).denominator())
    return result

  def get_point_in_plane(self):
    result = [0,0,0]
    for i in xrange(3):
      if (self.n[i] != 0):
        result[i] = -rational.int(self.c) / self.n[i]
        return result
    raise RuntimeError("cut_plane normal vector is the null vector.")

  def as_float_cut_plane(self):
    from cctbx.crystal.direct_space_asu import float_cut_plane
    result = float_cut_plane(n=[float(e) for e in self.n], c=float(self.c))
    return result
