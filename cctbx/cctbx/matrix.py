import math

class rec:

  def __init__(self, elems, n):
    assert len(n) == 2
    assert len(elems) == n[0] * n[1]
    self.elems = tuple(elems)
    self.n = tuple(n)

  def n_rows(self):
    return self.n[0]

  def n_columns(self):
    return self.n[1]

  def __neg__(self):
    return rec([-e for e in self.elems], self.n)

  def __add__(self, other):
    assert self.n == other.n
    a = self.elems
    b = other.elems
    return rec([a[i] + b[i] for i in xrange(len(a))], self.n)

  def __sub__(self, other):
    assert self.n == other.n
    a = self.elems
    b = other.elems
    return rec([a[i] - b[i] for i in xrange(len(a))], self.n)

  def __mul__(self, other):
    if (not hasattr(other, "elems")):
      return rec([x * other for x in self.elems], self.n)
    a = self.elems
    ar = self.n_rows()
    ac = self.n_columns()
    b = other.elems
    assert other.n_rows() == ac, "Incompatible matrices."
    bc = other.n_columns()
    result = []
    for i in xrange(ar):
      for k in xrange(bc):
        s = 0
        for j in xrange(ac):
          s += a[i * ac + j] * b[j * bc + k]
        result.append(s)
    if (ar == bc):
      return sqr(result)
    return rec(result, (ar, bc))

  def __div__(self, other):
    return rec([e/other for e in self.elems], self.n)

  def __call__(self, ir, ic):
    return self.elems[ir * self.n_columns() + ic]

  def norm(self):
    assert self.n_rows() == 1 or self.n_columns() == 1
    result = 0
    for e in self.elems:
      result += e*e
    return result

  def __abs__(self):
    return math.sqrt(self.norm())

  def dot(self, other):
    assert self.n_rows() == 1 or self.n_columns() == 1
    assert other.n_rows() == 1 or other.n_columns() == 1
    result = 0
    for i in xrange(len(self.elems)):
      result += self.elems[i] * other.elems[i]
    return result

  def cross(self, other):
    assert self.n in ((3,1), (1,3))
    assert self.n == other.n
    a = self.elems
    b = other.elems
    return col((
      a[1] * b[2] - b[1] * a[2],
      a[2] * b[0] - b[2] * a[0],
      a[0] * b[1] - b[0] * a[1]))

  def determinant(self):
    assert self.n == (3,3)
    m = self.elems
    return   m[0] * (m[4] * m[8] - m[5] * m[7]) \
           - m[1] * (m[3] * m[8] - m[5] * m[6]) \
           + m[2] * (m[3] * m[7] - m[4] * m[6])

  def co_factor_matrix_transposed(self):
    assert self.n == (3,3)
    m = self.elems
    return sqr((
       m[4] * m[8] - m[5] * m[7],
      -m[1] * m[8] + m[2] * m[7],
       m[1] * m[5] - m[2] * m[4],
      -m[3] * m[8] + m[5] * m[6],
       m[0] * m[8] - m[2] * m[6],
      -m[0] * m[5] + m[2] * m[3],
       m[3] * m[7] - m[4] * m[6],
      -m[0] * m[7] + m[1] * m[6],
       m[0] * m[4] - m[1] * m[3]))

  def inverse(self):
    determinant = self.determinant()
    assert determinant != 0
    return self.co_factor_matrix_transposed() / determinant

  def transpose(self):
    elems = []
    for j in xrange(self.n_columns()):
      for i in xrange(self.n_rows()):
        elems.append(self(i,j))
    return rec(elems, (self.n_columns(), self.n_rows()))

  def mathematica_form(self, label=""):
    s = ""
    if (label): s = label + "="
    s += "{"
    for ir in xrange(self.n_rows()):
      s += "{"
      for ic in xrange(self.n_columns()):
        s += str(self(ir, ic))
        s += ", "
      s = s[:-2] + "}, "
    return s[:-2] + "}"

class row(rec):

  def __init__(self, elems):
    rec.__init__(self, elems, (1, len(elems)))

class col(rec):

  def __init__(self, elems):
    rec.__init__(self, elems, (len(elems), 1))

class sqr(rec):

  def __init__(self, elems):
    l = len(elems)
    n = int(l**(.5) + 0.5)
    assert l == n * n
    rec.__init__(self, elems, (n,n))

class sym(rec):

  def __init__(self, elems):
    l = len(elems)
    n = int((-1 + (1+8*l)**(.5))/2. + 0.5)
    assert 2 * l == (n**2 + n), "Wrong number of unique elements."
    assert n == 3, "Not implemented."
    rec.__init__(self, (elems[0], elems[3], elems[4],
                        elems[3], elems[1], elems[5],
                        elems[4], elems[5], elems[2]), (n,n))

class rt:

  def __init__(self, tuple_r_t):
    if (hasattr(tuple_r_t[0], "elems")):
      self.r = sqr(tuple_r_t[0].elems)
    else:
      self.r = sqr(tuple_r_t[0])
    if (hasattr(tuple_r_t[1], "elems")):
      self.t = col(tuple_r_t[1].elems)
    else:
      self.t = col(tuple_r_t[1])
    assert self.r.n_rows() == self.t.n_rows()

  def __add__(self, other):
    return rt((self.r, self.t + other))

  def __sub__(self, other):
    return rt((self.r, self.t - other))

  def __mul__(self, other):
    try: return rt((self.r * other.r, self.r * other.t + self.t))
    except: pass
    try: return rt((self.r * other, self.t))
    except: pass
    try: return self.r * other + self.t
    except: pass
    if (len(other) == 9):
      return rt((self.r * sqr(other), self.t))
    return self.r * col(other) + self.t

if (__name__ == "__main__"):
  a = rec(range(1,7), (3,2))
  b = rec(range(1,7), (2,3))
  c = a * b
  d = rt((c, (1,2,3)))
  assert (-a).mathematica_form() == "{{-1, -2}, {-3, -4}, {-5, -6}}"
  assert d.r.mathematica_form() == "{{9, 12, 15}, {19, 26, 33}, {29, 40, 51}}"
  assert d.t.mathematica_form() == "{{1}, {2}, {3}}"
  e = d + col((3,5,6))
  assert e.r.mathematica_form() == "{{9, 12, 15}, {19, 26, 33}, {29, 40, 51}}"
  assert e.t.mathematica_form() == "{{4}, {7}, {9}}"
  f = e - col((1,2,3))
  assert f.r.mathematica_form() == "{{9, 12, 15}, {19, 26, 33}, {29, 40, 51}}"
  assert f.t.mathematica_form() == "{{3}, {5}, {6}}"
  ar = range(1,10)
  at = range(1,4)
  br = range(11,20)
  bt = range(4,7)
  g = rt((ar,at)) * rt((br,bt))
  assert g.r.mathematica_form() == \
    "{{90, 96, 102}, {216, 231, 246}, {342, 366, 390}}"
  assert g.t.mathematica_form() == "{{33}, {79}, {125}}"
  rt = g.r.transpose()
  assert rt.mathematica_form() == \
    "{{90, 216, 342}, {96, 231, 366}, {102, 246, 390}}"
  rtt = rt.transpose()
  assert rtt.mathematica_form() == \
    "{{90, 96, 102}, {216, 231, 246}, {342, 366, 390}}"
  tt = g.t.transpose()
  assert tt.mathematica_form() == "{{33, 79, 125}}"
  ttt = tt.transpose()
  assert ttt.mathematica_form() == "{{33}, {79}, {125}}"
  m = sqr((7, 7, -4, 3, 1, -1, 15, 16, -9))
  mi = m.inverse()
  assert mi.mathematica_form() == "{{7, -1, -3}, {12, -3, -5}, {33, -7, -14}}"
  assert (m*mi).mathematica_form() == "{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}"
  assert (mi*m).mathematica_form() == "{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}"
  print "OK"
