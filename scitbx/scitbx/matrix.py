try:
  from stdlib import math
except:
  import math

class rec:

  def __init__(self, elems, n):
    assert len(n) == 2
    if (not isinstance(elems, tuple)):
      elems = tuple(elems)
    assert len(elems) == n[0] * n[1]
    self.elems = elems
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

  def __rmul__(self, other):
    "scalar * matrix"
    return self * other

  def __div__(self, other):
    return rec([e/other for e in self.elems], self.n)

  def __call__(self, ir, ic):
    return self.elems[ir * self.n_columns() + ic]

  def __len__(self):
    return len(self.elems)

  def __getitem__(self, i):
    return self.elems[i]

  def __float__(self):
    return rec([float(e) for e in self.elems], self.n)

  def each_abs(self):
    return rec([abs(e) for e in self.elems], self.n)

  def min(self):
    result = None
    for e in self.elems:
      if (result is None or result > e):
        result = e
    return result

  def max(self):
    result = None
    for e in self.elems:
      if (result is None or result < e):
        result = e
    return result

  def min_index(self):
    result = None
    for i in xrange(len(self.elems)):
      if (result is None or self.elems[result] > self.elems[i]):
        result = i
    return result

  def max_index(self):
    result = None
    for i in xrange(len(self.elems)):
      if (result is None or self.elems[result] < self.elems[i]):
        result = i
    return result

  def sum(self):
    result = 0
    for e in self.elems:
      result += e
    return result

  def product(self):
    result = 1
    for e in self.elems:
      result *= e
    return result

  def norm(self):
    assert self.n_rows() == 1 or self.n_columns() == 1
    result = 0
    for e in self.elems:
      result += e*e
    return result

  def __abs__(self):
    return math.sqrt(self.norm())

  def normalize(self):
    return self / abs(self)

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

  def cos_angle(self, other, value_if_undefined=None):
    self_norm = self.norm()
    if (self_norm == 0): return value_if_undefined
    other_norm = other.norm()
    if (other_norm == 0): return value_if_undefined
    d = self_norm * other_norm
    if (d == 0): return value_if_undefined
    return self.dot(other) / math.sqrt(d)

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

  def mathematica_form(self, label="", one_row_per_line=00000, format=None):
    s = ""
    if (label): s = label + "="
    s += "{"
    for ir in xrange(self.n_rows()):
      s += "{"
      for ic in xrange(self.n_columns()):
        if (format is None):
          s += str(self(ir, ic))
        else:
          s += format % self(ir, ic)
        s += ", "
      s = s[:-2] + "},"
      if (one_row_per_line): s += "\n  "
      else: s += " "
    return s[:-2] + "}"

  def as_sym_mat3(self):
    assert self.n == (3,3)
    m = self.elems
    return (m[0],m[4],m[8],
            (m[1]+m[3])/2.,
            (m[2]+m[6])/2.,
            (m[5]+m[7])/2.)

  def extract_block(self, stop, start=(0,0), step=(1,1)):
    assert 0 <= stop[0] <= self.n[0]
    assert 0 <= stop[1] <= self.n[1]
    i_rows = range(start[0], stop[0], step[0])
    i_colums = range(start[1], stop[1], step[1])
    result = []
    for ir in i_rows:
      for ic in i_colums:
        result.append(self(ir,ic))
    return rec(result, (len(i_rows),len(i_colums)))

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

class diag(rec):

  def __init__(self, diag_elems):
    n = len(diag_elems)
    elems = [0 for i in xrange(n*n)]
    for i in xrange(n):
      elems[i*(n+1)] = diag_elems[i]
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
    if (isinstance(other, rt)):
      return rt((self.r + other.r, self.t + other.t))
    else:
      return rt((self.r, self.t + other))

  def __sub__(self, other):
    if (isinstance(other, rt)):
      return rt((self.r - other.r, self.t - other.t))
    else:
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

  def inverse(self):
    r_inv = self.r.inverse()
    return rt((r_inv, -(r_inv*self.t)))

  def __float__(self):
    return rt((float(self.r), float(self.t)))

  def as_augmented_matrix(self):
    assert self.r.n_rows() == self.r.n_columns()
    n = self.r.n_rows()
    result = []
    for i_row in xrange(n):
      result.extend(self.r.elems[i_row*n:(i_row+1)*n])
      result.append(self.t[i_row])
    result.extend([0]*n)
    result.append(1)
    return rec(result, (n+1,n+1))

if (__name__ == "__main__"):
  from libtbx.test_utils import approx_equal
  from boost import rational
  a = rec(range(1,7), (3,2))
  assert len(a) == 6
  assert a[1] == 2
  assert (a*3).mathematica_form() == "{{3, 6}, {9, 12}, {15, 18}}"
  assert (-2*a).mathematica_form() == "{{-2, -4}, {-6, -8}, {-10, -12}}"
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
  e = e + f
  assert e.r.mathematica_form() \
      == "{{18, 24, 30}, {38, 52, 66}, {58, 80, 102}}"
  assert e.t.mathematica_form() == "{{7}, {12}, {15}}"
  f = f - e
  assert f.r.mathematica_form() \
      == "{{-9, -12, -15}, {-19, -26, -33}, {-29, -40, -51}}"
  assert f.t.mathematica_form() == "{{-4}, {-7}, {-9}}"
  e = float(f)*.5
  assert e.r.mathematica_form() \
      == "{{-4.5, -6.0, -7.5}, {-9.5, -13.0, -16.5}, {-14.5, -20.0, -25.5}}"
  assert e.t.mathematica_form() == "{{-4.0}, {-7.0}, {-9.0}}"
  a = f.as_augmented_matrix()
  assert a.mathematica_form() == "{{-9, -12, -15, -4}, {-19, -26, -33, -7}," \
                               + " {-29, -40, -51, -9}, {0, 0, 0, 1}}"
  assert a.extract_block(stop=(1,1)).mathematica_form() \
      == "{{-9}}"
  assert a.extract_block(stop=(2,2)).mathematica_form() \
      == "{{-9, -12}, {-19, -26}}"
  assert a.extract_block(stop=(3,3)).mathematica_form() \
      == "{{-9, -12, -15}, {-19, -26, -33}, {-29, -40, -51}}"
  assert a.extract_block(stop=(4,4)).mathematica_form() \
      == a.mathematica_form()
  assert a.extract_block(stop=(4,4),step=(2,2)).mathematica_form() \
      == "{{-9, -15}, {-29, -51}}"
  assert a.extract_block(start=(1,1),stop=(4,4),step=(2,2)).mathematica_form()\
      == "{{-26, -7}, {0, 1}}"
  assert a.extract_block(start=(1,0),stop=(4,3),step=(2,1)).mathematica_form()\
      == "{{-19, -26, -33}, {0, 0, 0}}"
  ar = range(1,10)
  at = range(1,4)
  br = range(11,20)
  bt = range(4,7)
  g = rt((ar,at)) * rt((br,bt))
  assert g.r.mathematica_form() == \
    "{{90, 96, 102}, {216, 231, 246}, {342, 366, 390}}"
  assert g.t.mathematica_form() == "{{33}, {79}, {125}}"
  grt = g.r.transpose()
  assert grt.mathematica_form() == \
    "{{90, 216, 342}, {96, 231, 366}, {102, 246, 390}}"
  grtt = grt.transpose()
  assert grtt.mathematica_form() == \
    "{{90, 96, 102}, {216, 231, 246}, {342, 366, 390}}"
  gtt = g.t.transpose()
  assert gtt.mathematica_form() == "{{33, 79, 125}}"
  gttt = gtt.transpose()
  assert gttt.mathematica_form() == "{{33}, {79}, {125}}"
  m = sqr((7, 7, -4, 3, 1, -1, 15, 16, -9))
  mi = m.inverse()
  assert mi.mathematica_form() == "{{7, -1, -3}, {12, -3, -5}, {33, -7, -14}}"
  assert (m*mi).mathematica_form() == "{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}"
  assert (mi*m).mathematica_form() == "{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}"
  s = rt((m, (1,-2,3)))
  si = s.inverse()
  assert si.r.mathematica_form() == mi.mathematica_form()
  assert (s*si).r.mathematica_form() == "{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}"
  assert (si*s).r.mathematica_form() == "{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}"
  assert si.t.mathematica_form() == "{{0}, {-3}, {-5}}"
  assert (s*si).t.mathematica_form() == "{{0}, {0}, {0}}"
  assert (si*s).t.mathematica_form() == "{{0}, {0}, {0}}"
  assert approx_equal(float(col((rational.int(3,4),2,1.5))).elems,(0.75,2,1.5))
  assert approx_equal(col((-2,3,-6)).normalize().elems, (-2/7.,3/7.,-6/7.))
  assert col((-1,2,-3)).each_abs().elems == (1,2,3)
  assert col((5,3,4)).min() == 3
  assert col((4,5,3)).max() == 5
  assert col((5,3,4)).min_index() == 1
  assert col((4,5,3)).max_index() == 1
  assert col((4,5,3)).sum() == 12
  assert col((2,3,4)).product() == 2*3*4
  assert diag((1,2,3)).mathematica_form() == \
    "{{1, 0, 0}, {0, 2, 0}, {0, 0, 3}}"
  assert approx_equal(col((1,0,0)).cos_angle(col((1,1,0)))**2, 0.5)
  print "OK"
