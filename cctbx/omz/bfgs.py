from __future__ import division

def h0_scaling(sk, yk):
  "Nocedal & Wright (1999) Eq. 8.20; same as Eq. 9.6"
  yty = yk.dot(yk)
  assert yty > 0
  return yk.dot(sk) / yty

def h_update(hk, sk, yk):
  """Nocedal & Wright (1999) Eq. 9.2
     Pre-condition: hk must be positive definite"""
  yts = yk.dot(sk)
  assert yts > 0
  rho = 1 / yts
  v = -rho * yk.matrix_outer_product(sk)
  v.matrix_diagonal_add_in_place(1)
  hl = v.matrix_transpose_multiply(hk).matrix_multiply(v)
  hl += rho * sk.matrix_outer_product(sk)
  return hl

def b_update(bk, sk, yk):
  """Nocedal & Wright (1999) Eq. 8.19
     Pre-condition: bk must be positive definite"""
  yts = yk.dot(sk)
  assert yts > 0
  bsstb = bk.matrix_multiply(sk.matrix_outer_product(sk)).matrix_multiply(bk)
  stbs = sk.dot(bk.matrix_multiply(sk))
  assert stbs > 0
  yyt = yk.matrix_outer_product(yk)
  bl = bk - bsstb / stbs + yyt / yts
  return bl

class memory_element(object):

  __slots__ = ["s", "y", "rho"]

  def __init__(O, s, y):
    O.s = s
    O.y = y
    yts = y.dot(s)
    if (yts > 0):
      O.rho = 1 / yts
    else:
      O.rho = None

  def get(O):
    return O.s, O.y, O.rho

def hg_two_loop_recursion(memory, hk0, gk):
  """Nocedal & Wright (1999) Algorithm 9.1
     Pre-condition: hk0 must be positive definite"""
  q = gk.deep_copy()
  m = len(memory)
  alpha = [None] * m
  for i in xrange(m-1,-1,-1):
    s, y, rho = memory[i].get()
    alpha[i] = a = rho * s.dot(q)
    q = q - a * y
  if (hk0.is_square_matrix()):
    r = hk0.matrix_multiply(q)
  elif (hk0.is_trivial_1d()):
    r = hk0 * q
  else:
    raise ValueError("Improper hk0")
  for i in xrange(m):
    s, y, rho = memory[i].get()
    beta = rho * y.dot(r)
    r = r + s * (alpha[i] - beta)
  return r
