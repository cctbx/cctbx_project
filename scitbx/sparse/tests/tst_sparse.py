from scitbx.array_family import flex
from scitbx import sparse
from libtbx.test_utils import approx_equal, Exception_expected
import random

def exercise_vector():
  v = sparse.vector(5)
  assert v.size == 5
  assert v.is_structurally_zero()
  v[1] = 2
  v[2] = 0
  v[3] = 6
  assert list(v) == [(1,2.), (2,0.), (3,6.)]
  assert list(v.compact()) == [(1,2.), (2,0.), (3,6.)]
  assert [ v[i] for i in xrange(5) ] == [0, 2, 0, 6, 0]
  p = flex.size_t([1,2,3,4,0])
  assert list(v.permute(p)) == [(2,2.), (3,0.), (4,6.)]

  v = sparse.vector(10)
  v[7] = -5
  v[1] = -1
  v[4] = 0
  v[1] = 2
  v[9] = 9.
  v[7] = 6
  v[4] = 1
  v[1] = 3
  v[4] = 0
  assert list(v.compact()) == [(1,3.), (4,0.), (7,6.), (9,9.)]

  v = sparse.vector(10)
  v[4] += 1
  v[5] += 2
  v[4] += 2
  v[5] = 1
  v[3] = 2
  v[5] += 3
  assert list(v.compact()) == [ (3,2.), (4,3.), (5,4.) ]

  v = sparse.vector(6)
  v[3] = 1
  v[2] = 1
  v[5] = 1
  assert v.size == 6
  v[7] = 1
  assert v[7] == 0
  assert v.size == 6

  u = flex.double((1, -1, 2, 0, -2))
  v = sparse.vector(5)
  v[0] = 10
  v[3] = 4
  v[4] = 5
  assert u*v == 0

  approx_equal = sparse.approx_equal(tolerance=0.1)

  u = sparse.vector(4)
  u[0] = 1.01
  v = sparse.vector(4)
  v[0] = 1.02
  v[3] = 0.001
  assert approx_equal(u,v)

  u = sparse.vector(4)
  v = sparse.vector(4)
  v[3] = 0.001
  assert approx_equal(u,v)

  u = sparse.vector(5, {3: 0.3, 1: 0.1})
  assert list(u.as_dense_vector()) == [ 0, 0.1, 0, 0.3, 0 ]

  try:
    sparse.vector(4, [1, 2, 3, 4])
    raise Exception_expected
  except Exception, e:
    assert e.__class__.__module__ == 'Boost.Python'
    assert e.__class__.__name__ == 'ArgumentError'

def exercise_matrix():
  a = sparse.matrix(10,7)
  assert a.n_rows == 10 and a.n_cols == 7
  for c in a.cols():
    assert c.is_structurally_zero()
  a[0,1] = 1.
  a[9,5] = 2.
  for i in xrange(10):
    for j in xrange(7):
      if (i,j) == (0,1): assert a[i,j] == 1.
      elif (i,j) == (9,5): assert a[i,j] == 2.
      else: assert a[i,j] == 0

  a = sparse.matrix(6, 3)
  assert a.n_rows == 6
  a[1,1] = 1.
  a[3,2] = 2.
  a[5,1] = 2.
  a[4,0] = 1.
  assert a.n_rows == 6
  a[7,0] = 1.
  assert a[7,0] == 0
  assert a.n_rows == 6

  a = sparse.matrix(4,3)
  a[0,1] = 1.01
  b = sparse.matrix(4,3)
  b[0,1] = 1.02
  b[3,2] = 0.001
  approx_equal = sparse.approx_equal(tolerance=0.1)
  assert approx_equal(a,b)

  m = 10
  a = sparse.matrix(m, 2)
  columns = ( sparse.vector(m, {1:0.1, 2:0.2}),
              sparse.vector(m, {4:0.4, 8:0.8}) )
  a[:,0], a[:, 1] = columns
  assert a[:,0], a[:,1] == columns

  a = sparse.matrix(10, 3,
                    elements_by_columns=[ { 1: 1, 4: 4, },
                                          { 0: -1, 8:8, },
                                          { 6: 6, 9: 9, } ])
  assert "\n%s" % a == """
{
{ 0, -1, 0 },
{ 1, 0, 0 },
{ 0, 0, 0 },
{ 0, 0, 0 },
{ 4, 0, 0 },
{ 0, 0, 0 },
{ 0, 0, 6 },
{ 0, 0, 0 },
{ 0, 8, 0 },
{ 0, 0, 9 }
}
"""


def random_sparse_vector(n):
  x = sparse.vector(n)
  seen = {}
  for k in xrange(random.randint(1,2*n//3)):
    while True:
      i = random.randint(0,n-1)
      if i not in seen: break
    seen[i] = True
    val = random.uniform(-2., 2.)
    x[i] = val
  return x

def random_sparse_matrix(m,n):
  a = sparse.matrix(m,n)
  seen = {}
  for k in xrange(random.randint(3,2*m*n//3)):
    while True:
      i = random.randint(0,m-1)
      j = random.randint(0,n-1)
      if (i,j) not in seen: break
    seen[(i,j)] = True
    val = random.uniform(-3., 3.)
    a[i,j] = val
  return a

def exercise_matrix_x_vector():
  for m,n in [(5,5), (3,5), (5,3)]:
    for n_test in xrange(50):
      a = random_sparse_matrix(m,n)
      x = random_sparse_vector(n)
      y = a*x
      aa = a.as_dense_matrix()
      xx = x.as_dense_vector()
      yy1 = y.as_dense_vector()
      yy2 = aa.matrix_multiply(xx)
      assert approx_equal(yy1,yy2)

def exercise_matrix_x_matrix():
  a,b = random_sparse_matrix(3,4), random_sparse_matrix(4,2)
  c = a*b
  aa, bb, cc = [ m.as_dense_matrix() for m in (a,b,c) ]
  cc1 = aa.matrix_multiply(bb)
  assert approx_equal(cc, cc1)

def exercise_a_tr_b_a():
  for m,n in [(5,5), (3,5), (5,3)]:
    for n_test in xrange(50):
      b = flex.random_double(m*(m+1)/2)
      a = random_sparse_matrix(m,n)
      c = a.self_transpose_times_symmetric_times_self(b)
      aa = a.as_dense_matrix()
      bb = b.matrix_packed_u_as_symmetric()
      cc = c.matrix_packed_u_as_symmetric()
      assert approx_equal(
        cc,
        aa.matrix_transpose().matrix_multiply(bb.matrix_multiply(aa)))

def exercise_row_vector_x_matrix():
  u = flex.double((1,2,3))
  a = sparse.matrix(3,5)
  a[1,0] = 1
  a[2,1] = 1
  a[0,2] = 1
  a[-1,3] = 1
  a[-2,4] = 1
  v = u*a
  assert list(v) == [ 2, 3, 1, -2, -6 ]

def run():
  exercise_vector()
  exercise_matrix()
  exercise_matrix_x_vector()
  exercise_matrix_x_matrix()
  exercise_a_tr_b_a()
  print 'OK'

if __name__ == '__main__':
  run()
