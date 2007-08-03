from scitbx.array_family import flex
from scitbx.sparse import vector, matrix
from libtbx.test_utils import approx_equal
from random import uniform, randint

def exercise_vector():
  v = vector(5)
  assert v.size == 5
  assert v.is_structurally_zero()
  v[1] = 2
  v[2] = 0
  v[3] = 6
  assert [ v[i] for i in xrange(5) ] == [0, 2, 0, 6, 0]
  assert list(v) == [(1,2.), (2,0.), (3,6.)]
  p = flex.size_t([1,2,3,4,0])
  assert list(v.sort_indices().permute(p)) == [(2,2.), (3,0.), (4,6.)]
  v = vector(10)
  v[7] = -5
  v[1] = -1
  v[4] = 0
  v[1] = 2
  v[9] = 9.
  v[7] = 6
  v[4] = 1
  v[1] = 3
  v[4] = 0
  assert list(v.sort_indices()) == [(1,3.), (4,0.), (7,6.), (9,9.)]

def exercise_matrix():
  a = matrix(10,7)
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

def exercise_matrix_x_vector():
  for m,n in [(5,5), (3,5), (5,3)]:
    for n_test in xrange(50):
      a = matrix(m,n)
      x = vector(n)
      seen = set()
      for k in xrange(randint(3,2*m*n//3)):
        while True:
          i = randint(0,m-1)
          j = randint(0,n-1)
          if (i,j) not in seen: break
        seen.add((i,j))
        val = uniform(-3., 3.)
        a[i,j] = val
      seen = set()
      for k in xrange(randint(1,2*n//3)):
        while True:
          i = randint(0,n-1)
          if i not in seen: break
        seen.add(i)
        val = uniform(-2., 2.)
        x[i] = val
      y = a*x
      aa = a.as_dense_matrix()
      xx = x.as_dense_vector()
      yy1 = y.as_dense_vector()
      yy2 = aa.matrix_multiply(xx)
      assert approx_equal(yy1,yy2)

def run():
  exercise_vector()
  exercise_matrix()
  exercise_matrix_x_vector()
  print 'OK'

if __name__ == '__main__':
  run()
