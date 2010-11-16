from scitbx import sparse
from libtbx.test_utils import approx_equal

class test(object):

  def __init__(self, lu_factorization_type):
    self.lu_factorization_type = lu_factorization_type

  def __call__(self, a, **kwds):
    lu = self.lu_factorization_type(a, **kwds)
    assert lu.l().is_unit_lower_triangular()
    assert lu.u().is_upper_triangular()
    a1 = lu.l() * lu.u()
    a2 = lu.factored().clone().permute_rows(lu.rows_permutation())
    assert approx_equal(a1.as_dense_matrix(), a2.as_dense_matrix())


def exercise_gilbert_peierls_lu_factorization():
  check = test(sparse.gilbert_peierls_lu_factorization)

  """ mathematica
        a := SparseArray[ {
              {3, j_}  -> 1.5 - j/5,
              {7, j_}  -> -0.8 + j/5,
              {_, 5}   -> 2.1,
              {i_, i_} -> i
              }, {8, 8} ]
  """
  a = sparse.matrix(8,8)
  for j,c in enumerate(a.cols()):
    j += 1
    for i in xrange(a.n_rows):
      i += 1
      if i == 3:
          c[i-1] = 1.5 - j/5.
      elif i == 7:
          c[i-1] = -0.8 + j/5.
      else:
          if    j == 5: c[i-1] = 2.1
          elif  i == j: c[i-1] = i
  check(a)

  """ rectangular matrix m x n with m < n """
  b = sparse.matrix(5,8)
  b[4,0] = 1.
  b[1,1] = -1.
  b[1,2] = 0.5
  b[2,1] = 1.8
  b[2,2] = -2.
  b[0,3] = 1.;
  b[2,4] = -1.
  b[2,5] = 1.
  b[2,6] = 0.5
  b[3,4] = 1.
  b[3,5] = 0.5
  b[3,6] = 1.
  b[0,7] = 0.1
  b[1,7] = 0.2
  check(b);

  """ rectangular matrix m x n with m > n """
  c = b.transpose()
  check(c)


def run():
  exercise_gilbert_peierls_lu_factorization()
  print 'OK'

if __name__ == '__main__':
  run()
