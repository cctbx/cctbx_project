import scitbx.math
import scitbx.linalg
from scitbx import matrix
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
import random

def exercise_cholesky_decomposition():
  from scitbx.examples import immoptibox_ports
  immoptibox_ports.py_cholesky_decomposition \
    = immoptibox_ports.cholesky_decomposition
  immoptibox_ports.cholesky_decomposition \
    = exercise_scitbx_cholesky_decomposition
  immoptibox_ports.tst_flex_counts = 0
  immoptibox_ports.exercise_cholesky()
  immoptibox_ports.cholesky_decomposition \
    = immoptibox_ports.py_cholesky_decomposition
  assert immoptibox_ports.tst_flex_counts == 299
  del immoptibox_ports.tst_flex_counts

def exercise_scitbx_cholesky_decomposition(a):
  from scitbx.examples import immoptibox_ports
  c = immoptibox_ports.py_cholesky_decomposition(a)
  al = a.matrix_symmetric_as_packed_l()
  chol = scitbx.linalg.l_l_transpose_cholesky_decomposition_in_place(al)
  cl = al
  if (c is None):
    assert chol.failure
  else:
    assert approx_equal(cl, c.matrix_lower_triangle_as_packed_l())
    for i_trial in xrange(10):
      b = flex.random_double(size=a.focus()[0], factor=2)-1
      x = chol.solve(b)
      assert approx_equal(a.matrix_multiply(x), b)
  immoptibox_ports.tst_flex_counts += 1
  return c


def exercise_gill_murray_wright_cholesky_decomposition():

  def p_as_mx(p):
    n = len(p)
    m = flex.double(flex.grid(n,n))
    m.matrix_diagonal_set_in_place(1)
    for i in xrange(n):
      if p[i] != i:
        m.matrix_swap_rows_in_place(i, p[i])
    return matrix.sqr(m)

  def core(a):
    c = flex.double(a)
    c.resize(flex.grid(a.n))
    u = c.matrix_upper_triangle_as_packed_u()
    gwm = scitbx.linalg.gill_murray_wright_cholesky_decomposition_in_place(
      u,
      epsilon=1.e-8)
    assert gwm.epsilon == 1.e-8
    u = c.matrix_upper_triangle_as_packed_u()
    gwm = scitbx.linalg.gill_murray_wright_cholesky_decomposition_in_place(u)
    assert gwm.epsilon == scitbx.math.floating_point_epsilon_double_get()
    assert gwm.packed_u.id() == u.id()
    p, e = gwm.pivots, gwm.e
    r = matrix.sqr(u.matrix_packed_u_as_upper_triangle())
    if a.n != (0,0):
      rtr = r.transpose() * r
      pm = p_as_mx(p)
      papt = pm * a * pm.transpose()
      paept = papt + matrix.diag(e)
      delta_decomposition = scitbx.linalg.matrix_equality_ratio(paept, rtr)
      assert delta_decomposition < 10, delta_decomposition
      b = flex.random_double(size=a.n[0], factor=2)-1
      x = gwm.solve(b=b)
      px = pm * matrix.col(x)
      pb = pm * matrix.col(b)
      if 0:
        eigen = scitbx.linalg.eigensystem.real_symmetric(
          paept.as_flex_double_matrix())
        lambda_ = eigen.values()
        print "condition number: ", lambda_[0]/lambda_[-1]
      delta_solve = scitbx.linalg.matrix_cholesky_test_ratio(
        a=paept.as_flex_double_matrix(),
        x=flex.double(px),
        b=flex.double(pb),
        epsilon=gwm.epsilon)
      assert delta_solve < 10, delta_solve
    return p, e, r
  # empty matrix
  a = matrix.sqr([])
  p, e, r = core(a)
  assert p.size() == 0
  assert e.size() == 0
  assert len(r) == 0
  n_max = 15
  n_trials_per_n = 10
  # identity matrices
  for n in xrange(1,n_max+1):
    a = matrix.diag([1]*n)
    p, e, r = core(a)
    assert list(p) == range(n)
    assert approx_equal(e, [0]*n)
    assert approx_equal(r, a)
  # null matrices
  for n in xrange(1,n_max+1):
    a = matrix.sqr([0]*n*n)
    p, e, r = core(a)
    assert list(p) == range(n)
    assert list(e) == [scitbx.math.floating_point_epsilon_double_get()]*n
    for i in xrange(n):
      for j in xrange(n):
        if (i != j): r(i,j) == 0
        else: r(i,j) == r(0,0)
  # random semi-positive diagonal matrices
  for n in xrange(1,n_max+1):
    for i_trial in xrange(n_trials_per_n):
      a = matrix.diag(flex.random_double(size=n))
      p, e, r = core(a)
      assert approx_equal(e, [0]*n)
      for i in xrange(n):
        for j in xrange(n):
          if (i != j): approx_equal(r(i,j), 0)
  # random diagonal matrices
  for n in xrange(1,n_max+1):
    for i_trial in xrange(n_trials_per_n):
      a = matrix.diag(flex.random_double(size=n, factor=2)-1)
      p, e, r = core(a)
      for i in xrange(n):
        for j in xrange(n):
          if (i != j): approx_equal(r(i,j), 0)
  # random semi-positive definite matrices
  for n in xrange(1,n_max+1):
    for i_trial in xrange(n_trials_per_n):
      m = matrix.sqr(flex.random_double(size=n*n, factor=2)-1)
      a = m.transpose_multiply()
      p, e, r = core(a)
      assert approx_equal(e, [0]*n)
  # random matrices
  for n in xrange(1,n_max+1):
    size = n*(n+1)//2
    for i_trial in xrange(n_trials_per_n):
      a = (flex.random_double(size=size, factor=2)-1) \
            .matrix_packed_u_as_symmetric()
      core(matrix.sqr(a))
      a.matrix_diagonal_set_in_place(0)
      core(matrix.sqr(a))
  # J. Nocedal and S. Wright:
  # Numerical Optimization.
  # Springer, New York, 1999, pp. 145-150.
  for i in xrange(3):
    for j in xrange(3):
      a = flex.double([[4,2,1],[2,6,3],[1,3,-0.004]])
      a.matrix_swap_rows_in_place(i=i, j=j)
      a.matrix_swap_columns_in_place(i=i, j=j)
      p, e, r = core(matrix.sqr(a))
      if (i == 0 and j == 0):
        assert list(p) == [1,1,2] # swap row 0 and 1 and nothing else
      assert approx_equal(e, [0.0, 0.0, 3.008])
      assert approx_equal(r,
        [2.4494897427831779, 0.81649658092772592, 1.2247448713915889,
         0.0, 1.8257418583505538, 0.0,
         0.0, 0.0, 1.2263767773404712])

def run():
  exercise_cholesky_decomposition()
  exercise_gill_murray_wright_cholesky_decomposition()
  print 'OK'

if __name__ == '__main__':
  run()
