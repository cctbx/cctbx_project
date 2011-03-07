def exercise_rational():
  from scitbx.matrix import row_echelon
  from scitbx import matrix
  from libtbx.utils import flat_list
  from boost import rational
  import random
  rng = random.Random(0)
  #
  m = [[0]]
  t = [0]
  free_vars = row_echelon.form_rational(m, t)
  assert m == [[0]]
  assert t == [0]
  assert free_vars == [0]
  sol = row_echelon.back_substitution_rational(m, t, free_vars, [1])
  assert sol == [1]
  sol = row_echelon.back_substitution_rational(m, None, free_vars, [2])
  assert sol == [2]
  #
  m = [[0]]
  t = [1]
  free_vars = row_echelon.form_rational(m, t)
  assert m == [[0]]
  assert t == [1]
  assert free_vars == [0]
  sol = row_echelon.back_substitution_rational(m, t, free_vars, [1])
  assert sol is None
  #
  m = [[1]]
  t = [2]
  free_vars = row_echelon.form_rational(m, t)
  assert m == [[1]]
  assert t == [2]
  assert free_vars == []
  sol = row_echelon.back_substitution_rational(m, t, free_vars, [1])
  assert sol == [2]
  #
  def rr():
    return rational.int(rng.randrange(-5,6), rng.randrange(1,10))
  #
  for i_trial in xrange(10):
    for nr in [1,2,3]:
      for nc in [1,2,3]:
        m = []
        for ir in xrange(nr):
          m.append([rr() for ic in xrange(nc)])
        m_orig = matrix.rec(flat_list(m), (nr,nc))
        sol_orig = [rr() for ic in xrange(nc)]
        t_orig = list(m_orig * matrix.col(sol_orig))
        t = list(t_orig)
        free_vars = row_echelon.form_rational(m, t)
        sol = [None] * nc
        for ic in free_vars:
          sol[ic] = sol_orig[ic]
        sol = row_echelon.back_substitution_rational(m, t, free_vars, sol)
        assert sol is not None
        assert sol.count(None) == 0
        assert sol == sol_orig
        sol = [1] * nc
        sol = row_echelon.back_substitution_rational(m, None, free_vars, sol)
        assert sol is not None
        assert (m_orig * matrix.col(sol)).dot() == 0
  #
  for i_trial in xrange(10):
    from itertools import count
    for i in count(10):
      a = matrix.col([rr(), rr(), rr()])
      b = matrix.col([rr(), rr(), rr()])
      if (a.cross(b).dot() != 0):
        break
    else:
      raise RuntimeError
    p = rng.randrange(-5,6)
    q = rng.randrange(-5,6)
    def check(a, b, c, expected_free_vars, expected_sol):
      m = []
      t = []
      for i in xrange(3):
        m.append([a[i], b[i]])
        t.append(c[i])
      m_orig = matrix.rec(flat_list(m), (3,2))
      t_orig = list(t)
      free_vars = row_echelon.form_rational(m, t)
      assert free_vars == expected_free_vars
      sol = row_echelon.back_substitution_rational(m, t, free_vars, [3, 11])
      assert sol == expected_sol
      if (sol is not None):
        assert list(m_orig * sol) == t_orig
    check(a, b, p*a+q*b, [], [p,q])
    check(a, b, a.cross(b), [], None)
    check(a, 5*a, -7*a, [1], [-62,11])
    check(a, 5*a, b, [1], None)
    check([0,0,0], [0,0,0], [0,0,0], [0,1], [3,11])

def run(args):
  assert len(args) == 0
  exercise_rational()
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
