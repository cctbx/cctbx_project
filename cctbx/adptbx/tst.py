from cctbx_boost import uctbx
from cctbx_boost import adptbx
from cctbx.development import debug_utils

def exercise_debye_waller():
  ucell = uctbx.UnitCell((5,7,9,80,100,130))
  assert abs(adptbx.DebyeWallerFactorUiso(ucell.Q((7, 8, 9)) / 4, 0.025)
             - 0.006625427) < 1.e-6
  assert abs(adptbx.DebyeWallerFactorUiso(ucell, (7, 8, 9), 0.025)
             - 0.006625427) < 1.e-6
  Uij = (0.1, 0.2, 0.3, -0.04, 0.05, -0.06)
  assert abs(adptbx.DebyeWallerFactorUcif(ucell, (1, 2, 3), Uij)
             - 0.335822877927) < 1.e-6

def check_eigenvalue(adp, x):
  r = -adp[0] - adp[1] - adp[2]
  s =   adp[0] * adp[1] + adp[0] * adp[2] + adp[1] * adp[2] \
      - adp[3] * adp[3] - adp[4] * adp[4] - adp[5] * adp[5]
  t =   adp[0] * adp[5] * adp[5] - adp[0] * adp[1] * adp[2] \
      + adp[2] * adp[3] * adp[3] + adp[1] * adp[4] * adp[4] \
      - 2. * adp[3] * adp[4] * adp[5]
  f = x**3 + r * x**2 + s * x + t
  assert abs(f) < 1.e-4, f

def exercise_eigen_core(diag):
  u = debug_utils.random_rotate_ellipsoid(diag + [0.,0.,0.])
  ev = list(adptbx.Eigenvalues(u))
  diag.sort()
  ev.sort()
  for i in xrange(3):
    check_eigenvalue(u, ev[i])
  for i in xrange(3):
    assert abs(diag[i] - ev[i]) < 1.e-4
  if (adptbx.isPositiveDefinite(ev)):
    es = adptbx.Eigensystem(u)
    ev = []
    for i in xrange(3): ev.append(es.values(i))
    ev.sort()
    for i in xrange(3):
      check_eigenvalue(u, ev[i])
    for i in xrange(3):
      assert abs(diag[i] - ev[i]) < 1.e-4

def exercise_eigen(n_trials=100):
  exercise_eigen_core([0,0,0])
  for n_equal in xrange(3):
    for i_trial in xrange(n_trials):
      if (n_equal == 0):
        diag = [debug_utils.random.random() for i in xrange(3)]
      if (n_equal == 1):
        i = debug_utils.random.randrange(3)
        diag[i] = debug_utils.random.random()
        diag[(i+1)%3] = debug_utils.random.random()
        diag[(i+2)%3] = diag[(i+1)%3]
      if (n_equal == 2):
        diag = [debug_utils.random.random()] * 3
      exercise_eigen_core(diag)

def run():
  exercise_debye_waller()
  exercise_eigen()
  print "OK"

if (__name__ == "__main__"):
  run()
