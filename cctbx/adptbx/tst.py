from cctbx_boost import uctbx
from cctbx_boost import adptbx

def check_eigenvalue(adp, x):
  r = -adp[0] - adp[1] - adp[2]
  s =   adp[0] * adp[1] + adp[0] * adp[2] + adp[1] * adp[2] \
      - adp[3] * adp[3] - adp[4] * adp[4] - adp[5] * adp[5]
  t =   adp[0] * adp[5] * adp[5] - adp[0] * adp[1] * adp[2] \
      + adp[2] * adp[3] * adp[3] + adp[1] * adp[4] * adp[4] \
      - 2. * adp[3] * adp[4] * adp[5]
  f = x**3 + r * x**2 + s * x + t
  assert abs(f) < 1.e-5, f

u = uctbx.UnitCell((5,7,9,80,100,130))
print adptbx.DebyeWallerFactorUiso(u.Q((7, 8, 9)) / 4, 0.025)
print adptbx.DebyeWallerFactorUiso(u, (7, 8, 9), 0.025)
Uij = (0.1, 0.2, 0.3, -0.04, 0.05, -0.06)
print adptbx.DebyeWallerFactorUcif(u, (7, 8, 9), Uij)
ev = adptbx.Eigenvalues(Uij)
print ev
# XXX for x in ev: check_eigenvalue(Uij, x)
ES = adptbx.Eigensystem(Uij)
for i in xrange(3):
  print ES.vectors(i), ES.values(i)
  check_eigenvalue(Uij, ES.values(i))
