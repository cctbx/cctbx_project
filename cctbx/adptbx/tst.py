from cctbx_boost import uctbx
from cctbx_boost import adptbx

u = uctbx.UnitCell((5,7,9,80,100,130))
print adptbx.DebyeWallerFactorUiso(u.Q((7, 8, 9)) / 4, 0.025)
print adptbx.DebyeWallerFactorUiso(u, (7, 8, 9), 0.025)
Uij = (0.1, 0.2, 0.3, -0.04, 0.05, -0.06)
print adptbx.DebyeWallerFactorUcif(u, (7, 8, 9), Uij)
print adptbx.Eigenvalues(Uij)
ES = adptbx.Eigensystem(Uij)
for i in xrange(3):
  print ES.vectors(i), ES.values(i)
