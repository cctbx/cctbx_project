# Simple example for the use of the adptbx.

from cctbx_boost import uctbx # unit cell toolbox
from cctbx_boost import sgtbx # space group toolbox
from cctbx_boost import adptbx # anisotropic displacement parameter toolbox

UnitCell = uctbx.UnitCell((10.67, 10.67, 4.68, 90, 90, 120))
SpaceGroup = sgtbx.SpaceGroup(sgtbx.SpaceGroupSymbols("P 3"))
Coordinates = ((0, 0, 0.236))
Ucif = ((0.17, 0.17, 0.19, 0.09, 0, 0))

spsp = sgtbx.SpecialPositionSnapParameters(UnitCell, SpaceGroup)
SiteSymmetry = sgtbx.SiteSymmetry(spsp, Coordinates)

print "Input Ucif:", Ucif
Ustar = adptbx.Ucif_as_Ustar(UnitCell, Ucif)
if (not SiteSymmetry.isCompatibleUstar(Ustar)):
  print "Warning: ADP tensor is incompatible with site symmetry."
Ustar = SiteSymmetry.AverageUstar(Ustar)
Ucif = adptbx.Ustar_as_Ucif(UnitCell, Ustar)
print "Averaged Ucif:", Ucif

Ucart = adptbx.Ustar_as_Ucart(UnitCell, Ustar)
Eigenvalues = adptbx.Eigenvalues(Ucart)
if (not adptbx.isPositiveDefinite(Eigenvalues)):
  print "ADP tensor is not positive definite."

print "Eigenvectors and values:"
Eigensystem = adptbx.Eigensystem(Ucart)
for i in xrange(3):
  print "  v=(%.5f %.5f %.5f) " % Eigensystem.vectors(i),
  print "lambda=%.4f" % (Eigensystem.values(i),)
