# Simple example for the use of the adptbx.

from cctbx import uctbx # unit cell toolbox
from cctbx import sgtbx # space group toolbox
from cctbx import adptbx # anisotropic displacement parameter toolbox

UnitCell = uctbx.UnitCell((10.67, 10.67, 4.68, 90, 90, 120))
SpaceGroup = sgtbx.SpaceGroup(sgtbx.SpaceGroupSymbols("P 3"))
Coordinates = ((0, 0, 0.236))
Uuvrs = ((0.17, 0.17, 0.19, 0.09, 0, 0))

spsp = sgtbx.SpecialPositionSnapParameters(UnitCell, SpaceGroup)
SiteSymmetry = sgtbx.SiteSymmetry(spsp, Coordinates)

print "Input Uuvrs:", Uuvrs
Ustar = adptbx.Uuvrs_as_Ustar(UnitCell, Uuvrs)
if (not SiteSymmetry.isCompatibleUstar(Ustar)):
  print "Warning: ADP tensor is incompatible with site symmetry."
Ustar = SiteSymmetry.AverageUstar(Ustar)
Uuvrs = adptbx.Ustar_as_Uuvrs(UnitCell, Ustar)
print "Averaged Uuvrs:", Uuvrs

Ucart = adptbx.Ustar_as_Ucart(UnitCell, Ustar)
Eigenvalues = adptbx.Eigenvalues(Ucart)
if (not adptbx.isPositiveDefinite(Eigenvalues)):
  print "ADP tensor is not positive definite."

print "Eigenvectors and values:"
Eigensystem = adptbx.Eigensystem(Ucart)
for i in xrange(3):
  print "  v=(%.5f %.5f %.5f) " % Eigensystem.vectors(i),
  print "lambda=%.4f" % (Eigensystem.values(i),)
