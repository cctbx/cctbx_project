import uctbx
import sgtbx
from eltbx.caasf_wk1995 import CAASF_WK1995
import sftbx
def OneCycle():
  UnitCell = uctbx.UnitCell((10.002,10.002,34.141,90.0,90.0,90.0))
  SgOps = sgtbx.SpaceGroup(sgtbx.SpaceGroupSymbols("P42/NCM:2"))
  MillerIndices = sftbx.BuildMillerIndices(UnitCell, SgOps, 5.0)
  Sites = sftbx.vector_of_XrayScatterer()
  for Label, Coordinates in (
  ("SI1",       (0.09714,   0.70886,   0.90221)),
  ("SI2",       (0.04299,   0.04299,   0.45595)),
  ("SI3",       (0.09257,   0.09257,   0.27500)),
  ("SI4",       (0.10979,   0.10979,   0.74879)),
  ("SI5",       (0.94191,   0.94191,   0.33573)),
  ("SI6",       (0.14346,   0.14346,   0.06303)),
  ("SI7",       (0.25000,   0.25000,   0.48773)),
  ("O1",        (0.00998,   0.82408,   0.92409)),
  ("O2",        (0.10591,   0.58184,   0.93180)),
  ("O3",        (0.01904,   0.66920,   0.86208)),
  ("O4",        (0.25000,   0.75000,   0.89153)),
  ("O5",        (0.15766,   0.15766,   0.45932)),
  ("O6",        (0.00000,   0.00000,   0.50000)),
  ("O7",        (0.05008,   0.05008,   0.31935)),
  ("O8",        (0.06170,   0.25056,   0.26863)),
  ("O9",        (0.00421,   0.00421,   0.24440)),
  ("O10",       (0.12932,   0.12932,   0.70176)),
  ("O11",       (0.15642,   0.15642,   0.01556)),
  ("O12",       (0.25000,   0.25000,   0.08135)),
  ):
    SF = CAASF_WK1995(Label)
    Site = sftbx.XrayScatterer(Label, SF, 0j, Coordinates, 1., 0.03)
    Site.DetermineMultiplicity(UnitCell, SgOps, 0.5)
    Sites.append(Site)
  for Site in Sites:
    print Site.Label(), Site.Coordinates()
  Fcalc = sftbx.StructureFactorVector(UnitCell, SgOps, MillerIndices, Sites)
  for i in xrange(len(MillerIndices)):
    print MillerIndices[i], Fcalc[i]

if (__name__ == "__main__"):
  import sys
  Endless = "--Endless" in sys.argv
  while 1:
    OneCycle()
    if (not Endless): break
