import iotbx
from iotbx.mtz import Mtz
from iotbx.PrettyPrint.tables import binTable

def mtzio(file):
  p = Mtz(file)

  unitcell = tuple(p.UnitCell(0))
  N = p.size()

  H = p.getColumn("H")
  K = p.getColumn("K")
  L = p.getColumn("L")
  Iz    = p.getColumn("I")
  SIGIz = p.getColumn("SIGI")
  XDETz = p.getColumn("XDET")
  YDETz = p.getColumn("YDET")
  BATCHz = p.getColumn("BATCH")

for x in xrange(10000):
  mtzio("test.mtz")
  print x
