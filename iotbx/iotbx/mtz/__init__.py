from cctbx.uctbx import unit_cell

from scitbx.python_utils.misc import import_regular_symbols
from iotbx_boost import mtz as ext
import_regular_symbols(globals(), ext.__dict__)
del import_regular_symbols

from iotbx.mtz import writer
from cctbx import sgtbx

column_type_legend_source = \
  "http://www.ccp4.ac.uk/dist/html/mtzlib.html#fileformat"
column_type_legend = {
  "H": "index h,k,l",
  "J": "intensity",
  "F": "structure factor amplitude",
  "D": "anomalous difference",
  "Q": "standard deviation",
  "G": "F(+) or F(-)",
  "L": "standard deviation",
  "K": "I(+) or I(-)",
  "M": "standard deviation",
  "P": "phase angle in degrees",
  "W": "weight (of some sort)",
  "A": "phase probability coefficients (Hendrickson/Lattmann)",
  "B": "BATCH number",
  "Y": "M/ISYM, packed partial/reject flag and symmetry number",
  "I": "integer",
  "R": "real",
}

class Mtz (ext.Mtz):
  def __init__(self,s):
    ext.Mtz.__init__(self,s)

  def __getattr__(self,s):
    assert type(s) == type(str())
    return self.getColumn(s)

  def label_to_crystal(self, label):
    assert label in self.columns()
    for i in xrange(self.ncrystals()):
      cryst = self.getCrystal(i)
      for j in xrange(cryst.ndatasets()):
        data = cryst.getDataset(j)
        for k in xrange(data.ncolumns()):
          if data.getColumn(k).label() == label:
            return cryst

  def get_space_group_info(self):
    return sgtbx.space_group_info(group=self.getSgtbxSpaceGroup())

MtzWriter.add_miller_array = writer.add_miller_array
MtzWriter.setSpaceGroup = writer.setSpaceGroup
