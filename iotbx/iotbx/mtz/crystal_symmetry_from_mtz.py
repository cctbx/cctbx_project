"Extracts crystal symmetry from MTZ file."

from iotbx import mtz
from cctbx import crystal

def extract_from(file_name):
  mtz_file = mtz.Mtz(file_name)
  assert mtz_file.nsym() > 0
  assert mtz_file.ncrystals() > 0
  cryst = mtz_file.getCrystal(0)
  return crystal.symmetry(
    unit_cell=cryst.UnitCell(),
    space_group_info=mtz_file.get_space_group_info())
