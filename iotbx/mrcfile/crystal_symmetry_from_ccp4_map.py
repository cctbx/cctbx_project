from __future__ import absolute_import, division, print_function
import iotbx.mrcfile
from cctbx import crystal

def extract_from(file_name):
  # XXX This is just to be parallel to ccp4_map.

  m = iotbx.mrcfile.map_reader(file_name=file_name,header_only=True)
  ucp = m.unit_cell().parameters()
  sgn = max(1, m.unit_cell_crystal_symmetry().space_group_number())
  return crystal.symmetry(ucp, sgn)
