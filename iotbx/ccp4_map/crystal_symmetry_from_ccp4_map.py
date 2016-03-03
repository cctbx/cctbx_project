from __future__ import division
import iotbx.ccp4_map
from cctbx import crystal

def extract_from(file_name):
  m = iotbx.ccp4_map.map_reader(file_name=file_name)
  ucp = m.unit_cell_parameters
  sgn = max(1, m.space_group_number)
  return crystal.symmetry(ucp, sgn)
