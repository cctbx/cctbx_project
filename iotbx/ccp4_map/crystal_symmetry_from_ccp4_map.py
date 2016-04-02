from __future__ import division
import iotbx.ccp4_map
from cctbx import crystal
import sys

def extract_from(file_name):
  # XXX This is to hide stdout from ccp4io code (ccp4_cmap_open) about input
  # XXX being non-ccp4 map. That's why we never have io in c/c++ code!
  import cStringIO
  save_stdout = sys.stdout
  sys.stdout = cStringIO.StringIO()
  m = iotbx.ccp4_map.map_reader(file_name=file_name)
  sys.stdout = save_stdout
  ucp = m.unit_cell_parameters
  sgn = max(1, m.space_group_number)
  return crystal.symmetry(ucp, sgn)
