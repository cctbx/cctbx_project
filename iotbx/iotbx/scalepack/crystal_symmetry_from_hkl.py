from iotbx.scalepack import merge
from cctbx import crystal

def extract_from(file):
  scalepack_file = merge.reader(file, header_only=0001)
  return crystal.symmetry(
    unit_cell=scalepack_file.unit_cell,
    space_group_info=scalepack_file.space_group_info)
