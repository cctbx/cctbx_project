from iotbx.scalepack import reader
from cctbx import crystal

def extract_from(file):
  scalepack_file = reader.scalepack_file(file, header_only=0001)
  assert scalepack_file.unit_cell is not None
  assert scalepack_file.space_group_info is not None
  return crystal.symmetry(
    unit_cell=scalepack_file.unit_cell,
    space_group_info=scalepack_file.space_group_info)
