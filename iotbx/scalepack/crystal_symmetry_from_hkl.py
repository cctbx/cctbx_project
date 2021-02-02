from __future__ import absolute_import, division, print_function
from iotbx.scalepack import merge
from iotbx.scalepack import no_merge_original_index
from cctbx import crystal

def extract_from_merge(file):
  return merge.reader(file, header_only=True).crystal_symmetry()

def extract_from_no_merge_original_index(file_name):
  scalepack_file = no_merge_original_index.reader(file_name, header_only=True)
  return crystal.symmetry(
    unit_cell=None,
    space_group_info=scalepack_file.space_group_info())

def extract_from(file_name):
  try:
    with open(file_name) as f:
      result = extract_from_merge(f)
    return result
  except merge.FormatError: pass
  return extract_from_no_merge_original_index(file_name)
