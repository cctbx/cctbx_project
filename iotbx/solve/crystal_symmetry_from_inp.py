from __future__ import absolute_import, division, print_function
from cctbx import crystal
from libtbx.utils import detect_binary_file

def extract_from(file_name=None, file=None, monitor_initial=None):
  assert [file_name, file].count(None) == 1
  if (file is None):
    file = open(file_name)
  lines = file.readlines()
  file.close()
  detect_binary = detect_binary_file(monitor_initial=monitor_initial)
  unit_cell = None
  space_group_symbol = None
  for line in lines:
    if (detect_binary is not None):
      is_binary = detect_binary.is_binary_file(block=line)
      if (is_binary is not None):
        if (is_binary): break
        detect_binary = None
    flds = line.split("!",1)[0].split()
    if (len(flds) > 0):
      keyword = flds[0].upper()
      if (keyword == "CELL"):
        assert len(flds) > 1
        unit_cell = " ".join(flds[1:])
        if (space_group_symbol is not None):
          break
      elif (keyword == "SYMFILE"):
        assert len(flds) > 1
        space_group_symbol = flds[1].replace("\\","/").split("/")[-1]
        if (space_group_symbol.lower().endswith(".sym")):
          space_group_symbol = space_group_symbol[:-4]
        if (unit_cell is not None):
          break
  assert [unit_cell, space_group_symbol].count(None) < 2
  return crystal.symmetry(
    unit_cell=unit_cell,
    space_group_symbol=space_group_symbol)
