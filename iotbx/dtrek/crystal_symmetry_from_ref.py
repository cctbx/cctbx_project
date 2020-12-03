from __future__ import absolute_import, division, print_function
from iotbx.dtrek import reflnlist_reader

def extract_from(file_name):
  with open(file_name) as f:
    cs = reflnlist_reader.reflnlist(f, header_only=True).crystal_symmetry()
  return cs
