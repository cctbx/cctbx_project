from __future__ import absolute_import, division, print_function
from iotbx.xds import read_ascii

def extract_from(file_name):
  with open(file_name) as f:
    cs = read_ascii.reader(f, header_only=True).crystal_symmetry()
  return cs
