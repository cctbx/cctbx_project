"""Format crystal symmetry"""
from __future__ import absolute_import, division, print_function

def crystal_symmetry(cs):
  if (cs.unit_cell() is None):
    u = "None"
  else:
    u = "(%.6g, %.6g, %.6g, %.6g, %.6g, %.6g)" % cs.unit_cell().parameters()
  if (cs.space_group_info() is None):
    s = "None"
  else:
    s = "'%s'" % str(cs.space_group_info()).replace("'", "\\'")
  print("""\
crystal.symmetry(
  unit_cell=%s,
  space_group_symbol=%s)""" % (u, s))
