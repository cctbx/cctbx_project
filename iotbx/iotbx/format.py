def crystal_symmetry(cs):
  print """\
crystal.symmetry(
  unit_cell=(%.6g, %.6g, %.6g, %.6g, %.6g, %.6g),
  space_group_symbol='%s')""" % (
    cs.unit_cell().parameters()
    + (str(cs.space_group_info()),))
