from cctbx import sgtbx

def cmp_orthorhombic_cell_parameters(lhs, rhs):
  for i in xrange(3):
    if (lhs[i] < rhs[i]): return -1
    if (lhs[i] > rhs[i]): return  1
  return 0

def cmp_monoclinic_cell_parameters(lhs, rhs, unique_axis, angular_tolerance):
  lhs_ang = lhs[unique_axis]
  rhs_ang = rhs[unique_axis]
  if (abs(lhs_ang - rhs_ang) < angular_tolerance):
    return cmp_orthorhombic_cell_parameters(lhs, rhs)
  lhs_ang_d90 = abs(lhs_ang - 90)
  rhs_ang_d90 = abs(rhs_ang - 90)
  if (abs(lhs_ang_d90 - rhs_ang_d90) > angular_tolerance):
    if (lhs_ang_d90 < rhs_ang_d90): return -1
    if (lhs_ang_d90 > rhs_ang_d90): return  1
  else:
    if (lhs_ang > 90 and rhs_ang < 90): return -1
    if (lhs_ang < 90 and rhs_ang > 90): return  1
  if (lhs_ang > rhs_ang): return -1
  if (lhs_ang < rhs_ang): return  1
  return 0

class find_best_cell:

  def __init__(self, input_symmetry, angular_tolerance):
    space_group_number = input_symmetry.space_group_info().type().number()
    if (space_group_number < 3 or space_group_number >= 75):
      self._cb_op = sgtbx.change_of_basis_op()
      self._symmetry = input_symmetry
      return
    standard_info = sgtbx.space_group_info(
      symbol=space_group_number,
      table_id="A1983")
    cb_op_inp_ref = input_symmetry.space_group_info().type().cb_op()
    cb_op_std_ref = standard_info.type().cb_op()
    cb_op_std_inp = cb_op_inp_ref.inverse() * cb_op_std_ref
    assert standard_info.group().change_basis(cb_op_std_inp) == input_symmetry.space_group()
    if (space_group_number <= 15):
      affine_hall_symbols = ("P 4 2 (y,z,x)", "P 6 2 (y,z,x)")
      two_fold_info = sgtbx.rot_mx_info(input_symmetry.space_group()(1).r())
      assert abs(two_fold_info.type()) == 2
      ev = list(two_fold_info.ev())
      assert ev.count(0) == 2
      unique_axis = ev.index(1) + 3
    else:
      assert not str(standard_info).endswith(":2")
      affine_hall_symbols = ("P 4 3*",)
      unique_axis = None
    best_cell_parameters = input_symmetry.unit_cell().parameters()
    best_cb_op = sgtbx.change_of_basis_op()
    best_symmetry = input_symmetry
    for affine_hall_symbol in affine_hall_symbols:
      affine_group = sgtbx.space_group(affine_hall_symbol).change_basis(
        cb_op_std_inp)
      for affine_s in affine_group:
        affine_cb_op = sgtbx.change_of_basis_op(affine_s) \
          .new_denominators(best_cb_op)
        alt_symmetry = input_symmetry.change_basis(affine_cb_op)
        if (alt_symmetry.space_group() == input_symmetry.space_group()):
          alt_cell_parameters = alt_symmetry.unit_cell().parameters()
          if (unique_axis is not None):
            cmp_result = cmp_monoclinic_cell_parameters(
              best_cell_parameters, alt_cell_parameters, unique_axis,
              angular_tolerance)
          else:
            cmp_result = cmp_orthorhombic_cell_parameters(
              best_cell_parameters, alt_cell_parameters)
          if (cmp_result > 0):
            best_cell_parameters = alt_cell_parameters
            best_cb_op = affine_cb_op
            best_symmetry = alt_symmetry
    self._cb_op = best_cb_op
    self._symmetry = best_symmetry

  def cb_op(self):
    return self._cb_op

  def symmetry(self):
    return self._symmetry

def exercise():
  from cctbx import crystal
  cb_op = sgtbx.change_of_basis_op("y,z,x")
  for space_group_number in xrange(3,76):
    sgi = sgtbx.space_group_info(symbol=space_group_number)
    uc = sgi.any_compatible_unit_cell(volume=1000)
    best = find_best_cell(
      crystal.symmetry(unit_cell=uc, space_group_info=sgi).change_basis(cb_op),
      angular_tolerance=3)
    best.symmetry().show_summary()
    print best.cb_op().c()
    print

if (__name__ == "__main__"):
  exercise()
