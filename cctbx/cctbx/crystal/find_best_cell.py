from cctbx import sgtbx
from cctbx.array_family import flex
from cctbx import matrix

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

  def __init__(self, input_symmetry, angular_tolerance=None):
    if (angular_tolerance is None):
      angular_tolerance = 3
    self._all_cells = []
    space_group_number = input_symmetry.space_group_info().type().number()
    if (space_group_number == 1):
      self._cb_op = input_symmetry.change_of_basis_op_to_niggli_cell()
      self._symmetry = input_symmetry.change_basis(self._cb_op)
      self._all_cells.append(self._symmetry)
      return
    if (space_group_number < 3 or space_group_number >= 75):
      self._cb_op = sgtbx.change_of_basis_op()
      self._symmetry = input_symmetry
      self._all_cells.append(self._symmetry)
      return
    standard_info = sgtbx.space_group_info(
      symbol=space_group_number,
      table_id="A1983")
    cb_op_inp_ref = input_symmetry.space_group_info().type().cb_op()
    cb_op_std_ref = standard_info.type().cb_op()
    cb_op_std_inp = cb_op_inp_ref.inverse() * cb_op_std_ref
    assert standard_info.group().change_basis(cb_op_std_inp) == input_symmetry.space_group()
    best_cell_parameters = input_symmetry.unit_cell().parameters()
    best_cb_op = sgtbx.change_of_basis_op()
    best_symmetry = input_symmetry
    if (space_group_number <= 15):
      two_fold_info = sgtbx.rot_mx_info(input_symmetry.space_group()(1).r())
      assert abs(two_fold_info.type()) == 2
      ev = list(two_fold_info.ev())
      assert ev.count(0) == 2
      unique_axis = ev.index(1) + 3
      for elems in flex.nested_loop([-2]*9,[2]*9,00000):
        m = matrix.rec(elems, (3,3))
        if (m.determinant() != 1): continue
        cb_op = sgtbx.change_of_basis_op(sgtbx.rt_mx(sgtbx.rot_mx(elems)))
        alt_symmetry = input_symmetry.change_basis(cb_op)
        if (alt_symmetry.space_group() == input_symmetry.space_group()):
          self._all_cells.append(alt_symmetry)
          alt_cell_parameters = alt_symmetry.unit_cell().parameters()
          cmp_result = cmp_monoclinic_cell_parameters(
            best_cell_parameters, alt_cell_parameters, unique_axis,
            angular_tolerance)
          if (cmp_result > 0):
            best_cell_parameters = alt_cell_parameters
            best_cb_op = cb_op
            best_symmetry = alt_symmetry
    else:
      assert not str(standard_info).endswith(":2")
      affine_group = sgtbx.space_group("P 4 3*").change_basis(
        cb_op_std_inp)
      for affine_s in affine_group:
        cb_op = sgtbx.change_of_basis_op(affine_s) \
          .new_denominators(best_cb_op)
        alt_symmetry = input_symmetry.change_basis(cb_op)
        if (alt_symmetry.space_group() == input_symmetry.space_group()):
          self._all_cells.append(alt_symmetry)
          alt_cell_parameters = alt_symmetry.unit_cell().parameters()
          cmp_result = cmp_orthorhombic_cell_parameters(
            best_cell_parameters, alt_cell_parameters)
          if (cmp_result > 0):
            best_cell_parameters = alt_cell_parameters
            best_cb_op = cb_op
            best_symmetry = alt_symmetry
    self._cb_op = best_cb_op
    self._symmetry = best_symmetry

  def cb_op(self):
    return self._cb_op

  def symmetry(self):
    return self._symmetry

  def all_cells(self):
    return self._all_cells

def exercise():
  from cctbx import crystal
  cb_op = sgtbx.change_of_basis_op("y,z,x")
  for space_group_number in [1] + range(3,76):
    sgi = sgtbx.space_group_info(symbol=space_group_number)
    uc = sgi.any_compatible_unit_cell(volume=1000)
    best = find_best_cell(
      crystal.symmetry(unit_cell=uc, space_group_info=sgi).change_basis(cb_op),
      angular_tolerance=3)
    best.symmetry().show_summary()
    print best.cb_op().c()
    print len(best.all_cells())
    print

if (__name__ == "__main__"):
  exercise()
