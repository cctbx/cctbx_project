from __future__ import absolute_import, division, print_function
from cctbx import crystal
from cctbx import sgtbx
from cctbx.sgtbx import subgroups
from cctbx.sgtbx import lattice_symmetry
from cctbx.sgtbx import bravais_types
from cctbx.array_family import flex
from libtbx.utils import format_cpu_times
from six.moves import cStringIO as StringIO
import sys
from six.moves import range
from six.moves import zip

def exercise_quick():
  for space_group_symbol in ("P-1",
                             "P2/m",
                             "C2/m",
                             "Pmmm",
                             "Cmmm",
                             "Fmmm",
                             "Immm",
                             "P4/mmm",
                             "I4/mmm",
                             "R-3m",
                             "P6/mmm",
                             "Pm-3m",
                             "Im-3m",
                             "Fm-3m"):
    parent_group_info = sgtbx.space_group_info(space_group_symbol)
    non_centric = sgtbx.space_group()
    for i_ltr in range(parent_group_info.group().n_ltr()):
      for i_smx in range(parent_group_info.group().n_smx()):
        s = parent_group_info.group()(i_ltr,0,i_smx)
        non_centric.expand_smx(s)
    assert non_centric.f_inv() == 1
    assert non_centric.order_z() * 2 == parent_group_info.group().order_z()
    non_centric_info = sgtbx.space_group_info(group=non_centric)
    unit_cell = non_centric_info.any_compatible_unit_cell(volume=1000)
    crystal_symmetry = crystal.symmetry(
      unit_cell=unit_cell,
      space_group_info=non_centric_info)
    minimum_symmetry = crystal_symmetry.minimum_cell()
    lattice_group = lattice_symmetry.group(
      minimum_symmetry.unit_cell(), max_delta=0.5)
    lattice_group_info = sgtbx.space_group_info(group=lattice_group)
    assert lattice_group_info.group() == minimum_symmetry.space_group()
    subgrs = subgroups.subgroups(lattice_group_info).groups_parent_setting()
    for group in subgrs:
      subsym = crystal.symmetry(
        unit_cell=minimum_symmetry.unit_cell(),
        space_group=group,
        assert_is_compatible_unit_cell=False)
      assert subsym.unit_cell().is_similar_to(minimum_symmetry.unit_cell())
      assert lattice_symmetry.find_max_delta(
        reduced_cell=minimum_symmetry.unit_cell(),
        space_group=group) < 0.6
  minimum_symmetry = crystal.symmetry(
    unit_cell="106.04, 181.78, 110.12, 90, 90, 90",
    space_group_symbol="P 1").minimum_cell()
  for max_delta in range(10,100,10):
    lattice_group = lattice_symmetry.group(
      minimum_symmetry.unit_cell(), max_delta=max_delta)
    lattice_group_info = sgtbx.space_group_info(group=lattice_group)
    assert str(lattice_group_info) == "P 4 2 2"

def exercise_comprehensive(args):
  if ("--verbose" in args):
    out = sys.stdout
  else:
    out = StringIO()
  if ("--paranoid" in args):
    cb_range = 2
  else:
    cb_range = 1
  for symbol in bravais_types.acentric:
    print("bravais type:", symbol)
    sym = sgtbx.space_group_info(symbol=symbol) \
      .any_compatible_crystal_symmetry(volume=1000) \
      .niggli_cell()
    abc = list(sym.unit_cell().parameters()[:3])
    abc.sort()
    for cb_elements in flex.nested_loop([-cb_range]*9,[cb_range+1]*9):
      r = sgtbx.rot_mx(cb_elements)
      if (r.determinant() != 1): continue
      cb_op = sgtbx.change_of_basis_op(sgtbx.rt_mx(r))
      sym_cb = sym.change_basis(cb_op)
      abc_cb = list(sym_cb.unit_cell().parameters()[:3])
      abc_cb.sort()
      for x,y in zip(abc, abc_cb):
        assert y-x > -1.e-6
        if (y-x > 1.e-4): break
      else:
        print("cb_ob:", cb_op.c(), cb_elements, file=out)
        assert min(cb_elements) >= -1
        assert max(cb_elements) <= 1
        for s in sym_cb.space_group():
          assert s.r().den() == 1
          r_num = s.r().num()
          print("r:", r_num, file=out)
          assert min(r_num) >= -1
          assert max(r_num) <= 1
        for enforce in [False, True]:
          lattice_group = sgtbx.lattice_symmetry.group(
            reduced_cell=sym_cb.unit_cell(),
            max_delta=1.4,
            enforce_max_delta_for_generated_two_folds=enforce)
        assert lattice_group == sym_cb.space_group()
        sys.stdout.flush()
    print(file=out)

def run(args):
  exercise_quick()
  exercise_comprehensive(args)
  print(format_cpu_times())

if (__name__ == "__main__"):
  run(sys.argv[1:])
