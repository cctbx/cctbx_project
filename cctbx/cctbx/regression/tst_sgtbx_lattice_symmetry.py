from cctbx import crystal
from cctbx import sgtbx
from cctbx.sgtbx import subgroups
from cctbx.sgtbx import lattice_symmetry

def run():
  assert lattice_symmetry.group.n_potential_axes() == 2391
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
    for i_ltr in xrange(parent_group_info.group().n_ltr()):
      for i_smx in xrange(parent_group_info.group().n_smx()):
        s = parent_group_info.group()(i_ltr,0,i_smx)
        non_centric.expand_smx(s)
    assert non_centric.f_inv() == 1
    assert non_centric.order_z() * 2 == parent_group_info.group().order_z()
    non_centric_info = sgtbx.space_group_info(group=non_centric)
    unit_cell = non_centric_info.any_compatible_unit_cell(volume=1000)
    crystal_symmetry = crystal.symmetry(
      unit_cell=unit_cell,
      space_group_info=non_centric_info)
    niggli_symmetry = crystal_symmetry.niggli_cell()
    lattice_group = lattice_symmetry.group(
      niggli_symmetry.unit_cell(), max_delta=0.5)
    lattice_group_info = sgtbx.space_group_info(group=lattice_group)
    assert lattice_group_info.group() == niggli_symmetry.space_group()
    subgrs = subgroups.subgroups(lattice_group_info).groups_parent_setting()
    for group in subgrs:
      subsym = crystal.symmetry(
        unit_cell=niggli_symmetry.unit_cell(),
        space_group=group,
        assert_is_compatible_unit_cell=00000)
      assert subsym.unit_cell().is_similar_to(niggli_symmetry.unit_cell())
      assert lattice_symmetry.find_max_delta(
        niggli_cell=niggli_symmetry.unit_cell(),
        group=group) < 0.6
  print "OK"

if (__name__ == "__main__"):
  run()
