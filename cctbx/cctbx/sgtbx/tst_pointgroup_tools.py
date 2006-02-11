import sys,os
from cctbx import uctbx
from cctbx import sgtbx
from cctbx.array_family import flex
from cctbx import crystal
from cctbx import miller
from cctbx.sgtbx import pointgroup_tools as pt

def tst_pgtools():
  unit_cell = uctbx.unit_cell('40, 40, 60, 90.0, 90.0, 90.0')
  mi = flex.miller_index(((2,4,6), (2,4,8)))
  xs = crystal.symmetry(unit_cell, "P 21")
  ms = miller.set(xs, mi)

  # Go to the minimum cell, for safety
  cob_min_cell = ms.change_of_basis_op_to_minimum_cell()

  ms_new = ms.change_basis( cob_min_cell )

  lattice_group = sgtbx.lattice_symmetry.group(
    ms_new.unit_cell(),
    max_delta=5.0)

  point_group_low = ms_new.space_group().build_derived_point_group()
  point_group_high = lattice_group.build_derived_point_group()

  pgtree = pt.point_group_graph(point_group_low,point_group_high)

  # find the possible routes from 'P 2' to 'P 4 2 2'
  atlas = pgtree.graph.find_all_paths( 'P 1 2 1', 'P 4 2 2')
  route_1 = ['P 1 2 1', 'P 4 2 2']
  route_2 = ['P 1 2 1', 'P 2 2 2', 'P 4 2 2']
  assert route_1 in atlas
  assert route_2 in atlas
  assert len(atlas)==2

  # Now lets 'disqualify' point group 'P 2 2 2'
  pgtree.remove_point_group_and_its_super_groups_from_graph( sgtbx.space_group_info(16).group()  )
  assert len(pgtree.graph.node_objects)==1
  assert pgtree.graph.node_objects.has_key ( 'P 1 2 1' )

  # test *VERY* briefly the space group generator
  tst = pt.find_compatible_space_groups( lattice_group,
                                         ms.space_group(),
                                         ms.unit_cell() )
  candidates = [ 'P 4 21 2', 'P 41 21 2', 'P 43 21 2', 'P 42 21 2', 'I 4 2 2', 'I 41 2 2' ]
  for sg in tst.likely_sg:
    assert  str( sgtbx.space_group_info(group=sg) ) in candidates


def run():
  tst_pgtools()
  print "OK"



run()
