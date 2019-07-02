from __future__ import absolute_import, division, print_function
import sys
from cctbx import uctbx
from cctbx import sgtbx
from cctbx.array_family import flex
from cctbx import crystal
from cctbx.crystal.find_best_cell import alternative_find_best_cell as fbc
from cctbx import miller
from cctbx.sgtbx import pointgroup_tools as pt
from libtbx.test_utils import approx_equal
from six.moves import cStringIO as StringIO
from six.moves import range

def tst_pgtools():
  unit_cell = uctbx.unit_cell('40, 40, 60, 90.0, 90.0, 90.0')
  mi = flex.miller_index(((2,4,6), (2,4,8)))
  xs = crystal.symmetry(unit_cell, "P 1 2 1")
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
  pgtree.remove_point_group_and_its_super_groups_from_graph(
    str(sgtbx.space_group_info(16)))
  assert len(pgtree.graph.node_objects)==1
  assert 'P 1 2 1' in pgtree.graph.node_objects


def tst_sg_tools():
  unit_cell = uctbx.unit_cell('40, 50, 60, 90.0, 90.0, 90.0')
  mi = flex.miller_index(((2,4,6), (2,4,8)))
  xs = crystal.symmetry(unit_cell, "P 1 1 21")
  ms = miller.set(xs, mi)

  tmp_choice = pt.space_group_graph_from_cell_and_sg(unit_cell,
                                                     xs.space_group())

  xs1=crystal.symmetry(uctbx.unit_cell('40.00 60.00 50.00 90.00 90.00 90.00'),
                       'P 1 21 1')

  xs2=crystal.symmetry(uctbx.unit_cell('40.00 50.00 60.00 90.00 90.00 90.00'),
                       'P 2 2 21')

  xs3=crystal.symmetry(uctbx.unit_cell('40.00 60.00 50.00 90.00 90.00 90.00'),
                       'P 21 21 2')

  xs4=crystal.symmetry(uctbx.unit_cell('50.00 60.00 40.00 90.00 90.00 90.00'),
                       'P 21 21 2')

  xs5=crystal.symmetry(uctbx.unit_cell('40.00 50.00 60.00 90.00 90.00 90.00'),
                       'P 21 21 21')

  p222_dict = { str( xs2.unit_cell().parameters())
                + " " + str( xs2.space_group_info() ): None,
                str( xs3.unit_cell().parameters())
                + " " + str( xs3.space_group_info() ):None,
                str( xs4.unit_cell().parameters())
                + " " + str( xs4.space_group_info() ):None,
                str( xs5.unit_cell().parameters())
                + " " + str( xs5.space_group_info() ):None
                }


  p112 = tmp_choice.pg_graph.graph.node_objects['P 1 1 2'].allowed_xtal_syms
  # check the cell parameters
  assert approx_equal(p112[0][0].unit_cell().parameters(),
                      xs1.unit_cell().parameters())
  # check the sg
  assert p112[0][0].space_group() == xs1.space_group()
  # check the change of basis operator
  assert approx_equal(xs.change_basis(p112[0][1]).unit_cell().parameters(),
                      xs1.unit_cell().parameters())

  p222 = tmp_choice.pg_graph.graph.node_objects['P 2 2 2'].allowed_xtal_syms

  for xs in p222:
    comp_string = (str(xs[0].unit_cell().parameters())
                   + " " + str(xs[0].space_group_info()))
    assert comp_string in p222_dict



def test_extensively( this_chunk ):
  n_chunks = 1
  i_chunk = int(this_chunk)

  assert i_chunk <= n_chunks
  for i_sg, sg in enumerate(sgtbx.space_group_symbol_iterator()):
    if i_sg % n_chunks != i_chunk:
      continue

    buffer = StringIO()

    group = sgtbx.space_group(sg.hall())

    #if group != sgtbx.space_group_info( symbol = 'I 21 21 21' ).group():
    #  continue

    unit_cell =sgtbx.space_group_info(group=group).any_compatible_unit_cell(
      volume=57*57*76)

    xs_in = crystal.symmetry(unit_cell=unit_cell, space_group=group)
    xs_in = xs_in.best_cell()

    # Just make sure we have the 'best cell' to start out with

    xs_minimum = xs_in.change_basis(xs_in.change_of_basis_op_to_niggli_cell())

    sg_min = xs_minimum.space_group()
    sg_min_info = sgtbx.space_group_info(group = sg_min)
    sg_min_to_ref = sg_min_info.change_of_basis_op_to_reference_setting()
    sg_in_to_min = xs_in.change_of_basis_op_to_niggli_cell()

    xs_ref = xs_minimum.change_basis(sg_min_to_ref)
    best_cell_finder = fbc(xs_ref.unit_cell(),
                           xs_ref.space_group())
    xs_ref = best_cell_finder.return_best_xs()

    xs_lattice_pg = sgtbx.lattice_symmetry.group(
      xs_minimum.unit_cell(),
      max_delta=5.0)
    sg_in_ref_setting = sg_min.change_basis(sg_min_to_ref)

    # -----------
    queue = []
    if sg_in_ref_setting.is_chiral():
      print("Testing : ",\
            sgtbx.space_group_info(group=group),\
            "( or ", sgtbx.space_group_info(group=sg_in_ref_setting),\
            " in reference setting)")
      for s in sg_min:
        new_sg = sgtbx.space_group()
        new_sg.expand_smx(s)
        if new_sg in queue:
          continue
        else:
          queue.append(new_sg)
        xs_cheat = crystal.symmetry(xs_minimum.unit_cell(),
                                    space_group=new_sg)
        #print
        #print "Using symop", s

        xs_cheat_min = xs_cheat.change_basis(
          xs_cheat.change_of_basis_op_to_niggli_cell() )

        #print "This results in space group: ", xs_cheat.space_group_info()

        sg_clues = pt.space_group_graph_from_cell_and_sg(
          xs_cheat_min.unit_cell(),
          xs_cheat_min.space_group() )

        #sg_clues.show()

        found_it = False
        #print
        #print " --- Spacegroups consistent with input parameters ---"
        for sg_and_uc in sg_clues.return_likely_sg_and_cell():
          check_sg = False
          check_uc = False
          xs = crystal.symmetry(sg_and_uc[1], space_group=sg_and_uc[0])
          if approx_equal(
            sg_and_uc[1].parameters(),
            xs_ref.unit_cell().parameters(),
            eps=0.001,
            out=buffer):
            if sg_and_uc[0] == sg_in_ref_setting:
              found_it = True

          #print sgtbx.space_group_info( group=sg_and_uc[0] ), \
          #      sg_and_uc[1].parameters()
        if not found_it:
          print("FAILURE: ", sg.hall())
          assert found_it



def test_reference_setting_choices():
  buffer = StringIO()

  for space_group_info in sgtbx.reference_space_group_infos():
    space_group = space_group_info.group()

    uc = space_group_info.any_compatible_unit_cell(volume=57*57*76)
    xs = crystal.symmetry(uc, space_group=space_group)

    cobs = pt.reference_setting_choices(space_group)

    if len(cobs)>1:
      tmp_array = []
      for cob in cobs:
        xs_new =crystal.symmetry(xs.change_basis( cob ).unit_cell(),
                                 space_group=xs.space_group() )
        best_cell_finder = fbc(xs_new.unit_cell(),
                               xs_new.space_group() )

        xs_new = best_cell_finder.return_best_xs()
        tmp_array.append(xs_new)
      count = 0
      for tmp_xs1 in range(len(tmp_array)):
        for tmp_xs2 in range(len(tmp_array)):

          if (tmp_xs1 != tmp_xs2):
            assert (tmp_array[tmp_xs1].space_group()
                    == tmp_array[tmp_xs2].space_group())
            assert not approx_equal(
              tmp_array[tmp_xs1].unit_cell().parameters(),
              tmp_array[tmp_xs2].unit_cell().parameters(),
              eps=0.001,
              out=buffer)

def exercise_compatible_symmetries():
  pg = sgtbx.space_group('P 2x')
  ops = [ op.as_xyz() for op in pt.compatible_symmetries(pg) ]
  assert ops == [ 'x,-y,-z', 'x+1/2,-y,-z' ]

  pg = sgtbx.space_group('P 2x 2y')
  ops = [ op.as_xyz() for op in pt.compatible_symmetries(pg) ]
  assert ops == [ 'x,-y,-z', 'x+1/2,-y,-z',
                  '-x,y,-z', '-x,y+1/2,-z',
                  '-x,-y,z', '-x,-y,z+1/2' ]

  pg = sgtbx.space_group('P 2x 3*')
  ops = [ op.as_xyz() for op in pt.compatible_symmetries(pg) ]
  assert 'x,-y,-z' in ops
  assert 'x+1/2,-y,-z' in ops

  assert '-x,y,-z' in ops
  assert '-x,y+1/2,-z' in ops

  assert '-x,-y,z' in ops
  assert '-x,-y,z+1/2' in ops

  assert 'z,x,y' in ops
  assert 'z+1/3,x+1/3,y+1/3' in ops

  assert 'y,z,x' in ops
  assert 'y+1/3,z+1/3,x+1/3' in ops

  assert 'z,-x,-y' in ops
  assert 'z+1/3,-x-1/3,-y+1/3' in ops

  assert '-y,z,-x' in ops
  assert '-y-1/3,z+1/3,-x+1/3' in ops

  assert '-y,-z,x' in ops
  assert '-y+1/3,-z-1/3,x+1/3' in ops

  assert 'y,-z,-x' in ops
  assert 'y-1/3,-z-1/3,-x+1/3' in ops

  assert '-z,-x,y' in ops
  assert '-z-1/3,-x+1/3,y+1/3' in ops

  assert '-z,x,-y' in ops
  assert '-z-1/3,x-1/3,-y+1/3' in ops

  pg = sgtbx.space_group('-P -2y')
  ops = [ op.as_xyz() for op in pt.compatible_symmetries(pg) ]
  assert ops == ['x,-y,z', 'x+1/2,-y,z', 'x,-y,z+1/2', 'x+1/2,-y,z+1/2',
                 '-x,-y,-z', '-x,y,-z', '-x,y+1/2,-z']

def run():
  exercise_compatible_symmetries()
  tst_pgtools()
  tst_sg_tools()
  test_reference_setting_choices()
  print("OK")


if (__name__ == "__main__"):
  if len(sys.argv)>1:
    if sys.argv[1]=='insane':
      run()
      test_extensively(0) # Should take not more then 45 minutes
                          #on a reasonable machine
      print("OK")
  else:
    run()
