from cctbx.crystal.find_best_cell import alternative_find_best_cell as fbc
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from libtbx.test_utils import approx_equal

def tst_find_best_cell():
  uc_array=[ uctbx.unit_cell( '40, 50, 60, 90, 90, 90' ),
             uctbx.unit_cell( '40, 60, 50, 90, 90, 90' ),
             uctbx.unit_cell( '50, 40, 60, 90, 90, 90' ),
             uctbx.unit_cell( '50, 60, 40, 90, 90, 90' ),
             uctbx.unit_cell( '60, 40, 50, 90, 90, 90' ),
             uctbx.unit_cell( '60, 50, 40, 90, 90, 90' ) ]

  uc_correct = [ uctbx.unit_cell( '40, 50, 60, 90, 90, 90' ),
                 uctbx.unit_cell( '40, 60, 50, 90, 90, 90' ),
                 uctbx.unit_cell( '40, 50, 60, 90, 90, 90' ),
                 uctbx.unit_cell( '50, 60, 40, 90, 90, 90' ),
                 uctbx.unit_cell( '40, 60, 50, 90, 90, 90' ),
                 uctbx.unit_cell( '50, 60, 40, 90, 90, 90' ) ]

  sg_info =  sgtbx.space_group_info( 'P 21 21 2' )
  sg_info_2 = sgtbx.space_group_info( 'I 21 21 21' )
  sg = sg_info.group()
  sg_2 = sg_info_2.group()

  for uc, correct in zip(uc_array,uc_correct):
    best_cell_finder = fbc( uc, sg )
    assert approx_equal( correct.parameters(),
                         best_cell_finder.return_best_cell().parameters() )

    cb_op = best_cell_finder.return_change_of_basis_op_to_best_cell()
    xs = crystal.symmetry( uc, space_group=sg)
    assert approx_equal( correct.parameters(),
                         xs.change_basis(cb_op).unit_cell().parameters() )


    best_cell_finder = fbc( uc, sg_2 )
    assert approx_equal( uc_array[0].parameters(),
                         best_cell_finder.return_best_cell().parameters() )

    xs = crystal.symmetry( uc, space_group=sg_2)
    cb_op = best_cell_finder.return_change_of_basis_op_to_best_cell()
    assert approx_equal( uc_array[0].parameters(),
                         xs.change_basis(cb_op).unit_cell().parameters() )



def run():
  tst_find_best_cell()
  print "OK"

run()
