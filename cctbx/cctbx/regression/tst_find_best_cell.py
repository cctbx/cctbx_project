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

  sg_info =  sgtbx.space_group_info( 'I 21 21 21' )
  sg = sg_info.group()

  for uc in uc_array:
    best_cell_finder = fbc( uc, sg )

    assert approx_equal( uc_array[0].parameters(),
                         best_cell_finder.return_best_cell().parameters() )
    best_cell_finder.return_best_xs().show_summary()

def run():
  tst_find_best_cell()
  print "OK"

run()
