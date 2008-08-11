# Change the hand of a set of coordinates (useful in heavy atom location).

from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.web import io_utils
from cctbx.web import cgi_utils

def interpret_form_data(form):
  inp = cgi_utils.inp_from_form(form,
    (("ucparams", "1 1 1 90 90 90"),
     ("sgsymbol", "P1"),
     ("convention", ""),
     ("coor_type", None),
     ("skip_columns", "0")))
  inp.coordinates = cgi_utils.coordinates_from_form(form)
  return inp

def run(server_info, inp, status):
  print "<pre>"

  io_utils.show_input_symbol(inp.sgsymbol, inp.convention)
  crystal_symmetry = crystal.symmetry(
      unit_cell=uctbx.unit_cell(inp.ucparams),
      space_group_info=sgtbx.space_group_info(
        symbol=inp.sgsymbol,
        table_id=inp.convention))
  crystal_symmetry.show_summary()
  print

  change_of_hand_op \
    = crystal_symmetry.space_group_info().type().change_of_hand_op()
  print "Change-of-hand matrix:", change_of_hand_op.c()
  print "              Inverse:", change_of_hand_op.c_inv()
  print

  print inp.coor_type, "coordinates:"
  print

  skip_columns = io_utils.interpret_skip_columns(inp.skip_columns)

  for line in inp.coordinates:
    skipped, coordinates = io_utils.interpret_coordinate_line(line,skip_columns)
    if (inp.coor_type != "Fractional"):
      coordinates = crystal_symmetry.unit_cell().fractionalize(coordinates)
    flipped_coordinates = change_of_hand_op(coordinates)
    if (inp.coor_type != "Fractional"):
      flipped_coordinates \
        = crystal_symmetry.unit_cell().orthogonalize(flipped_coordinates)
    print skipped, "%.6g %.6g %.6g" % tuple(flipped_coordinates)

  print "</pre>"
