# Change the hand of a set of coordinates (useful in heavy atom location).

from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.web import utils

class empty: pass

def interpret_form_data(form):
  inp = empty()
  for key in (("ucparams", "1 1 1 90 90 90"),
              ("sgsymbol", "P1"),
              ("convention", ""),
              ("coor_type", None),
              ("skip_columns", "0")):
    if (form.has_key(key[0])):
      inp.__dict__[key[0]] = form[key[0]].value.strip()
    else:
      inp.__dict__[key[0]] = key[1]
  inp.coordinates = []
  if (form.has_key("coordinates")):
    lines = form["coordinates"].value.split("\015\012")
    for l in lines:
      s = l.strip()
      if (len(s) != 0): inp.coordinates.append(s)
  return inp

def run(server_info, inp, status):
  print "<pre>"

  utils.show_input_symbol(inp.sgsymbol, inp.convention)
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

  skip_columns = utils.interpret_skip_columns(inp.skip_columns)

  for line in inp.coordinates:
    skipped, coordinates = utils.interpret_coordinate_line(line,skip_columns)
    if (inp.coor_type != "Fractional"):
      coordinates = crystal_symmetry.unit_cell().fractionalize(coordinates)
    flipped_coordinates = change_of_hand_op(coordinates)
    if (inp.coor_type != "Fractional"):
      flipped_coordinates \
        = crystal_symmetry.unit_cell().orthogonalize(flipped_coordinates)
    print skipped, "%.6g %.6g %.6g" % tuple(flipped_coordinates)

  print "</pre>"
