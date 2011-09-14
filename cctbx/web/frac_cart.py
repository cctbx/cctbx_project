# Convert coordinates from/to Fractional or Orthogonal.

from cctbx import uctbx
from cctbx.web import io_utils
from cctbx.web import cgi_utils

def interpret_form_data(form):
  inp = cgi_utils.inp_from_form(form,
    (("ucparams", "1 1 1 90 90 90"),
     ("coor_type", None),
     ("skip_columns", "0")))
  inp.coordinates = cgi_utils.coordinates_from_form(form)
  return inp

def run(server_info, inp, status):
  print "<pre>"
  unit_cell = uctbx.unit_cell(inp.ucparams)
  unit_cell.show_parameters()
  print

  if (inp.coor_type == "Fractional"):
    print "Cartesian coordinates:"
  else:
    print "Fractional coordinates:"
  print

  skip_columns = io_utils.interpret_skip_columns(inp.skip_columns)

  for line in inp.coordinates:
    skipped, coordinates = io_utils.interpret_coordinate_line(line,skip_columns)
    if (inp.coor_type == "Fractional"):
      c = unit_cell.orthogonalize(coordinates)
    else:
      c = unit_cell.fractionalize(coordinates)
    print skipped, "%.6g %.6g %.6g" % tuple(c)

  print "</pre>"
