# Convert coordinates from/to Fractional or Orthogonal.

from cctbx import uctbx
from cctbx.web import utils
import sys
import traceback

class empty: pass

def interpret_form_data(form):
  inp = empty()
  for key in (("ucparams", "1 1 1 90 90 90"),
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

def run(cctbx_url, inp):
  print "Content-type: text/html"
  print

  print "<pre>"
  try:
    unit_cell = uctbx.unit_cell(inp.ucparams)
    uctbx.show_parameters(unit_cell)
    print

    if (inp.coor_type == "Fractional"):
      print "Cartesian coordinates:"
    else:
      print "Fractional coordinates:"
    print

    skip_columns = int(inp.skip_columns)
    if (skip_columns < 0):
      raise ValueError, "Negative number for columns to skip."

    for line in inp.coordinates:
      skipped, coordinates = utils.interpret_coordinate_line(line,skip_columns)
      if (inp.coor_type == "Fractional"):
        c = unit_cell.orthogonalize(coordinates)
      else:
        c = unit_cell.fractionalize(coordinates)
      print skipped, "%.6g %.6g %.6g" % tuple(c)

  except RuntimeError, e:
    print e
  except (AssertionError, ValueError):
    ei = sys.exc_info()
    print traceback.format_exception_only(ei[0], ei[1])[0]

  print "</pre>"
