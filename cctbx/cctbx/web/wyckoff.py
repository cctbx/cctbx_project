# When studying a crystal structure it can be helpful to know the Wyckoff
# positions of the atoms in the structure. This script can be used to
# assign Wyckoff letters.

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
              ("min_distance_sym_equiv", "0.5"),
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
  special_position_settings = crystal.special_position_settings(
    crystal.symmetry(
      unit_cell=uctbx.unit_cell(inp.ucparams),
      space_group_info=sgtbx.space_group_info(
        symbol=inp.sgsymbol,
        table_id=inp.convention)),
    min_distance_sym_equiv=float(inp.min_distance_sym_equiv))
  special_position_settings.show_summary()
  print "Minimum distance between symmetrically equivalent sites:",
  print special_position_settings.min_distance_sym_equiv()
  print

  skip_columns = utils.interpret_skip_columns(inp.skip_columns)

  wyckoff_table=special_position_settings.space_group_info().wyckoff_table()
  unit_cell = special_position_settings.unit_cell()
  print "</pre><table border=2 cellpadding=2>"
  status.in_table = True
  print "<tr>"
  if (skip_columns): print "<th>"
  print "<th colspan=3>" + inp.coor_type + " coordinates"
  print "<th>Multiplicity"
  print "<th>Wyckoff letter"
  print "<th>Site symmetry<br>point group type"
  print "<th>Special position operator"
  print "</tr>"
  for line in inp.coordinates:
    skipped, coordinates = utils.interpret_coordinate_line(line,skip_columns)
    if (inp.coor_type != "Fractional"):
      coordinates = unit_cell.fractionalize(coordinates)
    site_symmetry = special_position_settings.site_symmetry(coordinates)
    exact_site = site_symmetry.exact_site()
    if (inp.coor_type != "Fractional"):
      exact_site = unit_cell.orthogonalize(exact_site)
    wyckoff_mapping = wyckoff_table.mapping(site_symmetry)
    print "<tr>"
    if (skip_columns): print "<td>", skipped
    for x in exact_site: print "<td><tt>%.6g</tt>" % (x,)
    print "<td align=center>", wyckoff_mapping.position().multiplicity()
    print "<td align=center>", wyckoff_mapping.position().letter()
    print "<td align=center>", site_symmetry.point_group_type()
    print "<td><tt>" + str(site_symmetry.special_op()) + "</tt>"
    print "</tr>"
  print "</table><pre>"
  status.in_table = False

  print "</pre>"
