from cctbx import xray
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.eltbx.caasf import wk1995
from cctbx import adptbx
from cctbx.web import cgi_utils

def show_input_symbol(sgsymbol, convention, label="Input"):
  if (sgsymbol != ""):
    print label, "space group symbol:", sgsymbol
    print "Convention:",
    if   (convention == "A1983"):
      print "International Tables for Crystallography, Volume A 1983"
    elif (convention == "I1952"):
      print "International Tables for Crystallography, Volume I 1952"
    elif (convention == "Hall"):
      print "Hall symbol"
    else:
      print "Default"
    print

def interpret_skip_columns(skip_columns):
  result = int(skip_columns)
  if (result < 0):
    raise ValueError, "Negative number for columns to skip."
  return result

def interpret_coordinate_line(line, skip_columns):
  flds = line.split()
  if (len(flds) < skip_columns + 3): raise FormatError, line
  coordinates = [0,0,0]
  for i in xrange(3):
    try: coordinates[i] = float(flds[skip_columns + i])
    except: raise FormatError, line
  return " ".join(flds[:skip_columns]), coordinates

def read_scatterer(flds, default_b_iso=3.0):
  scatterer = xray.scatterer("const")
  # Label [ScatFact] x y z [Occ [Biso]]
  try:
    scatterer.label = flds[0]
    try:
      float(flds[1])
    except:
      offs = 2
      scatterer.caasf = wk1995(flds[1], 0001)
    else:
      offs = 1
      scatterer.caasf = wk1995(flds[0], 00000)
    site = flds[offs : offs + 3]
    for i in xrange(3):
      site[i] = float(site[i])
    scatterer.site = site
    scatterer.occupancy = 1.
    scatterer.anisotropic_flag = 00000
    scatterer.u_iso = adptbx.b_as_u(default_b_iso)
    if (len(flds) >= offs + 4):
      scatterer.occupancy = float(flds[offs + 3])
      if (len(flds) == offs + 5):
        scatterer.u_iso = adptbx.b_as_u(float(flds[offs + 4]))
      else:
        assert (len(flds) < offs + 5)
  except:
    raise cgi_utils.FormatError, flds
  return scatterer

def special_position_settings_from_inp(inp):
  return crystal.special_position_settings(
    crystal.symmetry(
      unit_cell=uctbx.unit_cell(inp.ucparams),
      space_group_info=sgtbx.space_group_info(
        symbol=inp.sgsymbol,
        table_id=inp.convention)),
    min_distance_sym_equiv=float(inp.min_distance_sym_equiv))

def structure_from_inp(inp, status, special_position_settings):
  wyckoff_table=special_position_settings.space_group_info().wyckoff_table()
  print "</pre><table border=2 cellpadding=2>"
  status.in_table = 0001
  print "<tr>"
  print "<th>Label"
  print "<th>Scattering<br>factor<br>label"
  print "<th>Multiplicty"
  print "<th>Wyckoff<br>position"
  print "<th>Site<br>symmetry"
  print "<th colspan=3>Fractional coordinates"
  print "<th>Occupancy<br>factor"
  print "<th>Biso"
  print "<tr>"
  structure = xray.structure(special_position_settings)
  print
  for line in inp.coordinates:
    scatterer = read_scatterer(line.split())
    if (inp.coor_type != "Fractional"):
      scatterer.site = structure.unit_cell().fractionalize(scatterer.site)
    structure.add_scatterer(scatterer)
    site_symmetry = structure.site_symmetry(scatterer.site)
    wyckoff_mapping = wyckoff_table.mapping(site_symmetry)
    wyckoff_position = wyckoff_mapping.position()
    print "<tr>"
    print (  "<td>%s<td>%s"
           + "<td align=center>%d<td align=center>%s<td align=center>%s"
           + "<td><tt>%.6g</tt><td><tt>%.6g</tt><td><tt>%.6g</tt>"
           + "<td align=center><tt>%.6g</tt>"
           + "<td align=center><tt>%.6g</tt>") % (
      (scatterer.label, scatterer.caasf.label(),
       wyckoff_position.multiplicity(), wyckoff_position.letter(),
       site_symmetry.point_group_type())
     + scatterer.site
     + (scatterer.occupancy, adptbx.u_as_b(scatterer.u_iso)))
  print "</table><pre>"
  status.in_table = 00000
  print
  return structure
