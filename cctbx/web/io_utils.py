import iotbx.pdb
from cctbx import xray
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
import cctbx.eltbx.xray_scattering
from cctbx import eltbx
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
    except Exception: raise FormatError, line
  return " ".join(flds[:skip_columns]), coordinates

def read_scatterer(flds, default_b_iso=3.0):
  scatterer = xray.scatterer(scattering_type="const")
  # Label [ScatFact] x y z [Occ [Biso]]
  try:
    scatterer.label = flds[0]
    try:
      float(flds[1])
    except Exception:
      offs = 2
      scatterer.scattering_type = eltbx.xray_scattering.get_standard_label(
        label=flds[1], exact=True)
    else:
      offs = 1
      scatterer.scattering_type = eltbx.xray_scattering.get_standard_label(
        label=flds[0], exact=False)
    site = flds[offs : offs + 3]
    for i in xrange(3):
      site[i] = float(site[i])
    scatterer.site = site
    scatterer.occupancy = 1.
    scatterer.set_use_u_iso_only()
    scatterer.u_iso = adptbx.b_as_u(default_b_iso)
    if (len(flds) >= offs + 4):
      scatterer.occupancy = float(flds[offs + 3])
      if (len(flds) == offs + 5):
        scatterer.u_iso = adptbx.b_as_u(float(flds[offs + 4]))
      else:
        assert (len(flds) < offs + 5)
  except Exception:
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
  status.in_table = True
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
      (scatterer.label, scatterer.scattering_type,
       wyckoff_position.multiplicity(), wyckoff_position.letter(),
       site_symmetry.point_group_type())
     + scatterer.site
     + (scatterer.occupancy, adptbx.u_as_b(scatterer.u_iso)))
  print "</table><pre>"
  status.in_table = False
  print
  return structure

def structure_from_inp_pdb(inp, status):
  pdb_file = inp.coordinates
  print "Input PDB file content:\n"
  for line in pdb_file:
    print line
  print
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_file)
  xray_structure = pdb_inp.xray_structure_simple(
    enable_scattering_type_unknown=True)
  u_tidy = xray_structure.scatterers().extract_u_cart(
    xray_structure.unit_cell())
  have_header = False
  for ut, atom in zip(u_tidy, pdb_inp.atoms_with_labels()):
    ui = atom.uij
    sui, sut = [(" %7.4f" * 6) % u for u in [ui, ut]]
    if (sui != sut):
      if (not have_header):
        print "Anisotropic displacement parameters (ANISOU) of the following"
        print "atoms were modified:"
        have_header = True
      print " ", atom.quote()
      print "    Ucart(input) =" + sui
      print "    Ucart(tidy)  =" + sut
  if (have_header):
    print
  return xray_structure
