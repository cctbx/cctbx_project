from cctbx import xray
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.eltbx.caasf import wk1995
from cctbx import adptbx
from cctbx.web import utils
from scitbx.python_utils import complex_math

class empty: pass

def interpret_form_data(form):
  inp = empty()
  for key in (("ucparams", "1 1 1 90 90 90"),
              ("sgsymbol", "P1"),
              ("convention", ""),
              ("d_min", "1"),
              ("min_distance_sym_equiv", "0.5"),
              ("coor_type", None)):
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

def read_scatterer(flds, default_b_iso=3.0):
  scatterer = xray.scatterer("const")
  # Label [ScatFact] x y z [Occ [Biso]]
  try:
    scatterer.label = flds[0]
    try:
      float(flds[1])
    except:
      offs = 2
      scatterer.caasf = wk1995(flds[1], True)
    else:
      offs = 1
      scatterer.caasf = wk1995(flds[0], False)
    site = flds[offs : offs + 3]
    for i in xrange(3):
      site[i] = float(site[i])
    scatterer.site = site
    scatterer.occupancy = 1.
    scatterer.anisotropic_flag = False
    scatterer.u_iso = adptbx.b_as_u(default_b_iso)
    if (len(flds) >= offs + 4):
      scatterer.occupancy = float(flds[offs + 3])
      if (len(flds) == offs + 5):
        scatterer.u_iso = adptbx.b_as_u(float(flds[offs + 4]))
      else:
        assert (len(flds) < offs + 5)
  except:
    raise utils.FormatError, flds
  return scatterer

def run(server_info, inp, status):
  print "<pre>"
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

  d_min = float(inp.d_min)
  print "Minimum d-spacing:", d_min
  if (d_min <= 0.):
    raise ValueError, "d-spacing must be greater than zero."
  print

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
      (scatterer.label, scatterer.caasf.label(),
       wyckoff_position.multiplicity(), wyckoff_position.letter(),
       site_symmetry.point_group_type())
     + scatterer.site
     + (scatterer.occupancy, adptbx.u_as_b(scatterer.u_iso)))
  print "</table><pre>"
  status.in_table = False
  print

  f_calc_array = xray.structure_factors(
    xray_structure=structure,
    miller_set=miller.build_set(
      crystal_symmetry=structure,
      anomalous_flag=False,
      d_min=d_min)).f_calc_array()
  print "Number of Miller indices:", f_calc_array.indices().size()
  print
  print "</pre><table border=2 cellpadding=2>"
  status.in_table = True
  print "<tr>"
  print "<th>hkl<th>Amplitude<th>Phase"
  for i,h in f_calc_array.indices().items():
    print "<tr>"
    print "<td>%3d %3d %3d<td>%.6g<td align=right>%.3f" % (
      h + complex_math.abs_arg(f_calc_array.data()[i], deg=True))
  print "</table><pre>"
  status.in_table = False
  print

  print "</pre>"
