from cctbx import adptbx
import sys

def format_float(value, zero_threshold=None):
  if (value is None): return "."
  if (zero_threshold is not None and abs(value) < zero_threshold):
    return "0"
  s = "%.6g" % value
  i = s.find("e")
  if (i < 0): return s
  j = i + 1
  exp_sign = s[j]
  assert exp_sign in ["-", "+"]
  j += 1
  n = len(s)
  while (j < n): # remove leading zeros from exponent literal
    if (s[j] != "0"): return s[:i] + exp_sign + s[j:]
    j += 1
  return s[:i]

def format_floats(values, zero_threshold=None):
  return " ".join([format_float(value=value, zero_threshold=zero_threshold)
    for value in values])

def quote_string(s):
  if (len(s) == 0): return "."
  if (s.find("'") < 0): return "'" + s + "'"
  if (s.find('"') < 0): return '"' + s + '"'
  raise RuntimeError(
    "string to be quoted for CIF output contains both a single and"
    " a double quote: sorry, this is currently not supported")

def as_cif_simple(self, out=None):
  if (out is None): out = sys.stdout
  print >> out, "data_global"
  sg = self.space_group()
  if (sg is not None):
    print >> out, "loop_"
    print >> out, "    _symmetry_equiv_pos_as_xyz"
    for s in sg:
      print >> out, "    '%s'" % s.as_xyz(
        decimal=False, t_first=False, symbol_letters="xyz", separator=', ')
  uc = self.unit_cell()
  if (uc is not None):
    a,b,c,alpha,beta,gamma = uc.parameters()
    print >> out, "_cell_length_a", format_float(value=a)
    print >> out, "_cell_length_b", format_float(value=b)
    print >> out, "_cell_length_c", format_float(value=c)
    print >> out, "_cell_angle_alpha", format_float(value=alpha)
    print >> out, "_cell_angle_beta", format_float(value=beta)
    print >> out, "_cell_angle_gamma", format_float(value=gamma)
  scatterers = self.scatterers()
  if (scatterers.size() != 0):
    print >> out, """\
loop_
    _atom_site_label
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    _atom_site_U_iso_or_equiv
    _atom_site_occupancy
    _atom_site_type_symbol"""
  aniso_scatterers = []
  for sc in scatterers:
    print >> out, "   ", \
      quote_string(sc.label), \
      format_floats(
        values=sc.site,
        zero_threshold=1.e-12), \
      format_float(
        value=sc.u_iso_or_equiv(unit_cell=uc),
        zero_threshold=1.e-12), \
      format_float(value=sc.occupancy), \
      quote_string(sc.scattering_type)
    if (sc.flags.use_u_aniso()):
      aniso_scatterers.append(sc)
  if (len(aniso_scatterers) != 0):
    print >> out, """\
loop_
    _atom_site_aniso_label
    _atom_site_aniso_U_11
    _atom_site_aniso_U_22
    _atom_site_aniso_U_33
    _atom_site_aniso_U_12
    _atom_site_aniso_U_13
    _atom_site_aniso_U_23"""
  for sc in aniso_scatterers:
    print >> out, "   ", \
      quote_string(sc.label), \
      format_floats(
        values=sc.u_cart_plus_u_iso(unit_cell=uc),
        zero_threshold=1.e-12)
