# This script reports a number of space group properties given a space
# group symbol or symmetry matrices, or a combination of the two.

from cctbx import sgtbx
from cctbx.sgtbx import harker
from cctbx.web import io_utils
from cctbx.web import cgi_utils

def interpret_form_data(form):
  inp = cgi_utils.inp_from_form(form,
    (("sgsymbol", ""),
     ("convention", "")))
  inp.shelx_latt = []
  inp.symxyz = []
  if (form.has_key("symxyz")):
    lines = form["symxyz"].value.split("\015\012")
    for l in lines:
      # Treat SHELX LATT & SYMM cards
      s = l.strip()
      CARD = s[:4].upper()
      if   (CARD == "LATT"):
        inp.shelx_latt.append(s[4:])
      elif (CARD == "SYMM"):
        inp.symxyz.append(s[4:].strip())
      else:
        # Plain symmetry operations
        for s in l.split(";"):
          s = s.strip()
          if (s != ""): inp.symxyz.append(s)
  return inp

def str_ev(ev):
  return "[%d,%d,%d]" % ev

def rt_mx_analysis_header():
  print "<tr>"
  print "<th>Matrix"
  print "<th>Rotation-part type"
  print "<th>Axis direction"
  print "<th>Screw/glide component"
  print "<th>Origin shift"
  print "</tr>"

def rt_mx_analysis(s):
  print "<tr>"
  print "<td><tt>" + str(s) + "</tt>"
  r_info = sgtbx.rot_mx_info(s.r())
  t_info = sgtbx.translation_part_info(s)
  if (r_info.type() == 1):
    print "<td>1<td>-<td>-<td>-"
  elif (r_info.type() == -1):
    print "<td>%d<td>-<td>-<td>%s" % (
      r_info.type(),
      str(t_info.origin_shift()))
  elif (abs(r_info.type()) == 2):
    print "<td>%d<td>%s<td>%s<td>%s" % (
      r_info.type(),
      str_ev(r_info.ev()),
      str(t_info.intrinsic_part()),
      str(t_info.origin_shift()))
  else:
    print "<td>%d^%d<td>%s<td>%s<td>%s" % (
       r_info.type(),
       r_info.sense(),
       str_ev(r_info.ev()),
       str(t_info.intrinsic_part()),
       str(t_info.origin_shift()))
  print "</tr>"

def show_group_generic(sg_type, status):
  sg = sg_type.group()
  print "Number of lattice translations:", sg.n_ltr()
  if (sg.is_centric()):
    print "Space group is centric."
  else:
    print "Space group is acentric."
  if (sg.is_chiral()):
    print "Space group is chiral."
  if (sg_type.is_enantiomorphic()):
    print "Space group is enantiomorphic."
  print "Number of representative symmetry operations:", sg.n_smx()
  print "Total number of symmetry operations:", sg.order_z()
  print
  print "Parallelepiped containing an asymmetric unit:"
  try: brick = sgtbx.brick(sg_type)
  except RuntimeError, e:
    print " ", e
  else:
    print " ", str(brick).replace("<", "&lt;").replace(">", "&gt;")
  print
  print "List of symmetry operations:"
  print "</pre><table border=2 cellpadding=2>"
  status.in_table = True
  rt_mx_analysis_header()
  for s in sg: rt_mx_analysis(s)
  print "</table><pre>"
  status.in_table = False
  print

def show_symbols(symbols):
  print "  Space group number:", symbols.number()
  print "  Schoenflies symbol:", symbols.schoenflies()
  print "  Hermann-Mauguin symbol:", symbols.hermann_mauguin()
  e = symbols.extension()
  if (e != "\0"):
    if (e in "12"):
      print "  Origin choice:", e
    elif (e == "H"):
      print "  Trigonal using hexagonal axes"
    elif (e == "R"):
      print "  Trigonal using rhombohedral axes"
    else:
      raise RuntimeError, "Internal error"
  q = symbols.qualifier()
  if (q != ""):
    if (symbols.number() < 16):
      if (q[-1] in "123"):
        unique_axis = q[:-1]
        cell_choice = q[-1]
      else:
        unique_axis = q
        cell_choice = ""
      print "  Unique axis:", unique_axis
      if (cell_choice != ""):
        print "  Cell choice:", cell_choice
    else:
      print "  Relation to standard setting:", q
  print "  Hall symbol:", symbols.hall().strip()

def expand_shelx_latt(sg, n_fld):
  z_dict = {
    "P": 1,
    "I": 2,
    "R": 3,
    "F": 4,
    "A": 5,
    "B": 6,
    "C": 7,
  }
  n_dict = {}
  for z in z_dict.keys(): n_dict[z_dict[z]] = z
  try:
    n = int(n_fld)
    z = n_dict[abs(n)]
  except Exception:
    raise RuntimeError, "Format Error: LATT " + str(n_fld)
  print "Addition of SHELX LATT " + str(n) + ":"
  if (n > 0):
    print "  Addition of centre of inversion at the origin."
    sg.expand_smx(sgtbx.rt_mx("-x,-y,-z"))
  print "  Addition of lattice translations for centring type " + str(z) + "."
  sg.expand_conventional_centring_type(z)
  print

def run(server_info, inp, status):
  print "<pre>"
  io_utils.show_input_symbol(inp.sgsymbol, inp.convention)

  symbols_inp = None
  lookup_symbol = inp.sgsymbol
  if (lookup_symbol == ""): lookup_symbol = "P 1"
  if (inp.convention == "Hall"):
    hall_symbol = lookup_symbol
  else:
    symbols_inp = sgtbx.space_group_symbols(lookup_symbol, inp.convention)
    hall_symbol = symbols_inp.hall()
    if (symbols_inp.number() == 0):
      symbols_inp = None
      inp.convention = "Hall"
    else:
      print "Result of symbol lookup:"
      show_symbols(symbols_inp)
      print

  try:
    ps = sgtbx.parse_string(hall_symbol)
    sg = sgtbx.space_group(ps)
  except RuntimeError, e:
    print "--&gt;" + ps.string() + "&lt;--"
    print ("-" * (ps.where() + 3)) + "^"
    raise

  if (len(inp.shelx_latt) != 0):
    for n_fld in inp.shelx_latt:
      expand_shelx_latt(sg, n_fld)

  if (len(inp.symxyz) != 0):
    print "Addition of symmetry operations:"
    print "</pre><table border=2 cellpadding=2>"
    status.in_table = True
    rt_mx_analysis_header()
    for s in inp.symxyz:
      ps = sgtbx.parse_string(s)
      try:
        s = sgtbx.rt_mx(ps)
      except RuntimeError, e:
        print "</table><pre>"
        status.in_table = False
        print "--&gt;" + ps.string() + "&lt;--"
        print ("-" * (ps.where() + 3)) + "^"
        raise
      rt_mx_analysis(s)
      sg.expand_smx(s)
    print "</table><pre>"
    status.in_table = False
    print

  sg_type = sg.type()
  show_group_generic(sg_type, status)

  if (   inp.convention == "Hall"
      or len(inp.shelx_latt) != 0
      or len(inp.symxyz) != 0):
    symbols_match = sg.match_tabulated_settings()
    if (symbols_match.number() != 0):
      if (   symbols_inp is None
          or    symbols_inp.universal_hermann_mauguin()
             != symbols_match.universal_hermann_mauguin()):
        print "Symmetry operations match:"
        show_symbols(symbols_match)
        print
      else:
        print "Additional symmetry operations are redundant."
        print
    else:
      print "Space group number:", sg_type.number()
      print "Conventional Hermann-Mauguin symbol:", \
        sgtbx.space_group_symbols(sg_type.number()) \
        .universal_hermann_mauguin()
      print "Universal    Hermann-Mauguin symbol:", \
        sg_type.universal_hermann_mauguin_symbol()
      print "Hall symbol:", sg_type.hall_symbol()
      print "Change-of-basis matrix:", sg_type.cb_op().c()
      print "               Inverse:", sg_type.cb_op().c_inv()
      print

  wyckoff_table = sgtbx.wyckoff_table(sg_type)
  print "List of Wyckoff positions:"
  print "</pre><table border=2 cellpadding=2>"
  status.in_table = True
  print "<tr>"
  print "<th>Wyckoff letter"
  print "<th>Multiplicity"
  print "<th>Site symmetry<br>point group type"
  print "<th>Representative special position operator"
  print "</tr>"
  for i_position in xrange(wyckoff_table.size()):
    position = wyckoff_table.position(i_position)
    print "<tr>"
    print "<td>%s<td>%d<td>%s<td><tt>%s</tt>" % (
      position.letter(),
      position.multiplicity(),
      position.point_group_type(),
      str(position.special_op_simplified()))
    print "</tr>"
  print "</table><pre>"
  status.in_table = False
  print

  print "Harker planes:"
  print "</pre><table border=2 cellpadding=2>"
  status.in_table = True
  print "<tr>"
  print "<th>Algebraic"
  print "<th>Normal vector"
  print "<th>A point in the plane"
  print "</tr>"
  planes = harker.planes_fractional(sg)
  for plane in planes.list:
    print "<tr>"
    print "<td><tt>%s</tt><td>%s<td>%s" % (
      plane.algebraic(), str_ev(plane.n), str(plane.p))
    print "</tr>"
  print "</table><pre>"
  status.in_table = False
  print

  print "Additional generators of Euclidean normalizer:"
  ss = sgtbx.structure_seminvariants(sg)
  ss_vm = ss.vectors_and_moduli()
  print "  Number of structure-seminvariant vectors and moduli:", len(ss_vm)
  if (len(ss_vm)):
    print "    Vector    Modulus"
    for vm in ss_vm: print "   ", vm.v, vm.m
  k2l = sg_type.addl_generators_of_euclidean_normalizer(True, False)
  l2n = sg_type.addl_generators_of_euclidean_normalizer(False, True)
  if (len(k2l)):
    print "  Inversion through a centre at:",
    assert len(k2l) == 1
    print sgtbx.translation_part_info(k2l[0]).origin_shift()
  if (len(l2n)):
    print "  Further generators:"
    print "</pre><table border=2 cellpadding=2>"
    status.in_table = True
    rt_mx_analysis_header()
    for s in l2n: rt_mx_analysis(s)
    print "</table><pre>"
    status.in_table = False
  print

  print "Grid factors implied by symmetries:"
  grid_sg = sg.gridding()
  grid_ss = ss.gridding()
  eucl_sg = sg_type.expand_addl_generators_of_euclidean_normalizer(True,True)
  grid_eucl = eucl_sg.refine_gridding(grid_ss)
  print "  Space group:", grid_sg
  print "  Structure-seminvariant vectors and moduli:", grid_ss
  print "  Euclidean normalizer:", grid_eucl
  print
  print "  All points of a grid over the unit cell are mapped"
  print "  exactly onto other grid points only if the factors"
  print "  shown above are factors of the grid."
  print

  print "</pre>"
