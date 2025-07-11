"""Given PDB file and two atom selections that allow to define two lines compute the angle between these two lines"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.angle

import sys, os
import iotbx.pdb
from libtbx.utils import Sorry
from scitbx.linalg import eigensystem
import math
import scitbx.math
from scitbx import matrix
from six.moves import zip

legend = """phenix.angle:
  Given PDB file and two atom selections that allow to define two lines compute
  the angle between these two lines. If atom selection defines two points then
  the line is defined uniquely and passes through these points. If atom
  selection defines more than two points then line coincides with the longest
  axis of the cloud of points.

How to run:
  phenix.angle model.pdb "chain A and (resseq 1 and name CA or resseq 2 and name CA)" \
    "chain B and (resseq 1 and name CA or resseq 2 and name CA)"
  phenix.angle model.pdb "chain A" "chain B"

Feedback:
  PAfonine@lbl.gov
  phenixbb@phenix-online.org"""

def process_args(args):
  pdb_file_name, line_selections = None,[]
  if(len(args) != 3):
    raise Sorry(
      "Three arguments expected: PDB file, two atom selections to define axes.")
  for arg in args:
    if(os.path.isfile(arg) and iotbx.pdb.is_pdb_file(file_name=arg)):
      pdb_file_name = arg
    else:
      line_selections.append(arg)
  if(pdb_file_name is None):
    raise Sorry("PDB file must be provided.")
  ph = iotbx.pdb.input(file_name = pdb_file_name).construct_hierarchy()
  asc = ph.atom_selection_cache()
  # If it is protein: leave only backbone
  if(len(list(ph.chains()))==1):
    chain = list(ph.chains())[0]
    if(chain.is_protein()):
      sel = asc.selection(
        string = "pepnames and (name ca or name n or name c or name o)")
      ph = ph.select(sel)
      asc = ph.atom_selection_cache()
  #
  if(len(line_selections) != 2):
    raise Sorry("Two atom selections to define two axes must be provided.")
  try:
    sel1 = asc.selection(string=line_selections[0])
  except Exception:
    raise Sorry("Invalid atom selection: %s"%line_selections[0])
  try:
    sel2 = asc.selection(string=line_selections[1])
  except Exception:
    raise Sorry("Invalid atom selection: %s"%line_selections[1])
  for sel, ls in zip([sel1, sel2], line_selections):
    if(sel.count(True)<2):
      raise Sorry(
        "Atom selection '%s' selects less than two points."%ls)
  return ph, asc, sel1, sel2

def vector_from_two_points(s1, s2):
  a = [s2[0]-s1[0], s2[1]-s1[1], s2[2]-s1[2]]
  norm = math.sqrt(a[0]**2 + a[1]**2 + a[2]**2)
  if(abs(norm)<1.e-9):
    raise Sorry("Two points defining axis coincide.")
  return matrix.col((a[0]/norm, a[1]/norm, a[2]/norm))

def get_axis_from_xrs(xrs):
  if xrs.scatterers().size() > 2:
    sites_cart_moving = xrs.sites_cart()-xrs.center_of_mass()
    es = scitbx.math.principal_axes_of_inertia(points=sites_cart_moving).eigensystem()
    vecs = es.vectors()
    axis = vecs[6], vecs[7], vecs[8]
    return matrix.col((axis[0], axis[1], axis[2]))
  elif xrs.scatterers().size() == 2:
    sites_cart = xrs.sites_cart()
    assert sites_cart.size() == 2
    s1,s2 = sites_cart[0], sites_cart[1]
    return vector_from_two_points(s1, s2)
  return None

def calculate_axes_and_angle(xrs1, xrs2):
  a1 = get_axis_from_xrs(xrs1)
  a2 = get_axis_from_xrs(xrs2)
  angle = a1.angle(a2)*180./math.pi
  return a1, a2, angle

def calculate_axes_and_angle_directional(xrs1, xrs2):
  """ The same as above, but check the direction of vectors assuming that
  atoms in xrs1 and xrs2 are ordered from N to C terminus. This will enable
  the function to produce angles > 90 degrees.
  Uses rough estimate of direction of xrs using the first and the last atom
  and inverting vectors a1 and a2 when necessary.
  Used in iotbx.pdb.secondary_structure:concatenate_consecutive_helices()"""
  a1 = get_axis_from_xrs(xrs1)
  a2 = get_axis_from_xrs(xrs2)
  v1 = vector_from_two_points(xrs1.sites_cart()[0], xrs1.sites_cart()[-1])
  v2 = vector_from_two_points(xrs2.sites_cart()[0], xrs2.sites_cart()[-1])
  if a1.angle(v1)*180./math.pi > 90:
    a1 = -a1
  if a2.angle(v2)*180./math.pi > 90:
    a2 = -a2
  angle = a1.angle(a2)*180./math.pi
  return a1, a2, angle

def run(args, log=sys.stdout):
  if(len(args)==0 or (len(args)==1 and
     ("-h" in args or "--h" in args or "-help" in args or "--help" in args))):
    print("-"*79, file=log)
    print(legend, file=log)
    print("-"*79, file=log)
    sys.exit(0)
  ph, asc, sel1, sel2 = process_args(args=args)
  sel12 = sel1 | sel2
  xrs = ph.extract_xray_structure()
  #
  xrs  = xrs.select(sel12)
  sel1 = sel1.select(sel12)
  sel2 = sel2.select(sel12)
  xrs  = xrs.orthorhombic_unit_cell_around_centered_scatterers(buffer_size = 3)
  #
  a1, a2, angle = calculate_axes_and_angle(xrs.select(sel1), xrs.select(sel2))
  print("Axis 1: %6.4f %6.4f %6.4f"%(a1[0], a1[1], a1[2]), file=log)
  print("Axis 2: %6.4f %6.4f %6.4f"%(a2[0], a2[1], a2[2]), file=log)
  print("Angle : %6.4f" % angle, file=log)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])

