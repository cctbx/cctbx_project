from __future__ import division
import sys, os
import time
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
import mmtbx.utils
from mmtbx import monomer_library
import iotbx.phil
from cctbx import geometry_restraints
from scitbx import matrix
#from libtbx.utils import Sorry
from libtbx import adopt_init_args
from libtbx import group_args
from libtbx.utils import null_out
from libtbx.utils import multi_out
from stdlib import math
from scitbx.math import dihedral_angle

import hydrogen_connectivity

legend = """\
phenix.hydrogen_parameterization:
Computes the parameters for the compuationf of H atom positions
based on their parent atoms

Inputs:
  - Model file in PDB format (other files are not supported)

Usage examples:
   phenix.hydrogen_parameterization model.pdb

Output:
   This code produces a dictionary
   H atom index is the key, pointing to an object which contains necessary
   parameters to build the H atom again (knowing the coordinates of
   the parent atoms)
   H (iseq) --> parameterization_info

Important:
  Hydrogen atoms need to be present in the file.
  H atoms can be generated with phenix.reduce
"""

class parameterization_info(object):
  def __init__(
    self,
    htype      = None,  # type of hydrogen environment
    a0         = None,  # parent atom index
    a1         = None,  # 1-3 neighbor index
    a2         = None,  # 1-3 or 1-4 neighbor index
    a3         = None,  # 1-3 or 1-4 neighbor index
    a          = None,  # coefficient for reconstruction
    b          = None,  # coefficient for reconstruction
    h          = None,  # coefficient for reconstruction
    phi        = None,  # angle
    alpha      = None,  # angle
    n          = None,  # parameter for propeller and H2 groups
    dist_h     = None): # measured or ideal distance
    adopt_init_args(self, locals())

# this function catches three cases:
# 1. planar geometry
# 2. two tetragonal CH2 geometry
# 3. H out of plane of its 3 neighbors (should be rare and not in AA)
def get_coefficients(ih, a0, a1, a2, ih2, idealize, sites_cart, typeh):
  rh = matrix.col(sites_cart[ih])
  r0 = matrix.col(sites_cart[a0.iseq])
  r1 = matrix.col(sites_cart[a1.iseq])
  r2 = matrix.col(sites_cart[a2.iseq])
  uh0 = (rh - r0).normalize()
  u10 = (r1 - r0).normalize()
  u20 = (r2 - r0).normalize()
  if idealize:
    alpha0 = math.radians(a0.angle_ideal)
    alpha1 = math.radians(a1.angle_ideal)
    alpha2 = math.radians(a2.angle_ideal)
    c0, c1, c2 = math.cos(alpha0), math.cos(alpha1), math.cos(alpha2)
  else:
    alpha0 = (u10).angle(u20)
    alpha1 = (u10).angle(uh0)
    alpha2 = (uh0).angle(u20)
    c0 = (u10).dot(u20)
    c1 = (u10).dot(uh0)
    c2 = (uh0).dot(u20)
  sumang = alpha0 + alpha1 + alpha2
  denom = (1.0-c0**2)
  if(denom==0):
    raise RuntimeError(
      "Denominator zero: (1-c0*c0) in get_h_parameterization.")
  a = (c1-c0*c2)/(1-c0*c0)
  b = (c2-c0*c1)/(1-c0*c0)
  root = None
  # check if H, A0, A1, A2 are in a plane
  if (sumang < (2*math.pi + 0.05) and (sumang > 2*math.pi - 0.05)):
    h = None
  elif (sumang > (2*math.pi + 0.05) and 1-c1*c1-c2*c2-c0*c0+2*c0*c1*c2 < 0):
    root = 1-c1*c1-c2*c2-c0*c0+2*c0*c1*c2
    h = None
    return sumang, a, b, h, root
  else:
    # two tetragonal geometry: e.g. CH2 group
    if (ih2 is not None):
      rh2 = matrix.col(sites_cart[ih2.iseq])
      uh02 = (rh2 - r0).normalize()
      if idealize:
        h = math.radians(ih2.angle_ideal) * 0.5
      else:
        h = (uh0).angle(uh02) * 0.5
      #test if vector v points to same 'side' as uh0
      if((u10.cross(u20)).dot(uh0) < 0):
        h =  -h
    else:
    # if H is out of plane, but not in tetrahedral geometry
      root = 1-c1*c1-c2*c2-c0*c0+2*c0*c1*c2
      if(root < 0):
        raise RuntimeError(
          "Expression in square root < 0 in get_h_parameterization.")
      denom = math.sin(alpha0)
      if(denom==0):
        raise RuntimeError(
          "Denominator zero: sin(alpha0)in get_h_parameterization.")
      cz = (math.sqrt(1-c1*c1-c2*c2-c0*c0+2*c0*c1*c2))/math.sin(alpha0)
      h = cz
      #test if vector v points to same 'side' as uh0
      if((u10.cross(u20)).dot(uh0) < 0):
        h = -h
  return sumang, a, b, h, root

# obtain coefficients for tetragonal H (such as HA) using Cramer's rule
def get_coefficients_alg3(rh, a0, a1, a2, a3, idealize, sites_cart):
  r0 = matrix.col(sites_cart[a0.iseq])
  r1 = matrix.col(sites_cart[a1.iseq])
  r2 = matrix.col(sites_cart[a2.iseq])
  r3 = matrix.col(sites_cart[a3.iseq])
  uh0 = (rh - r0).normalize()
  u10 = (r1 - r0).normalize()
  u20 = (r2 - r0).normalize()
  u30 = (r3 - r0).normalize()
  if idealize:
    alpha0 = math.radians(a1.angle_ideal)
    alpha1 = math.radians(a2.angle_ideal)
    alpha2 = math.radians(a3.angle_ideal)
    c1, c2, c3 = math.cos(alpha0), math.cos(alpha1), math.cos(alpha2)
    omega0 = math.radians(a0.angle_ideal[0])
    omega1 = math.radians(a0.angle_ideal[1])
    omega2 = math.radians(a0.angle_ideal[2])
    w12, w23, w13 = math.cos(omega0), math.cos(omega1), math.cos(omega2)
  else:
    c1 = (uh0).dot(u10)
    c2 = (uh0).dot(u20)
    c3 = (uh0).dot(u30)
    w12 = (u10).dot(u20)
    w23 = (u20).dot(u30)
    w13 = (u10).dot(u30)
  matrix_d = matrix.sqr([
    1,   w12, w13,
    w12, 1,   w23,
    w13, w23, 1   ])
  #
  matrix_x = matrix.sqr([
    c1, w12, w13,
    c2, 1,   w23,
    c3, w23, 1   ])
  #
  matrix_y = matrix.sqr([
    1,   c1,  w13,
    w12, c2,  w23,
    w13, c3,  1   ])
  #
  matrix_z = matrix.sqr([
    1,   w12,  c1,
    w12, 1,    c2,
    w13, w23,  c3 ])
  if(matrix_d.determinant()==0):
    raise RuntimeError(
      "Denominator zero: matrix_d in get_h_parameterization.")
  a = matrix_x.determinant()/matrix_d.determinant()
  b = matrix_y.determinant()/matrix_d.determinant()
  c = matrix_z.determinant()/matrix_d.determinant()
  return a, b, c


# for every H atom, determine the type of bond
# for debugging
#def get_h_parameterization(connectivity, sites_cart, idealize, atoms_list, names):
def get_h_parameterization(connectivity, sites_cart, idealize):
  h_parameterization = {}
  for ih in connectivity.keys():
   #for debugging
    #print 'atom:', names[ih]+' ('+str(ih)+ ') residue:', \
    #  atoms_list[ih].resseq, 'chain', atoms_list[ih].chain_id
    # if entry exists already, skip it
    if ih in h_parameterization:
      continue
    a0 = connectivity[ih][0]
    count_H  = a0.count_H
    reduced_neighbs = connectivity[ih][1]
    n_red_neigbs = len(reduced_neighbs)
    rh = matrix.col(sites_cart[ih])
    r0 = matrix.col(sites_cart[a0.iseq])
    if idealize:
      dist_h = a0.dist_ideal
    else:
      dist_h = (r0 - rh).length()
  # alg2a, 2tetra, 2neigbs
    if(n_red_neigbs == 2):
      a1, a2  = reduced_neighbs[0], reduced_neighbs[1]
      # if H is second neighbor, gets its index
      if (count_H == 1):
        hlist = connectivity[ih][2]
        if hlist:
          ih2 = (hlist[0])
          i_h2 = (hlist[0]).iseq
      else:
        ih2 = None
      sumang, a, b, h, root = get_coefficients(
        ih         = ih,
        a0         = a0,
        a1         = a1,
        a2         = a2,
        ih2        = ih2,
        idealize   = idealize,
        sites_cart = sites_cart,
        typeh      = 'alg2')
      h_parameterization[ih] = parameterization_info(
        a0     = a0.iseq,
        a1     = a1.iseq,
        a2     = a2.iseq,
        a      = a,
        b      = b,
        dist_h = dist_h)
      # alg2a
      if (sumang > (2*math.pi + 0.05) and root < 0):
        h_parameterization[ih].htype = 'unk_ideal'
      elif (sumang < (2*math.pi + 0.05) and (sumang > 2*math.pi - 0.05)):
        h_parameterization[ih].htype = 'flat_2neigbs'
      else:
        if (count_H == 1):
        # 2 tetragonal geometry
          h_parameterization[ih].htype = '2tetra'
          h_parameterization[ih].alpha = h
          h_parameterization[i_h2] = parameterization_info(
            a0     = a0.iseq,
            a1     = a1.iseq,
            a2     = a2.iseq,
            a      = a,
            b      = b,
            alpha  = -h,
            dist_h = dist_h,
            htype  = '2tetra')
        else:
          # 2neigbs
          h_parameterization[ih].h = h
          h_parameterization[ih].htype = '2neigbs'
    # tetragonal geometry: 3neigbs
    elif (n_red_neigbs == 3 and count_H == 0):
      a1, a2, a3  = reduced_neighbs[0], reduced_neighbs[1], reduced_neighbs[2]
      a, b, h = get_coefficients_alg3(
        rh         = rh,
        a0         = a0,
        a1         = a1,
        a2         = a2,
        a3         = a3,
        idealize   = idealize,
        sites_cart = sites_cart)
      h_parameterization[ih] = parameterization_info(
        a0     = a0.iseq,
        a1     = a1.iseq,
        a2     = a2.iseq,
        a3     = a3.iseq,
        a      = a,
        b      = b,
        h      = h,
        dist_h = dist_h,
        htype  = '3neigbs')
    # alg1a: X-H2 planar groups, such as in ARG, ASN, GLN
    # requires that dihedral angle restraint exists
    elif(n_red_neigbs == 1 and count_H == 1 and len(connectivity[ih])==4):
      a1 = reduced_neighbs[0]
      r1 = matrix.col(sites_cart[a1.iseq])
      hlist = connectivity[ih][2]
      ih_2 = hlist[0].iseq
      #if(len(connectivity[ih])!=4):
      #  continue
      b1 = (connectivity[ih][3])[0]
      if (b1.dihedral_ideal == None):
        continue
      #iseq_b1 = b1.iseq
      dihedral = dihedral_angle(
        sites=[sites_cart[b1.iseq], sites_cart[a1.iseq],
        sites_cart[a0.iseq],sites_cart[ih]])
      #print 'dihedrals', a0.dihedral, dihedral, a0.dihedral_ideal
      rb1 = matrix.col(sites_cart[b1.iseq])
      uh0 = (rh - r0).normalize()
      u10 = (r1 - r0).normalize()
      if idealize:
        alpha = math.radians(a1.angle_ideal)
        phi = math.radians(b1.dihedral_ideal)
        #phi = dihedral
      else:
        alpha = (u10).angle(uh0)
        #phi = a0.dihedral
        phi = dihedral
      u1 = (r0 - r1).normalize()
      rb10 = rb1 - r1
      u2 = (rb10 - ((rb10).dot(u10)) * u10).normalize()
      u3 = u1.cross(u2)
      h_parameterization[ih] = parameterization_info(
        htype  = 'alg1a',
        a0     = a0.iseq,
        a1     = a1.iseq,
        a2     = b1.iseq,
        phi    = phi,
        n      = 0,
        alpha  = alpha,
        dist_h = dist_h)
      h_parameterization[ih_2] = parameterization_info(
        htype  = 'alg1a',
        a0     = a0.iseq,
        a1     = a1.iseq,
        a2     = b1.iseq,
        phi    = phi+math.pi,
        n      = 0,
        alpha  = alpha,
        dist_h = dist_h)
    # case 1b
# a0.dihedral = dihedral angle between angle ideal and actual position
    elif(n_red_neigbs == 1 and (count_H == 0 or count_H ==2)):
      #if(count_H == 0 and len(connectivity[ih])!=4):
      #  print 'the culprit is ', ih, names[ih], atoms_list[ih].resseq
      if (len(connectivity[ih])!=4):
        continue
      a1 = reduced_neighbs[0]
      b1 = (connectivity[ih][3])[0]
      r1 = matrix.col(sites_cart[a1.iseq])
      rb1 = matrix.col(sites_cart[b1.iseq])
      uh0 = (rh - r0).normalize()
      u10 = (r1 - r0).normalize()
      dihedral = dihedral_angle(
        sites=[sites_cart[ih], sites_cart[a0.iseq],
        sites_cart[a1.iseq],sites_cart[b1.iseq]])
      if idealize:
        alpha = math.radians(a1.angle_ideal)
        #phi = math.radians(b1.dihedral_ideal)
        #allow for rotation even for idealize = True
        phi = dihedral
      else:
        alpha = (u10).angle(uh0)
        phi = dihedral
      u1 = (r0 - r1).normalize()
      rb10 = rb1 - r1
      u2 = (rb10 - ((rb10).dot(u1)) * u1).normalize()
      u3 = u1.cross(u2)
      h_parameterization[ih] = parameterization_info(
        htype  = 'alg1b',
        a0     = a0.iseq,
        a1     = a1.iseq,
        a2     = b1.iseq,
        phi    = phi,
        n      = 0,
        alpha  = alpha,
        dist_h = dist_h)
      if (count_H == 2):
        h_parameterization[ih].htype = 'prop'
        hlist = connectivity[ih][2]
        # TO DO: Can the order be reversed? To be kept in mind!!
        ih_2, ih_3 = hlist[0].iseq, hlist[1].iseq
        h_parameterization[ih_2] = parameterization_info(
          htype  = 'prop',
          a0     = a0.iseq,
          a1     = a1.iseq,
          a2     = b1.iseq,
          phi    = phi,
          n      = 1,
          alpha  = alpha,
          dist_h = dist_h)
        ih_2_coord = generate_H_positions(
          sites_cart        = sites_cart,
          ih                = ih_2,
          para_info         = h_parameterization[ih_2]).rH_gen
        h_parameterization[ih_3] = parameterization_info(
          htype  = 'prop',
          a0     = a0.iseq,
          a1     = a1.iseq,
          a2     = b1.iseq,
          phi    = phi,
          n      = 2,
          alpha  = alpha,
          dist_h = dist_h)
        if ((ih_2_coord - matrix.col(sites_cart[ih_3])).length() <
          (ih_2_coord - matrix.col(sites_cart[ih_2])).length() ):
          h_parameterization[ih_2].n = 2
          h_parameterization[ih_3].n = 1
    else:
      a1 = reduced_neighbs[0]
      h_parameterization[ih] = parameterization_info(
        htype  = 'unk',
        a0     = a0.iseq,
        a1     = a1.iseq)
  return h_parameterization


def generate_H_positions(sites_cart, ih, para_info):
  r0 = matrix.col(sites_cart[para_info.a0])
  r1 = matrix.col(sites_cart[para_info.a1])
  dh = para_info.dist_h
  a, b, h = para_info.a, para_info.b, para_info.h
  # alg2a
  if (para_info.htype == 'flat_2neigbs'):
    a, b = para_info.a, para_info.b
    r2 = matrix.col(sites_cart[para_info.a2])
    u10, u20 = (r1 - r0).normalize(), (r2 - r0).normalize()
    length = math.sqrt(a*a + b*b + 2*a*b*(u10).dot(u20))
    if(length==0):
      raise RuntimeError("Denominator zero: length in generate_H_positions")
    uh0 = (a * u10 + b * u20)/length
    rH_gen = r0 + dh * uh0
    deltaH = (rH_gen - matrix.col(sites_cart[ih])).length()
  # 2 neigbs
  elif (para_info.htype == '2neigbs'):
    a, b, h = para_info.a, para_info.b, para_info.h
    r2 = matrix.col(sites_cart[para_info.a2])
    u10, u20 = (r1 - r0).normalize(), (r2 - r0).normalize()
    v0 = (u10.cross(u20)).normalize()
    rh0 = (a * u10 + b * u20 + h * v0)
    length = math.sqrt(rh0.dot(rh0))
    if(length==0):
      raise RuntimeError("Denominator zero: length in generate_H_positions")
    uh0 = rh0/length
    rH_gen = r0 + dh * uh0
    deltaH = (rH_gen - matrix.col(sites_cart[ih])).length()
  # 2tetrahedral
  elif (para_info.htype == '2tetra'):
    a, b, delta = para_info.a, para_info.b, para_info.alpha
    r2 = matrix.col(sites_cart[para_info.a2])
    u10, u20 = (r1 - r0).normalize(), (r2 - r0).normalize()
    v0 = (u10.cross(u20)).normalize()
    d0 = (a * u10 + b * u20).normalize()
    rH_gen = r0 + dh * (math.cos(delta) * d0 + math.sin(delta) * v0)
    deltaH = (rH_gen - matrix.col(sites_cart[ih])).length()
  # tetragonal alg3
  elif (para_info.htype == '3neigbs'):
    a, b, h = para_info.a, para_info.b, para_info.h
    r2 = matrix.col(sites_cart[para_info.a2])
    r3 = matrix.col(sites_cart[para_info.a3])
    u10 = (r1 - r0).normalize()
    u20 = (r2 - r0).normalize()
    u30 = (r3 - r0).normalize()
    rh0 = (a*u10 + b*u20 + h*u30)
    length = math.sqrt(rh0.dot(rh0))
    if(length==0):
      raise RuntimeError("Denominator zero: length in generate_H_positions")
    uh0 = rh0/length
    rH_gen = r0 + dh * uh0
    deltaH = (rH_gen - matrix.col(sites_cart[ih])).length()
# alg1b or alg1a or propeller group
  elif (para_info.htype in ['alg1b', 'alg1a', 'prop']):
    rb1 = matrix.col(sites_cart[para_info.a2])
    n = para_info.n
    phi = para_info.phi + n*2*math.pi/3
    #phi = para_info.phi
    alpha = para_info.alpha
    salpha = math.sin(alpha)
    calpha = math.cos(alpha)
    sphi = math.sin(phi)
    cphi = math.cos(phi)
    u1 = (r0 - r1).normalize()
    rb10 = rb1 - r1
    u2 = (rb10 - ((rb10).dot(u1)) * u1).normalize()
    u3 = u1.cross(u2)
    rH_gen = r0 + dh * (salpha*(cphi*u2 + sphi*u3) - calpha*u1)
    deltaH = (rH_gen - matrix.col(sites_cart[ih])).length()
  else:
    deltaH, rH_gen = None, None
  return group_args(
    distance = deltaH,
    rH_gen   = rH_gen)

def run(args, out=sys.stdout):
  log = multi_out()
  log.register("stdout", out)
  log_file_name = "hydrogen_parameterization.log"
  logfile = open(log_file_name, "w")
  log.register("logfile", logfile)
  if (len(args) == 0):
    print >>log, legend
    return

  print >> log, "phenix.hydrogen_parameterization is running..."
  print >> log, "input parameters:\n", args

# parse through params --> switch off CDL or not
  params_line = grand_master_phil_str
  params = iotbx.phil.parse(
      input_string=params_line, process_includes=True).extract()
  params.pdb_interpretation.restraints_library.cdl=False

  processed_args = mmtbx.utils.process_command_line_args(
    args=args, log=null_out())
  pdb_filename = processed_args.pdb_file_names[0]
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  if pdb_filename is not None:
    processed_pdb_file = monomer_library.pdb_interpretation.process(
      params         = params.pdb_interpretation,
      mon_lib_srv    = mon_lib_srv,
      ener_lib       = ener_lib,
      file_name      = pdb_filename,
      force_symmetry = True)

  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  xray_structure = processed_pdb_file.xray_structure()
  if pdb_filename is not None:
    pdb_str = pdb_hierarchy.as_pdb_string()

  geometry_restraints = processed_pdb_file.geometry_restraints_manager(
    show_energies = False)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry_restraints,
    normalization = False)

  bond_proxies_simple, asu = restraints_manager.geometry.get_all_bond_proxies(
    sites_cart = xray_structure.sites_cart())
  angle_proxies = restraints_manager.geometry.get_all_angle_proxies()
  dihedral_proxies = restraints_manager.geometry.dihedral_proxies
  hd_selection = xray_structure.hd_selection()
  names = list(pdb_hierarchy.atoms().extract_name())
  sites_cart = xray_structure.sites_cart()
  atoms = pdb_hierarchy.atoms()
  names = list(atoms.extract_name())
  #scatterers = xray_structure.scatterers()
  atoms_list = list(pdb_hierarchy.atoms_with_labels())

  idealize = True

  print >>log, '\nNow determining connectivity table for H atoms...'
  connectivity = hydrogen_connectivity.determine_H_neighbors(
    geometry_restraints = geometry_restraints,
    bond_proxies        = bond_proxies_simple,
    angle_proxies       = angle_proxies,
    dihedral_proxies    = dihedral_proxies,
    hd_selection        = hd_selection,
    sites_cart          = sites_cart,
    atoms               = atoms)

  print >>log, '\nNow determining the parameterization for H atoms...'
  #h_parameterization = get_h_parameterization(
  #  connectivity   = connectivity,
  #  sites_cart     = sites_cart,
  #  idealize       = idealize)
  # for debugging
  h_parameterization = get_h_parameterization(
    connectivity   = connectivity,
    sites_cart     = sites_cart,
    idealize       = idealize,
    atoms_list     = atoms_list,
    names          = names)

  print >>log, '\nNow reconstructing H atoms...'
  # H atoms for which distance compared to input model is large
  long_distance_list = []
  # H atoms with unknown parameterization
  unk_list = []
  # H atoms with nonsensical ideal angles
  unk_ideal_list = []
  for ih in h_parameterization.keys():
    residue = atoms_list[ih].resseq
    hp = h_parameterization[ih]
    if (hp.htype == 'unk'):
      unk_list.append(ih)
    elif (hp.htype == 'unk_ideal'):
      unk_ideal_list.append(ih)
    else:
      h_obj = generate_H_positions(
        sites_cart        = sites_cart,
        ih                = ih,
        para_info         = hp)
      if(h_obj.distance is not None):
        print >> log, hp.htype, 'atom:', names[ih]+' ('+str(ih)+ ') residue:', \
          residue, 'distance:', h_obj.distance
        if(h_obj.distance > 0.03):
          long_distance_list.append(ih)
  print >>log, '*'*79
  # list residues with unknown parameterization
  if unk_list:
    print >>log, 'Warning: The following atoms where not assigned an H type'
    for ih in unk_list:
      residue = atoms_list[ih].resseq
      hp = h_parameterization[ih]
      print >> log, 'atom:', names[ih], 'residue:', residue, \
        'chain', atoms_list[ih].chain_id
    print >>log, '*'*79

  # list residues with nonsensical angles
  if unk_ideal_list:
    print >>log, 'Warning: The following atoms have nonsensical ideal angles.'
    print >>log, 'Check the geo file for ideal angles involving these H atoms.'
    for ih in unk_ideal_list:
      residue = atoms_list[ih].resseq
      hp = h_parameterization[ih]
      print >> log, 'atom:', names[ih], 'residue:', residue, \
        'chain', atoms_list[ih].chain_id
    print >>log, '*'*79

  # list output for residues where position is not reproduced
  if long_distance_list:
    print >>log, 'Warning: The position of the following H atoms was not reproduced'
    for ih in long_distance_list:
      residue = atoms_list[ih].resseq
      hp = h_parameterization[ih]
      h_obj = generate_H_positions(
        sites_cart        = sites_cart,
        ih                = ih,
        para_info         = hp)
      if(h_obj.distance is not None and h_obj.distance > 0.05):
        print >> log, hp.htype, 'atom:', names[ih]+' ('+str(ih)+ ') residue:', \
          residue, 'chain', atoms_list[ih].chain_id, 'distance:', h_obj.distance
      sites_cart[ih] = h_obj.rH_gen
  xray_structure.set_sites_cart(sites_cart)
  pdb_hierarchy.adopt_xray_structure(xray_structure)

  if pdb_filename is not None:
    pdb_basename = os.path.basename(pdb_filename.split(".")[0])
    pdb_hierarchy.write_pdb_file(
      file_name        = pdb_basename+"_H.pdb",
      crystal_symmetry = xray_structure.crystal_symmetry())
  print >>log, '*'*79

  #for ih in h_parameterization.keys():
  #  hp = h_parameterization[ih]
  #  #if(ih != 220):
  #  #  continue
  #  print 'htype = ', hp.htype, 'a0 = ', hp.a0, 'a1 = ', hp.a1, 'a2 = ', hp.a2, \
  #    'a = ', hp.a, 'b = ', hp.b, 'h = ', hp.h, 'phi = ', hp.phi, \
  #    'alpha = ', hp.alpha, 'dist_h =', hp.dist_h

if (__name__ == "__main__"):
  t0 = time.time()
  run(sys.argv[1:])
  print "Time:", round(time.time()-t0, 2)
