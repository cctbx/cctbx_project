from __future__ import division
import sys
import time
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
import mmtbx.utils
from mmtbx import monomer_library
import iotbx.phil
#from mmtbx import hydrogens
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
    htype      = None,  # type of bond
    a0         = None,  # parent atom index
    a1         = None,  # 1-3 neighbor index
    a2         = None,  # 1-3 or 1-4 neighbor index
    a          = None,  # coefficient for reconstruction
    b          = None,  # coefficient for reconstruction
    h          = None,  # coefficient for reconstruction
    phi        = None,  # angle
    alpha      = None,  # angle
    dist_h     = None): # measured or ideal distance
    adopt_init_args(self, locals())


# for every H atom, determine type of bond
def get_h_parameterization(connectivity, sites_cart, idealize):
  h_parameterization = {}
  for ih in connectivity.keys():
    a0 = connectivity[ih][0]
    count_H, reduced_neighbs = a0.count_H, a0.reduced_neighbs
    n_red_neigbs = len(reduced_neighbs)
    rh = matrix.col(sites_cart[ih])
    r0 = matrix.col(sites_cart[a0.iseq])
    if idealize:
      dist_h = a0.dist_ideal
    else:
      dist_h = (r0 - rh).length()
  # case 2a and 2b, case 3
    if(n_red_neigbs == 2 or
      (n_red_neigbs == 3 and count_H == 0)):
      a1, a2  = reduced_neighbs[0], reduced_neighbs[1]
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
      h_parameterization[ih] = parameterization_info(
        a0     = a0.iseq,
        a1     = a1.iseq,
        a2     = a2.iseq,
        a      = a,
        b      = b,
        dist_h = dist_h)
      #if ((sumang < 361 and sumang > 359) and idealize == True ):
      #if (0):
      #if (sumang < 361 and sumang > 359):
      if (sumang < (2*math.pi + 0.05) and (sumang > 2*math.pi - 0.05)):
        h_parameterization[ih].htype = 'flat_2neigbs'
      else:
        root = 1-c1*c1-c2*c2-c0*c0+2*c0*c1*c2
        if(root < 0):
          raise RuntimeError(
            "Expression in square root < 0 in get_h_parameterization.")
        denom = math.sin(alpha0)
        if(denom==0):
          raise RuntimeError(
            "Denominator zero: sin(alpha0)in get_h_parameterization.")
        cz = (math.sqrt(1-c1*c1-c2*c2-c0*c0+2*c0*c1*c2))/math.sin(alpha0)
        h = cz/math.sin(alpha0)
        #test if vector v points to same 'side' as uh0
        if((u10.cross(u20)).dot(uh0) < 0):
          h = -h
        h_parameterization[ih].h = h
        if (n_red_neigbs == 2): # case 2b
          h_parameterization[ih].htype = '2neigbs'
        elif (n_red_neigbs == 3): # case 3
          h_parameterization[ih].htype = '3neigbs'
        #if(count_H == 1):
        #  print h_parameterization[ih].htype
    # case 1a
    elif(n_red_neigbs == 1 and count_H == 1 and len(connectivity[ih][2])==2):
      neigbs_14 = connectivity[ih][2]
      a1 = reduced_neighbs[0]
      b1, b2 = neigbs_14[0], neigbs_14[1]
      r1 = matrix.col(sites_cart[a1.iseq])
      rb1 = matrix.col(sites_cart[b1.iseq])
      rb2 = matrix.col(sites_cart[b2.iseq])
      # chose 1-4 neighbor which is closer - important!
      if((rh-rb2).length() < (rh-rb1).length()):
        neigbs_14[0], neigbs_14[1] = neigbs_14[1], neigbs_14[0]
      rb1 = matrix.col(sites_cart[neigbs_14[0].iseq])
      r2 = r0 + (r1 - rb1)
      uh0 = (rh - r0).normalize()
      u10 = (r1 - r0).normalize()
      u20 = (r2 - r0).normalize()
      alpha0 = (u10).angle(u20)
      alpha1 = (u10).angle(uh0)
      alpha2 = (uh0).angle(u20)
      c0, c1, c2 = math.cos(alpha0), math.cos(alpha1), math.cos(alpha2)
      denom = (1-c0*c0)
      if(denom==0):
        raise RuntimeError(
          "Denominator zero: (1-c0*c0) in get_h_parameterization.")
      a = (c1-c0*c2)/(1-c0*c0)
      b = (c2-c0*c1)/(1-c0*c0)
      root = 1-c1*c1-c2*c2-c0*c0+2*c0*c1*c2
      if(root < 0):
        raise RuntimeError(
          "Expression in square root < 0 in get_h_parameterization.")
      denom = math.sin(alpha0)
      if(denom==0):
        raise RuntimeError(
          "Denominator zero: sin(alpha0)in get_h_parameterization.")
      cz = (math.sqrt(1-c1*c1-c2*c2-c0*c0+2*c0*c1*c2))/math.sin(alpha0)
      h = cz/math.sin(alpha0)
      # test if vector v points to same 'side' as uh0
      if((u10.cross(u20)).dot(uh0) < 0):
        h = -h
      h_parameterization[ih] = parameterization_info(
        htype  = 'alg1a',
        a0     = a0.iseq,
        a1     = a1.iseq,
        a2     = neigbs_14[0].iseq,
        a      = a,
        b      = b,
        h      = h,
        dist_h = dist_h)
    # case 1b
    elif(n_red_neigbs == 1 and (count_H == 0 or count_H ==2)):
      a1 = reduced_neighbs[0]
      sec_neigbs = connectivity[ih][2]
      b1 = sec_neigbs[0]
      r1 = matrix.col(sites_cart[a1.iseq])
      rb1 = matrix.col(sites_cart[b1.iseq])
      uh0 = (rh - r0).normalize()
      u10 = (r1 - r0).normalize()
      if idealize:
        alpha = math.radians(a1.angle_ideal)
      else:
        alpha = (u10).angle(uh0)
      phi = dihedral_angle(
        sites=[sites_cart[ih], sites_cart[a0.iseq],
        sites_cart[a1.iseq],sites_cart[b1.iseq]])
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
        alpha  = alpha,
        dist_h = dist_h)
    else:
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
    uh0 = (a * u10 + b * u20)
    rH_gen = r0 + dh * uh0
    deltaH = (rH_gen - matrix.col(sites_cart[ih])).length()
  # alg2b and alg3
  elif (para_info.htype == '2neigbs' or para_info.htype == '3neigbs'):
    a, b, h = para_info.a, para_info.b, para_info.h
    r2 = matrix.col(sites_cart[para_info.a2])
    u10, u20 = (r1 - r0).normalize(), (r2 - r0).normalize()
    v = u10.cross(u20)
    uh0 = (a * u10 + b * u20 + h * v)
    rH_gen = r0 + dh * uh0
    deltaH = (rH_gen - matrix.col(sites_cart[ih])).length()
  # alg1a
  elif (para_info.htype == 'alg1a'):
    a, b, h = para_info.a, para_info.b, para_info.h
    rb1 = matrix.col(sites_cart[para_info.a2])
    r2 = r0 + (r1 - rb1)
    u10, u20 = (r1 - r0).normalize(), (r2 - r0).normalize()
    v = u10.cross(u20)
    uh0 = (a * u10 + b * u20 + h * v)
    rH_gen = r0 + dh * uh0
    deltaH = (rH_gen - matrix.col(sites_cart[ih])).length()
#
  elif (para_info.htype == 'alg1b'):
    rb1 = matrix.col(sites_cart[para_info.a2])
    phi = para_info.phi
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
  hd_selection = xray_structure.hd_selection()
  names = list(pdb_hierarchy.atoms().extract_name())
  sites_cart = xray_structure.sites_cart()
  #scatterers = xray_structure.scatterers()

  atoms_list = list(pdb_hierarchy.atoms_with_labels())

  print >>log, '\nNow determining connectivity table for H atoms...'
  connectivity = hydrogen_connectivity.determine_H_neighbors(
    geometry_restraints   = geometry_restraints,
    bond_proxies          = bond_proxies_simple,
    angle_proxies         = angle_proxies,
    hd_selection          = hd_selection,
    sites_cart            = sites_cart)

  print >>log, '\nNow determining the parameterization for H atoms...'
  h_parameterization = get_h_parameterization(
    connectivity   = connectivity,
    sites_cart     = sites_cart,
    idealize       = True)

  print >>log, '\nNow reconstructing H atoms...'
  long_distance_list = []
  unk_list = []
  for ih in h_parameterization.keys():
    residue = atoms_list[ih].resseq
    hp = h_parameterization[ih]
    h_obj = generate_H_positions(
      sites_cart        = sites_cart,
      ih                = ih,
      para_info         = hp)
    if(h_obj.distance is not None):
      print >> log, hp.htype, 'atom:', names[ih]+' ('+str(ih)+ ') residue:', \
        residue, 'distance:', h_obj.distance
      if(h_obj.distance > 0.05):
        long_distance_list.append(ih)
    else:
      print >> log, hp.htype, 'atom:', names[ih]+' ('+str(ih)+ ') residue:', residue
      unk_list.append(ih)

  # some informative output for residues with unknown algorithm
  if unk_list:
    print >>log, '*'*79
    print >>log, 'Warning: The following atoms where not assigned an H type'
    for ih in unk_list:
      residue = atoms_list[ih].resseq
      hp = h_parameterization[ih]
      print >> log, 'atom:', names[ih], 'residue:', residue, \
        'chain', atoms_list[ih].chain_id
    print >>log, '*'*79

  # some informative output for residues where position is NOT reproduced
  # -> wronlgy assigned
  if long_distance_list:
    print >>log, '*'*79
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
  print >>log, '*'*79

  #for ih in h_parameterization.keys():
  #  hp = h_parameterization[ih]
  #  print 'htype = ', hp.htype, 'a0 = ', hp.a0, 'a1 = ', hp.a1, 'a2 = ', hp.a2, \
  #    'a = ', hp.a, 'b = ', hp.b, 'h = ', hp.h, 'phi = ', hp.phi, \
  #    'alpha = ', hp.alpha, 'dist_h =', hp.dist_h

if (__name__ == "__main__"):
  t0 = time.time()
  run(sys.argv[1:])
  print "Time:", round(time.time()-t0, 2)
