from __future__ import division
from scitbx import matrix
#from scitbx.array_family import flex
#from libtbx.utils import Sorry
from libtbx import adopt_init_args
from libtbx import group_args
from stdlib import math
from scitbx.math import dihedral_angle

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

# this function determines parameters for three cases:
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
    alpha0 = math.radians(a0.angle_ideal[0])
    alpha1 = math.radians(a1.angle_ideal)
    alpha2 = math.radians(a2.angle_ideal)
    c0, c1, c2 = math.cos(alpha0), math.cos(alpha1), math.cos(alpha2)
  else:
    alpha0 = (u10).angle(u20)
    alpha0 = math.acos(u10.dot(u20))
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
def get_h_parameterization(h_connectivity, sites_cart, idealize):
  h_parameterization = {}
  n_atoms = len(sites_cart)
  for ih in h_connectivity.keys():
    #if (ih != 22):
    #  continue
   #for debugging
    #print 'atom:', names[ih]+' ('+str(ih)+ ') residue:', \
    #  atoms_list[ih].resseq, 'chain', atoms_list[ih].chain_id
    # if entry exists already, skip it
    if ih in h_parameterization.keys():
      continue
    a0 = h_connectivity[ih][0]
    count_H  = a0.count_H
    reduced_neighbs = h_connectivity[ih][1]
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
        hlist = h_connectivity[ih][2]
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
    # requires that dihedral angle restraint exists for at least one H atom
    elif(n_red_neigbs == 1 and count_H == 1 and len(h_connectivity[ih])==4):
      b1 = None
      a1 = reduced_neighbs[0]
      r1 = matrix.col(sites_cart[a1.iseq])
      hlist = h_connectivity[ih][2]
      ih_2 = hlist[0].iseq
      for b1_test in h_connectivity[ih][3]:
        if (b1_test.dihedral_ideal is None):
          continue
        else:
          b1 = b1_test
      if (b1 == None):
        for b1_test in h_connectivity[ih_2][3]:
          if (b1_test.dihedral_ideal is None):
            continue
          else:
            b1 = b1_test
            ih, ih_2 = ih_2, ih
      if (b1 == None):
        h_parameterization[ih] = parameterization_info(
          htype  = 'unk',
          a0     = a0.iseq)
        continue
      # check if angle is typical for propeller
      # catches case of missing propeller atom
      if (hlist[0].angle_ideal > 107 and
          hlist[0].angle_ideal < 111):
        h_parameterization[ih] = parameterization_info(
          htype  = 'unk',
          a0     = a0.iseq)
        h_parameterization[ih_2] = parameterization_info(
          htype  = 'unk',
          a0     = a0.iseq)
        continue
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
      #print names[ih], names[b1.iseq], names[a1.iseq]
      if ih not in h_parameterization:
        h_parameterization[ih] = parameterization_info(
          htype  = 'alg1a',
          a0     = a0.iseq,
          a1     = a1.iseq,
          a2     = b1.iseq,
          phi    = phi,
          n      = 0,
          alpha  = alpha,
          dist_h = dist_h)
      if ih_2  not in h_parameterization:
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
      if (len(h_connectivity[ih])!=4):
        continue
      a1 = reduced_neighbs[0]
      b1 = (h_connectivity[ih][3])[0]
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
        hlist = h_connectivity[ih][2]
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
        # check if order is reversed
        ih_2_coord = compute_H_position(
          sites_cart = sites_cart,
          ih         = ih_2,
          hp         = h_parameterization[ih_2])
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
      h_parameterization[ih] = parameterization_info(
        htype  = 'unk',
        a0     = a0.iseq)
  return h_parameterization

def compute_H_position(ih, sites_cart, hp):
  r0 = matrix.col(sites_cart[hp.a0])
  r1 = matrix.col(sites_cart[hp.a1])
  dh = hp.dist_h
  a, b, h = hp.a, hp.b, hp.h
  # alg2a
  if (hp.htype == 'flat_2neigbs'):
    a, b = hp.a, hp.b
    r2 = matrix.col(sites_cart[hp.a2])
    u10, u20 = (r1 - r0).normalize(), (r2 - r0).normalize()
    length = math.sqrt(a*a + b*b + 2*a*b*(u10).dot(u20))
    if(length==0):
      raise RuntimeError("Denominator zero: length in generate_H_positions")
    uh0 = (a * u10 + b * u20)/length
    rh_calc = r0 + dh * uh0
    h_distance = (rh_calc - matrix.col(sites_cart[ih])).length()
  # 2 neigbs
  elif (hp.htype == '2neigbs'):
    a, b, h = hp.a, hp.b, hp.h
    r2 = matrix.col(sites_cart[hp.a2])
    u10, u20 = (r1 - r0).normalize(), (r2 - r0).normalize()
    v0 = (u10.cross(u20)).normalize()
    rh0 = (a * u10 + b * u20 + h * v0)
    length = math.sqrt(rh0.dot(rh0))
    if(length==0):
      raise RuntimeError("Denominator zero: length in generate_H_positions")
    uh0 = rh0/length
    rh_calc = r0 + dh * uh0
    h_distance = (rh_calc - matrix.col(sites_cart[ih])).length()
  # 2tetrahedral
  elif (hp.htype == '2tetra'):
    a, b, delta = hp.a, hp.b, hp.alpha
    r2 = matrix.col(sites_cart[hp.a2])
    u10, u20 = (r1 - r0).normalize(), (r2 - r0).normalize()
    v0 = (u10.cross(u20)).normalize()
    d0 = (a * u10 + b * u20).normalize()
    rh_calc = r0 + dh * (math.cos(delta) * d0 + math.sin(delta) * v0)
    h_distance = (rh_calc - matrix.col(sites_cart[ih])).length()
  # tetragonal alg3
  elif (hp.htype == '3neigbs'):
    a, b, h = hp.a, hp.b, hp.h
    r2 = matrix.col(sites_cart[hp.a2])
    r3 = matrix.col(sites_cart[hp.a3])
    u10 = (r1 - r0).normalize()
    u20 = (r2 - r0).normalize()
    u30 = (r3 - r0).normalize()
    rh0 = (a*u10 + b*u20 + h*u30)
    length = math.sqrt(rh0.dot(rh0))
    if(length==0):
      raise RuntimeError("Denominator zero: length in generate_H_positions")
    uh0 = rh0/length
    rh_calc = r0 + dh * uh0
    h_distance = (rh_calc - matrix.col(sites_cart[ih])).length()
# alg1b or alg1a or propeller group
  elif (hp.htype in ['alg1b', 'alg1a', 'prop']):
    rb1 = matrix.col(sites_cart[hp.a2])
    n = hp.n
    phi = hp.phi + n*2*math.pi/3
    alpha = hp.alpha
    salpha = math.sin(alpha)
    calpha = math.cos(alpha)
    sphi = math.sin(phi)
    cphi = math.cos(phi)
    u1 = (r0 - r1).normalize()
    rb10 = rb1 - r1
    u2 = (rb10 - ((rb10).dot(u1)) * u1).normalize()
    u3 = u1.cross(u2)
    rh_calc = r0 + dh * (salpha*(cphi*u2 + sphi*u3) - calpha*u1)
    h_distance = (rh_calc - matrix.col(sites_cart[ih])).length()
  else:
    rh_calc = sites_cart[ih]
  return rh_calc

def diagnostics_parameterization(connectivity_obj,
  h_parameterization, sites_cart, threshold):
  number_h = connectivity_obj.number_h
  double_H = connectivity_obj.double_H
  h_connectivity = connectivity_obj.h_connectivity
  h_distances = {}
  unk_list = []
  unk_ideal_list = []
  long_distance_list = []
  for ih in sorted(h_parameterization.keys()):
    hp = h_parameterization[ih]
    rh = matrix.col(sites_cart[ih])
    if (hp.htype == 'unk'):
      h_distance = None
      unk_list.append(ih)
    elif (hp.htype == 'unk_ideal'):
      h_distance = None
      unk_ideal_list.append(ih)
    else:
      rh_calc = compute_H_position(
        ih         = ih,
        sites_cart = sites_cart,
        hp         = hp)
      if (rh_calc is not None):
        h_distance = (rh_calc - rh).length()
    if (h_distance is not None):
      h_distances[ih] = h_distance
      if (h_distance > threshold):
        long_distance_list.append(ih)
  set_temp = set(list(h_parameterization.keys()))
  slipped = [x for x in list(h_connectivity.keys()) if x not in set_temp]
  return group_args(
    number_h           = number_h,
    double_H           = double_H,
    h_distances        = h_distances,
    unk_list           = unk_list,
    unk_ideal_list     = unk_ideal_list,
    long_distance_list = long_distance_list,
    n_connect          = len(h_connectivity.keys()),
    slipped            = slipped,
    threshold          = threshold)

