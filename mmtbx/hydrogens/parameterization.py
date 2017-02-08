from __future__ import division
from scitbx import matrix
#from scitbx.array_family import flex
#from libtbx.utils import Sorry
from libtbx import group_args
from stdlib import math
from scitbx.math import dihedral_angle

class parameterization_info(object):
  def __init__(
    self,
    htype   = None,  # type of hydrogen geometry
    a0      = None,  # parent atom index
    a1      = None,  # 1-3 neighbor index
    a2      = None,  # 1-3 or 1-4 neighbor index
    a3      = None,  # 1-3 or 1-4 neighbor index
    a       = None,  # coefficient for reconstruction
    b       = None,  # coefficient for reconstruction
    h       = None,  # coefficient for reconstruction
    n       = None,  # parameter for propeller and H2 groups: integer
    dist_h  = None): # measured or ideal distance
    self.htype  = htype
    self.a0     = a0
    self.a1     = a1
    self.a2     = a2
    self.a3     = a3
    self.a      = a
    self.b      = b
    self.h      = h
    self.n      = n
    self.dist_h = dist_h

class manager(object):
  def __init__(self,
      h_connectivity,
      sites_cart,
      use_ideal_bonds_angles):
    self.h_connectivity = h_connectivity
    self.sites_cart = sites_cart
    self.use_ideal_bonds_angles = use_ideal_bonds_angles

  def test_print(self, ih, neighbors):
    print 'now at atom %s' % (ih)

# for every H atom, determine the type of bond
  def determine_parameterization(self):
    self.h_parameterization = {}
    for neighbors in self.h_connectivity:
      if (neighbors is None): continue
      ih = neighbors.ih
      #if ih in h_parameterization.keys():
      if ih in self.h_parameterization:
        continue
      number_h_neighbors = neighbors.number_h_neighbors
      number_non_h_neighbors = neighbors.number_non_h_neighbors
      # alg2a, 2tetra, 2neigbs
      if (number_non_h_neighbors == 2):
        self.process_2_neighbors(neighbors = neighbors)
      # tetragonal geometry: 3neigbs
      elif (number_non_h_neighbors == 3 and number_h_neighbors == 0):
        self.process_3_neighbors(neighbors = neighbors)
      # Free rotation and propeller groups
      elif(number_non_h_neighbors == 1 and
        (number_h_neighbors == 0 or number_h_neighbors ==2)):
        self.process_1_neighbor(neighbors = neighbors)
      # planar Y-X-H2 groups such as in ARG head
      elif(number_non_h_neighbors == 1 and number_h_neighbors == 1):
        self.process_1_neighbor_type_arg(neighbors = neighbors)
      else:
        self.h_parameterization[ih] = parameterization_info(
          htype = 'unk',
          a0    = neighbors.a0['iseq'])
    return self.h_parameterization

  def process_1_neighbor(self, neighbors):
    ih = neighbors.ih
    i_a0 = neighbors.a0['iseq']
    rh = matrix.col(self.sites_cart[ih])
    r0 = matrix.col(self.sites_cart[i_a0])
    if self.use_ideal_bonds_angles:
      dist_h = neighbors.a0['dist_ideal']
    else:
      dist_h = (r0 - rh).length()
    i_a1 = neighbors.a1['iseq']
    i_b1 = neighbors.b1['iseq']
    r1 = matrix.col(self.sites_cart[i_a1])
    rb1 = matrix.col(self.sites_cart[i_b1])
    uh0 = (rh - r0).normalize()
    u10 = (r1 - r0).normalize()
    dihedral = dihedral_angle(
      sites=[self.sites_cart[ih], self.sites_cart[i_a0],
      self.sites_cart[i_a1],self.sites_cart[i_b1]])
    if self.use_ideal_bonds_angles:
      alpha = math.radians(neighbors.a1['angle_ideal'])
      #allow for rotation even for idealize = True
      #phi = math.radians(b1.dihedral_ideal)
      phi = dihedral
    else:
      alpha = (u10).angle(uh0)
      phi = dihedral
    u1 = (r0 - r1).normalize()
    rb10 = rb1 - r1
    u2 = (rb10 - ((rb10).dot(u1)) * u1).normalize()
    u3 = u1.cross(u2)
    self.h_parameterization[ih] = parameterization_info(
      htype  = 'alg1b',
      a0     = i_a0,
      a1     = i_a1,
      a2     = i_b1,
      a3     = 0,
      a      = alpha,
      b      = phi,
      h      = 0,
      n      = 0,
      dist_h = dist_h)
    if (neighbors.number_h_neighbors == 2):
      self.h_parameterization[ih].htype = 'prop'
      i_h1, i_h2 = neighbors.h1['iseq'], neighbors.h2['iseq']
      self.h_parameterization[i_h1] = parameterization_info(
        htype  = 'prop',
        a0     = i_a0,
        a1     = i_a1,
        a2     = i_b1,
        a3     = 0,
        a      = alpha,
        n      = 1,
        b      = phi,
        h      = 0,
        dist_h = dist_h)
      # check if order is reversed
      # this can maybe be done earlier, in connectivity?
      i_h1_coord = compute_H_position(
        sites_cart = self.sites_cart,
        ih         = i_h1,
        hp         = self.h_parameterization[i_h1])
      self.h_parameterization[i_h2] = parameterization_info(
        htype  = 'prop',
        a0     = i_a0,
        a1     = i_a1,
        a2     = i_b1,
        a3     = 0,
        a      = alpha,
        n      = 2,
        b      = phi,
        h      = 0,
        dist_h = dist_h)
      if ((i_h1_coord - matrix.col(self.sites_cart[i_h2])).length() <
        (i_h1_coord - matrix.col(self.sites_cart[i_h1])).length()):
        self.h_parameterization[i_h1].n = 2
        self.h_parameterization[i_h2].n = 1
#info: a0.dihedral = dihedral angle between angle ideal and actual position

#    # alg1a: X-H2 planar groups, such as in ARG, ASN, GLN
#    # requires that dihedral angle restraint exists for at least one H atom
  def process_1_neighbor_type_arg(self, neighbors):
    ih = neighbors.ih
    i_h1 = neighbors.h1['iseq']
    i_a0 = neighbors.a0['iseq']
    rh = matrix.col(self.sites_cart[ih])
    r0 = matrix.col(self.sites_cart[i_a0])
    if self.use_ideal_bonds_angles:
      dist_h = neighbors.a0['dist_ideal']
    else:
      dist_h = (r0 - rh).length()
    i_a1 = neighbors.a1['iseq']
    r1 = matrix.col(self.sites_cart[i_a1])
    if ('dihedral_ideal' in neighbors.b1):
      ih_dihedral = ih
      ih_no_dihedral = i_h1
    else:
      if ('dihedral_ideal' in self.h_connectivity[i_h1].b1):
        ih_dihedral = i_h1
        ih_no_dihedral = ih
      else:
        self.h_parameterization[ih] = parameterization_info(
          htype  = 'unk',
          a0     = i_a0)
        self.h_parameterization[i_h1] = parameterization_info(
          htype  = 'unk',
          a0     = i_a0)
        return
    i_b1 = self.h_connectivity[ih_dihedral].b1['iseq']
    rb1 = matrix.col(self.sites_cart[i_b1])
    # check if angle is typical for propeller
    # catches case of missing propeller atom
    if (neighbors.h1['angle_ideal'] >107 and neighbors.h1['angle_ideal'] <111):
      self.h_parameterization[ih] = parameterization_info(
        htype  = 'unk',
        a0     = a0.iseq)
      self.h_parameterization[i_h1] = parameterization_info(
        htype  = 'unk',
        a0     = a0.iseq)
    else:
      dihedral = dihedral_angle(
        sites=[self.sites_cart[i_b1], self.sites_cart[i_a1],
        self.sites_cart[i_a0], self.sites_cart[ih_dihedral]])
      uh0 = (rh - r0).normalize()
      u10 = (r1 - r0).normalize()
      if self.use_ideal_bonds_angles:
        alpha = math.radians(neighbors.a1['angle_ideal'])
        phi = math.radians(self.h_connectivity[ih_dihedral].b1['dihedral_ideal'])
      else:
        alpha = (u10).angle(uh0)
        phi = dihedral
      u1 = (r0 - r1).normalize()
      rb10 = rb1 - r1
      u2 = (rb10 - ((rb10).dot(u10)) * u10).normalize()
      u3 = u1.cross(u2)
      if ih_dihedral not in self.h_parameterization:
        self.h_parameterization[ih_dihedral] = parameterization_info(
          htype  = 'alg1a',
          a0     = i_a0,
          a1     = i_a1,
          a2     = i_b1,
          a3     = 0,
          a      = alpha,
          n      = 0,
          b      = phi,
          h      = 0,
          dist_h = dist_h)
      if ih_no_dihedral not in self.h_parameterization:
        self.h_parameterization[ih_no_dihedral] = parameterization_info(
          htype  = 'alg1a',
          a0     = i_a0,
          a1     = i_a1,
          a2     = i_b1,
          a3     = 0,
          a      = alpha,
          b      = phi+math.pi,
          n      = 0,
          h      = 0,
          dist_h = dist_h)

  # alg2a, 2tetra, 2neigbs
  def process_2_neighbors(self, neighbors):
    ih = neighbors.ih
    i_a0 = neighbors.a0['iseq']
    rh = matrix.col(self.sites_cart[ih])
    r0 = matrix.col(self.sites_cart[i_a0])
    if self.use_ideal_bonds_angles:
      dist_h = neighbors.a0['dist_ideal']
    else:
      dist_h = (r0 - rh).length()
    # if H is second neighbor, get its index
    if (neighbors.number_h_neighbors == 1):
      i_h1 = neighbors.h1['iseq']
    else:
      i_h1 = None
    sumang, a, b, h, root = self.get_coefficients(
      ih                     = ih,
      use_ideal_bonds_angles = self.use_ideal_bonds_angles,
      sites_cart             = self.sites_cart)
    self.h_parameterization[ih] = parameterization_info(
      a0     = neighbors.a0['iseq'],
      a1     = neighbors.a1['iseq'],
      a2     = neighbors.a2['iseq'],
      a3     = 0,
      a      = a,
      b      = b,
      h      = 0,
      n      = 0,
      dist_h = dist_h)
    # alg2a
    if (sumang > (2*math.pi + 0.05) and root < 0):
      self.h_parameterization[ih].htype = 'unk_ideal'
    elif (sumang < (2*math.pi + 0.05) and (sumang > 2*math.pi - 0.05)):
      self.h_parameterization[ih].htype = 'flat_2neigbs'
    else:
      if (neighbors.number_h_neighbors == 1):
      # 2 tetragonal geometry
        self.h_parameterization[ih].htype = '2tetra'
        self.h_parameterization[ih].h = h
        self.h_parameterization[ih].n = 0
        self.h_parameterization[i_h1] = parameterization_info(
          a0     = neighbors.a0['iseq'],
          a1     = neighbors.a1['iseq'],
          a2     = neighbors.a2['iseq'],
          a3     = 0,
          a      = a,
          b      = b,
          h      = -h,
          n      = 0,
          dist_h = dist_h,
          htype  = '2tetra')
      else:
        # 2neigbs
        self.h_parameterization[ih].h = h
        self.h_parameterization[ih].htype = '2neigbs'
        self.h_parameterization[ih].n = 0

  def process_3_neighbors(self, neighbors):
    ih = neighbors.ih
    i_a0 = neighbors.a0['iseq']
    rh = matrix.col(self.sites_cart[ih])
    r0 = matrix.col(self.sites_cart[i_a0])
    if self.use_ideal_bonds_angles:
      dist_h = neighbors.a0['dist_ideal']
    else:
      dist_h = (r0 - rh).length()
    a, b, h = self.get_coefficients_alg3(
      neighbors              = neighbors,
      use_ideal_bonds_angles = self.use_ideal_bonds_angles,
      sites_cart             = self.sites_cart)
    self.h_parameterization[ih] = parameterization_info(
      a0     = neighbors.a0['iseq'],
      a1     = neighbors.a1['iseq'],
      a2     = neighbors.a2['iseq'],
      a3     = neighbors.a3['iseq'],
      a      = a,
      b      = b,
      h      = h,
      n      = 0,
      dist_h = dist_h,
      htype  = '3neigbs')

# this function determines parameters for three cases:
# 1. planar geometry
# 2. two tetragonal CH2 geometry
# 3. H out of plane of its 3 neighbors (should be rare and not in AA)
  def get_coefficients(self, ih, use_ideal_bonds_angles, sites_cart):
    neighbors = self.h_connectivity[ih]
    ih = neighbors.ih
    if (neighbors.number_h_neighbors == 1):
      i_h1 = neighbors.h1['iseq']
    else:
      i_h1 = None
    i_a0 = neighbors.a0['iseq']
    i_a1 = neighbors.a1['iseq']
    i_a2 = neighbors.a2['iseq']
    rh = matrix.col(sites_cart[ih])
    r0 = matrix.col(sites_cart[i_a0])
    r1 = matrix.col(sites_cart[i_a1])
    r2 = matrix.col(sites_cart[i_a2])
    uh0 = (rh - r0).normalize()
    u10 = (r1 - r0).normalize()
    u20 = (r2 - r0).normalize()
    if use_ideal_bonds_angles:
      alpha0 = math.radians(neighbors.a0['angle_a1a0a2'])
      alpha1 = math.radians(neighbors.a1['angle_ideal'])
      alpha2 = math.radians(neighbors.a2['angle_ideal'])
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
  #  # check if H, A0, A1, A2 are in a plane
    if (sumang < (2*math.pi + 0.05) and (sumang > 2*math.pi - 0.05)):
      h = None
    elif (sumang > (2*math.pi + 0.05) and 1-c1*c1-c2*c2-c0*c0+2*c0*c1*c2 < 0):
      root = 1-c1*c1-c2*c2-c0*c0+2*c0*c1*c2
      h = None
      return sumang, a, b, h, root
    else:
      # two tetragonal geometry: e.g. CH2 group
      if (i_h1 is not None):
        rh2 = matrix.col(sites_cart[neighbors.h1['iseq']])
        uh02 = (rh2 - r0).normalize()
        if use_ideal_bonds_angles:
          h = math.radians(neighbors.h1['angle_ideal']) * 0.5
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

#
# obtain coefficients for tetragonal H (such as HA) using Cramer's rule
  def get_coefficients_alg3(self,neighbors, use_ideal_bonds_angles, sites_cart):
    ih = neighbors.ih
    i_a0 = neighbors.a0['iseq']
    i_a1 = neighbors.a1['iseq']
    i_a2 = neighbors.a2['iseq']
    i_a3 = neighbors.a3['iseq']
    rh = matrix.col(sites_cart[ih])
    r0 = matrix.col(sites_cart[i_a0])
    r1 = matrix.col(sites_cart[i_a1])
    r2 = matrix.col(sites_cart[i_a2])
    r3 = matrix.col(sites_cart[i_a3])
    uh0 = (rh - r0).normalize()
    u10 = (r1 - r0).normalize()
    u20 = (r2 - r0).normalize()
    u30 = (r3 - r0).normalize()
    if use_ideal_bonds_angles:
      alpha0 = math.radians(neighbors.a1['angle_ideal'])
      alpha1 = math.radians(neighbors.a2['angle_ideal'])
      alpha2 = math.radians(neighbors.a3['angle_ideal'])
      c1, c2, c3 = math.cos(alpha0), math.cos(alpha1), math.cos(alpha2)
      omega0 = math.radians(neighbors.a0['angle_a1a0a2'])
      omega1 = math.radians(neighbors.a0['angle_a2a0a3'])
      omega2 = math.radians(neighbors.a0['angle_a3a0a1'])
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
    a, b, delta = hp.a, hp.b, hp.h
    r2 = matrix.col(sites_cart[hp.a2])
    u10, u20 = (r1 - r0).normalize(), (r2 - r0).normalize()
    v0 = (u10.cross(u20)).normalize()
    d0 = (a * u10 + b * u20).normalize()
    rh_calc = r0 + dh * (math.cos(delta) * d0 + math.sin(delta) * v0)
    h_distance = (rh_calc - matrix.col(sites_cart[ih])).length()
  # tetrahedral coordination
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
    phi = hp.b + n*2*math.pi/3
    alpha = hp.a
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

def count_h(h_connectivity):
  number_h = 0
  for item in h_connectivity:
    if item: number_h = number_h + 1
  return number_h

def diagnostics(sites_cart, threshold, h_parameterization, h_connectivity):
  number_h = count_h(h_connectivity = h_connectivity)
  #double_H = self.h_connectivity.double_H
  h_distances = {}
  unk_list = []
  unk_ideal_list = []
  long_distance_list = []
  for ih in h_parameterization:
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
      else:
        h_distance = None
    if (h_distance is not None):
      h_distances[ih] = h_distance
      if (h_distance > threshold):
        long_distance_list.append(ih)
  set_temp = set(list(h_parameterization.keys()))
  slipped = [x for x in h_connectivity if x not in set_temp]
  return group_args(
    number_h           = number_h,
    #double_H           = double_H,
    h_distances        = h_distances,
    unk_list           = unk_list,
    unk_ideal_list     = unk_ideal_list,
    long_distance_list = long_distance_list,
    n_connect          = len(h_connectivity),
    slipped            = slipped,
    threshold          = threshold)

