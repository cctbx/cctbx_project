from __future__ import absolute_import, division, print_function
from scitbx import matrix
import math
from scitbx.math import dihedral_angle
from mmtbx_hydrogens_ext import *
from libtbx.utils import Sorry
from six.moves import zip

class manager(object):
  def __init__(self,
      h_connectivity,
      sites_cart,
      use_ideal_bonds_angles,
      site_labels,
      use_ideal_dihedral = False,
      ignore_h_with_dof  = False):
    self.h_connectivity = h_connectivity
    self.sites_cart = sites_cart
    self.use_ideal_bonds_angles = use_ideal_bonds_angles
    self.site_labels = site_labels
    self.use_ideal_dihedral = use_ideal_dihedral
    self.ignore_h_with_dof = ignore_h_with_dof
    self.determine_parameterization()

#-------------------------------------------------------------------------------

  def determine_parameterization(self):
    """
    For every H atom, determine the type of geometry
    """
    self.unk_list, self.unk_ideal_list = [], []
    self.h_parameterization = [None]*len(self.h_connectivity)
    for neighbors in self.h_connectivity:
      if (neighbors is None): continue
      ih = neighbors.ih
      if self.h_parameterization[ih] is not None:
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
        (number_h_neighbors == 0 or number_h_neighbors == 2) and
        not self.ignore_h_with_dof):
        self.process_1_neighbor(neighbors = neighbors)
      # planar Y-X-H2 groups such as in ARG head
      elif(number_non_h_neighbors == 1 and number_h_neighbors == 1):
        self.process_1_neighbor_type_arg(neighbors = neighbors)
      else:
        self.unk_list.append(ih)

#-------------------------------------------------------------------------------

  def process_1_neighbor(self, neighbors):
    ih = neighbors.ih
    # if used for hydrogenate, make sure that first we use the H with dihedral angle
    # However, this needs some tweaking for neutron H/D situations
    if (neighbors.number_h_neighbors == 2):
      i_h1, i_h2 = neighbors.h1['iseq'], neighbors.h2['iseq']
      if ('dihedral_ideal' in neighbors.b1):
        neighbors = self.h_connectivity[ih]
      elif ('dihedral_ideal' in self.h_connectivity[i_h1].b1):
        if self.h_parameterization[i_h1] is None:
          neighbors = self.h_connectivity[i_h1]
      elif ('dihedral_ideal' in self.h_connectivity[i_h2].b1):
        if self.h_parameterization[i_h2] is None:
          neighbors = self.h_connectivity[i_h2]
    ih = neighbors.ih
    #print(self.site_labels[ih])
    i_a0 = neighbors.a0['iseq']
    rh = matrix.col(self.sites_cart[ih])
    r0 = matrix.col(self.sites_cart[i_a0])
    if self.use_ideal_bonds_angles:
      disth = neighbors.a0['dist_ideal']
    else:
      disth = (r0 - rh).length()
    if (not neighbors.a1 or not neighbors.b1):
      self.unk_list.append(ih)
      return
    i_a1 = neighbors.a1['iseq']
    i_b1 = neighbors.b1['iseq']
    r1 = matrix.col(self.sites_cart[i_a1])
    rb1 = matrix.col(self.sites_cart[i_b1])
    self.check_if_atoms_superposed(rh, r0, ih, i_a0)
    self.check_if_atoms_superposed(r1, r0, i_a1, i_a0)
    uh0 = (rh - r0).normalize()
    u10 = (r1 - r0).normalize()
    dihedral = dihedral_angle(
      sites=[self.sites_cart[ih], self.sites_cart[i_a0],
      self.sites_cart[i_a1],self.sites_cart[i_b1]])
    if self.use_ideal_bonds_angles:
      alpha = math.radians(neighbors.a1['angle_ideal'])
      #allow for rotation even for idealize = True
      phi = dihedral
      if self.use_ideal_dihedral:
        #phi = math.radians(b1.dihedral_ideal)
        if 'dihedral_ideal' in neighbors.b1:
          phi = math.radians(neighbors.b1['dihedral_ideal'])
    else:
      alpha = (u10).angle(uh0)
      phi = dihedral
    #print(math.degrees(phi))
    u1 = (r0 - r1).normalize()
    rb10 = rb1 - r1
    # TODO check needed?
    u2 = (rb10 - ((rb10).dot(u1)) * u1).normalize()
    u3 = u1.cross(u2)
    if (neighbors.number_h_neighbors == 0):
      self.h_parameterization[ih] = riding_coefficients(
        htype  = 'alg1b',
        ih     = ih,
        a0     = i_a0,
        a1     = i_a1,
        a2     = i_b1,
        a3     = -1,
        a      = alpha,
        b      = phi,
        h      = 0,
        n      = 0,
        disth = disth)
    if (neighbors.number_h_neighbors == 2):
      i_h1, i_h2 = neighbors.h1['iseq'], neighbors.h2['iseq']
      i_h1, i_h2 = self.check_propeller_order(
        i_a0 = i_a0,
        i_a1 = i_a1,
        ih   = ih,
        i_h1 = i_h1,
        i_h2 = i_h2)
      for nprop, hprop in zip([0,1,2],[ih,i_h1,i_h2]):
        self.h_parameterization[hprop] = riding_coefficients(
          htype  = 'prop',
          ih     = hprop,
          a0     = i_a0,
          a1     = i_a1,
          a2     = i_b1,
          a3     = -1,
          a      = alpha,
          n      = nprop,
          b      = phi,
          h      = 0,
          disth = disth)
#a0.dihedral : dihedral angle between angle ideal and actual position

#-------------------------------------------------------------------------------

  def process_1_neighbor_type_arg(self, neighbors):
    """
      alg1a: X-H2 planar groups, such as in ARG, ASN, GLN
      requires that dihedral angle restraint exists for at least one H atom
    """
    ih = neighbors.ih
    i_h1 = neighbors.h1['iseq']
    i_a0 = neighbors.a0['iseq']
    rh = matrix.col(self.sites_cart[ih])
    r0 = matrix.col(self.sites_cart[i_a0])
    if self.use_ideal_bonds_angles:
      disth = neighbors.a0['dist_ideal']
    else:
      disth = (r0 - rh).length()
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
        self.unk_list.append(ih)
        return
    i_b1 = self.h_connectivity[ih_dihedral].b1['iseq']
    rb1 = matrix.col(self.sites_cart[i_b1])
    frame = None
    if self.is_pyramidal_nh2(neighbors, i_h1):
      frame = self.nh2_frame(ih_dihedral, ih_no_dihedral, i_b1)
    if frame is not None:
      # a rotor, like NH3+: refined individually against neutron data
      if self.ignore_h_with_dof:
        self.unk_list.append(ih)
        return
      self.process_pyramidal_nh2(ih_dihedral, ih_no_dihedral, *frame)
      return
    # check if angle is typical for propeller
    # catches case of missing propeller atom
    if (neighbors.h1['angle_ideal'] >107 and neighbors.h1['angle_ideal'] <111):
      self.unk_list.append(ih)
    else:
      dihedral = dihedral_angle(
        sites=[self.sites_cart[i_b1], self.sites_cart[i_a1],
        self.sites_cart[i_a0], self.sites_cart[ih_dihedral]])
      self.check_if_atoms_superposed(rh, r0, ih, i_a0)
      self.check_if_atoms_superposed(r1, r0, i_a1, i_a0)
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
      # TODO check needed?
      u2 = (rb10 - ((rb10).dot(u10)) * u10).normalize()
      u3 = u1.cross(u2)
      for ih_alg1a, phi_alg1a in zip(
        [ih_dihedral,ih_no_dihedral],[phi, phi+math.pi]):
        if self.h_parameterization[ih_alg1a] is None:
          self.h_parameterization[ih_alg1a] = riding_coefficients(
            htype  = 'alg1a',
            ih     = ih_alg1a,
            a0     = i_a0,
            a1     = i_a1,
            a2     = i_b1,
            a3     = -1,
            a      = alpha,
            b      = phi_alg1a,
            n      = 0,
            h      = 0,
            disth = disth)

  def is_pyramidal_nh2(self, neighbors, i_h1):
    """
    Amine NH2 (cif: N with 2 H and one sp3 C, N or O neighbour X; see
    connectivity._dictionary_pyramidal_nh2), neither H in a plane restraint
    """
    other = self.h_connectivity[i_h1]
    if other is None: return False
    if not (getattr(neighbors, 'pyramidal_nh2', False) and
            getattr(other, 'pyramidal_nh2', False)):
      return False
    for nb in (neighbors, other):
      if not nb.a1 or 'angle_ideal' not in nb.a1: return False
      if not nb.a0 or 'dist_ideal' not in nb.a0: return False
    if not neighbors.h1 or 'angle_ideal' not in neighbors.h1: return False
    return self.nh2_dphi(math.radians(neighbors.a1['angle_ideal']),
      math.radians(other.a1['angle_ideal']),
      math.radians(neighbors.h1['angle_ideal'])) is not None

  @staticmethod
  def nh2_dphi(a1, a2, hnh):
    """
    Torsion difference of the two H about X-N from the X-N-H angles a1, a2
    and H-N-H (radians); pi for planar angles (sum 360). None if no
    torsion gives hnh: H-N-H outside [|a1 - a2|, min(a1 + a2, 360 - a1 - a2)]
    (above the upper bound: angle sum above 360).
    """
    c = (math.cos(hnh) - math.cos(a1) * math.cos(a2)) / (
      math.sin(a1) * math.sin(a2))
    if abs(c) > 1 + 1.e-6: return None
    return math.acos(max(-1., min(1., c)))

  def nh2_frame(self, ih_dihedral, ih_no_dihedral, i_b1):
    """
    (frame atom, heavy-third torsions per H, placement phi without
    torsions) of an amine NH2.
    Frame: i_b1 if heavy (phi: connectivity's staggered default), else the
    third atom of a heavy-third torsion, else a heavy atom bonded to X
    (phi: anti, as the torsion minimum nearest 180). With no heavy atom on X
    (methylamine whose torsions are H-C-N-H only) the frame is the H i_b1:
    during placement it is still at its seed, so the NH2 is pyramidal at the
    dictionary angles with the first H anti to that H, and its orientation
    is arbitrary.
    """
    nb = {ih_dihedral: self.h_connectivity[ih_dihedral],
          ih_no_dihedral: self.h_connectivity[ih_no_dihedral]}
    # connectivity keeps torsions with a heavy third atom only (CCD fallback
    # torsions H-C-N-H have an H, still at its seed during placement)
    torsions = dict((i, list(getattr(nb[i], 'torsions', []))) for i in nb)
    if not getattr(nb[ih_dihedral], 'b1_is_hd', True):
      return i_b1, torsions, math.radians(
        nb[ih_dihedral].b1['dihedral_ideal'])
    thirds = ([t[0] for i in nb for t in torsions[i]] +
      list(getattr(nb[ih_dihedral], 'heavy_thirds', [])))
    if thirds: return thirds[0], torsions, math.pi
    return i_b1, torsions, math.pi

  def process_pyramidal_nh2(self, ih_dihedral, ih_no_dihedral, i_b1,
      torsions, phi_default):
    """
    Amine NH2 as two alg1a H, each with its own X-N-H and N-H, H-N-H from
    the dictionary (dphi). With use_ideal_dihedral (placement): one H at
    its torsion restraint minimum nearest 180 (anti), the other at
    phi +/- dphi on the side nearer its own torsion restraint, + without
    one; independent of where the H were seeded. Otherwise the first H
    keeps its current torsion (the NH2 is a rotor) and the other takes the
    side nearer its current site.
    """
    nb = {ih_dihedral: self.h_connectivity[ih_dihedral],
          ih_no_dihedral: self.h_connectivity[ih_no_dihedral]}
    i_a0 = nb[ih_dihedral].a0['iseq']
    i_a1 = nb[ih_dihedral].a1['iseq']
    r0 = matrix.col(self.sites_cart[i_a0])
    r1 = matrix.col(self.sites_cart[i_a1])
    rb1 = matrix.col(self.sites_cart[i_b1])
    for i in nb:
      self.check_if_atoms_superposed(matrix.col(self.sites_cart[i]), r0, i, i_a0)
    self.check_if_atoms_superposed(r1, r0, i_a1, i_a0)
    # the alg1a frame (hydrogens.h): phi of an atom = its torsion about
    # a1-a0 from i_b1
    u1 = (r0 - r1).normalize()
    rb10 = rb1 - r1
    u2 = (rb10 - (rb10.dot(u1)) * u1).normalize()
    u3 = u1.cross(u2)
    def site(alpha, phi, d):
      return r0 + d * (math.sin(alpha) * (math.cos(phi) * u2 +
        math.sin(phi) * u3) - math.cos(alpha) * u1)
    def frame_phi(r, origin):
      v = r - origin
      return math.atan2(v.dot(u3), v.dot(u2))
    def adiff(x, y):
      return abs((x - y + math.pi) % (2 * math.pi) - math.pi)
    def minima(i):
      """phi of the minima of the H's torsion restraints"""
      result = []
      for i_third, ideal, p in torsions[i]:
        off = frame_phi(matrix.col(self.sites_cart[i_third]), r1)
        n = p if p > 0 else 1
        result.extend(off + math.radians(ideal + k * 360. / n)
          for k in range(n))
      return result
    alpha, phi, d = {}, {}, {}
    for i in nb:
      if self.use_ideal_bonds_angles:
        alpha[i] = math.radians(nb[i].a1['angle_ideal'])
        d[i] = nb[i].a0['dist_ideal']
      else:
        v = matrix.col(self.sites_cart[i]) - r0
        alpha[i] = (r1 - r0).angle(v)
        d[i] = v.length()
    if not self.use_ideal_bonds_angles:
      for i in nb:
        phi[i] = frame_phi(matrix.col(self.sites_cart[i]), r0)
    else:
      first, second = ih_dihedral, ih_no_dihedral
      if self.use_ideal_dihedral:
        if not minima(first) and minima(second):
          first, second = second, first
        m1 = minima(first)
        if m1:
          phi1 = min(m1, key=lambda t: (adiff(t, math.pi), -t))
        else:
          phi1 = phi_default
      else:
        phi1 = frame_phi(matrix.col(self.sites_cart[first]), r0)
      dphi = self.nh2_dphi(alpha[first], alpha[second],
        math.radians(nb[first].h1['angle_ideal']))
      candidates = [phi1 + dphi, phi1 - dphi]
      if self.use_ideal_dihedral:
        m2 = minima(second)
        if m2:
          score = [min(adiff(c, t) for t in m2) for c in candidates]
        else:
          score = [0, 1]
      else:
        rh = matrix.col(self.sites_cart[second])
        score = [(site(alpha[second], c, d[second]) - rh).length()
          for c in candidates]
      phi[first] = phi1
      phi[second] = candidates[0 if score[0] <= score[1] else 1]
    for i in (ih_dihedral, ih_no_dihedral):
      if self.h_parameterization[i] is None:
        self.h_parameterization[i] = riding_coefficients(
          htype  = 'alg1a',
          ih     = i,
          a0     = i_a0,
          a1     = i_a1,
          a2     = i_b1,
          a3     = -1,
          a      = alpha[i],
          b      = phi[i],
          n      = 0,
          h      = 0,
          disth  = d[i])

  # alg2a, 2tetra, 2neigbs
  def process_2_neighbors(self, neighbors):
    ih = neighbors.ih
    i_a0 = neighbors.a0['iseq']
    rh = matrix.col(self.sites_cart[ih])
    r0 = matrix.col(self.sites_cart[i_a0])
    if self.use_ideal_bonds_angles:
      disth = neighbors.a0['dist_ideal']
    else:
      disth = (r0 - rh).length()
    # if H is second neighbor, get its index
    if (neighbors.number_h_neighbors == 1):
      i_h1 = neighbors.h1['iseq']
    else:
      i_h1 = None
    sumang, a, b, h, root = self.get_coefficients(ih = ih)
    # alg2a
    if (sumang > (2*math.pi + 0.05) and root < 0):
      self.unk_ideal_list.append(ih)
      return
    elif (sumang < (2*math.pi + 0.05) and (sumang > 2*math.pi - 0.05)):
      htype = 'flat_2neigbs'
    else:
      if (neighbors.number_h_neighbors == 1):
      # 2 tetragonal geometry
        htype = '2tetra'
        self.h_parameterization[i_h1] = riding_coefficients(
          ih     = i_h1,
          a0     = neighbors.a0['iseq'],
          a1     = neighbors.a1['iseq'],
          a2     = neighbors.a2['iseq'],
          a3     = -1,
          a      = a,
          b      = b,
          h      = -h,
          n      = 0,
          disth  = disth,
          htype  = '2tetra')
      else:
        # 2neigbs
        htype = '2neigbs'
        # A single H whose parent is restrained planar (e.g. an amide N-H such
        # as the lactam N of PCA) must lie in the plane of its two neighbours,
        # even when the ideal angles around the parent do not sum to 360. Place
        # it flat instead of pyramidal.
        if getattr(neighbors, 'is_in_plane', False):
          htype = 'flat_2neigbs'
    if (h is None):
      h = 0
    self.h_parameterization[ih] = riding_coefficients(
      htype = htype,
      ih    = ih,
      a0    = neighbors.a0['iseq'],
      a1    = neighbors.a1['iseq'],
      a2    = neighbors.a2['iseq'],
      a3    = -1,
      a     = a,
      b     = b,
      h     = h,
      n     = 0,
      disth = disth)

#-------------------------------------------------------------------------------

  def process_3_neighbors(self, neighbors):
    ih = neighbors.ih
    i_a0 = neighbors.a0['iseq']
    rh = matrix.col(self.sites_cart[ih])
    r0 = matrix.col(self.sites_cart[i_a0])
    if self.use_ideal_bonds_angles:
      disth = neighbors.a0['dist_ideal']
    else:
      disth = (r0 - rh).length()
    a, b, h = self.get_coefficients_alg3(ih = ih)
    self.h_parameterization[ih] = riding_coefficients(
      ih     = ih,
      a0     = neighbors.a0['iseq'],
      a1     = neighbors.a1['iseq'],
      a2     = neighbors.a2['iseq'],
      a3     = neighbors.a3['iseq'],
      a      = a,
      b      = b,
      h      = h,
      n      = 0,
      disth = disth,
      htype  = '3neigbs')

#-------------------------------------------------------------------------------

  def get_coefficients(self, ih):
    """
    This function determines parameters for three cases:
    1. planar geometry
    2. two tetragonal CH2 geometry
    3. H out of plane of its 3 neighbors (should be rare and not in AA)
    """
    neighbors = self.h_connectivity[ih]
    if (neighbors.number_h_neighbors == 1):
      i_h1 = neighbors.h1['iseq']
    else:
      i_h1 = None
    i_a0 = neighbors.a0['iseq']
    i_a1 = neighbors.a1['iseq']
    i_a2 = neighbors.a2['iseq']
    rh = matrix.col(self.sites_cart[ih])
    r0 = matrix.col(self.sites_cart[i_a0])
    r1 = matrix.col(self.sites_cart[i_a1])
    r2 = matrix.col(self.sites_cart[i_a2])
    self.check_if_atoms_superposed(rh, r0, ih, i_a0)
    self.check_if_atoms_superposed(r1, r0, i_a1, i_a0)
    self.check_if_atoms_superposed(r2, r0, i_a2, i_a0)
    uh0 = (rh - r0).normalize()
    u10 = (r1 - r0).normalize()
    u20 = (r2 - r0).normalize()
    if self.use_ideal_bonds_angles:
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
      self.broadcast_problem(ih = ih, i_a0 = i_a0)
      #raise RuntimeError(
      #  "Denominator zero: (1-c0*c0) in get_h_parameterization.")
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
        rh2 = matrix.col(self.sites_cart[neighbors.h1['iseq']])
        #print(i_h1, neighbors.h1['iseq'], i_a0)
        self.check_if_atoms_superposed(rh2, r0, neighbors.h1['iseq'], i_a0)
        uh02 = (rh2 - r0).normalize()
        if self.use_ideal_bonds_angles:
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
          self.broadcast_problem(ih = ih, i_a0 = i_a0,
            msg='(Square root of zero in get_coefficients: H out of plane)')
        denom = math.sin(alpha0)
        if(denom==0):
          self.broadcast_problem(ih = ih, i_a0 = i_a0,
            msg='(Denominator zero in get_coefficients: H out of plane)')
        cz = (math.sqrt(1-c1*c1-c2*c2-c0*c0+2*c0*c1*c2))/math.sin(alpha0)
        h = cz
        #test if vector v points to same 'side' as uh0
        if((u10.cross(u20)).dot(uh0) < 0):
          h = -h
    return sumang, a, b, h, root

#-------------------------------------------------------------------------------

  def get_coefficients_alg3(self, ih):
    """
    Obtain coefficients for tetragonal H (such as HA) using Cramer's rule
    """
    neighbors = self.h_connectivity[ih]
    i_a0 = neighbors.a0['iseq']
    i_a1 = neighbors.a1['iseq']
    i_a2 = neighbors.a2['iseq']
    i_a3 = neighbors.a3['iseq']
    rh = matrix.col(self.sites_cart[ih])
    r0 = matrix.col(self.sites_cart[i_a0])
    r1 = matrix.col(self.sites_cart[i_a1])
    r2 = matrix.col(self.sites_cart[i_a2])
    r3 = matrix.col(self.sites_cart[i_a3])
    self.check_if_atoms_superposed(rh, r0, ih, i_a0)
    self.check_if_atoms_superposed(r1, r0, i_a1, i_a0)
    self.check_if_atoms_superposed(r2, r0, i_a2, i_a0)
    self.check_if_atoms_superposed(r3, r0, i_a3, i_a0)
    uh0 = (rh - r0).normalize()
    u10 = (r1 - r0).normalize()
    u20 = (r2 - r0).normalize()
    u30 = (r3 - r0).normalize()
    if self.use_ideal_bonds_angles:
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
      self.broadcast_problem(ih = ih, i_a0 = i_a0)
      #raise RuntimeError(
      #  "Denominator zero: matrix_d in get_h_parameterization.")
    a = matrix_x.determinant()/matrix_d.determinant()
    b = matrix_y.determinant()/matrix_d.determinant()
    c = matrix_z.determinant()/matrix_d.determinant()
    return a, b, c

#-------------------------------------------------------------------------------

  def check_propeller_order(self, i_a0, i_a1, ih, i_h1, i_h2):
    rh = matrix.col(self.sites_cart[ih])
    rh_2 = matrix.col(self.sites_cart[i_h2])
    r0 = matrix.col(self.sites_cart[i_a0])
    r1 = matrix.col(self.sites_cart[i_a1])
    if (((rh-r0).cross(rh_2-r0)).dot(r1-r0) >= 0):
      return i_h1, i_h2
    else:
      return i_h2, i_h1

#-------------------------------------------------------------------------------

  def check_if_atoms_superposed(self,r1, r2, i1, i2):
    if abs((r1-r2).length()) < 0.001:
      sorry_str = '''Atoms %s and %s are superposed or very close.
Fix your model before proceeding.''' % (self.site_labels[i1], self.site_labels[i2])
      raise Sorry(sorry_str)

#-------------------------------------------------------------------------------

  def broadcast_problem(self, ih, i_a0, msg=None,):
    default_msg = '''Please double check atom %s, bound to %s
as well as nearest neighbors. The input geometry is most likely wrong.
Solution: Fix the geometry or delete the H atom.'''
    #if msg is not None:
    #  default_msg = default_msg + '\n' + msg
    raise Sorry(default_msg % (self.site_labels[ih], self.site_labels[i_a0]))
