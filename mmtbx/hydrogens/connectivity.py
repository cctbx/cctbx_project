from __future__ import absolute_import, division, print_function
import math
from cctbx import geometry_restraints
from libtbx.utils import Sorry
from libtbx import group_args
from scitbx.array_family import flex
from scitbx.math import dihedral_angle
from six.moves import zip,range

class neighbors(object):
  def __init__(self,
      ih = None,
      a0 = {},
      a1 = {},
      a2 = {},
      a3 = {},
      h1 = {},
      h2 = {},
      b1 = {},
      number_h_neighbors = None,
      number_non_h_neighbors = None):
    self.ih = ih
    self.a0 = a0
    self.a1 = a1
    self.a2 = a2
    self.a3 = a3
    self.h1 = h1
    self.h2 = h2
    self.b1 = b1
    self.number_h_neighbors = number_h_neighbors
    self.number_non_h_neighbors = number_non_h_neighbors

class determine_connectivity(object):
  """ Obtain information about the necessary number of neighbors to reconstruct
  the position of every H atom ("connectivity") and to determine the geometry
  of the H atom. Store also information about ideal angles involving H atom,
  and non-H atoms as well as dihedral angle
  :returns: an array of "neighbors" objects (defined above)
  :rtype: list[]"""
  def __init__(self,
      pdb_hierarchy,
      geometry_restraints):
    geometry = geometry_restraints
    self.atoms = pdb_hierarchy.atoms()
    self.sites_cart = self.atoms.extract_xyz()
    bond_proxies_simple, asu = \
      geometry.get_all_bond_proxies(sites_cart = self.sites_cart)
    angle_proxies = geometry.get_all_angle_proxies()
    dihedral_proxies = geometry.dihedral_proxies # this should be function in GRM, like previous
    planarity_proxies = geometry.planarity_proxies
    fsc0 = geometry.shell_sym_tables[0].full_simple_connectivity()
    self.n_atoms = pdb_hierarchy.atoms_size()
    self.hd_sel = self.hd_selection()
    self.h_connectivity = [None]*self.n_atoms
    # 1. Find parent atoms and ideal A0-H bond distances
    self.find_first_neighbors(
      bond_proxies_simple = bond_proxies_simple,
      fsc0                = fsc0)
    # Check if number H in connectivity and number H in model is the same
    self.count_H()
    # 2. find preliminary list of second neighbors
    self.find_second_neighbors_raw(angle_proxies = angle_proxies)

    # 3. Get plane proxies --> useful for NH2 groups without dihedrals
    self.process_plane_proxies(planarity_proxies = planarity_proxies)

    # 4. process preliminary list to eliminate atoms in double conformation
    self.process_second_neighbors()

    # 5. Find third neighbors via dihedral proxies
    self.find_third_neighbors(dihedral_proxies = dihedral_proxies)

    # 6. Find angles involving a0 and covalently bound non-H atoms and
    #    find preliminary list of third neighbors in cases where no dihedral
    #    proxy is present
    self.determine_a0_angles_and_third_neighbors_without_dihedral(
      angle_proxies = angle_proxies)

    # 7. Assign the angles found previously and process preliminary list of
    #    third neighbors
    self.process_a0_angles_and_third_neighbors_without_dihedral()

    # Add slipped H atoms to h_connectivity
    self.add_slipped()


  def find_first_neighbors(self, bond_proxies_simple, fsc0):
    """ Find first neighbors by looping through bond proxies
    Fills in dictionary of 'a0' in object 'neighbors'.
    Keys: i_seq, dist_ideal"""
    self.double_H = {}
    self.parents = set()
    for bproxy in bond_proxies_simple:
      i_seq, j_seq = bproxy.i_seqs
      is_i_hd = self.hd_sel[i_seq]
      is_j_hd = self.hd_sel[j_seq]
      if(not is_i_hd and not is_j_hd): continue
      elif(is_i_hd and is_j_hd):       assert 0
      else:
        if  (is_i_hd): ih, i_parent = i_seq, j_seq
        elif(is_j_hd): ih, i_parent = j_seq, i_seq
        # else case should not happen...
        #else:
        #  raise Sorry("Something went wrong in bond proxies")
      # if neighbor exists, use only first one found
      if self.h_connectivity[ih] is not None: continue
      # find H atoms bound to two parent atoms, store info in list 'double_H'
      if (fsc0[ih].size() > 1):
        self.double_H[ih] = list(fsc0[ih])
      #else:
      self.h_connectivity[ih] = neighbors(
        ih = ih,
        a0 = {'iseq':i_parent, 'dist_ideal': bproxy.distance_ideal,
          'angles':[]})
      self.parents.add(i_parent) #parent atoms stored in a set

  def find_second_neighbors_raw(self, angle_proxies):
    """Get a an array listing all second neighbors for every H atom.
    :returns: a list of lists: [[iseq1, iseq2, ...], ...]
    :rtype: [[],[].[]]
    """
    self.second_neighbors_raw = [[] for i in range(self.n_atoms)]
    self.angle_dict = {}
    for ap in angle_proxies:
      for i_test in ap.i_seqs:
        if (self.h_connectivity[i_test] is None and not self.hd_sel[i_test]):
          continue
        ih = i_test
        i_parent = ap.i_seqs[1]
        # if one H bound to two parents (case double_H)
        if (i_parent not in self.parents): continue
        bonded = [ih, i_parent]
        i_second = [x for x in ap.i_seqs if x not in bonded][0]
        if (self.h_connectivity[ih] is None): continue
        #assert(self.h_connectivity[ih].a0['iseq'] == i_parent)
        if (self.h_connectivity[ih].a0['iseq'] != i_parent):
          if ih in self.double_H:
          #if self.double_H[ih]:
            continue
          raise Sorry  (
            "It looks like angle and bond restraints are conflicting.\n\
             Bond proxies:  H atom %s is bound to %s. \n\
             Angle proxies: H atom %s is bound to %s" % \
    (self.atoms[ih].id_str(),
     self.atoms[self.h_connectivity[ih].a0['iseq']].id_str(),
     self.atoms[ih].id_str(), self.atoms[i_parent].id_str()  )   )
        self.second_neighbors_raw[ih].append(i_second)
        self.angle_dict[(ih, i_parent, i_second)] = ap.angle_ideal


  def process_second_neighbors(self):
    """Once candidates for second neighbors are determined, they are further
    processed, mainly to avoid alternative conformations of the same atom"""
    self.a0a1_dict = {}
    self.a1_atoms = set()
    for neighbor_obj in self.h_connectivity:
      if (neighbor_obj is None): continue
      ih = neighbor_obj.ih
      i_parent = neighbor_obj.a0['iseq']
      second_neighbors_reduced = []
      alt_conf_neighbors = []
      for i_second in self.second_neighbors_raw[ih]:
        altloc_i_second = self.atoms[i_second].parent().altloc
        if (altloc_i_second == ''):
          second_neighbors_reduced.append(i_second)
        else:
          alt_conf_neighbors.append(i_second)
      second_neighbors_reduced.extend(self.process_alternate_neighbors(
        alt_conf_neighbors = alt_conf_neighbors))
      second_neighbors_H = self.determine_second_neighbors_H(
        second_neighbors_reduced = second_neighbors_reduced)
      second_neighbors_non_H = list(
        set(second_neighbors_reduced) - set(second_neighbors_H))
      self.assign_second_neighbors(
        ih               = ih,
        i_parent         = i_parent,
        neighbors_list   = second_neighbors_non_H,
        neighbors_list_H = second_neighbors_H)
      if (neighbor_obj.number_non_h_neighbors == 1):
        i_a1 = self.h_connectivity[ih].a1['iseq']
        self.a1_atoms.add(i_a1)
        if i_parent in self.a0a1_dict.keys():
          self.a0a1_dict[i_parent].append(i_a1)
        else:
          self.a0a1_dict[i_parent]=[i_a1]
          #self.a0a1_dict[i_parent] = i_a1


  def determine_a0_angles_and_third_neighbors_without_dihedral(self, angle_proxies):
    """Loop through angle proxies to find angles involving a0 and second
    neighbors. Find raw list of third neighbors, which don't have dihedral
    proxies."""
    self.parent_angles = [{} for i in range(self.n_atoms)]
    self.third_neighbors_raw = [[] for i in range(self.n_atoms)]
    for ap in angle_proxies:
      ix, iy, iz = ap.i_seqs
      is_hd_ix = self.hd_sel[ix]
      is_hd_iz = self.hd_sel[iz]
      # get all X1-A0-X2 angles if A0 is parent atom
      if (iy in self.parents and not is_hd_ix and not is_hd_iz):
          self.parent_angles[iy][(ix, iz)] = ap.angle_ideal
      # for third neighbors, a1 atom is central
      if (is_hd_ix or is_hd_iz): continue
      i_third = None
      if (iy in self.a1_atoms):
        if (ix in self.parents and ix in self.a0a1_dict):
          if (iy in self.a0a1_dict[ix] and not is_hd_iz):
            i_parent = ix
            i_third = iz
          elif (i_third == None):
            raise  Sorry(
  "It looks like angle restraints involving an H atom are missing.\n\
  Check H atoms bound to %s and with second neighbor %s" % \
    (self.atoms[ix].id_str(), self.atoms[iy].id_str()))
        elif (iz in self.parents and iz in self.a0a1_dict):
          if (iy in self.a0a1_dict[iz] and not is_hd_ix):
            i_parent = iz
            i_third = ix
          elif (i_third == None):
            raise  Sorry(
  "It looks like angle restraints involving an H atom are missing.\n\
  Check H atoms bound to %s and with second neighbor %s" % \
    (self.atoms[iz].id_str(), self.atoms[iy].id_str()))
        else:
          continue
        if (i_third not in self.third_neighbors_raw[i_parent]):
          self.third_neighbors_raw[i_parent].append(i_third)
          if (i_third in self.parents):
            self.third_neighbors_raw[i_third].append(i_parent)

  def process_a0_angles_and_third_neighbors_without_dihedral(self):
    """Process raw list of third neighbors withouth ideal dihedral proxy."""
    for neighbor_obj in self.h_connectivity:
      if (neighbor_obj is None): continue
      ih = neighbor_obj.ih
      i_parent = neighbor_obj.a0['iseq']
      self.assign_a0_angles(ih = ih)
      if (neighbor_obj.number_non_h_neighbors != 1 or 'iseq' in neighbor_obj.b1):
        continue
      third_neighbors = self.third_neighbors_raw[i_parent]
      third_neighbors_reduced = []
      alt_conf_neighbors = []
      for i_third in third_neighbors:
        altloc_i_third = self.atoms[i_third].parent().altloc
        if (altloc_i_third == ''):
          third_neighbors_reduced.append(i_third)
        else:
          alt_conf_neighbors.append(i_third)
      third_neighbors_reduced.extend(self.process_alternate_neighbors(
        alt_conf_neighbors = alt_conf_neighbors))
      # If there is no dihedral ideal angle, use randomly first atom
      # in list of third neighbors
      if third_neighbors_reduced:
        if (neighbor_obj.number_h_neighbors == 2):
          self.h_connectivity[ih].b1 = {'iseq': third_neighbors_reduced[0]}
        else:
          self.h_connectivity[ih] = neighbors(
            ih = ih,
            number_non_h_neighbors = 0)
        #self.h_connectivity[ih].b1 = {'iseq': third_neighbors_reduced[0]}
        #self.check_for_plane_proxy(ih)

  def process_plane_proxies(self, planarity_proxies):
    self.plane_h = {}
    for pp in planarity_proxies:
      hlist = []
      for i_test in pp.i_seqs:
        if self.hd_sel[i_test]:
          hlist.append(i_test)
      if hlist:
        self.plane_h[hlist[0]]=hlist[1:]

  def check_for_plane_proxy(self, ih):
    neighbors = self.h_connectivity[ih]
    a0 = neighbors.a0['iseq']
    a1 = neighbors.a1['iseq']
    b1 = neighbors.b1['iseq']
    if (neighbors.h1 and not neighbors.h2): # if there is only 1 H atom as second neighbor
      ih2 = neighbors.h1['iseq']
      if ('dihedral_ideal' not in neighbors.b1 or
          'dihedral_ideal' not in self.h_connectivity[ih2].b1):
        if ih in self.plane_h:
          dihedral = dihedral_angle(
                sites = [self.sites_cart[ih], self.sites_cart[a0],
                self.sites_cart[a1],self.sites_cart[b1]])
          neighbors.b1['dihedral_ideal'] = math.degrees(dihedral)

  def assign_a0_angles(self, ih):
    """ Having a list of dictionaries for the angles involving atom a0,
    assign the angles to the correct set of three atoms. """
    partners = self.h_connectivity[ih]
    if (partners.number_non_h_neighbors > 1):
      a0_iseq = partners.a0['iseq']
      a1_iseq = partners.a1['iseq']
      a2_iseq = partners.a2['iseq']
      if (a1_iseq, a2_iseq) in self.parent_angles[a0_iseq]:
        self.h_connectivity[ih].a0['angle_a1a0a2'] = \
          self.parent_angles[a0_iseq][(a1_iseq, a2_iseq)]
      elif (a2_iseq, a1_iseq) in self.parent_angles[a0_iseq]:
        self.h_connectivity[ih].a0['angle_a1a0a2'] = \
          self.parent_angles[a0_iseq][(a2_iseq, a1_iseq)]
      if ('angle_a1a0a2' not in partners.a0):
        self.h_connectivity[ih] = neighbors(
          ih = ih,
          number_non_h_neighbors = 0)
        return
      if (partners.number_non_h_neighbors == 3):
        a3_iseq = partners.a3['iseq']
        if (a2_iseq, a3_iseq) in self.parent_angles[a0_iseq]:
          self.h_connectivity[ih].a0['angle_a2a0a3'] = \
            self.parent_angles[a0_iseq][(a2_iseq, a3_iseq)]
        elif (a3_iseq, a2_iseq) in self.parent_angles[a0_iseq]:
          self.h_connectivity[ih].a0['angle_a2a0a3'] = \
            self.parent_angles[a0_iseq][(a3_iseq, a2_iseq)]
        if (a3_iseq, a1_iseq) in self.parent_angles[a0_iseq]:
          self.h_connectivity[ih].a0['angle_a3a0a1'] = \
            self.parent_angles[a0_iseq][(a3_iseq, a1_iseq)]
        elif (a1_iseq, a3_iseq) in self.parent_angles[a0_iseq]:
          self.h_connectivity[ih].a0['angle_a3a0a1'] = \
            self.parent_angles[a0_iseq][(a1_iseq, a3_iseq)]
        if ('angle_a2a0a3' not in partners.a0 or 'angle_a3a0a1' not in partners.a0):
          self.h_connectivity[ih] = neighbors(
            ih = ih,
            number_non_h_neighbors = 0)


  def find_third_neighbors(self, dihedral_proxies):
    """ Loop through dihedral angle proxies to find third neighbor
    Fill in neighbors.b1 with iseq and angle proxy"""
    for dp in dihedral_proxies:
      for i_test in dp.i_seqs:
        if (self.h_connectivity[i_test] is None and not self.hd_sel[i_test]):
          continue
        ih = i_test
        i1, i2, i3, i4 = dp.i_seqs
        if (ih == i1):
          i_third = i4
        if (ih == i4):
          i_third = i1
        dihedral = dihedral_angle(
              sites = [self.sites_cart[i1], self.sites_cart[i2],
              self.sites_cart[i3],self.sites_cart[i4]])
        if dihedral is None:
          return
        dihedral_id = dp.angle_ideal
        delta = geometry_restraints.angle_delta_deg(
          angle_1 = math.degrees(dihedral),
          angle_2 = dihedral_id,
          periodicity = dp.periodicity)
        dihedral_ideal = math.degrees(dihedral) + delta
        b1 = {'iseq': i_third, 'dihedral_ideal': dihedral_ideal}
        self.h_connectivity[ih].b1 = b1
        self.assign_b1_for_H_atom_groups(ih = ih, i_third = i_third)

  def assign_b1_for_H_atom_groups(self, ih, i_third):
    """ For atom groups (such as propeller), only one H atom has dihedral proxy
    For the other H atoms of such groups, fill in only iseq"""
    number_h_neighbors = self.h_connectivity[ih].number_h_neighbors
    if (number_h_neighbors > 0):
      i_h1 = self.h_connectivity[ih].h1['iseq']
      if (self.h_connectivity[i_h1] is not None):
        self.h_connectivity[i_h1].b1 = {'iseq': i_third}
      if (number_h_neighbors == 2):
        i_h2 = self.h_connectivity[ih].h2['iseq']
        if (self.h_connectivity[i_h2] is not None):
          self.h_connectivity[i_h2].b1 = {'iseq': i_third}

  def determine_second_neighbors_H(self, second_neighbors_reduced):
    """Determine if there are H atoms among second neighbors and store them
    in a list 'second_neighbors_H' """
    second_neighbors_H = []
    for iseq in second_neighbors_reduced:
      if (self.hd_sel[iseq]):
        second_neighbors_H.append(iseq)
    return second_neighbors_H

  def assign_second_neighbors(
          self, ih, i_parent, neighbors_list, neighbors_list_H):
    """With the information of second neighbors, fill in the dictionaries for
    each atom (a1, a2, a3, according to which is present)."""
    number_h_neighbors = 0
    number_non_h_neighbors = 0
    for iseq, n in zip(neighbors_list, [1,2,3]):
      number_non_h_neighbors = number_non_h_neighbors + 1
      a = self.make_neighbor_dict(iseq = iseq, ih = ih, i_parent = i_parent)
      if (n == 1): self.h_connectivity[ih].a1 = a
      if (n == 2): self.h_connectivity[ih].a2 = a
      if (n == 3): self.h_connectivity[ih].a3 = a
    for iseqh, nh in zip(neighbors_list_H, [1,2]):
      number_h_neighbors = number_h_neighbors + 1
      ah = self.make_neighbor_dict(iseq = iseqh, ih = ih, i_parent = i_parent)
      if (nh == 1): self.h_connectivity[ih].h1 = ah
      if (nh == 2): self.h_connectivity[ih].h2 = ah
    self.h_connectivity[ih].number_h_neighbors = number_h_neighbors
    self.h_connectivity[ih].number_non_h_neighbors = number_non_h_neighbors

  def make_neighbor_dict(self, iseq, ih, i_parent):
    angle = self.angle_dict[(ih, i_parent, iseq)]
    neighbor = {'iseq': iseq, 'angle_ideal': angle}
    return neighbor

  def process_alternate_neighbors(self, alt_conf_neighbors):
    """ For a list of atoms in alternative conformations, retain singles and
    those with higher occupancy
    Example: [CA-A, N-A, N-B] --> [CA-A, N-B] if occ(N-B) > occ(N-A)"""
    alt_conf_neighbors_temp = []
    alt_conf_neighbors_reduced = []
    # if only one atom in alt conf list, not necessary to search for other atoms
    if (len(alt_conf_neighbors) == 1):
      alt_conf_neighbors_reduced.append(alt_conf_neighbors[0])
    # Go through each atom in dc, get the name and make temporary dictionary
    # with iseq as key which points to occupancy: {i_1:occA, i_2:occB, i_3:occC}
    # Then choose atom with maximum occupancy
    else:
      for i_second in alt_conf_neighbors:
        altloc_dict_temp = {}
        if (i_second in alt_conf_neighbors_temp): continue
        name_i_second = self.atoms[i_second].name
        altloc_i_second = self.atoms[i_second].parent().altloc
        # check if all neighbors are in the same alt conf (otherwise no angle proxy)
        if alt_conf_neighbors_reduced:
          i_previous = alt_conf_neighbors_reduced[0]
          altloc_previous = self.atoms[i_previous].parent().altloc
          if (altloc_i_second != altloc_previous):
            alt_conf_neighbors_temp.append(i_second)
            continue
        altloc_dict_temp[i_second] = self.atoms[i_second].occ
        for ag in self.atoms[i_second].parent().parent().atom_groups():
          for ag_atom in ag.atoms():
            if (ag_atom.i_seq == altloc_i_second): continue
            if (ag_atom.i_seq in alt_conf_neighbors_temp): continue
            if (ag_atom.name == name_i_second and
                ag_atom.i_seq in alt_conf_neighbors):
              altloc_dict_temp[ag_atom.i_seq] = ag_atom.occ
        i_second_max = max(altloc_dict_temp, key=lambda k: altloc_dict_temp[k])
        alt_conf_neighbors_reduced.append(i_second_max)
        # Store "used" atoms in temp list --> avoids going through atoms twice
        for index in altloc_dict_temp.keys():
          alt_conf_neighbors_temp.append(index)
    return alt_conf_neighbors_reduced

  def count_H(self):
    """ Check if number H/D atoms in the in the model and in h_connectivity
    are the same"""
    self.connectivity_slipped = []
    number_h_input_model = self.hd_sel.count(True)
    number_h_connectivity = \
      len(self.h_connectivity) - self.h_connectivity.count(None)
    if (number_h_input_model != number_h_connectivity):
      self.find_mismatch()

  def find_mismatch(self):
    list_H_connect = []
    list_H = []
    for item in self.h_connectivity:
      if item: list_H_connect.append(item.ih)
    for atom in self.atoms:
      if (atom.element_is_hydrogen()):
        list_H.append(atom.i_seq)
    set_list_H_connect = set(list_H_connect)
    slipped = [x for x in list_H if x not in set_list_H_connect]
    self.connectivity_slipped = slipped

  def add_slipped(self):
    for ih in self.connectivity_slipped:
      self.h_connectivity[ih] = neighbors(
        ih = ih,
        number_non_h_neighbors = 0)

  def get_diagnostics(self):
    h_in_connectivity = []
    for neighbor_obj in self.h_connectivity:
      if (neighbor_obj is not None):
        h_in_connectivity.append(neighbor_obj.ih)
    return group_args(
      double_H = self.double_H,
      connectivity_slipped = self.connectivity_slipped,
      h_in_connectivity = h_in_connectivity)

  def hd_selection(self):
    """Get a selector array for all hydrogen and deuterium scatterers of the structure.
    :returns: an array to select all H and D scatterers of the structure
    :rtype: boolean[]
    """
    result = flex.bool()
    self.list_H = []
    for atom in self.atoms:
      result.append(atom.element_is_hydrogen())
    return result

