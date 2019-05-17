from __future__ import division
#from libtbx.utils import Sorry
#from scitbx.array_family import flex
#from libtbx import easy_run
#import iotbx.pdb
#import math
#import sys
from libtbx import group_args
from scitbx import matrix


def check_if_1_5_interaction(
      i_seq,
      j_seq,
      hd_sel,
      full_connectivity_table):
  """
  Checks if there is 1-5 interaction between a hydrogen (H) and heavy atom (X): H-A-A-A-X

  Parameters:
    i_seq (int): atom i_seq
    j_seq (int): atom i_seq
    hd_sel (bool array)        hd_sel[i] returns True of False if atom i is H or not
    full_connectivity_table (dict of lists of int):  dictionary with a list
                                             of all atoms connected to atom i

  Returns:
    bool: True/False if there is 1-5 hydrogen and heavy atom interaction
  """
  # check if we have hydrogen - heavy atom interaction
  xor = lambda a,b: (a or b) and not (a and b)
  if xor(hd_sel[i_seq],hd_sel[j_seq]):
    # starting with hydrogen will make process shorter
    if not hd_sel[i_seq]:
      i_seq,j_seq = j_seq,i_seq
    # build connection table of i_seq, 4 steps deep
    atoms_numbers = dict([(i_seq, 0)]) # i_seq is in zero distance
    used_connections = {i_seq}
    new_connections = {i_seq}
    for i in range(2,6):
      connections = set()
      for key in new_connections:
        # add all new connections for the current step
        connections = connections.union(set(full_connectivity_table[key]))
      # Remove the connection that were already used
      new_connections = connections - used_connections
      # Add the new connection to the used once
      used_connections = used_connections.union(connections)
      # Add the new atoms with their distance from key
      for new_atom in new_connections:
        atoms_numbers[new_atom] = i
    # return true if j_seq in the is 1-5 connection
    return (j_seq in atoms_numbers) and (atoms_numbers[j_seq] == 5)
  else:
    return False



def cos_vec(u,v,w):
  """(tuple,tuple) -> float

  Calculate the cosine to evaluate whether clashing atoms are inline
  A1 clashes with A2 and A3. Find out if A2 and A3 are inline.
  A1 ~~~ A2
  A1 ~~~ A3

  Parameters:
  u: vector of clashing atom (A2 or A3, order does not matter)
  v: vector of clashing atom (A2 or A3, order does not matter)
  w: vector of common clashing atom (A1)

  Returns:
    float (cos_angle): the cosine of the angle between center of the common atom
      and the mid point between the other two atoms.
  """
  u = matrix.col(u)
  v = matrix.col(v)
  w = matrix.col(w)

  vec1 = w - (u/2 + v/2)
  vec2 = u - v

  try:
    cos_angle = abs(vec1.normalize().dot(vec2.normalize()))
  except ZeroDivisionError:
    cos_angle = 1
  return cos_angle




class clashes(object):
  """
  Class for clashes
  """
  def __init__(self, clashes_dict, model):
    """
    clashes_dict  {(iseq, jseq):(distance, sum_vdw_radii)}
    iseq          atom i
    jseq          atom j
    distance      distance between atom i and atom j
    sum_vdw_radii sum of vdW radii
    """
    self._clashes_dict = clashes_dict
    self.model = model

  def show(self):
    pass

  def is_clashing(self, iseq):
    pass

  def sort_clashes(self):
    pass

  def get_results(self):
    n_clashes = len(self._clashes_dict)
    n_atoms = self.model.size()
    clashscore = n_clashes * 1000 / n_atoms
    #print(dir(self.model))
    return group_args(
      n_clashes = n_clashes,
      clashscore = clashscore)



class hbonds(object):
  """
  Class for hbonds
  """
  def __init__(self, hbonds_dict):
    """
    hbonds_dict = {(iseq, jseq, kseq):(H_A_distance, X_A_distance, X_H_A_angle)}
    hydrogen bond: X-H...A
    iseq          atom X (donor heavy atom)
    jseq          atom H (donor H atom)
    kseq          atom A (acceptor atom)
    H_A_distance
    X_A_distance
    X_H_A_angle
    """
    self.hbonds_dict = hbonds_dict

  def show(self):
    pass

  def forms_hbond(self, iseq):
    pass

  def sort_hbonds(self, sort_distances = True, sort_angles = False):
    pass


class manager():

  def __init__(self,
               model):
    self.model = model

    #
    self._clashes = None
    self._hbonds  = None

    # in manager or do we enfore that input model has H?
    # self._add_H_atoms() ????

  def get_clashes(self):
    """
    Accessor for clashes object
    """
    if self._clashes is None:
      self._process_nonbonded_proxies(find_clashes = True)
      return self._clashes
    else:
      return self._clashes

  def get_hbonds(self):
    """
    Accessor for hbonds object
    """
    if not self._hbonds:
      self._process_nonbonded_proxies(find_hbonds = True)
    else:
      return self._hbonds

  def has_hbonds(self):
    # necessary?
    pass

  def has_clashes(self):
    # necessary?
    pass

  def show(self):
    """
    Print information
    """
    if self.has_clashes():
      self._clashes.show()
    if self.has_hbonds():
      self._hbonds.show()

  def _process_nonbonded_proxies(self,
                                 find_clashes = True,
                                 find_hbonds = False):
    """
    Here is where the calculations are done
    Either all is done at once (clashes, hbonds, other?)
    or it will be modular (use find_clashes and find_hbodns parameters)
    """
    grm = self.model.get_restraints_manager().geometry
    xrs = self.model.get_xray_structure()
    sites_cart = self.model.get_sites_cart()
    site_labels = xrs.scatterers().extract_labels()
    hd_sel = self.model.get_hd_selection()

    pair_proxies = grm.pair_proxies(
                        sites_cart  = sites_cart,
                        site_labels = site_labels)
    proxies_info_nonbonded = pair_proxies.nonbonded_proxies.get_sorted(
      by_value    = "delta",
      sites_cart  = sites_cart,
      site_labels = site_labels)

    if proxies_info_nonbonded is not None:
      nonbonded_list = proxies_info_nonbonded[0]
    else:
      nonbonded_list = []
      # create 'empty' instance of results class
      self._clashes = clashes(clashes_dict = dict())
      self._hbonds = hbonds(hbonds_dict = dict())
      return

    fsc0=grm.shell_sym_tables[0].full_simple_connectivity()
    fsc2=grm.shell_sym_tables[2].full_simple_connectivity()

    self._clashes_dict = dict()
    self._hbonds_dict  = dict()

    # loop over nonbonded proxies do stuff and fill in the dicts:
    # self._clashes_dict[(iseq, jseq)] = (relevant info)
    # self._hbonds_dict[(iseq, jseq)] = (relevant info)
    for item in nonbonded_list:
      i_seq          = item[1]
      j_seq          = item[2]
      model_distance = item[3]
      vdw_sum        = item[4]
      symop_str      = item[5] # TODO probably not necessary
      symop          = item[6]

      # check for overlap
      delta = model_distance - vdw_sum
      if (delta < -0.40):
        is_clash = self._is_clash(
                    i_seq = i_seq,
                    j_seq = j_seq,
                    hd_sel = hd_sel,
                    fsc0 = fsc0,
                    sites_cart = sites_cart)
        if is_clash:
          self._clashes_dict[(i_seq, j_seq)] = [model_distance, vdw_sum, symop_str, symop]

    self._clashes = clashes(
                      clashes_dict = self._clashes_dict,
                      model        = self.model)
    print(self._clashes)
    print(self._clashes_dict.keys())

  def _is_clash(self,
                i_seq,
                j_seq,
                hd_sel,
                fsc0,
                sites_cart):
    is_clash = False
    # Check if there is 1-5 interaction
    is_1_5_interaction = check_if_1_5_interaction(
             i_seq = i_seq,
             j_seq = j_seq,
             hd_sel = hd_sel,
             full_connectivity_table = fsc0)
    if not is_1_5_interaction:
      if i_seq > j_seq:
        i_seq, j_seq = j_seq, i_seq
      iseq_tuple = (i_seq, j_seq)
      # Check if either of the atoms is already involved in another clash
      for i in iseq_tuple:
        multiples = [item for item in self._clashes_dict.keys() if i in item] # is this slow?
        # for atoms that overlap more than once, check for inline overlaps
        if multiples:
          for multiple in multiples:
            multiple_atoms = [x for x in list(multiple + iseq_tuple) if i != x]
            # other way:
            # multiple_atoms = list(set(iseq_tuple + multiple) - (set(iseq_tuple) & set(multiple)))
            # ignore overlaps that are cause by symmetry operation -->
            # TODO: not sure why these should be ignored
            if (len(multiple_atoms) == 2): # should be always 2, per definitionem, probably assert?
              multiple_1 = multiple_atoms[0]
              multiple_2 = multiple_atoms[1]
            # test inline only if the two atoms that overlap with the common atom are connected
              if multiple_1 in fsc0[multiple_2]:
                atom_1_xyz = sites_cart[multiple_1]
                atom_2_xyz = sites_cart[multiple_2]
                atom_i_xyz = sites_cart[i]
                cos_angle = cos_vec(atom_1_xyz, atom_2_xyz, atom_i_xyz)
                # atoms are inline if cosine of the angle between vectors > 0.707 (45 degrees)
                # TODO where does that number come from?
                if abs(cos_angle) > 0.707 and (atom_1_xyz != atom_2_xyz):
                  if self._clashes_dict[multiple][0] > model_distance:
                    del self._clashes_dict[multiple]
                    is_clash = True
                  else:
                    continue
              else:
                is_clash = True
        else:
          is_clash = True
    #
    return is_clash

    # create class:
    # self._clashes = clashes(clashes_dict = clashes_dict)
    # self._hbonds = hbonds(hbonds_dict = hbonds_dict)

