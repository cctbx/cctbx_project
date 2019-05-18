from __future__ import division, print_function
#from libtbx.utils import Sorry
#from scitbx.array_family import flex
#from libtbx import easy_run
#import iotbx.pdb
#import math
#import sys
from libtbx import group_args
from scitbx import matrix
from libtbx.str_utils import make_sub_header
import sys


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
  # check if there is hydrogen - heavy atom interaction
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
      # Add the new connection to the used ones
      used_connections = used_connections.union(connections)
      # Add the new atoms with their distance from key
      for new_atom in new_connections:
        atoms_numbers[new_atom] = i
    # return true if j_seq in the is 1-5 connection
    return (j_seq in atoms_numbers) and (atoms_numbers[j_seq] == 5)
  else:
    return False


def cos_vec(u, v, w):
  """
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
    '''
    clashes_dict  {(iseq, jseq):(distance, sum_vdw_radii)}
    iseq          atom i
    jseq          atom j
    distance      distance between atom i and atom j
    sum_vdw_radii sum of vdW radii
    '''
    self._clashes_dict = clashes_dict
    self.model = model


# TODO: what's the log?
  def show(self):
    '''
    Print all clashes in a table.
    '''
    # TODO : make sure self._clashes_dict is not empty
    make_sub_header(' Nonbonded overlaps', out=sys.stdout)
    # General information
    results = self.get_results()
    result_str = '{:<18} : {:5.2f}'
    print(result_str.format(' Number of clashes', results.n_clashes))
    print(result_str.format(' Clashscore', results.clashscore) + '\n')
    # print table with all overlaps
    labels =  ["Overlapping residues info","model distance","overlap",
               "symmetry"]
    lbl_str = '{:^33}|{:^16}|{:^11}|{:^10}'
    table_str = '{:>16}|{:>16}|{:^16.2f}|{:^11.2}|{:^10}|'
    print(lbl_str.format(*labels))
    print('-'*74)
    atoms = self.model.get_atoms()
    for iseq_tuple, record in self._clashes_dict.iteritems():
      i_seq, j_seq = iseq_tuple
      delta = record[1] - record[0]
      if record[3] is not None:
        symop = record[3]
      else: symop = ''
      i_id_str = atoms[i_seq].id_str().replace('pdb=','').replace('"','')
      j_id_str = atoms[j_seq].id_str().replace('pdb=','').replace('"','')
      line = [i_id_str, j_id_str,round(record[0], 2),round(delta, 2), symop]
      print(table_str.format(*line))
    print('-'*74)


  def is_clashing(self, iseq):
    pass


  def sort_clashes(self):
    pass


  def _obtain_symmetry_clashes(self):
    self._symmetry_clashes_dict = dict()
    for iseq_tuple, record in self._clashes_dict.iteritems():
      if record[3] is not None:
        self._symmetry_clashes_dict[iseq_tuple] = record

  def _obtain_macro_mol_clashes(self):
    self._macro_mol_clashes_dict = dict()
    proxies = self.model.all_chain_proxies
    cache = proxies.pdb_hierarchy.atom_selection_cache()
    macro_mol_sel = proxies.selection(
      cache  = cache,
      string = 'protein or dna or rna')
    for iseq_tuple, record in self._clashes_dict.iteritems():
      if macro_mol_sel[iseq_tuple[0]] and macro_mol_sel[iseq_tuple[1]] and record[3] is None:
        self._macro_mol_clashes_dict[iseq_tuple] = record


  def get_results(self):
    # overall
    n_clashes = len(self._clashes_dict)
    n_atoms = self.model.size()
    clashscore = n_clashes * 1000 / n_atoms
    # due to symmetry
    self._obtain_symmetry_clashes()
    if self._symmetry_clashes_dict:
      n_clashes_sym = len(self._symmetry_clashes_dict)
      # Does this number actually make sense?
      clashscore_sym = n_clashes_sym * 1000 / n_atoms
    else:
      # None or 0?
      n_clashes_sym = 0
      clashscore_sym = 0
    self._obtain_macro_mol_clashes()
    if self._macro_mol_clashes_dict:
      n_clashes_macro_mol = len(self._macro_mol_clashes_dict)
      clashscore_macro_mol = n_clashes_macro_mol * 1000 / n_atoms
    else:
      n_clashes_macro_mol = 0
      clashscore_macro_mol = 0

    return group_args(
             n_clashes      = n_clashes,
             clashscore     = clashscore,
             n_clashes_sym  = n_clashes_sym,
             clashscore_sym = clashscore_sym,
             n_clashes_macro_mol  = n_clashes_macro_mol,
             clashscore_macro_mol = clashscore_macro_mol)


class hbonds(object):
  '''
  Class for hbonds
  '''
  def __init__(self, hbonds_dict):
    '''
    hbonds_dict = {(iseq, jseq, kseq):(H_A_distance, X_A_distance, X_H_A_angle)}
    hydrogen bond: X-H...A
    iseq          atom X (donor heavy atom)
    jseq          atom H (donor H atom)
    kseq          atom A (acceptor atom)
    H_A_distance
    X_A_distance
    X_H_A_angle
    '''
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

    # add H in manager or do we enfore that input model has H?
    # self._add_H_atoms() ????

  def get_clashes(self):
    '''
    Accessor for clashes object
    '''
    if self._clashes is None:
      self._process_nonbonded_proxies(find_clashes = True)
      return self._clashes
    else:
      return self._clashes

  def get_hbonds(self):
    '''
    Accessor for hbonds object
    '''
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
    '''
    Print information
    '''
    if self.has_clashes():
      self._clashes.show()
    if self.has_hbonds():
      self._hbonds.show()


  def _process_nonbonded_proxies(self,
                                 find_clashes = True,
                                 find_hbonds = False):
    '''
    Here is where the calculations are done
    Either all is done at once (clashes, hbonds, other?)
    or it will be modular (use find_clashes and find_hbodns parameters)
    '''
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
      self._hbonds  = hbonds(hbonds_dict = dict())
      return

    fsc0 = grm.shell_sym_tables[0].full_simple_connectivity()
    fsc2 = grm.shell_sym_tables[2].full_simple_connectivity()

    self._clashes_dict = dict()
    self._hbonds_dict  = dict()
    self._mult_clash_dict = dict()

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

      # Find all clashes
      delta = model_distance - vdw_sum
      if (delta < -0.40):
        is_clash = self._is_clash(
                    i_seq = i_seq,
                    j_seq = j_seq,
                    hd_sel = hd_sel,
                    fsc0 = fsc0,
                    model_distance = model_distance)
        if is_clash:
          self._clashes_dict[(i_seq, j_seq)] = [model_distance, vdw_sum, symop_str, symop]

    # Remove clashes involving common atoms (cannot be done in first loop!)
    self._process_clashes(sites_cart = sites_cart, fsc0 = fsc0)
    # Create clashes class
    self._clashes = clashes(
                      clashes_dict = self._clashes_dict,
                      model        = self.model)


  def _is_clash(self,
                i_seq,
                j_seq,
                hd_sel,
                fsc0,
                model_distance):
    '''
    Determine if a nonbonded proxy is a clash.

    Parameters:
      i_seq (int): atom i_seq
      j_seq (int): atom i_seq
      hd_sel (bool array)        hd_sel[i] returns True of False if atom i is H or not
      fsc0 (dict of lists of int): dictionary with a list of all atoms connected to an atom
      model_distance (float): distance between atom i and atom j

    Returns:
      bool (is_clash): if a nonbonded proxy is a clash
    '''
    is_clash = False
    is_1_5_interaction = check_if_1_5_interaction(
             i_seq = i_seq,
             j_seq = j_seq,
             hd_sel = hd_sel,
             full_connectivity_table = fsc0)
    if not is_1_5_interaction:
      if i_seq > j_seq:
        i_seq, j_seq = j_seq, i_seq
      # Check to prevent that ymmetry clashes are counted twice
      if (i_seq, j_seq) not in self._clashes_dict.keys():
        is_clash = True
        if (i_seq not in self._mult_clash_dict): self._mult_clash_dict[i_seq] = list()
        if (j_seq not in self._mult_clash_dict): self._mult_clash_dict[j_seq] = list()
        self._mult_clash_dict[i_seq].append(j_seq)
        self._mult_clash_dict[j_seq].append(i_seq)

    return is_clash


  def _process_clashes(self, sites_cart, fsc0):
    clashes_to_be_removed = list()
    for i_seq, j_seq_list in self._mult_clash_dict.iteritems():
      n_multiples = len(j_seq_list)
      if n_multiples <= 1: continue
      for i in range(n_multiples-1):
        for j in range(i+1, n_multiples):
          multiple_1 = j_seq_list[i]
          multiple_2 = j_seq_list[j]
        # test inline only if the two atoms that overlap with the common atom are connected
          if multiple_1 in fsc0[multiple_2]:
            atom_1_xyz = sites_cart[multiple_1]
            atom_2_xyz = sites_cart[multiple_2]
            atom_i_xyz = sites_cart[i_seq]
            cos_angle = cos_vec(atom_1_xyz, atom_2_xyz, atom_i_xyz)

            if abs(cos_angle) > 0.707 and (atom_1_xyz != atom_2_xyz):
              tuple1 = [i_seq, multiple_1]
              tuple2 = [i_seq, multiple_2]
              tuple1.sort()
              tuple2.sort()
              tuple1 = tuple(tuple1)
              tuple2 = tuple(tuple2)
              if tuple1 in self._clashes_dict.keys() and tuple2 in self._clashes_dict.keys():
                if self._clashes_dict[tuple1][0] < self._clashes_dict[tuple2][0]:
                  clashes_to_be_removed.append(tuple2)
                else:
                  clashes_to_be_removed.append(tuple1)
    for clash_tuple in clashes_to_be_removed:
      if clash_tuple in self._clashes_dict.keys():
        del self._clashes_dict[clash_tuple]
#
