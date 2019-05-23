from __future__ import division, print_function
from libtbx import group_args
from scitbx import matrix
from libtbx.str_utils import make_sub_header
from libtbx.utils import null_out


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


def unknown_pairs_present(model):
  """
  Test if PDB file contains unknown type pairs

  Parameters:
    model (obj): model object

  Returns:
    (bool): True if PDB file contains unknown type pairs
  """
  grm = model.get_restraints_manager()
  sites_cart = model.get_sites_cart()
  site_labels = model.get_xray_structure.scatterers().extract_labels()
  pp= grm.pair_proxies(sites_cart=sites_cart,site_labels=site_labels)
  return (pp.nonbonded_proxies.n_unknown_nonbonded_type_pairs != 0)


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
    #
    self.sort_clashes()


  def show(self, log=null_out(), show_clashscore=True):
    """
    Print all clashes in a table.
    """
    make_sub_header(' Nonbonded overlaps', out=log)
    if self._clashes_dict:
      # General information
      results = self.get_results()
      result_str = '{:<18} : {:5d}'
      print(result_str.format(' Number of clashes', results.n_clashes), file=log)
      result_str = '{:<18} : {:5.2f}'
      if show_clashscore:
        print(result_str.format(' Clashscore', results.clashscore), file=log)
      # print table with all overlaps
      labels =  ['\n' + "Overlapping residues info","model distance","overlap",
                 "symmetry"]
      lbl_str = '{:^33}|{:^16}|{:^11}|{:^15}'
      table_str = '{:>16}|{:>16}|{:^16.2f}|{:^11.2}|{:^15}|'
      print(lbl_str.format(*labels), file=log)
      print('-'*78, file=log)
      atoms = self.model.get_atoms()
      for iseq_tuple, record in self._clashes_dict.iteritems():
        i_seq, j_seq = iseq_tuple
        overlap = record[2]
        if record[4] is not None:
          symop = record[4]
        else: symop = ''
        i_id_str = atoms[i_seq].id_str().replace('pdb=','').replace('"','')
        j_id_str = atoms[j_seq].id_str().replace('pdb=','').replace('"','')
        line = [i_id_str, j_id_str,round(record[0], 2),round(overlap, 2), symop]
        print(table_str.format(*line), file=log)
      print('-'*78, file=log)
    else:
      print('No clashes found', file=log)


  def is_clashing(self, iseq):
    """
    Test if a particular atom is involved in a clash.

    Parameters:
      iseq (int): atom i_seq number

    Returns:
      (bool): True if the atom is involved in a clash
    """
    is_clashing = False
    if self._clashes_dict:
      i_seqs, j_seqs = zip(*self._clashes_dict)
      if iseq in i_seqs or iseq in j_seqs:
        is_clashing = True
    return is_clashing


  def sort_clashes(self,
                   sort_vdW            = False,
                   sort_model_distance = False,
                   sort_overlap        = False,
                   sort_symmetry       = False):
    """
    Sort clashes according to vdW distance, model distance, overlap or symmetry
    """
    from collections import OrderedDict
    options = [sort_vdW, sort_model_distance, sort_overlap, sort_symmetry]
    if (options.count(True) == 0):
      sort_overlap = True
    elif (options.count(True) > 1):
      raise Sorry('Can only sort by one value.')
    if sort_model_distance: key = 0
    if sort_vdW: key = 1
    if sort_overlap: key = 2
    if sort_symmetry: key = 4
    self._clashes_dict = OrderedDict(
      sorted(self._clashes_dict.items(), key=lambda x: x[1][key]))


  def _obtain_symmetry_clashes(self):
    """
    Get clashes due to symmetry
    """
    self._symmetry_clashes_dict = dict()
    n_clashes_sym = 0
    clashscore_sym = 0
    for iseq_tuple, record in self._clashes_dict.iteritems():
      if record[4] is not None:
        self._symmetry_clashes_dict[iseq_tuple] = record
    if self._symmetry_clashes_dict:
      n_clashes_sym = len(self._symmetry_clashes_dict)
      # Does clashscore_sym actually make sense?
      n_atoms = self.model.size()
      clashscore_sym = n_clashes_sym * 1000 / n_atoms
    return n_clashes_sym, clashscore_sym


  def _obtain_macro_mol_clashes(self):
    """
    Get clashes involving macro-mol atoms only
    """
    self._macro_mol_clashes_dict = dict()
    n_clashes_macro_mol = 0
    clashscore_macro_mol = 0
    macro_mol_sel = self.model.selection(string = 'protein')
    for iseq_tuple, record in self._clashes_dict.iteritems():
      if (macro_mol_sel[iseq_tuple[0]] and macro_mol_sel[iseq_tuple[1]]
          and record[4] is None):
        self._macro_mol_clashes_dict[iseq_tuple] = record
    if self._macro_mol_clashes_dict:
      n_clashes_macro_mol = len(self._macro_mol_clashes_dict)
      clashscore_macro_mol = n_clashes_macro_mol * 1000 / self.model.select(macro_mol_sel).size()
    return n_clashes_macro_mol, clashscore_macro_mol


  def get_results(self):
    """
    Accessor for results
    """
    # overall
    n_clashes = len(self._clashes_dict)
    n_atoms = self.model.size()
    clashscore = n_clashes * 1000 / n_atoms
    # due to symmetry
    n_clashes_sym, clashscore_sym = self._obtain_symmetry_clashes()
    # macromolecule ('protein')
    n_clashes_macro_mol, clashscore_macro_mol = self._obtain_macro_mol_clashes()

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
    """
    Print all hbonds in a table.
    """
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
    """
    True/False if any clashes were found.
    """
    has_clashes = False
    if not self._clashes:
      self._process_nonbonded_proxies(find_clashes = True)
      if self._clashes_dict:
        has_clashes = True
    return has_clashes


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
    Process nonbonded_proxies to find bonds, interactions and clashes.
    Clashes code refactored from Youval Dar's code for nonbonded_overlaps (LBNL 2013)
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
      self._hbonds  = hbonds(hbonds_dict = dict())
      return

    fsc0 = grm.shell_sym_tables[0].full_simple_connectivity()
    fsc2 = grm.shell_sym_tables[2].full_simple_connectivity()

    self._clashes_dict = dict()
    self._hbonds_dict  = dict()
    self._mult_clash_dict = dict()

    # loop over nonbonded proxies do stuff and fill in the dicts:
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
          self._clashes_dict[(i_seq, j_seq)] = \
            [model_distance, vdw_sum, abs(delta), symop_str, symop]

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
    """
    Determine if a nonbonded proxy is a clash.

    Parameters:
      i_seq (int): atom i_seq
      j_seq (int): atom i_seq
      hd_sel (bool array)        hd_sel[i] returns True of False if atom i is H or not
      fsc0 (dict of lists of int): dictionary with a list of all atoms connected to an atom
      model_distance (float): distance between atom i and atom j

    Returns:
      bool (is_clash): if a nonbonded proxy is a clash
    """
    is_clash = False
    # Exclude 1-5 interaction of H atom and heavy atom
    is_1_5_interaction = check_if_1_5_interaction(
             i_seq = i_seq,
             j_seq = j_seq,
             hd_sel = hd_sel,
             full_connectivity_table = fsc0)
    if not is_1_5_interaction:
      #if i_seq > j_seq:
      #  i_seq, j_seq = j_seq, i_seq
      # Check to prevent that symmetry clashes are counted twice
      #if (i_seq, j_seq) not in self._clashes_dict.keys():
      is_clash = True
      if (i_seq not in self._mult_clash_dict): self._mult_clash_dict[i_seq] = list()
      if (j_seq not in self._mult_clash_dict): self._mult_clash_dict[j_seq] = list()
      self._mult_clash_dict[i_seq].append(j_seq)
      self._mult_clash_dict[j_seq].append(i_seq)

    return is_clash


  def _process_clashes(self, sites_cart, fsc0):
    """
    Process clashes from previous loop through nonbonded_proxies.

    This step is necessary to filter out clashes with common atoms.
    X-H ~~~ Y might produce two clashes, one between X and Y, the other
    between H and Y. This step filters the raw results and keeps the shorter
    of the two clashes (if an angular cutoff is above a limit)
    """
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
            # check if atoms are inline
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
