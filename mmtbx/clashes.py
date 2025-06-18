"""Obtain information on clashing atoms"""
from __future__ import absolute_import, division, print_function
import mmtbx.monomer_library
import mmtbx.monomer_library.pdb_interpretation
import sys
from cctbx.geometry_restraints.linking_class import linking_class
from six.moves import zip
origin_ids = linking_class()
from libtbx import adopt_init_args
import mmtbx.model


# XXX Need to support model_manager object input

# MARKED_FOR_DELETION_OLEG
# REASON: bad practices for obtaining GRM. remove_clashes should be used instead.
def from_xray_structure(
      clash_threshold,
      xray_structure,
      geometry_restraints_manager):
  # Find covalently bonded atoms to exclude them from non-bonded clash list
  bond_proxies_simple, asu = geometry_restraints_manager.get_all_bond_proxies(
      sites_cart = xray_structure.sites_cart())
  bonded_i_seqs = []
  for bp in bond_proxies_simple:
    bonded_i_seqs.append(bp.i_seqs)
  # Find all pairs within clash_threshold
  pair_asu_table = xray_structure.pair_asu_table(
    distance_cutoff=clash_threshold)
  pair_sym_table = pair_asu_table.extract_pair_sym_table()
  atom_pairs_i_seqs = pair_sym_table.simple_edge_list()
  # Get only non-bonded clashes
  return list(set(atom_pairs_i_seqs).difference(set(bonded_i_seqs)))

class from_pdb(object):

  def __init__(self, pdb_str,
      clash_threshold, remove_clash_guard=False,
      remove_max_reasonable_bond_distance=False,
      allow_polymer_cross_special_position=False):
    params = mmtbx.monomer_library.pdb_interpretation.master_params.extract()
    if remove_clash_guard:
      params.clash_guard.nonbonded_distance_threshold=None
      params.proceed_with_excessive_length_bonds=True
    if remove_max_reasonable_bond_distance:
      params.max_reasonable_bond_distance=None
    if allow_polymer_cross_special_position:
      params.allow_polymer_cross_special_position=True

    ppf = mmtbx.monomer_library.pdb_interpretation.process(
      mon_lib_srv              = mmtbx.monomer_library.server.server(),
      ener_lib                 = mmtbx.monomer_library.server.ener_lib(),
      raw_records              = pdb_str,
      params                   = params,
      strict_conflict_handling = True,
      force_symmetry           = True,
      log                      = None)
    self.atoms = list(ppf.all_chain_proxies.pdb_hierarchy.atoms())
    self._pairs = from_xray_structure(
      clash_threshold             = clash_threshold,
      xray_structure              = ppf.xray_structure(),
      geometry_restraints_manager = ppf.geometry_restraints_manager())
    self._ppf=ppf  #so that parents are remembered

  def format_clash_string(self, pair):
    a1 = self.atoms[pair[0]]
    a2 = self.atoms[pair[1]]
    l1,l2 = a1.pdb_label_columns(), a2.pdb_label_columns()
    dist = a1.distance(a2)
    fmt = "'%s' clash with '%s' distance %4.2f A"
    return fmt%(l1,l2,dist)

  def clashing_pairs(self):
    return self._pairs

  def show(self, log=None):
    if(log is None): log = sys.stdout
    for pair in self._pairs:
      print(self.format_clash_string(pair=pair), file=log)

# END_MARKED_FOR_DELETION_OLEG

main_chain_plus_cb=[
       "N","CA","C","O","CB",
       "O5'","C5'","C4'","O4'",
              "C1'","C2'","C3'",
              "O3'","P","OP1","OP2","OP3","O1P","O2P","O3P",
       "O2'" ]

class remove_clashes(object):
  # Remove side-chains/entire residues that are causing clashes or have
  #  bad bonds

  def __init__(self,
      model=None,
      non_bonded_deviation_threshold=0.4,
      bonded_deviation_threshold=0.2,
      max_fraction=None):

    adopt_init_args(self, locals())
    self.side_chains_removed=0
    self.residues_removed=0

    self.get_labels_and_proxies()
    if max_fraction:
      max_items=int(max_fraction*self.model.get_number_of_atoms())
    else:
      max_items=None
    self.worst_non_bonded_pairs_and_delta=self.get_bad_nonbonded(
       max_items=max_items)
    self.worst_bonded_pairs_and_delta=self.get_bad_bonded(
       max_items=max_items)

    self.get_remove_selection_string()
    # redo model after removing residues/side chains specified
    self.pare_model()


  def pare_model(self):
    if not self.remove_selection_string:
      self.new_model=self.model.deep_copy() # nothing to do
    else:
      sel1 = self.model.selection(string = " NOT (%s) " %(
         self.remove_selection_string))
      self.new_model=self.model.select(sel1)

  def get_remove_selection_string(self):
    # Now choose the worst side-chains or main chain to remove.
    # prioritize on shorter fragments and side chains first
    remove_residue_list=[]
    remove_side_chain_list=[]
    for [i,j,delta],threshold in zip(
       self.worst_non_bonded_pairs_and_delta+self.worst_bonded_pairs_and_delta,
        len(self.worst_non_bonded_pairs_and_delta)*
            [self.non_bonded_deviation_threshold]+
        len(self.worst_bonded_pairs_and_delta)*
            [self.bonded_deviation_threshold]):
      if delta > -threshold: continue

      i_is_shorter_fragment = \
         (self.segment_lengths[i] < self.segment_lengths[j])
      i_is_side_chain = not (self.is_main_chain_plus_cb(self.get_atomname(i)))
      j_is_side_chain = not (self.is_main_chain_plus_cb(self.get_atomname(j)))
      key_i=self.get_key(i)
      key_j=self.get_key(j)

      if (i_is_side_chain or j_is_side_chain):
        if not i_is_side_chain:  # take j
          if not key_j in remove_side_chain_list:remove_side_chain_list.append(
              key_j)
        elif not j_is_side_chain: #  take i
          if not key_i in remove_side_chain_list:remove_side_chain_list.append(
              key_i)
        elif i_is_shorter_fragment: # take i as shorter side chain
          if not key_i in remove_side_chain_list:remove_side_chain_list.append(
              key_i)
        else:
          if not key_j in remove_side_chain_list:remove_side_chain_list.append(
              key_j)
      else:  # neither is a side chain. Take shorter fragment
        if i_is_shorter_fragment: # take i whole residue
          if not key_i in remove_residue_list:remove_residue_list.append(key_i)
        else:
          if not key_j in remove_residue_list:remove_residue_list.append(key_j)

    # Take anything that is in remove_residue_list out of remove_side_chain_list
    new_remove_side_chain_list=[]
    for x in remove_side_chain_list:
      if not x in remove_residue_list:
        new_remove_side_chain_list.append(x)
    remove_side_chain_list=new_remove_side_chain_list
    self.side_chains_removed=len(remove_side_chain_list)
    self.residues_removed=len(remove_residue_list)

    remove_residue_selection_list=[]
    for x in remove_residue_list:
      chain_id,resseq=x.split(":")
      remove_residue_selection_list.append("(chain %s and resseq %s)" % (
         chain_id,resseq))
    if remove_residue_selection_list:
      remove_residue_selection_string=" ( %s ) " %( " OR ".join(
          remove_residue_selection_list))
    else:
      remove_residue_selection_string=""

    remove_side_chain_selection_list=[]
    for x in remove_side_chain_list:
      chain_id,resseq=x.split(":")
      remove_side_chain_selection_list.append(
        "chain %s and resseq %s" % (chain_id,resseq))
    if remove_side_chain_selection_list:
      side_chain_string=" ( NOT ( NAME %s ) )" % (
         " OR NAME ".join(main_chain_plus_cb))
      remove_side_chain_selection_string=" (%s AND ( %s  ) ) " %(
         side_chain_string, " OR ".join(remove_side_chain_selection_list))
    else:
      remove_side_chain_selection_string=""

    if remove_residue_selection_string and remove_side_chain_selection_string:
      self.remove_selection_string="%s OR %s" %(
        remove_residue_selection_string,remove_side_chain_selection_string)
    elif remove_residue_selection_string:
       self.remove_selection_string=remove_residue_selection_string
    else:
       self.remove_selection_string=remove_side_chain_selection_string

  def get_atomname(self,i):
    return self.model.get_atoms()[i].name

  def get_key(self,i):
      atom=self.model.get_atoms()[i]
      ag=atom.parent()
      rg=ag.parent()
      resseq_i=rg.resseq
      chain_id_i=rg.parent().id
      key_i="%s:%s" %(chain_id_i,resseq_i)
      return key_i

  def is_main_chain_plus_cb(self, atomname):
    if atomname.strip().upper() in main_chain_plus_cb:
      return True
    else:
      return False

  def get_bad_nonbonded(self,max_items=None):

    sites_cart=self.model.get_sites_cart()
    bad_nonbonded = []
    sorted_table, n_not_shown = self.pair_proxies.nonbonded_proxies.get_sorted(
        by_value="delta",
        sites_cart=sites_cart,
        site_labels=self.site_labels,
        max_items=max_items)
    for info in sorted_table:
      labels, i_seq, j_seq, delta, vdw_distance, sym_op_j, rt_mx = info
      if delta - vdw_distance  <= -self.non_bonded_deviation_threshold:
        bad_nonbonded.append([i_seq, j_seq, delta - vdw_distance])
    return bad_nonbonded

  def get_bad_bonded(self,max_items=None):

    # guess max items as percentage of sites_cart
    sites_cart=self.model.get_sites_cart()
    bad_bonded = []
    sorted_table, n_not_shown = self.pair_proxies.bond_proxies.get_sorted(
        by_value="residual",
        sites_cart=sites_cart,
        site_labels=self.site_labels,
        max_items=max_items,
        origin_id=origin_ids.get_origin_id('covalent geometry'))
    if sorted_table is not None:
      for restraint_info in sorted_table:
        (i_seq,j_seq,
            labels, distance_ideal, distance_model, slack, delta, sigma, weight,
            residual, sym_op_j, rt_mx) = restraint_info
        bad_bonded.append([i_seq,j_seq,delta])
    return bad_bonded


  def get_labels_and_proxies(self):
    # These are not really necessary
    self.site_labels = self.model.get_site_labels()

    self.pair_proxies = \
      self.model.get_restraints_manager().geometry.pair_proxies(
        sites_cart=self.model.get_sites_cart(),
        site_labels=self.site_labels)

    # Identify length of segments so we can prioritize on shorter segment
    self.get_segment_lengths()

  def get_segment_lengths(self):
    # Length of segment (residues) that each atom is in
    last_ca=None
    count_residues=0
    count_atoms=0
    self.segment_lengths=[]
    for atom in self.model.get_hierarchy().atoms():
      count_atoms+=1
      if atom.name.strip().upper() in ['CA','P']:
        count_residues+=1
        if last_ca is None or atom.distance(last_ca)< 4.5: #
           last_ca=atom
        else:
           self.segment_lengths+=count_atoms*[count_residues]  # count
           last_ca=None
           count_residues=0
           count_atoms=0
    if count_atoms:
      self.segment_lengths+=count_atoms*[count_residues]  # count
    assert len(self.segment_lengths) == self.model.get_number_of_atoms()
