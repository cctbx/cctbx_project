from __future__ import division
from __future__ import print_function
from builtins import object
import mmtbx.monomer_library
import mmtbx.monomer_library.pdb_interpretation
import sys

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

  def __init__(self, pdb_str, clash_threshold, remove_clash_guard=False,
      remove_max_reasonable_bond_distance=False):
    params = mmtbx.monomer_library.pdb_interpretation.master_params.extract()
    if remove_clash_guard:
      params.clash_guard.nonbonded_distance_threshold=None
      params.proceed_with_excessive_length_bonds=True
    if remove_max_reasonable_bond_distance:
      params.max_reasonable_bond_distance=None

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
