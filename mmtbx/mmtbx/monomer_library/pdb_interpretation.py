from __future__ import generators
from iotbx import pdb
import iotbx.pdb.interpretation
from mmtbx.monomer_library import server
from mmtbx.monomer_library import cif_types
from cctbx import geometry_restraints
import cctbx.geometry_restraints.manager
from cctbx import crystal
import cctbx.crystal.coordination_sequences
from cctbx.array_family import flex
from scitbx.python_utils import dicts
from scitbx.python_utils.misc import user_plus_sys_time, plural_s
from libtbx.utils import flat_list, Sorry
from cStringIO import StringIO
import string
import sys

def flush_log(log):
  if (log is not None):
    flush = getattr(log, "flush", None)
    if (flush is not None): flush()

class counters(object):

  def __init__(self, label):
    self.label = label
    self.corrupt_monomer_library_definitions = 0
    self.already_assigned_to_first_conformer = 0
    self.unresolved_non_hydrogen = 0
    self.unresolved_hydrogen = 0
    self.undefined = 0
    self.resolved = 0
    self.discarded_because_of_special_positions = 0

def involves_special_positions(special_position_indices, i_seqs):
  if (special_position_indices is None): return False
  for i_seq in i_seqs:
    if (i_seq in special_position_indices):
      return True
  return False

class source_info_server(object):

  def __init__(self, m_i, m_j):
    self.m_i, self.m_j = m_i, m_j

  def labels(self):
    if (self.m_j is None):
      return "residue: %s" % self.m_i.residue_conformer_label()
    return "residues: %s + %s" % (
      self.m_i.residue_conformer_label(),
      self.m_j.residue_conformer_label())

  def n_expected_atoms(self):
    if (self.m_j is None):
      return len(self.m_i.expected_atom_i_seqs)
    return len(self.m_i.expected_atom_i_seqs) \
         + len(self.m_j.expected_atom_i_seqs)

def format_exception_message(
      m_i,
      m_j,
      i_seqs,
      base_message,
      source_labels=None,
      show_residue_names=True,
      lines=[]):
  s = StringIO()
  print >> s, base_message
  for line in lines:
    print >> s, " ", line
  if (source_labels is not None):
    for i,label in enumerate(source_labels):
      print >> s, "  %d. definition from: %s" % (i+1, label)
  if (show_residue_names):
    print >> s, "  " + source_info_server(m_i, m_j).labels()
  print >> s, "  atom%s:" % plural_s(len(i_seqs))[1]
  stage_1 = m_i.pdb_residue.chain.conformer.model.stage_1
  stage_1.show_atom_labels(i_seqs=i_seqs, f=s, prefix="    ", max_lines=10)
  return s.getvalue()[:-1]

def discard_conflicting_pdb_element_column(mm):
  trusted_library_definitions = [
    "GLY", "VAL", "ALA", "LEU", "ILE", "PRO", "MET", "PHE", "TRP", "SER",
    "THR", "TYR", "CYS", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS",
    "MSE"]
  if (mm.residue_name in trusted_library_definitions): return True
  if (mm.residue_name[3:] == "%COO"
      and mm.residue_name[:3] in trusted_library_definitions):
    return True
  return False

class type_symbol_registry_base(object):

  def __init__(self, type_label, symbols, strict_conflict_handling):
    assert type_label in ["scattering", "nonbonded energy"]
    self.type_label = type_label
    self.symbols = symbols
    self.strict_conflict_handling = strict_conflict_handling
    self.n_resolved_conflicts = 0
    self.source_labels = flex.std_string(symbols.size())
    self.source_n_expected_atoms = flex.int(symbols.size(), -1)

  def discard_tables(self):
    self.source_labels = None
    self.source_n_expected_atoms = None

  def assign_from_monomer_mapping(self, conformer_label, mm):
    atom_dict = mm.monomer_atom_dict
    for atom_id,i_seq in mm.expected_atom_i_seqs.items():
      if (self.type_label == "scattering"):
        symbol = atom_dict[atom_id].type_symbol
        if (symbol == "H" and self.symbols[i_seq] == "D"):
          symbol = "D"
      else:
        symbol = atom_dict[atom_id].type_energy
      if (symbol is None): continue
      source_label = mm.residue_conformer_label()
      source_n_expected_atoms = len(mm.expected_atom_i_seqs)
      prev_symbol = self.symbols[i_seq]
      prev_source_label = self.source_labels[i_seq]
      prev_source_n_expected_atoms = self.source_n_expected_atoms[i_seq]
      assign = False
      raise_conflict = False
      if (prev_symbol == ""):
        assign = True
      elif (prev_symbol != symbol):
        if (self.strict_conflict_handling):
          raise_conflict = True
        elif (self.type_label == "nonbonded energy"):
          if (prev_source_n_expected_atoms == source_n_expected_atoms):
            raise_conflict = True
          else:
            self.n_resolved_conflicts += 1
            if (prev_source_n_expected_atoms < source_n_expected_atoms):
              assign = True
        elif (self.type_label == "scattering"):
          if (prev_source_label == ""
              and discard_conflicting_pdb_element_column(mm)):
            self.n_resolved_conflicts += 1
            assign = True
          else:
            raise_conflict = True
      assert not (assign and raise_conflict)
      if (assign):
        self.symbols[i_seq] = symbol
        self.source_labels[i_seq] = source_label
        self.source_n_expected_atoms[i_seq] = source_n_expected_atoms
      elif (raise_conflict):
        source = "with residue name %s" % source_label
        if (prev_source_label == ""):
          assert self.type_label == "scattering"
          prev_source = "from pdb element column"
        else:
          if (prev_source_label == source_label):
            source = "also " + source
          prev_source = "with residue name %s" % prev_source_label
        raise Sorry(format_exception_message(
          m_i=mm,
          m_j=None,
          i_seqs=[i_seq],
          base_message="Conflicting %s type symbols:" % self.type_label,
          show_residue_names=False,
          lines=['initial symbol: "%s" (%s)' % (prev_symbol, prev_source),
                 '    new symbol: "%s" (%s)' % (symbol, source)]))

  def n_unknown_type_symbols(self):
    return self.symbols.count("")

  def report(self, stage_1, log, prefix, max_lines=10):
    n_unknown = self.n_unknown_type_symbols()
    if (n_unknown > 0):
      print >> log, "%sNumber of atoms with unknown %s type symbols: %d" % (
        prefix, self.type_label, n_unknown)
      i_seqs = (self.symbols == "").iselection()
      stage_1.show_atom_labels(
        i_seqs=i_seqs, f=log, prefix=prefix+"  ", max_lines=max_lines)
    if (self.n_resolved_conflicts > 0):
      print >> log, "%sNumber of resolved %s type symbol conflicts: %d" % (
        prefix, self.type_label, self.n_resolved_conflicts)

class scattering_type_registry(type_symbol_registry_base):

  def __init__(self, scattering_types, strict_conflict_handling):
    type_symbol_registry_base.__init__(self,
      type_label="scattering",
      symbols=scattering_types,
      strict_conflict_handling=strict_conflict_handling)

class nonbonded_energy_type_registry(type_symbol_registry_base):

  def __init__(self, n_seq, strict_conflict_handling):
    type_symbol_registry_base.__init__(self,
      type_label="nonbonded energy",
      symbols=flex.std_string(n_seq),
      strict_conflict_handling=strict_conflict_handling)

class monomer_mapping_summary(object):

  __slots__ = [
    "conformer_label",
    "residue_name",
    "expected_atom_i_seqs",
    "unexpected_atom_i_seqs",
    "duplicate_atom_i_seqs",
    "ignored_atom_i_seqs",
    "classification",
    "incomplete_info",
    "is_terminus",
    "is_unusual"]

  def __init__(self, **keyword_args):
    for key in monomer_mapping_summary.__slots__:
      setattr(self, key, keyword_args.get(key, None))

  def all_associated_i_seqs(self):
    return flex.size_t(
        self.expected_atom_i_seqs
      + self.unexpected_atom_i_seqs
      + self.duplicate_atom_i_seqs
      + self.ignored_atom_i_seqs)

class monomer_mapping(object):

  def __init__(self,
        mon_lib_srv,
        i_conformer,
        conformer_label,
        pdb_residue):
    self.i_conformer = i_conformer
    self.conformer_label = conformer_label
    self.pdb_residue = pdb_residue
    stage_1 = self.pdb_residue.chain.conformer.model.stage_1
    atom = stage_1.atom_attributes_list[self.pdb_residue.iselection[0]]
    self.residue_name = atom.resName
    self.monomer = mon_lib_srv.get_comp_comp_id(self.residue_name)
    if (self.monomer is not None):
      self._get_mappings(mon_lib_srv=mon_lib_srv)
      if (self.monomer.chem_comp.group == "rna_dna_placeholder"):
        if (self.missing_non_hydrogen_atoms.has_key("O2*")):
          second_character = "d"
        else:
          second_character = "r"
        self.monomer = mon_lib_srv.get_comp_comp_id(
          self.monomer.chem_comp.id[0]+second_character)
        self._get_mappings(mon_lib_srv=mon_lib_srv)
      self._set_incomplete_info()
      self.chem_mod_ids = []
      self.is_terminus = None
      self.monomer.set_classification()
      if (self.incomplete_info is None):
        self._resolve_unexpected(mon_lib_srv=mon_lib_srv)

  def _get_mappings(self, mon_lib_srv):
    self.monomer_atom_dict = atom_dict = self.monomer.atom_dict()
    processed_atom_names = {}
    self.expected_atom_i_seqs = {}
    self.unexpected_atom_i_seqs = {}
    self.ignored_atom_i_seqs = {}
    self.duplicate_atom_i_seqs = {}
    is_rna_dna = (self.monomer.chem_comp.group == "rna_dna_placeholder"
                  or self.monomer.is_rna_dna())
    stage_1 = self.pdb_residue.chain.conformer.model.stage_1
    for i_seq in self.pdb_residue.iselection:
      atom = stage_1.atom_attributes_list[i_seq]
      if (atom.element.strip() == "Q"):
        self.ignored_atom_i_seqs.setdefault(atom.name, []).append(i_seq)
        continue
      atom_name_given = atom.name.replace(" ","")
      atom_name = atom_name_given
      if (atom_name[0] in string.digits):
        atom_name = atom_name[1:] + atom_name[0] # AD HOC manipulation
      if (is_rna_dna):
        atom_name = atom_name.replace("'", "*")
      if (not atom_dict.has_key(atom_name)):
        atom_name = mon_lib_srv.comp_synonym_atom_list_dict.get(
          self.monomer.chem_comp.id, {}).get(atom_name, atom_name)
      i_seq_prev = processed_atom_names.get(atom_name, None)
      if (i_seq_prev is None):
        processed_atom_names[atom_name] = i_seq
        if (atom_dict.has_key(atom_name)):
          self.expected_atom_i_seqs[atom_name] = i_seq
        else:
          self.unexpected_atom_i_seqs[atom_name] = i_seq
      else:
        self.duplicate_atom_i_seqs.setdefault(atom_name, []).append(i_seq)
    if (self.monomer.is_peptide()):
      self._rename_ot1_ot2("OXT" in atom_dict)
    self._set_missing_atoms()

  def _rename_ot1_ot2(self, oxt_in_atom_dict):
    if (not self.expected_atom_i_seqs.has_key("O")):
      i_seq = self.unexpected_atom_i_seqs.get("OT1", None)
      if (i_seq is not None):
        self.expected_atom_i_seqs["O"] = i_seq
        del self.unexpected_atom_i_seqs["OT1"]
    if (oxt_in_atom_dict):
      oxt_dict = self.expected_atom_i_seqs
    else:
      oxt_dict = self.unexpected_atom_i_seqs
    if (not oxt_dict.has_key("OXT")):
      i_seq = self.unexpected_atom_i_seqs.get("OT2", None)
      if (i_seq is not None):
        oxt_dict["OXT"] = i_seq
        del self.unexpected_atom_i_seqs["OT2"]

  def _set_missing_atoms(self):
    self.missing_non_hydrogen_atoms = {}
    self.missing_hydrogen_atoms = {}
    for atom in self.monomer.atom_list:
      if (not self.expected_atom_i_seqs.has_key(atom.atom_id)):
        if (atom.type_symbol != "H"):
          self.missing_non_hydrogen_atoms[atom.atom_id] = atom
        else:
          self.missing_hydrogen_atoms[atom.atom_id] = atom

  def _get_incomplete_info(self):
    if (    len(self.unexpected_atom_i_seqs) == 0
        and len(self.missing_non_hydrogen_atoms) > 0):
      if (self.monomer.is_peptide()):
        atom_ids = self.expected_atom_i_seqs.keys()
        atom_ids.sort()
        atom_ids = " ".join(atom_ids)
        if (atom_ids == "CA"): return "c_alpha_only"
        if (atom_ids == "C CA N"): return "n_c_alpha_c_only"
        if (atom_ids == "C CA N O"): return "backbone_only"
        if (atom_ids == "C CA CB N O"): return "truncation_to_alanine"
      elif (self.monomer.is_rna_dna()):
        atom_ids = " ".join(self.expected_atom_i_seqs.keys())
        if (atom_ids == "P"): return "p_only"
    return None

  def _set_incomplete_info(self):
    self.incomplete_info = self._get_incomplete_info()

  def _resolve_unexpected(self, mon_lib_srv):
    mod_mod_ids = []
    if (self.monomer.classification == "peptide"):
      if ("OXT" in self.unexpected_atom_i_seqs):
        mod_mod_ids.append(mon_lib_srv.mod_mod_id_dict["COO"])
      if (      "H1" in self.unexpected_atom_i_seqs
            and "H2" in self.unexpected_atom_i_seqs
            and "H3" in self.unexpected_atom_i_seqs):
        mod_mod_ids.append(mon_lib_srv.mod_mod_id_dict["NH3"])
      elif (    "HN1" in self.unexpected_atom_i_seqs
            and "HN2" in self.unexpected_atom_i_seqs):
        mod_mod_ids.append(mon_lib_srv.mod_mod_id_dict["NH2"])
      elif (    "HN" in self.unexpected_atom_i_seqs):
        mod_mod_ids.append(mon_lib_srv.mod_mod_id_dict["NH1"])
      elif (    "H1" in self.unexpected_atom_i_seqs
            or  "H2" in self.unexpected_atom_i_seqs
            or  "H3" in self.unexpected_atom_i_seqs):
        mod_mod_ids.append(mon_lib_srv.mod_mod_id_dict["NH3"])
      elif (    "HN1" in self.unexpected_atom_i_seqs
            or  "HN2" in self.unexpected_atom_i_seqs):
        mod_mod_ids.append(mon_lib_srv.mod_mod_id_dict["NH2"])
    elif (self.monomer.classification in ["RNA", "DNA"]):
      if (    not "P"   in self.expected_atom_i_seqs
          and not "OP1" in self.expected_atom_i_seqs
          and not "OP2" in self.expected_atom_i_seqs):
        mod_mod_ids.append(mon_lib_srv.mod_mod_id_dict["5*END"])
      if ("HO3*" in self.unexpected_atom_i_seqs):
        mod_mod_ids.append(mon_lib_srv.mod_mod_id_dict["3*END"])
    for mod_mod_id in mod_mod_ids:
      try:
        mod_mon = self.monomer.apply_mod(mod_mod_id)
      except RuntimeError:
        return
      self.chem_mod_ids.append(mod_mod_id.chem_mod.id)
      self.residue_name += "%" + mod_mod_id.chem_mod.id
      mod_mon.classification = self.monomer.classification
      self.monomer = mod_mon
      if (mod_mod_id.chem_mod.name.lower().find("terminus") >= 0):
        self.is_terminus = True # AD HOC manipulation
      self._get_mappings(mon_lib_srv=mon_lib_srv)

  def residue_conformer_label(self):
    result = self.residue_name
    if (self.conformer_label != " "):
      result += ", conformer %s" % self.conformer_label
    return result

  def is_unusual(self):
    if (self.monomer.is_peptide()): return False
    if (self.monomer.is_rna_dna()): return False
    if (self.monomer.is_water()): return False
    return True

  def summary(self):
    return monomer_mapping_summary(
      conformer_label=self.conformer_label,
      residue_name=self.residue_name,
      expected_atom_i_seqs=self.expected_atom_i_seqs.values(),
      unexpected_atom_i_seqs=self.unexpected_atom_i_seqs.values(),
      duplicate_atom_i_seqs=flat_list(self.duplicate_atom_i_seqs.values()),
      ignored_atom_i_seqs=self.ignored_atom_i_seqs.values(),
      classification=self.monomer.classification,
      incomplete_info=self.incomplete_info,
      is_terminus=self.is_terminus,
      is_unusual=self.is_unusual())

  def add_bond_proxies(self, bond_simple_proxy_registry):
    self.bond_counters = add_bond_proxies(
      counters=counters(label="bond"),
      m_i=self,
      m_j=self,
      bond_list=self.monomer.bond_list,
      bond_simple_proxy_registry=bond_simple_proxy_registry).counters

  def add_angle_proxies(self, special_position_indices, angle_proxy_registry):
    self.angle_counters = add_angle_proxies(
      counters=counters(label="angle"),
      m_i=self,
      m_j=None,
      angle_list=self.monomer.angle_list,
      angle_proxy_registry=angle_proxy_registry,
      special_position_indices=special_position_indices).counters

  def add_dihedral_proxies(self, special_position_indices,
                                 dihedral_proxy_registry):
    self.dihedral_counters = add_dihedral_proxies(
      counters=counters(label="dihedral"),
      m_i=self,
      m_j=None,
      tor_list=self.monomer.tor_list,
      dihedral_proxy_registry=dihedral_proxy_registry,
      special_position_indices=special_position_indices).counters

  def add_chirality_proxies(self, special_position_indices,
                                  chirality_proxy_registry):
    self.chirality_counters = add_chirality_proxies(
      counters=counters(label="chirality"),
      m_i=self,
      m_j=None,
      chir_list=self.monomer.chir_list,
      chirality_proxy_registry=chirality_proxy_registry,
      special_position_indices=special_position_indices).counters

  def add_planarity_proxies(self, special_position_indices,
                                  planarity_proxy_registry):
    self.planarity_counters = add_planarity_proxies(
      counters=counters(label="planarity"),
      m_i=self,
      m_j=None,
      plane_list=self.monomer.get_planes(),
      planarity_proxy_registry=planarity_proxy_registry,
      special_position_indices=special_position_indices).counters

class link_match_one(object):

  def __init__(self, chem_link_comp_id, chem_link_group_comp,
                     comp_id, comp_group):
    if (comp_group in cif_types.peptide_comp_groups):
      comp_group = "peptide"
    elif (comp_group in cif_types.dna_rna_comp_groups):
      comp_group = "DNA/RNA"
    if (   chem_link_comp_id == ""
        or comp_id.lower() == chem_link_comp_id.lower()):
      self.is_comp_id_match = True
      self.len_comp_id_match = len(chem_link_comp_id)
    else:
      self.is_comp_id_match = False
      self.len_comp_id_match = 0
    if (   chem_link_group_comp == ""
        or comp_group.lower() == chem_link_group_comp.lower()):
      self.is_group_match = True
      self.len_group_match = len(chem_link_group_comp)
    else:
      self.is_group_match = False
      self.len_group_match = 0

  def is_match(self):
    return self.is_comp_id_match and self.is_group_match

class link_match(object):

  def __init__(self, link_link_id, comp_id_1, comp_group_1,
                                   comp_id_2, comp_group_2):
    self.link_link_id = None
    chem_link = link_link_id.chem_link
    match_1 = link_match_one(
      chem_link.comp_id_1, chem_link.group_comp_1,
      comp_id_1, comp_group_1)
    if (not match_1.is_match()): return
    match_2 = link_match_one(
      chem_link.comp_id_2, chem_link.group_comp_2,
      comp_id_2, comp_group_2)
    if (not match_2.is_match()): return
    self.link_link_id = link_link_id
    self.len_comp_id_match_1 = match_1.len_comp_id_match
    self.len_group_match_1 = match_1.len_group_match
    self.len_comp_id_match_2 = match_2.len_comp_id_match
    self.len_group_match_2 = match_2.len_group_match

  def __cmp__(self, other):
    if (self.n_unresolved_bonds < other.n_unresolved_bonds): return -1
    if (self.n_unresolved_bonds > other.n_unresolved_bonds): return  1
    if (self.n_unresolved_angles < other.n_unresolved_angles): return -1
    if (self.n_unresolved_angles > other.n_unresolved_angles): return  1
    if (self.len_comp_id_match_1 > other.len_comp_id_match_1): return -1
    if (self.len_comp_id_match_1 < other.len_comp_id_match_1): return  1
    if (self.len_comp_id_match_2 > other.len_comp_id_match_2): return -1
    if (self.len_comp_id_match_2 < other.len_comp_id_match_2): return  1
    if (self.len_group_match_1 > other.len_group_match_1): return -1
    if (self.len_group_match_1 < other.len_group_match_1): return  1
    if (self.len_group_match_2 > other.len_group_match_2): return -1
    if (self.len_group_match_2 < other.len_group_match_2): return  1
    return 0

def get_lib_link_peptide(mon_lib_srv, m_i, m_j):
  link_id = "TRANS"
  if (m_j.expected_atom_i_seqs.get("CN", None) is not None):
    link_id = "NM" + link_id
  elif (m_j.monomer.chem_comp.id == "PRO"):
    link_id = "P" + link_id
  return mon_lib_srv.link_link_id_dict[link_id]

def get_lib_link(mon_lib_srv, m_i, m_j):
  if (m_i.monomer.is_peptide() and m_j.monomer.is_peptide()):
    return get_lib_link_peptide(mon_lib_srv, m_i, m_j)
  elif (m_i.monomer.is_rna_dna() and m_j.monomer.is_rna_dna()):
    return mon_lib_srv.link_link_id_dict["p"]
  if (m_i.monomer.is_water() or m_j.monomer.is_water()): return None
  comp_id_1 = m_i.monomer.chem_comp.id
  comp_id_2 = m_j.monomer.chem_comp.id
  comp_1 = mon_lib_srv.get_comp_comp_id(comp_id_1)
  comp_2 = mon_lib_srv.get_comp_comp_id(comp_id_2)
  group_1 = comp_1.chem_comp.group
  group_2 = comp_2.chem_comp.group
  atom_dicts = [None, m_i.monomer_atom_dict, m_j.monomer_atom_dict]
  matches = []
  for link_link_id in mon_lib_srv.link_link_id_list:
    chem_link = link_link_id.chem_link
    if (chem_link.name in cif_types.non_chain_links): continue
    if (    chem_link.comp_id_1 == ""
        and chem_link.mod_id_1 == ""
        and chem_link.group_comp_1 == ""
        and chem_link.comp_id_2 == ""
        and chem_link.mod_id_2 == ""
        and chem_link.group_comp_2 == ""): continue
    match = link_match(link_link_id, comp_id_1, group_1, comp_id_2, group_2)
    if (match.link_link_id is not None):
      match.n_unresolved_bonds = 0
      for bond in match.link_link_id.bond_list:
        atoms = [
          atom_dicts[bond.atom_1_comp_id].get(bond.atom_id_1, None),
          atom_dicts[bond.atom_2_comp_id].get(bond.atom_id_2, None)]
        if (None in atoms):
          match.n_unresolved_bonds += 1
      match.n_unresolved_angles = 0
      for angle in match.link_link_id.angle_list:
        atoms = [
          atom_dicts[angle.atom_1_comp_id].get(angle.atom_id_1, None),
          atom_dicts[angle.atom_2_comp_id].get(angle.atom_id_2, None),
          atom_dicts[angle.atom_3_comp_id].get(angle.atom_id_3, None)]
        if (None in atoms):
          match.n_unresolved_angles += 1
      matches.append(match)
  if (len(matches) == 0): return None
  matches.sort()
  best_matches = []
  for m in matches:
    if (cmp(m, matches[0]) != 0): break
    best_matches.append(m)
  return best_matches[0].link_link_id

def evaluate_registry_process_result(
      proxy_label,
      m_i, m_j, i_seqs,
      registry_process_result,
      lines=[]):
  if (registry_process_result.is_conflicting):
    raise AssertionError(format_exception_message(
      m_i=m_i,
      m_j=m_j,
      i_seqs=i_seqs,
      base_message="Conflicting %s restraints:" % proxy_label,
      source_labels=registry_process_result.conflict_source_labels,
      show_residue_names=False,
      lines=lines))
  if (not registry_process_result.is_new
      and (   m_i.i_conformer == 0
           or not m_i.pdb_residue.chain.conformer.model
                    .stage_1.are_all_blank_altLocs(i_seqs))):
    raise AssertionError(format_exception_message(
      m_i=m_i,
      m_j=m_j,
      i_seqs=i_seqs,
      base_message="Duplicate %s restraints:" % proxy_label,
      lines=lines))

class add_bond_proxies(object):

  def __init__(self,
        counters,
        m_i,
        m_j,
        bond_list,
        bond_simple_proxy_registry,
        sites_cart=None,
        distance_cutoff=None):
    assert m_i.i_conformer == m_j.i_conformer
    self.counters = counters
    self.n_breaks = 0
    stage_1 = m_i.pdb_residue.chain.conformer.model.stage_1
    for bond in bond_list:
      if (   not m_i.monomer_atom_dict.has_key(bond.atom_id_1)
          or not m_j.monomer_atom_dict.has_key(bond.atom_id_2)):
        counters.corrupt_monomer_library_definitions += 1
        continue
      i_seqs = (m_i.expected_atom_i_seqs.get(bond.atom_id_1, None),
                m_j.expected_atom_i_seqs.get(bond.atom_id_2, None))
      if (None in i_seqs):
        if (   m_i.monomer_atom_dict[bond.atom_id_1].type_symbol == "H"
            or m_j.monomer_atom_dict[bond.atom_id_2].type_symbol == "H"):
          counters.unresolved_hydrogen += 1
        else:
          counters.unresolved_non_hydrogen += 1
      elif (   bond.value_dist is None
            or bond.value_dist_esd in [None, 0]):
        counters.undefined += 1
      else:
        counters.resolved += 1
        proxy = geometry_restraints.bond_simple_proxy(
          i_seqs=i_seqs,
          distance_ideal=bond.value_dist,
          weight=1/bond.value_dist_esd**2)
        is_large_distance = False
        if (sites_cart is not None):
          r = geometry_restraints.bond(sites_cart=sites_cart, proxy=proxy)
          if (r.distance_model > distance_cutoff):
            is_large_distance = True
            self.n_breaks += 1
        if (not is_large_distance):
          if (    m_i.i_conformer > 0
              and stage_1.are_all_blank_altLocs(i_seqs)):
            counters.already_assigned_to_first_conformer += 1
          else:
            registry_process_result = bond_simple_proxy_registry.process(
              source_info=source_info_server(m_i=m_i, m_j=m_j),
              proxy=proxy)
            evaluate_registry_process_result(
              proxy_label="bond_simple", m_i=m_i, m_j=m_j, i_seqs=i_seqs,
              registry_process_result=registry_process_result)

class add_angle_proxies(object):

  def __init__(self,
        counters,
        m_i,
        m_j,
        angle_list,
        angle_proxy_registry,
        special_position_indices):
    self.counters = counters
    if (m_j is None):
      m_1,m_2,m_3 = m_i,m_i,m_i
    else:
      assert m_i.i_conformer == m_j.i_conformer
    stage_1 = m_i.pdb_residue.chain.conformer.model.stage_1
    for angle in angle_list:
      if (m_j is not None):
        m_1,m_2,m_3 = [(m_i, m_j)[comp_id-1] for comp_id in (
          angle.atom_1_comp_id, angle.atom_2_comp_id, angle.atom_3_comp_id)]
      if (   not m_1.monomer_atom_dict.has_key(angle.atom_id_1)
          or not m_2.monomer_atom_dict.has_key(angle.atom_id_2)
          or not m_3.monomer_atom_dict.has_key(angle.atom_id_3)):
        counters.corrupt_monomer_library_definitions += 1
        continue
      i_seqs = (m_1.expected_atom_i_seqs.get(angle.atom_id_1, None),
                m_2.expected_atom_i_seqs.get(angle.atom_id_2, None),
                m_3.expected_atom_i_seqs.get(angle.atom_id_3, None))
      if (None in i_seqs):
        if (   m_1.monomer_atom_dict[angle.atom_id_1].type_symbol == "H"
            or m_2.monomer_atom_dict[angle.atom_id_2].type_symbol == "H"
            or m_3.monomer_atom_dict[angle.atom_id_3].type_symbol == "H"):
          counters.unresolved_hydrogen += 1
        else:
          counters.unresolved_non_hydrogen += 1
      elif (   angle.value_angle is None
            or angle.value_angle_esd in [None, 0]):
        counters.undefined += 1
      else:
        counters.resolved += 1
        if (involves_special_positions(special_position_indices, i_seqs)):
          counters.discarded_because_of_special_positions += 1
        else:
          registry_process_result = angle_proxy_registry.process(
            source_info=source_info_server(m_i=m_i, m_j=m_j),
            proxy=geometry_restraints.angle_proxy(
              i_seqs=i_seqs,
              angle_ideal=angle.value_angle,
              weight=1/angle.value_angle_esd**2))
          evaluate_registry_process_result(
            proxy_label="angle", m_i=m_i, m_j=m_j, i_seqs=i_seqs,
            registry_process_result=registry_process_result)

class add_dihedral_proxies(object):

  def __init__(self,
        counters,
        m_i,
        m_j,
        tor_list,
        dihedral_proxy_registry,
        special_position_indices,
        sites_cart=None,
        chem_link_id=None,
        discard_peptide_psi_phi=True,
        cis_threshold=45):
    self.counters = counters
    self.chem_link_id = chem_link_id
    if (chem_link_id not in ["TRANS", "PTRANS", "NMTRANS"]):
      sites_cart = None
    if (m_j is None):
      m_1,m_2,m_3,m_4 = m_i,m_i,m_i,m_i
    else:
      assert m_i.i_conformer == m_j.i_conformer
    stage_1 = m_i.pdb_residue.chain.conformer.model.stage_1
    for tor in tor_list:
      if (m_j is not None):
        m_1,m_2,m_3,m_4 = [(m_i, m_j)[comp_id-1] for comp_id in (
          tor.atom_1_comp_id,
          tor.atom_2_comp_id,
          tor.atom_3_comp_id,
          tor.atom_4_comp_id)]
      if (   not m_1.monomer_atom_dict.has_key(tor.atom_id_1)
          or not m_2.monomer_atom_dict.has_key(tor.atom_id_2)
          or not m_3.monomer_atom_dict.has_key(tor.atom_id_3)
          or not m_4.monomer_atom_dict.has_key(tor.atom_id_4)):
        counters.corrupt_monomer_library_definitions += 1
        continue
      i_seqs = (m_1.expected_atom_i_seqs.get(tor.atom_id_1, None),
                m_2.expected_atom_i_seqs.get(tor.atom_id_2, None),
                m_3.expected_atom_i_seqs.get(tor.atom_id_3, None),
                m_4.expected_atom_i_seqs.get(tor.atom_id_4, None))
      if (None in i_seqs):
        if (   m_1.monomer_atom_dict[tor.atom_id_1].type_symbol == "H"
            or m_2.monomer_atom_dict[tor.atom_id_2].type_symbol == "H"
            or m_3.monomer_atom_dict[tor.atom_id_3].type_symbol == "H"
            or m_4.monomer_atom_dict[tor.atom_id_4].type_symbol == "H"):
          counters.unresolved_hydrogen += 1
        else:
          counters.unresolved_non_hydrogen += 1
      elif (   tor.value_angle is None
            or tor.value_angle_esd in [None, 0]):
        counters.undefined += 1
      else:
        counters.resolved += 1
        if (involves_special_positions(special_position_indices, i_seqs)):
          counters.discarded_because_of_special_positions += 1
        elif (discard_peptide_psi_phi
              and tor.id in ["psi", "phi"]
              and self.chem_link_id in ["TRANS", "PTRANS", "NMTRANS",
                                        "CIS",   "PCIS",   "NMCIS"]):
          pass
        else:
          proxy = geometry_restraints.dihedral_proxy(
            i_seqs=i_seqs,
            angle_ideal=tor.value_angle,
            weight=1/tor.value_angle_esd**2,
            periodicity=max(1,tor.period))
          if (sites_cart is not None and tor.id == "omega"):
            assert abs(tor.value_angle - 180) < 1.e-6
            r = geometry_restraints.dihedral(
              sites_cart=sites_cart,
              proxy=proxy)
            if (abs(r.delta) > 180-cis_threshold):
              self.chem_link_id = self.chem_link_id.replace("TRANS", "CIS")
              proxy.angle_ideal = 0
          registry_process_result = dihedral_proxy_registry.process(
            source_info=source_info_server(m_i=m_i, m_j=m_j),
            proxy=proxy)
          evaluate_registry_process_result(
            proxy_label="dihedral", m_i=m_i, m_j=m_j, i_seqs=i_seqs,
            registry_process_result=registry_process_result,
            lines=["tor id: " + str(tor.id)])

class add_chirality_proxies(object):

  def __init__(self,
        counters,
        m_i,
        m_j,
        chir_list,
        chirality_proxy_registry,
        special_position_indices,
        lib_link=None,
        chir_volume_esd=0.2):
    self.counters = counters
    self.counters.unsupported_volume_sign = dicts.with_default_value(0)
    if (m_j is None):
      m_c,m_1,m_2,m_3 = m_i,m_i,m_i,m_i
    else:
      assert m_i.i_conformer == m_j.i_conformer
    stage_1 = m_i.pdb_residue.chain.conformer.model.stage_1
    for chir in chir_list:
      if (m_j is not None):
        m_c,m_1,m_2,m_3 = [(m_i, m_j)[comp_id-1] for comp_id in (
          chir.atom_centre_comp_id,
          chir.atom_1_comp_id,
          chir.atom_2_comp_id,
          chir.atom_3_comp_id)]
      volume_sign = chir.volume_sign
      if (volume_sign is not None):
        volume_sign = volume_sign[:4].lower()
      if (volume_sign not in ["posi", "nega", "both"]):
        counters.unsupported_volume_sign[volume_sign] += 1
        continue
      if (   not m_c.monomer_atom_dict.has_key(chir.atom_id_centre)
          or not m_1.monomer_atom_dict.has_key(chir.atom_id_1)
          or not m_2.monomer_atom_dict.has_key(chir.atom_id_2)
          or not m_3.monomer_atom_dict.has_key(chir.atom_id_3)):
        counters.corrupt_monomer_library_definitions += 1
        continue
      i_seqs = (m_c.expected_atom_i_seqs.get(chir.atom_id_centre, None),
                m_1.expected_atom_i_seqs.get(chir.atom_id_1, None),
                m_2.expected_atom_i_seqs.get(chir.atom_id_2, None),
                m_3.expected_atom_i_seqs.get(chir.atom_id_3, None))
      if (None in i_seqs):
        if (   m_c.monomer_atom_dict[chir.atom_id_centre].type_symbol == "H"
            or m_1.monomer_atom_dict[chir.atom_id_1].type_symbol == "H"
            or m_2.monomer_atom_dict[chir.atom_id_2].type_symbol == "H"
            or m_3.monomer_atom_dict[chir.atom_id_3].type_symbol == "H"):
          counters.unresolved_hydrogen += 1
        else:
          counters.unresolved_non_hydrogen += 1
      elif (   volume_sign is None
            or chir_volume_esd in [None, 0]):
        counters.undefined += 1
      else:
        if (m_j is None):
          volume_ideal = m_i.monomer.get_chir_volume_ideal(chir)
        else:
          volume_ideal = lib_link.get_chir_volume_ideal(
            m_i.monomer, m_j.monomer, chir)
        if (volume_ideal is None):
          counters.undefined += 1
        else:
          counters.resolved += 1
          if (involves_special_positions(special_position_indices, i_seqs)):
            counters.discarded_because_of_special_positions += 1
          else:
            registry_process_result = chirality_proxy_registry.process(
              source_info=source_info_server(m_i=m_i, m_j=m_j),
              proxy=geometry_restraints.chirality_proxy(
                i_seqs=i_seqs,
                volume_ideal=volume_ideal,
                both_signs=(volume_sign == "both"),
                weight=1/chir_volume_esd**2))
            evaluate_registry_process_result(
              proxy_label="chirality", m_i=m_i, m_j=m_j, i_seqs=i_seqs,
              registry_process_result=registry_process_result)

class add_planarity_proxies(object):

  def __init__(self,
        counters,
        m_i,
        m_j,
        plane_list,
        planarity_proxy_registry,
        special_position_indices):
    self.counters = counters
    self.counters.less_than_four_sites = dicts.with_default_value(0)
    if (m_j is not None):
      assert m_i.i_conformer == m_j.i_conformer
    stage_1 = m_i.pdb_residue.chain.conformer.model.stage_1
    for plane in plane_list:
      this_plane_has_unresolved_non_hydrogen = False
      i_seqs = []
      weights = []
      for plane_atom in plane.plane_atoms:
        if (m_j is None):
          m_x = m_i
        else:
          assert plane_atom.atom_comp_id in (1,2)
          m_x = (m_i, m_j)[plane_atom.atom_comp_id-1]
        if (not m_x.monomer_atom_dict.has_key(plane_atom.atom_id)):
          counters.corrupt_monomer_library_definitions += 1
          continue
        i_seq = m_x.expected_atom_i_seqs.get(plane_atom.atom_id, None)
        if (i_seq is None):
          if (m_x.monomer_atom_dict[plane_atom.atom_id].type_symbol == "H"):
            counters.unresolved_hydrogen += 1
          else:
            counters.unresolved_non_hydrogen += 1
            this_plane_has_unresolved_non_hydrogen = True
        elif (plane_atom.dist_esd in [None, 0]):
          counters.undefined += 1
        else:
          counters.resolved += 1
          if (special_position_indices is not None
              and i_seq in special_position_indices):
            counters.discarded_because_of_special_positions += 1
          else:
            i_seqs.append(i_seq)
            weights.append(1./plane_atom.dist_esd**2)
      if (len(i_seqs) < 4):
        if (this_plane_has_unresolved_non_hydrogen):
          counters.less_than_four_sites[plane.plane_id] += 1
      else:
        registry_process_result = planarity_proxy_registry.process(
          source_info=source_info_server(m_i=m_i, m_j=m_j),
          proxy=geometry_restraints.planarity_proxy(
            i_seqs=flex.size_t(i_seqs),
            weights=flex.double(weights)))
        evaluate_registry_process_result(
          proxy_label="planarity", m_i=m_i, m_j=m_j, i_seqs=i_seqs,
          registry_process_result=registry_process_result,
          lines=["plane id: " + str(plane.plane_id)])

# XXX TODO synonymes
def ener_lib_as_nonbonded_params(
      ener_lib,
      factor_1_4_interactions,
      default_distance,
      minimum_distance):
  params = geometry_restraints.nonbonded_params(
    factor_1_4_interactions=factor_1_4_interactions,
    const_shrink_1_4_interactions=0,
    default_distance=default_distance,
    minimum_distance=minimum_distance)
  for vdw in ener_lib.lib_vdw:
    entry = params.distance_table.setdefault(
      vdw.atom_type_1)[vdw.atom_type_2] = vdw.radius_min
  for atom_type,energy_lib_atom in ener_lib.lib_atom.items():
    if (len(atom_type) == 0): continue
    if (energy_lib_atom.vdw_radius is not None):
      params.radius_table[atom_type] = energy_lib_atom.vdw_radius
  return params

def is_same_model_as_before(model_type_indices, i_model, models):
  stage_1 = models[0].stage_1
  selection_cache = stage_1.selection_cache()
  atoms = stage_1.atom_attributes_list
  sel_i = selection_cache.get_MODELserial(models[i_model].serial)[0]
  for j_model in xrange(0, i_model):
    if (model_type_indices[j_model] != j_model): continue
    sel_j = selection_cache.get_MODELserial(models[j_model].serial)[0]
    if (sel_j.size() == sel_i.size()):
      is_same = True
      for i,j in zip(sel_i, sel_j):
        if (not atoms[i].is_label_equivalent(atoms[j])):
          is_same = False
          break
      if (is_same):
        model_type_indices[i_model] = j_model
        return True
  model_type_indices[i_model] = i_model
  return False

class processed_model(pdb.interpretation.model):

  def __init__(self, stage_1, serial):
    pdb.interpretation.model.__init__(self,
      stage_1=stage_1, serial=serial)

  def monomer_mappings(self):
    for conformer in self.conformers:
      for mm in conformer.monomer_mappings():
        yield mm

  def monomer_mapping_summaries(self):
    for conformer in self.conformers:
      for mm in conformer.monomer_mapping_summaries():
        yield mm

class processed_conformer(pdb.interpretation.conformer_base):

  def __init__(self, model, altLoc, iselection):
    pdb.interpretation.conformer_base.__init__(self,
      model=model, altLoc=altLoc, iselection=iselection)
    self.chains = []

  def add_chain(self, chain):
    self.chains.append(chain)

  def monomer_mappings(self):
    for chain in self.chains:
      assert chain.monomer_mappings is not None
      for mm in chain.monomer_mappings:
        yield mm

  def monomer_mapping_summaries(self):
    for chain in self.chains:
      for mm in chain.monomer_mapping_summaries():
        yield mm

class processed_chain(object):

  def __init__(self,
        model_index,
        conformer_index,
        chainID,
        segID,
        keep_monomer_mappings):
    self.model_index = model_index
    self.conformer_index = conformer_index
    self.chainID = chainID
    self.segID = segID
    if (keep_monomer_mappings):
      self.monomer_mappings = []
      self.monomer_mapping_summaries_ = None
    else:
      self.monomer_mappings = None
      self.monomer_mapping_summaries_ = []

  def add_residue(self, monomer_mapping):
    if (self.monomer_mappings is not None):
      self.monomer_mappings.append(monomer_mapping)
    else:
      self.monomer_mapping_summaries_.append(monomer_mapping.summary())

  def monomer_mapping_summaries(self):
    if (self.monomer_mappings is not None):
      for mm in self.monomer_mappings:
        yield mm.summary()
    else:
      for mm in self.monomer_mapping_summaries_:
        yield mm

  def residue_names(self):
    result = []
    for mm in self.monomer_mapping_summaries():
      result.append(mm.residue_name.split("%")[0])
    return result

class build_chain_proxies(object):

  def __init__(self,
        mon_lib_srv,
        link_distance_cutoff,
        stage_1,
        sites_cart,
        special_position_indices,
        keep_monomer_mappings,
        scattering_type_registry,
        nonbonded_energy_type_registry,
        geometry_proxy_registries,
        cystein_sulphur_i_seqs,
        cystein_monomer_mappings,
        i_model,
        is_unique_model,
        i_conformer,
        conformer,
        chain,
        log):
    self.i_model = i_model
    self.i_conformer = i_conformer
    self.processed = processed_chain(
      model_index=i_model,
      conformer_index=i_conformer,
      chainID=chain.chainID,
      segID=chain.segID,
      keep_monomer_mappings=keep_monomer_mappings)
    unknown_residues = dicts.with_default_value(0)
    unusual_residues = dicts.with_default_value(0)
    inner_chain_residues_flagged_as_termini = []
    n_expected_atoms = 0
    unexpected_atoms = dicts.with_default_value(0)
    ignored_atoms = dicts.with_default_value(0)
    duplicate_atoms = dicts.with_default_value(0)
    classifications = dicts.with_default_value(0)
    modifications_used = dicts.with_default_value(0)
    incomplete_infos = dicts.with_default_value(0)
    link_ids = dicts.with_default_value(0)
    n_unresolved_chain_links = 0
    n_chain_breaks = 0
    n_unresolved_chain_link_angles = 0
    n_unresolved_chain_link_dihedrals = 0
    n_unresolved_chain_link_chiralities = 0
    n_unresolved_chain_link_planarities = 0
    corrupt_monomer_library_definitions = dicts.with_default_value(0)
    n_bond_proxies_already_assigned_to_first_conformer = 0
    n_unresolved_non_hydrogen_bonds = 0
    n_unresolved_non_hydrogen_angles = 0
    n_angles_discarded_because_of_special_positions = 0
    n_unresolved_non_hydrogen_dihedrals = 0
    n_dihedrals_discarded_because_of_special_positions = 0
    unsupported_chir_volume_sign = dicts.with_default_value(0)
    n_unresolved_non_hydrogen_chiralities = 0
    n_chiralities_discarded_because_of_special_positions = 0
    planarities_with_less_than_four_sites = dicts.with_default_value(0)
    n_unresolved_non_hydrogen_planarities = 0
    n_planarities_discarded_because_of_special_positions = 0
    break_block_identifiers = stage_1.get_break_block_identifiers()
    prev_break_block_identifier = None
    prev_mm = None
    for i_residue,residue in enumerate(chain.residues):
      break_block_identifier = break_block_identifiers[residue.iselection[0]]
      mm = monomer_mapping(
        mon_lib_srv=mon_lib_srv,
        i_conformer=i_conformer,
        conformer_label=conformer.altLoc,
        pdb_residue=residue)
      if (mm.monomer is None):
        unknown_residues[mm.residue_name] += 1
        n_chain_breaks += 1
      elif (break_block_identifier != prev_break_block_identifier
            and prev_break_block_identifier is not None):
        n_chain_breaks += 1
      else:
        if (mm.is_unusual()):
          unusual_residues[mm.residue_name] += 1
        if (    mm.is_terminus == True
            and i_residue > 0
            and i_residue < len(chain.residues)-1):
          inner_chain_residues_flagged_as_termini.append(
            stage_1.atom_attributes_list[residue.iselection[0]]
              .residue_labels())
        n_expected_atoms += len(mm.expected_atom_i_seqs)
        for atom_name in mm.unexpected_atom_i_seqs.keys():
          unexpected_atoms[mm.residue_name+","+atom_name] += 1
        for atom_name,i_seqs in mm.ignored_atom_i_seqs.items():
          ignored_atoms[mm.residue_name+","+atom_name] += len(i_seqs)
        for atom_name,i_seqs in mm.duplicate_atom_i_seqs.items():
          duplicate_atoms[mm.residue_name+","+atom_name] += len(i_seqs)
        if (mm.incomplete_info is not None):
          incomplete_infos[mm.incomplete_info] += 1
        if (mm.monomer.classification is not None):
          classifications[mm.monomer.classification] += 1
        for chem_mod_id in mm.chem_mod_ids:
          modifications_used[chem_mod_id] += 1
        if (prev_mm is not None and prev_mm.monomer is not None):
          prev_mm.lib_link = get_lib_link(
            mon_lib_srv=mon_lib_srv,
            m_i=prev_mm,
            m_j=mm)
          if (prev_mm.lib_link is None):
            link_ids[None] += 1
          else:
            link_resolution = add_bond_proxies(
              counters=counters(label="link_bond"),
              m_i=prev_mm,
              m_j=mm,
              bond_list=prev_mm.lib_link.bond_list,
              bond_simple_proxy_registry=geometry_proxy_registries.bond_simple,
              sites_cart=sites_cart,
              distance_cutoff=link_distance_cutoff)
            n_bond_proxies_already_assigned_to_first_conformer += \
              link_resolution.counters.already_assigned_to_first_conformer
            n_unresolved_chain_links \
              += link_resolution.counters.unresolved_non_hydrogen
            n_chain_breaks += link_resolution.n_breaks
            link_resolution = add_angle_proxies(
              counters=counters(label="link_angle"),
              m_i=prev_mm,
              m_j=mm,
              angle_list=prev_mm.lib_link.angle_list,
              angle_proxy_registry=geometry_proxy_registries.angle,
              special_position_indices=special_position_indices)
            n_unresolved_chain_link_angles \
              += link_resolution.counters.unresolved_non_hydrogen
            link_resolution = add_dihedral_proxies(
              counters=counters(label="link_dihedral"),
              m_i=prev_mm,
              m_j=mm,
              tor_list=prev_mm.lib_link.tor_list,
              dihedral_proxy_registry=geometry_proxy_registries.dihedral,
              special_position_indices=special_position_indices,
              sites_cart=sites_cart,
              chem_link_id=prev_mm.lib_link.chem_link.id)
            n_unresolved_chain_link_dihedrals \
              += link_resolution.counters.unresolved_non_hydrogen
            link_ids[link_resolution.chem_link_id] += 1
            link_resolution = add_chirality_proxies(
              counters=counters(label="link_chirality"),
              m_i=prev_mm,
              m_j=mm,
              chir_list=prev_mm.lib_link.chir_list,
              chirality_proxy_registry=geometry_proxy_registries.chirality,
              special_position_indices=special_position_indices,
              lib_link=prev_mm.lib_link)
            n_unresolved_chain_link_chiralities \
              += link_resolution.counters.unresolved_non_hydrogen
            link_resolution = add_planarity_proxies(
              counters=counters(label="link_planarity"),
              m_i=prev_mm,
              m_j=mm,
              plane_list=prev_mm.lib_link.get_planes(),
              planarity_proxy_registry=geometry_proxy_registries.planarity,
              special_position_indices=special_position_indices)
            n_unresolved_chain_link_planarities \
              += link_resolution.counters.unresolved_non_hydrogen
        scattering_type_registry.assign_from_monomer_mapping(
          conformer_label=conformer.altLoc, mm=mm)
        nonbonded_energy_type_registry.assign_from_monomer_mapping(
          conformer_label=conformer.altLoc, mm=mm)
        mm.add_bond_proxies(
          bond_simple_proxy_registry=geometry_proxy_registries.bond_simple)
        n_bond_proxies_already_assigned_to_first_conformer += \
          mm.bond_counters.already_assigned_to_first_conformer
        if (mm.bond_counters.corrupt_monomer_library_definitions > 0):
          corrupt_monomer_library_definitions[mm.residue_name] \
            += mm.bond_counters.corrupt_monomer_library_definitions
        n_unresolved_non_hydrogen_bonds \
          += mm.bond_counters.unresolved_non_hydrogen
        mm.add_angle_proxies(
          special_position_indices=special_position_indices,
          angle_proxy_registry=geometry_proxy_registries.angle)
        if (mm.angle_counters.corrupt_monomer_library_definitions > 0):
          corrupt_monomer_library_definitions[mm.residue_name] \
            += mm.angle_counters.corrupt_monomer_library_definitions
        n_unresolved_non_hydrogen_angles \
          += mm.angle_counters.unresolved_non_hydrogen
        n_angles_discarded_because_of_special_positions \
          += mm.angle_counters.discarded_because_of_special_positions
        mm.add_dihedral_proxies(
          special_position_indices=special_position_indices,
          dihedral_proxy_registry=geometry_proxy_registries.dihedral)
        if (mm.dihedral_counters.corrupt_monomer_library_definitions > 0):
          corrupt_monomer_library_definitions[mm.residue_name] \
            += mm.dihedral_counters.corrupt_monomer_library_definitions
        n_unresolved_non_hydrogen_dihedrals \
          += mm.dihedral_counters.unresolved_non_hydrogen
        n_dihedrals_discarded_because_of_special_positions \
          += mm.dihedral_counters.discarded_because_of_special_positions
        mm.add_chirality_proxies(
          special_position_indices=special_position_indices,
          chirality_proxy_registry=geometry_proxy_registries.chirality)
        if (mm.chirality_counters.corrupt_monomer_library_definitions > 0):
          corrupt_monomer_library_definitions[mm.residue_name] \
            += mm.chirality_counters.corrupt_monomer_library_definitions
        for s,n in mm.chirality_counters.unsupported_volume_sign.items():
          unsupported_chir_volume_sign[s] += n
        n_unresolved_non_hydrogen_chiralities \
          += mm.chirality_counters.unresolved_non_hydrogen
        n_chiralities_discarded_because_of_special_positions \
          += mm.chirality_counters.discarded_because_of_special_positions
        mm.add_planarity_proxies(
          special_position_indices=special_position_indices,
          planarity_proxy_registry=geometry_proxy_registries.planarity)
        if (mm.planarity_counters.corrupt_monomer_library_definitions > 0):
          corrupt_monomer_library_definitions[mm.residue_name] \
            += mm.planarity_counters.corrupt_monomer_library_definitions
        for p,n in mm.planarity_counters.less_than_four_sites.items():
          planarities_with_less_than_four_sites[mm.residue_name+":"+p] += n
        n_unresolved_non_hydrogen_planarities \
          += mm.planarity_counters.unresolved_non_hydrogen
        n_planarities_discarded_because_of_special_positions \
          += mm.planarity_counters.discarded_because_of_special_positions
        if (mm.monomer.chem_comp.id == "CYS"):
          sulphur_i_seq = mm.expected_atom_i_seqs.get("SG", None)
          # XXX keep track of weights
          if (sulphur_i_seq is not None
              and sulphur_i_seq not in cystein_sulphur_i_seqs):
            cystein_sulphur_i_seqs.append(sulphur_i_seq)
            cystein_monomer_mappings.append(mm)
        self.processed.add_residue(monomer_mapping=mm)
      prev_break_block_identifier = break_block_identifier
      prev_mm = mm
      prev_mm.lib_link = None
    if (is_unique_model and log is not None):
      print >> log, "        Number of residues, atoms: %d, %d" % (
        len(chain.residues),
        n_expected_atoms + flex.sum(flex.long(unexpected_atoms.values())))
      if (len(unknown_residues) > 0):
        print >> log, "          Unknown residues:", unknown_residues
      if (len(unusual_residues) > 0):
        print >> log, "          Unusual residues:", unusual_residues
      if (len(inner_chain_residues_flagged_as_termini) > 0):
        print >> log, "          Inner-chain residues flagged as termini:", \
          inner_chain_residues_flagged_as_termini
      if (len(unexpected_atoms) > 0):
        print >> log, "          Unexpected atoms:", unexpected_atoms
      if (len(ignored_atoms) > 0):
        print >> log, "          Ignored atoms:", ignored_atoms
      if (len(duplicate_atoms) > 0):
        print >> log, "          Duplicate atoms:", duplicate_atoms
      if (len(classifications) > 0):
        print >> log, "          Classifications:", classifications
      if (len(modifications_used) > 0):
        print >> log, "          Modifications used:", modifications_used
      if (len(incomplete_infos) > 0):
        print >> log, "          Incomplete info:", incomplete_infos
    if (log is not None):
      if (len(link_ids) > 0):
        print >> log, "          Link IDs:", link_ids
    if (is_unique_model and log is not None):
      if (n_unresolved_chain_links > 0):
        print >> log, "          Unresolved chain links:", \
          n_unresolved_chain_links
    if (log is not None):
      if (n_chain_breaks > 0):
        print >> log, "          Chain breaks:", n_chain_breaks
    if (is_unique_model and log is not None):
      if (n_unresolved_chain_link_angles > 0):
        print >> log, "          Unresolved chain link angles:", \
          n_unresolved_chain_link_angles
      if (n_unresolved_chain_link_dihedrals > 0):
        print >> log, "          Unresolved chain link dihedrals:", \
          n_unresolved_chain_link_dihedrals
      if (n_unresolved_chain_link_chiralities > 0):
        print >> log, "          Unresolved chain link chiralities:", \
          n_unresolved_chain_link_chiralities
      if (n_unresolved_chain_link_planarities > 0):
        print >> log, "          Unresolved chain link planarities:", \
          n_unresolved_chain_link_planarities
      if (len(corrupt_monomer_library_definitions) > 0):
        print >> log, "          Corrupt monomer library definitions:", \
          corrupt_monomer_library_definitions
      if (n_unresolved_non_hydrogen_bonds > 0):
        print >> log, "          Unresolved non-hydrogen bonds:", \
          n_unresolved_non_hydrogen_bonds
      if (n_unresolved_non_hydrogen_angles > 0):
        print >> log, "          Unresolved non-hydrogen angles:", \
          n_unresolved_non_hydrogen_angles
    if (log is not None):
      if (n_angles_discarded_because_of_special_positions > 0):
        print >> log, \
          "          Angles discarded because of special positions:", \
          n_angles_discarded_because_of_special_positions
    if (is_unique_model and log is not None):
      if (n_unresolved_non_hydrogen_dihedrals > 0):
        print >> log, "          Unresolved non-hydrogen dihedrals:", \
          n_unresolved_non_hydrogen_dihedrals
    if (log is not None):
      if (n_dihedrals_discarded_because_of_special_positions > 0):
        print >> log, \
          "          Dihedrals discarded because of special positions:",\
          n_dihedrals_discarded_because_of_special_positions
    if (is_unique_model and log is not None):
      if (len(unsupported_chir_volume_sign) > 0):
        print >> log, "          Unsupported chir.volume_sign:", \
          unsupported_chir_volume_sign
      if (n_unresolved_non_hydrogen_chiralities > 0):
        print >> log, "          Unresolved non-hydrogen chiralities:", \
          n_unresolved_non_hydrogen_chiralities
    if (log is not None):
      if (n_chiralities_discarded_because_of_special_positions > 0):
        print >> log, \
          "          Chiralities discarded because of special positions:", \
          n_chiralities_discarded_because_of_special_positions
    if (is_unique_model and log is not None):
      if (len(planarities_with_less_than_four_sites) > 0):
        print >> log, "          Planarities with less than four sites:", \
          planarities_with_less_than_four_sites
      if (n_unresolved_non_hydrogen_planarities > 0):
        print >> log, "          Unresolved non-hydrogen planarities:", \
          n_unresolved_non_hydrogen_planarities
    if (log is not None):
      if (n_planarities_discarded_because_of_special_positions > 0):
        print >> log, \
          "          planarities discarded because of special positions:", \
          n_planarities_discarded_because_of_special_positions
      if (n_bond_proxies_already_assigned_to_first_conformer > 0):
        print >> log, \
          "          bond proxies already assigned to first conformer:", \
          n_bond_proxies_already_assigned_to_first_conformer

class geometry_restraints_proxy_registries(object):

  def __init__(self, n_seq, strict_conflict_handling):
    self.bond_simple = geometry_restraints.bond_simple_proxy_registry(
      n_seq=n_seq, strict_conflict_handling=strict_conflict_handling)
    self.angle = geometry_restraints.angle_proxy_registry(
      strict_conflict_handling=strict_conflict_handling)
    self.dihedral = geometry_restraints.dihedral_proxy_registry(
      strict_conflict_handling=strict_conflict_handling)
    self.chirality = geometry_restraints.chirality_proxy_registry(
      strict_conflict_handling=strict_conflict_handling)
    self.planarity = geometry_restraints.planarity_proxy_registry(
      strict_conflict_handling=strict_conflict_handling)

  # XXX TODO use counts to modify weights

  def initialize_tables(self):
    self.bond_simple.initialize_table()
    self.angle.initialize_table()
    self.dihedral.initialize_table()
    self.chirality.initialize_table()
    self.planarity.initialize_table()

  def discard_tables(self):
    self.bond_simple.discard_table()
    self.angle.discard_table()
    self.dihedral.discard_table()
    self.chirality.discard_table()
    self.planarity.discard_table()

  def report(self, prefix, log):
    if (self.bond_simple.n_resolved_conflicts > 0):
      print >> log, prefix + (
        "Number of resolved bond restraint conflicts: %d"
          % self.bond_simple.n_resolved_conflicts)
    if (self.angle.n_resolved_conflicts > 0):
      print >> log, prefix + (
        "Number of resolved angle restraint conflicts: %d"
          % self.angle.n_resolved_conflicts)
    if (self.dihedral.n_resolved_conflicts > 0):
      print >> log, prefix + (
        "Number of resolved dihedral restraint conflicts: %d"
          % self.dihedral.n_resolved_conflicts)
    if (self.chirality.n_resolved_conflicts > 0):
      print >> log, prefix + (
        "Number of resolved chirality restraint conflicts: %d"
          % self.chirality.n_resolved_conflicts)
    if (self.planarity.n_resolved_conflicts > 0):
      print >> log, prefix + (
        "Number of resolved planarity restraint conflicts: %d"
          % self.planarity.n_resolved_conflicts)

class build_all_chain_proxies(object):

  def __init__(self,
        mon_lib_srv,
        ener_lib,
        file_name=None,
        raw_records=None,
        special_position_settings=None,
        crystal_symmetry=None,
        force_symmetry=False,
        substitute_non_crystallographic_unit_cell_if_necessary=False,
        link_distance_cutoff=3,
        strict_conflict_handling=True,
        keep_monomer_mappings=False,
        max_atoms=None,
        log=None):
    timer = user_plus_sys_time()
    self.time_building_chain_proxies = None
    if (log is not None and file_name is not None):
      print >> log, file_name
    self.stage_1 = pdb.interpretation.stage_1(
      file_name=file_name,
      raw_records=raw_records)
    if (log is not None):
      print >> log, "  Total number of atoms:", \
        len(self.stage_1.atom_attributes_list)
    self.special_position_settings = None
    self._site_symmetry_table = None
    self.sites_cart = None
    self._sites_cart_exact = None
    if (max_atoms is not None
        and len(self.stage_1.atom_attributes_list) > max_atoms):
      if (log is not None):
        print >> log, "  More than %d atoms: no processing." % max_atoms
        return
    self.sites_cart = self.stage_1.get_sites_cart()
    self.special_position_settings \
      = self.stage_1.get_special_position_settings(
          special_position_settings=special_position_settings,
          crystal_symmetry=crystal_symmetry,
          force_symmetry=force_symmetry)
    if (self.special_position_settings is None
        and substitute_non_crystallographic_unit_cell_if_necessary):
      self.special_position_settings = crystal.non_crystallographic_symmetry(
        sites_cart=self.sites_cart).special_position_settings()
    if (self.special_position_settings is None):
      special_position_indices = None
    else:
      special_position_indices = \
        self.site_symmetry_table().special_position_indices()
    models = self.stage_1.get_models_and_conformers()
    if (log is not None):
      print >> log, "  Number of models:", len(models)
    model_type_indices = [-1] * len(models)
    self.scattering_type_registry = scattering_type_registry(
      scattering_types=self.stage_1.get_element_symbols(strip_symbols=True),
      strict_conflict_handling=strict_conflict_handling)
    self.nonbonded_energy_type_registry = nonbonded_energy_type_registry(
      n_seq=len(self.stage_1.atom_attributes_list),
      strict_conflict_handling=strict_conflict_handling)
    self.geometry_proxy_registries = geometry_restraints_proxy_registries(
      n_seq=len(self.stage_1.atom_attributes_list),
      strict_conflict_handling=strict_conflict_handling)
    self.cystein_sulphur_i_seqs = flex.size_t()
    self.cystein_monomer_mappings = []
    self.processed_models = []
    n_unique_models = 0
    for i_model,model in enumerate(models):
      if (log is not None):
        print >> log, "  Model:", model.serial
      is_unique_model = not is_same_model_as_before(
        model_type_indices, i_model, models)
      if (is_unique_model):
        n_unique_models += 1
      elif (log is not None):
        print >> log, "    Same as model", \
          models[model_type_indices[i_model]].serial
      processed_model_ = processed_model(
        stage_1=model.stage_1, serial=model.serial)
      self.processed_models.append(processed_model_)
      if (is_unique_model and log is not None):
        print >> log, "    Number of conformers:", len(model.conformers)
      self.geometry_proxy_registries.initialize_tables()
      for i_conformer,conformer in enumerate(model.conformers):
        if (is_unique_model and log is not None):
          print >> log, '    Conformer: "%s"' % conformer.altLoc
          print >> log, "      Number of atoms:", conformer.iselection.size()
        for j in xrange(i_conformer):
          other = model.conformers[j]
          n = conformer.iselection_common_atoms(other=other).size()
          if (is_unique_model and log is not None):
            print >> log, '      Common with "%s":' % other.altLoc, n
        chains = conformer.get_chains()
        if (is_unique_model and log is not None):
          print >> log, "      Number of chains:", len(chains)
        processed_conformer_ = processed_conformer(
          model=processed_model_,
          altLoc=conformer.altLoc,
          iselection=conformer.iselection)
        processed_model_.add_conformer(processed_conformer_)
        for chain in chains:
          chain_proxies = build_chain_proxies(
            mon_lib_srv=mon_lib_srv,
            link_distance_cutoff=link_distance_cutoff,
            stage_1=self.stage_1,
            sites_cart=self.sites_cart,
            special_position_indices=special_position_indices,
            keep_monomer_mappings=keep_monomer_mappings,
            scattering_type_registry=self.scattering_type_registry,
            nonbonded_energy_type_registry=self.nonbonded_energy_type_registry,
            geometry_proxy_registries=self.geometry_proxy_registries,
            cystein_sulphur_i_seqs=self.cystein_sulphur_i_seqs,
            cystein_monomer_mappings=self.cystein_monomer_mappings,
            i_model=i_model,
            is_unique_model=is_unique_model,
            i_conformer=i_conformer,
            conformer=conformer,
            chain=chain,
            log=log)
          processed_conformer_.add_chain(chain_proxies.processed)
          del chain_proxies
          flush_log(log)
    self.geometry_proxy_registries.discard_tables()
    self.scattering_type_registry.discard_tables()
    self.nonbonded_energy_type_registry.discard_tables()
    if (log is not None):
      if (n_unique_models != 1):
        print >> log, "  Number of unique models:", n_unique_models
      self.geometry_proxy_registries.report(log=log, prefix="  ")
      self.scattering_type_registry.report(
        stage_1=self.stage_1, log=log, prefix="  ")
      self.nonbonded_energy_type_registry.report(
        stage_1=self.stage_1, log=log, prefix="  ")
    self.time_building_chain_proxies = timer.elapsed()

  def monomer_mappings(self):
    for model in self.processed_models:
      for mm in model.monomer_mappings():
        yield mm

  def monomer_mapping_summaries(self):
    for model in self.processed_models:
      for mm in model.monomer_mapping_summaries():
        yield mm

  def site_symmetry_table(self):
    if (self._site_symmetry_table is None):
      assert self.special_position_settings is not None
      assert self.sites_cart is not None
      self._site_symmetry_table = \
        self.special_position_settings.site_symmetry_table(
          sites_cart=self.sites_cart)
    return self._site_symmetry_table

  def sites_cart_exact(self):
    if (self._sites_cart_exact is None):
      self._sites_cart_exact = self.site_symmetry_table().apply_symmetry_sites(
        unit_cell=self.special_position_settings.unit_cell(),
        sites_cart=self.sites_cart)
    return self._sites_cart_exact

  def sel_classification(self, classification):
    result = flex.bool(len(self.stage_1.atom_attributes_list), False)
    for summary in self.monomer_mapping_summaries():
      if (summary.classification == classification):
        result.set_selected(summary.all_associated_i_seqs(), True)
    return result

  def sel_backbone_or_sidechain(self, backbone_flag, sidechain_flag):
    result = flex.bool(len(self.stage_1.atom_attributes_list), False)
    atoms = self.stage_1.atom_attributes_list
    for summary in self.monomer_mapping_summaries():
      if (summary.classification == "peptide"):
        for i_seq in summary.expected_atom_i_seqs:
          if (atoms[i_seq].name.strip() in ["N", "CA", "C", "O"]):
            result[i_seq] = backbone_flag
          else:
            result[i_seq] = sidechain_flag
      elif (summary.classification in ["RNA", "DNA"]):
        for i_seq in summary.expected_atom_i_seqs:
          if (atoms[i_seq].name.strip()
                in ["P", "O1P", "O2P", "O3'", "O5'",
                                       "O3*", "O5*",
                    "O4'", "C1'", "C2'", "C3'", "C4'", "C5'",
                    "O4*", "C1*", "C2*", "C3*", "C4*", "C5*"]):
            result[i_seq] = backbone_flag
          else:
            result[i_seq] = sidechain_flag
    return result

  def sel_backbone(self):
    return self.sel_backbone_or_sidechain(True, False)

  def sel_sidechain(self):
    return self.sel_backbone_or_sidechain(False, True)

  def sel_phosphate(self):
    result = flex.bool(len(self.stage_1.atom_attributes_list), False)
    atoms = self.stage_1.atom_attributes_list
    for summary in self.monomer_mapping_summaries():
      if (summary.classification in ["RNA", "DNA"]):
        for i_seq in summary.expected_atom_i_seqs:
          if (atoms[i_seq].name.strip()
                in ["P", "O1P", "O2P", "O3'", "O5'",
                                       "O3*", "O5*"]):
            result[i_seq] = True
    return result

  def sel_ribose(self):
    result = flex.bool(len(self.stage_1.atom_attributes_list), False)
    atoms = self.stage_1.atom_attributes_list
    for summary in self.monomer_mapping_summaries():
      if (summary.classification in ["RNA", "DNA"]):
        for i_seq in summary.expected_atom_i_seqs:
          if (atoms[i_seq].name.strip()
                in ["O4'", "C1'", "C2'", "C3'", "C4'", "C5'",
                    "O4*", "C1*", "C2*", "C3*", "C4*", "C5*"]):
            result[i_seq] = True
    return result

  def sel_within(self, radius, primary_selection):
    assert radius > 0
    assert self.special_position_settings is not None
    return crystal.neighbors_fast_pair_generator(
      asu_mappings=self.special_position_settings.asu_mappings(
        buffer_thickness=radius,
        sites_cart=self.sites_cart),
      distance_cutoff=radius).neighbors_of(
        primary_selection=primary_selection)

  def _selection_callback(self, word, word_iterator, result_stack):
    lword = word.value.lower()
    if (lword in ["peptide", "protein"]):
      result_stack.append(self.sel_classification(classification="peptide"))
    elif (lword == "rna"):
      result_stack.append(self.sel_classification(classification="RNA"))
    elif (lword == "dna"):
      result_stack.append(self.sel_classification(classification="DNA"))
    elif (lword == "water"):
      result_stack.append(self.sel_classification(classification="water"))
    elif (lword in ["nucleotide", "nuc"]):
      result_stack.append(
          self.sel_classification(classification="RNA")
        | self.sel_classification(classification="DNA"))
    elif (lword == "backbone"):
      result_stack.append(self.sel_backbone())
    elif (lword == "sidechain"):
      result_stack.append(self.sel_sidechain())
    elif (lword == "phosphate"):
      result_stack.append(self.sel_phosphate())
    elif (lword == "ribose"):
      result_stack.append(self.sel_ribose())
    elif (lword == "within"):
      assert word_iterator.pop().value == "("
      radius = float(word_iterator.pop().value)
      assert word_iterator.pop().value == ","
      sel = self.stage_1.selection_cache().selection_parser(
        word_iterator=word_iterator,
        callback=self._selection_callback,
        expect_nonmatching_closing_parenthesis=True)
      result_stack.append(self.sel_within(radius=radius,primary_selection=sel))
    else:
      return False
    return True

  def selection(self, string):
    return self.stage_1.selection_cache().selection(
      string=string,
      callback=self._selection_callback)

  def iselection(self, string):
    return self.selection(string=string).iselection()

  def create_disulfides(self,
        disulfide_distance_cutoff,
        model_indices,
        conformer_indices,
        log=None):
    if (model_indices is not None):
      model_indices = model_indices.select(self.cystein_sulphur_i_seqs)
    conformer_indices = conformer_indices.select(self.cystein_sulphur_i_seqs)
    asu_mappings = self.special_position_settings.asu_mappings(
      buffer_thickness=disulfide_distance_cutoff)
    sulphur_sites_cart = self.sites_cart.select(self.cystein_sulphur_i_seqs)
    asu_mappings.process_sites_cart(
      original_sites=sulphur_sites_cart,
      site_symmetry_table=self.site_symmetry_table().select(
        self.cystein_sulphur_i_seqs))
    pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
    nonbonded_proxies = geometry_restraints.nonbonded_sorted_asu_proxies(
      model_indices=model_indices,
      conformer_indices=conformer_indices,
      nonbonded_params=geometry_restraints.nonbonded_params(
        default_distance=1),
      nonbonded_types=flex.std_string(conformer_indices.size()),
      nonbonded_distance_cutoff_plus_buffer=disulfide_distance_cutoff,
      shell_asu_tables=[pair_asu_table])
    labels = [self.stage_1.atom_attributes_list[i_seq].pdb_format()
      for i_seq in self.cystein_sulphur_i_seqs]
    for proxy in nonbonded_proxies.simple:
      pair_asu_table.add_pair(proxy.i_seqs)
    for proxy in nonbonded_proxies.asu:
      pair_asu_table.add_pair(proxy)
    pair_sym_table = pair_asu_table.extract_pair_sym_table()
    if (log is not None):
      n_simple, n_symmetry = 0, 0
      for sym_pair in pair_sym_table.iterator():
        if (sym_pair.rt_mx_ji.is_unit_mx()): n_simple += 1
        else:                                n_symmetry += 1
      print >> log, "  Number of disulfides: simple=%d, symmetry=%d" % (
        n_simple, n_symmetry)
      if (n_symmetry == 0):
         blanks = ""
      else:
         blanks = "  "
    max_distance_model = 0
    frac = asu_mappings.unit_cell().fractionalize
    orth = asu_mappings.unit_cell().orthogonalize
    for sym_pair in pair_sym_table.iterator():
      distance_model = geometry_restraints.bond(
        [sulphur_sites_cart[sym_pair.i_seq],
         orth(sym_pair.rt_mx_ji*frac(sulphur_sites_cart[sym_pair.j_seq]))],
        distance_ideal=0,
        weight=1).distance_model
      max_distance_model = max(max_distance_model, distance_model)
      if (log is not None):
        if (sym_pair.rt_mx_ji.is_unit_mx()):
          disulfide_type = "Simple disulfide:%s" % blanks
        else:
          disulfide_type = "Symmetry disulfide:"
        print >> log, "    %s %s - %s distance=%.2f" % tuple(
          [disulfide_type]
          + [labels[i_seq] for i_seq in sym_pair.i_seqs()]
          + [distance_model]),
        if (not sym_pair.rt_mx_ji.is_unit_mx()):
          print >> log, sym_pair.rt_mx_ji,
        print >> log
    return pair_sym_table, max_distance_model

  def construct_geometry_restraints_manager(self,
        ener_lib,
        disulfide_link,
        disulfide_distance_cutoff=3,
        nonbonded_distance_cutoff=None,
        default_vdw_distance=1,
        min_vdw_distance=1,
        nonbonded_buffer=1,
        vdw_1_4_factor=2/3.,
        plain_pairs_radius=None,
        log=None):
    assert self.special_position_settings is not None
    timer = user_plus_sys_time()
    sel_cache = self.stage_1.selection_cache()
    model_indices = sel_cache.get_model_indices()
    conformer_indices = sel_cache.get_conformer_indices()
    bond_params_table = geometry_restraints.extract_bond_params(
      n_seq=sel_cache.n_seq,
      bond_simple_proxies=self.geometry_proxy_registries.bond_simple.proxies)
    bond_distances_model = geometry_restraints.bond_distances_model(
      sites_cart=self.sites_cart,
      proxies=self.geometry_proxy_registries.bond_simple.proxies)
    if (bond_distances_model.size() > 0):
      excessive_bonds = (
        bond_distances_model > self.special_position_settings.unit_cell()
          .shortest_vector_sq()**.5).iselection()
      if (excessive_bonds.size() > 0):
        atoms = self.stage_1.atom_attributes_list
        proxies = self.geometry_proxy_registries.bond_simple.proxies
        print >> log, "  Bonds with excessive lengths:"
        for i_proxy in excessive_bonds:
          proxy = proxies[i_proxy]
          bond = geometry_restraints.bond(
            sites_cart=self.sites_cart, proxy=proxy)
          print >> log, "    Distance model: %.6g (ideal: %.6g)" % (
            bond.distance_model, bond.distance_ideal)
          for i,i_seq in enumerate(proxy.i_seqs):
            print >> log, "      atom %d: %s" % (i+1,atoms[i_seq].pdb_format())
        raise Sorry("Number of bonds with excessive lengths: %d" %
          excessive_bonds.size())
    disulfide_sym_table, max_disulfide_bond_distance = \
      self.create_disulfides(
        disulfide_distance_cutoff=disulfide_distance_cutoff,
        model_indices=model_indices,
        conformer_indices=conformer_indices,
        log=log)
    max_bond_distance = max_disulfide_bond_distance
    if (bond_distances_model.size() > 0):
      max_bond_distance = max(max_bond_distance,
        flex.max(bond_distances_model))
    asu_mappings = self.special_position_settings.asu_mappings(
      buffer_thickness=max_bond_distance*3)
    asu_mappings.process_sites_cart(
      original_sites=self.sites_cart,
      site_symmetry_table=self.site_symmetry_table())
    bond_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
    geometry_restraints.add_pairs(
      bond_asu_table, self.geometry_proxy_registries.bond_simple.proxies)
    self.geometry_proxy_registries.bond_simple.proxies = None # free memory
    #
    disulfide_bond = disulfide_link.bond_list[0]
    assert disulfide_bond.value_dist is not None
    assert disulfide_bond.value_dist > 0
    assert disulfide_bond.value_dist_esd is not None
    assert disulfide_bond.value_dist_esd > 0
    for sym_pair in disulfide_sym_table.iterator():
      i_seq = self.cystein_sulphur_i_seqs[sym_pair.i_seq]
      j_seq = self.cystein_sulphur_i_seqs[sym_pair.j_seq]
      params = geometry_restraints.bond_params(
        distance_ideal=disulfide_bond.value_dist,
        weight=1/disulfide_bond.value_dist_esd**2)
      if (i_seq <= j_seq): bond_params_table[i_seq][j_seq] = params
      else:                bond_params_table[j_seq][i_seq] = params
      bond_asu_table.add_pair(
        i_seq=i_seq,
        j_seq=j_seq,
        rt_mx_ji=sym_pair.rt_mx_ji)
    #
    shell_asu_tables = crystal.coordination_sequences.shell_asu_tables(
      pair_asu_table=bond_asu_table,
      max_shell=3)
    shell_sym_tables = [shell_asu_table.extract_pair_sym_table()
      for shell_asu_table in shell_asu_tables]
    #
    nonbonded_params = ener_lib_as_nonbonded_params(
      ener_lib=ener_lib,
      factor_1_4_interactions=vdw_1_4_factor,
      default_distance=default_vdw_distance,
      minimum_distance=min_vdw_distance)
    result = geometry_restraints.manager.manager(
      crystal_symmetry=self.special_position_settings,
      model_indices=model_indices,
      conformer_indices=conformer_indices,
      site_symmetry_table=self.site_symmetry_table(),
      bond_params_table=bond_params_table,
      shell_sym_tables=shell_sym_tables,
      nonbonded_params=nonbonded_params,
      nonbonded_types=self.nonbonded_energy_type_registry.symbols,
      nonbonded_function=geometry_restraints.prolsq_repulsion_function(),
      nonbonded_distance_cutoff=nonbonded_distance_cutoff,
      nonbonded_buffer=1,
      angle_proxies=self.geometry_proxy_registries.angle.proxies,
      dihedral_proxies=self.geometry_proxy_registries.dihedral.proxies,
      chirality_proxies=self.geometry_proxy_registries.chirality.proxies,
      planarity_proxies=self.geometry_proxy_registries.planarity.proxies,
      plain_pairs_radius=plain_pairs_radius)
    self.time_building_geometry_restraints_manager = timer.elapsed()
    return result

  def extract_xray_structure(self):
    assert self.special_position_settings is not None
    return self.stage_1.extract_xray_structure(
      special_position_settings=self.special_position_settings,
      sites_cart=self.sites_cart_exact(),
      site_symmetry_table=self.site_symmetry_table(),
      scattering_types=self.scattering_type_registry.symbols,
      unknown_scattering_type_substitute="?")

class process(object):

  def __init__(self,
        mon_lib_srv,
        ener_lib,
        file_name=None,
        raw_records=None,
        strict_conflict_handling=False,
        special_position_settings=None,
        crystal_symmetry=None,
        force_symmetry=False,
        keep_monomer_mappings=False,
        max_atoms=None,
        log=None):
    self.mon_lib_srv = mon_lib_srv
    self.ener_lib = ener_lib
    self.log = log
    self.all_chain_proxies = build_all_chain_proxies(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=file_name,
      raw_records=raw_records,
      special_position_settings=special_position_settings,
      crystal_symmetry=crystal_symmetry,
      force_symmetry=force_symmetry,
      strict_conflict_handling=strict_conflict_handling,
      keep_monomer_mappings=keep_monomer_mappings,
      max_atoms=max_atoms,
      log=log)
    if (log is not None
        and self.all_chain_proxies.time_building_chain_proxies is not None):
      print >> log, \
        "  Time building chain proxies: %.2f, per 1000 atoms: %.2f" % (
          self.all_chain_proxies.time_building_chain_proxies,
          self.all_chain_proxies.time_building_chain_proxies * 1000
            / len(self.all_chain_proxies.stage_1.atom_attributes_list))
    self._geometry_restraints_manager = None
    self._xray_structure = None

  def geometry_restraints_manager(self,
        plain_pairs_radius=None,
        show_energies=True):
    if (    self.all_chain_proxies.sites_cart is not None
        and self.all_chain_proxies.special_position_settings is not None
        and self._geometry_restraints_manager is None):
      self._geometry_restraints_manager \
        = self.all_chain_proxies.construct_geometry_restraints_manager(
            ener_lib=self.ener_lib,
            disulfide_link=self.mon_lib_srv.link_link_id_dict["SS"],
            plain_pairs_radius=plain_pairs_radius,
            log=self.log)
      if (self.log is not None):
        print >> self.log, \
          "  Time building geometry restraints manager: %.2f seconds" % (
            self.all_chain_proxies.time_building_geometry_restraints_manager)
        flush_log(self.log)
        labels = [atom.pdb_format()
          for atom in self.all_chain_proxies.stage_1.atom_attributes_list]
        pair_proxies = self._geometry_restraints_manager.pair_proxies(
          sites_cart=self.all_chain_proxies.sites_cart_exact())
        pair_proxies.bond_proxies.show_histogram_of_model_distances(
          sites_cart=self.all_chain_proxies.sites_cart_exact(),
          f=self.log,
          prefix="  ")
        pair_proxies.bond_proxies.show_sorted_by_residual(
          sites_cart=self.all_chain_proxies.sites_cart_exact(),
          labels=labels,
          f=self.log,
          prefix="  ",
          max_lines=5)
        self._geometry_restraints_manager.dihedral_proxies \
          .show_histogram_of_deltas(
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            f=self.log,
            prefix="  ")
        self._geometry_restraints_manager.dihedral_proxies \
          .show_sorted_by_residual(
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            labels=labels,
            f=self.log,
            prefix="  ",
            max_lines=3)
        flush_log(self.log)
        if (show_energies):
          timer = user_plus_sys_time()
          energies = self._geometry_restraints_manager.energies_sites(
            sites_cart=self.all_chain_proxies.sites_cart_exact())
          energies.show(f=self.log, prefix="  ")
          print >> self.log, "  Time first energy calculation" \
                             " (mainly nonbonded setup): %.2f" % (
            timer.elapsed())
          flush_log(self.log)
    return self._geometry_restraints_manager

  def xray_structure(self):
    if (    self.all_chain_proxies.sites_cart is not None
        and self.all_chain_proxies.special_position_settings is not None
        and self._xray_structure is None):
      self._xray_structure = self.all_chain_proxies.extract_xray_structure()
      self._xray_structure.scattering_type_registry(
        types_without_a_scattering_contribution=["?"])
      if (self.log is not None):
        self._xray_structure.show_summary(f=self.log, prefix="  ")
        self._xray_structure.show_special_position_shifts(
          sites_cart_original=self.all_chain_proxies.sites_cart,
          out=self.log, prefix="  ")
        self._xray_structure.scattering_type_registry().show(
          show_gaussians=False, out=self.log, prefix="  ")
        flush_log(self.log)
    return self._xray_structure

def run(args, strict_conflict_handling=True, max_atoms=None):
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  for file_name in args:
    processed_pdb_file = process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=file_name,
      strict_conflict_handling=strict_conflict_handling,
      max_atoms=max_atoms,
      log=sys.stdout)
    processed_pdb_file.geometry_restraints_manager()
    processed_pdb_file.xray_structure()
  return processed_pdb_file

if (__name__ == "__main__"):
  run(sys.argv[1:])
