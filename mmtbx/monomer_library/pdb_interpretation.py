from __future__ import division
from iotbx import pdb
import iotbx.phil
from mmtbx.monomer_library import server
from mmtbx.monomer_library import cif_types
from mmtbx.monomer_library import rna_sugar_pucker_analysis
from mmtbx.monomer_library import conformation_dependent_restraints
from cctbx import geometry_restraints
import cctbx.geometry_restraints.manager
from cctbx import crystal
import cctbx.crystal.coordination_sequences
from cctbx import sgtbx
from cctbx.array_family import flex
from scitbx.python_utils import dicts
from libtbx.str_utils import show_string
from libtbx.utils import flat_list, Sorry, user_plus_sys_time, plural_s
from libtbx.utils import format_exception
from libtbx import Auto, group_args, slots_getstate_setstate
from cStringIO import StringIO
import string
import sys, os

# see iotbx/pdb/common_residue_names.h
ad_hoc_single_atom_residue_element_types = """\
ZN CA MG CL NA MN K FE CU CD HG NI CO BR XE SR CS PT BA TL PB SM AU RB YB LI
KR MO LU CR OS GD TB LA F AR AG HO GA CE W SE RU RE PR IR EU AL V TE SB PD
""".split()

class ad_hoc_single_atom_residue(object):

  def __init__(self, residue_name, atom_name, atom_element):
    atom_element = atom_element.strip().upper()
    if (atom_element in ad_hoc_single_atom_residue_element_types):
      self.scattering_type = atom_element
      self.energy_type = atom_element
      return
    atom_name = atom_name.strip().upper()
    if (    len(atom_element) == 0
        and atom_name == residue_name.strip()
        and atom_name in ad_hoc_single_atom_residue_element_types):
      self.scattering_type = atom_name
      self.energy_type = atom_name
      return
    if (    residue_name == "NH3"
        and atom_element == "N"
        or (len(atom_element) == 0 and atom_name.startswith("N"))):
      self.scattering_type = "N"
      self.energy_type = "N"
      return
    if (    residue_name == "CH4"
        and atom_element == "C"
        or (len(atom_element) == 0 and atom_name.startswith("C"))):
      self.scattering_type = "C"
      self.energy_type = "C"
      return
    self.scattering_type = None
    self.energy_type = None

dihedral_function_type_params_str = """\
  dihedral_function_type = *determined_by_sign_of_periodicity \
                            all_sinusoidal \
                            all_harmonic
    .type=choice
    .optional=False
"""

clash_guard_params_str = """\
  clash_guard
    .short_caption = Clash guard
    .style = noauto box auto_align
    .expert_level=2
  {
    nonbonded_distance_threshold = 0.5
      .type = float
    max_number_of_distances_below_threshold = 100
      .type = int
    max_fraction_of_distances_below_threshold = 0.1
      .type = float
  }
"""

master_params_str = """\
  apply_cif_modification
    .optional = True
    .multiple = True
    .short_caption = Modify CIF
    .style = noauto auto_align
  {
    data_mod = None
      .type = str
    residue_selection = None
      .type = str
      .style = selection
  }
  apply_cif_link
    .optional = True
    .multiple = True
    .short_caption = Add link to CIF
    .style = noauto auto_align
  {
    data_link = None
      .type = str
    residue_selection_1 = None
      .type = str
      .style = selection
    residue_selection_2 = None
      .type = str
      .style = selection
  }
  link_distance_cutoff = 3
    .type=float
    .optional=False
  disulfide_distance_cutoff = 3
    .type=float
    .optional=False
  %(dihedral_function_type_params_str)s
  chir_volume_esd = 0.2
    .type=float
    .optional=False
  peptide_link
    .short_caption = Peptide link settings
    .style = box auto_align noauto
  {
    cis_threshold = 45
      .type = float
      .optional = False
    discard_psi_phi = True
      .type = bool
      .optional = False
    omega_esd_override_value = None
      .type = float
  }
  max_reasonable_bond_distance = 50.0
    .type=float
  nonbonded_distance_cutoff = None
    .type=float
  default_vdw_distance = 1
    .type=float
    .optional=False
  min_vdw_distance = 1
    .type=float
    .optional=False
  nonbonded_buffer = 1
    .type=float
    .optional=False
  vdw_1_4_factor = 0.8
    .type=float
    .optional=False
  translate_cns_dna_rna_residue_names = None
    .type=bool
    .optional=False
  proceed_with_excessive_length_bonds = False
    .type=bool
  rna_sugar_pucker_analysis
    .short_caption = RNA sugar pucker analysis
    .style = box noauto auto_align menu_item parent_submenu:advanced
  {
    include scope mmtbx.monomer_library.rna_sugar_pucker_analysis.master_phil
  }
  show_histogram_slots
    .style = box auto_align noauto
    .expert_level = 2
  {
    bond_lengths = 5
      .type=int
    nonbonded_interaction_distances = 5
      .type=int
    bond_angle_deviations_from_ideal = 5
      .type=int
    dihedral_angle_deviations_from_ideal = 5
      .type=int
    chiral_volume_deviations_from_ideal = 5
      .type=int
  }
  show_max_items
    .expert_level = 2
    .style = box auto_align noauto
  {
    not_linked = 5
      .type=int
    bond_restraints_sorted_by_residual = 5
      .type=int
    nonbonded_interactions_sorted_by_model_distance = 5
      .type=int
    bond_angle_restraints_sorted_by_residual = 5
      .type=int
    dihedral_angle_restraints_sorted_by_residual = 3
      .type=int
    chirality_restraints_sorted_by_residual = 3
      .type=int
    planarity_restraints_sorted_by_residual = 3
      .type=int
    residues_with_excluded_nonbonded_symmetry_interactions = 12
      .type=int
  }
  %(clash_guard_params_str)s
""" % vars()

master_params = iotbx.phil.parse(
  input_string=master_params_str,
  process_includes=True)

geometry_restraints_edits_str = """\
excessive_bond_distance_limit = 10
  .type = float
bond
  .optional = True
  .multiple = True
  .short_caption = Bond
  .style = auto_align
{
  action = *add delete change
    .type = choice
  atom_selection_1 = None
    .type = str
  atom_selection_2 = None
    .type = str
  symmetry_operation = None
    .help = "The bond is between atom_1 and symmetry_operation * atom_2,"
            " with atom_1 and atom_2 given in fractional coordinates."
            " Example: symmetry_operation = -x-1,-y,z"
    .type = str
  distance_ideal = None
    .type = float
  sigma = None
    .type = float
  slack = None
    .type = float
}
angle
  .optional = True
  .multiple = True
  .short_caption = Angle
  .style = auto_align
{
  action = *add delete change
    .type = choice
  atom_selection_1 = None
    .type = str
  atom_selection_2 = None
    .type = str
  atom_selection_3 = None
    .type = str
  angle_ideal = None
    .type = float
  sigma = None
    .type = float
}
planarity
  .optional = True
  .multiple = True
  .short_caption = Planarity
  .style = auto_align
{
  action = *add delete change
    .type = choice
  atom_selection = None
    .type = str
    .multiple = True
  sigma = None
    .type = float
}

"""

geometry_restraints_remove_str = """\
angles=None
  .optional=True
  .type=str
  .multiple=True
  .input_size=400
dihedrals=None
  .optional=True
  .type=str
  .multiple=True
  .input_size=400
chiralities=None
  .optional=True
  .type=str
  .multiple=True
  .input_size=400
planarities=None
  .optional=True
  .type=str
  .multiple=True
  .input_size=400
"""

grand_master_phil_str = """\
pdb_interpretation {
  %(master_params_str)s
}
geometry_restraints.edits {
  %(geometry_restraints_edits_str)s
}
geometry_restraints.remove {
  %(geometry_restraints_remove_str)s
}
""" % vars()

def flush_log(log):
  if (log is not None):
    flush = getattr(log, "flush", None)
    if (flush is not None): flush()

def all_atoms_are_in_main_conf(atoms):
  for atom in atoms:
    if (atom.parent().altloc != ""): return False
  return True

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

def involves_broken_bonds(broken_bond_i_seq_pairs, i_seqs):
  if (broken_bond_i_seq_pairs is None): return False
  i_seqs = sorted(i_seqs)
  for i in xrange(len(i_seqs)-1):
    for j in xrange(i+1,len(i_seqs)):
      if ((i_seqs[i],i_seqs[j]) in broken_bond_i_seq_pairs):
        return True
  return False

class source_info_server(object):

  def __init__(self, m_i, m_j):
    self.m_i, self.m_j = m_i, m_j

  def labels(self):
    if (self.m_j is None):
      return "residue: %s" % self.m_i.residue_altloc()
    return "residues: %s + %s" % (
      self.m_i.residue_altloc(),
      self.m_j.residue_altloc())

  def n_expected_atoms(self):
    if (self.m_j is None):
      return len(self.m_i.expected_atoms)
    return len(self.m_i.expected_atoms) \
         + len(self.m_j.expected_atoms)

def _show_atom_labels(pdb_atoms, i_seqs, out=None, prefix="", max_lines=None):
  if (out is None): out = sys.stdout
  for i_line,i_seq in enumerate(i_seqs):
    if (i_line == max_lines and len(i_seqs) > max_lines+1):
      print >> out, prefix + "... (remaining %d not shown)" % (
        len(i_seqs)-max_lines)
      break
    print >> out, prefix + pdb_atoms[i_seq].quote()

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
  _show_atom_labels(
    pdb_atoms=m_i.pdb_atoms, i_seqs=i_seqs, out=s, prefix="    ", max_lines=10)
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

  def assign_directly(self, i_seq, symbol):
    self.symbols[i_seq] = symbol

  def assign_from_monomer_mapping(self, conf_altloc, mm):
    atom_dict = mm.monomer_atom_dict
    for atom_id,atom in mm.expected_atoms.items():
      i_seq = atom.i_seq
      if (self.type_label == "scattering"):
        symbol = atom_dict[atom_id].type_symbol
        if (symbol == "H" and self.symbols[i_seq] == "D"):
          symbol = "D"
      else:
        symbol = atom_dict[atom_id].type_energy
      if (symbol is None): continue
      source_label = mm.residue_altloc()
      source_n_expected_atoms = len(mm.expected_atoms)
      prev_symbol = self.symbols[i_seq]
      prev_source_label = self.source_labels[i_seq]
      prev_source_n_expected_atoms = self.source_n_expected_atoms[i_seq]
      assign = False
      raise_conflict = False
      if (prev_symbol == ""
          or (self.type_label == "scattering") and prev_symbol in ["Q","X"]):
        assign = True
      elif (prev_symbol.upper() != symbol.upper()):
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

  def report_unknown_message(self):
    return "Number of atoms with unknown %s type symbols" % self.type_label

  def report(self, pdb_atoms, log, prefix, max_lines=10):
    n_unknown = self.n_unknown_type_symbols()
    if (n_unknown > 0):
      print >> log, "%s%s: %d" % (
        prefix, self.report_unknown_message(), n_unknown)
      i_seqs = (self.symbols == "").iselection()
      _show_atom_labels(
        pdb_atoms=pdb_atoms, i_seqs=i_seqs,
        out=log, prefix=prefix+"  ", max_lines=max_lines)
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

class monomer_mapping_summary(slots_getstate_setstate):

  __slots__ = [
    "conf_altloc",
    "residue_name",
    "expected_atoms",
    "unexpected_atoms",
    "duplicate_atoms",
    "ignored_atoms",
    "classification",
    "incomplete_info",
    "is_terminus",
    "is_unusual"]

  def __init__(self, **keyword_args):
    for key in monomer_mapping_summary.__slots__:
      setattr(self, key, keyword_args.get(key, None))

  def all_associated_i_seqs(self):
    return flex.size_t(
        [a.i_seq for a in self.expected_atoms]
      + [a.i_seq for a in self.unexpected_atoms]
      + [a.i_seq for a in self.duplicate_atoms]
      + [a.i_seq for a in self.ignored_atoms])

  def summary(self):
    return self

class monomer_mapping(slots_getstate_setstate):

  __slots__ = [
    "active_atoms",
    "angle_counters",
    "atom_name_interpretation",
    "atom_names_given",
    "bond_counters",
    "chem_mod_ids",
    "chirality_counters",
    "classification",
    "conf_altloc",
    "dihedral_counters",
    "duplicate_atoms",
    "expected_atoms",
    "i_conformer",
    "ignored_atoms",
    "incomplete_info",
    "is_first_conformer_in_chain",
    "is_rna2p",
    "is_rna_dna",
    "is_terminus",
    "is_unusual",
    "lib_link",
    "missing_hydrogen_atoms",
    "missing_non_hydrogen_atoms",
    "mon_lib_names",
    "mon_lib_srv",
    "monomer",
    "monomer_atom_dict",
    "pdb_atoms",
    "pdb_residue",
    "pdb_residue_id_str",
    "planarity_counters",
    "residue_name",
    "unexpected_atoms"]

  def __init__(self,
        pdb_atoms,
        mon_lib_srv,
        translate_cns_dna_rna_residue_names,
        rna_sugar_pucker_analysis_params,
        apply_cif_modifications,
        apply_cif_links_mm_pdbres_dict,
        i_model,
        i_conformer,
        is_first_conformer_in_chain,
        conf_altloc,
        pdb_residue,
        next_pdb_residue):
    self.pdb_atoms = pdb_atoms
    self.mon_lib_srv = mon_lib_srv
    self.i_conformer = i_conformer
    self.is_first_conformer_in_chain = is_first_conformer_in_chain
    self.conf_altloc = conf_altloc
    self.pdb_residue = pdb_residue
    self.pdb_residue_id_str = pdb_residue.id_str(suppress_segid=-1)
    self.residue_name = pdb_residue.resname
    atom_id_str_pdbres_list = self._collect_atom_names()
    self.monomer, self.atom_name_interpretation \
      = self.mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(
          residue_name=self.residue_name,
          atom_names=self.atom_names_given,
          translate_cns_dna_rna_residue_names
            =translate_cns_dna_rna_residue_names)
    if (self.atom_name_interpretation is None):
      self.mon_lib_names = None
    else:
      self.mon_lib_names = self.atom_name_interpretation.mon_lib_names()
    if (self.monomer is None):
      self.expected_atoms = {}
      self.unexpected_atoms = {}
      self.duplicate_atoms = {}
      for atom,atom_name in zip(self.active_atoms, self.atom_names_given):
        self.unexpected_atoms[atom_name] = atom
      self.incomplete_info = None
      self.is_terminus = None
    else:
      self.chem_mod_ids = set()
      self._rna_sugar_pucker_analysis(
        params=rna_sugar_pucker_analysis_params,
        next_pdb_residue=next_pdb_residue)
      self._get_mappings()
      for id_str in atom_id_str_pdbres_list:
        for apply_data_mod in apply_cif_modifications.get(id_str, []):
          self.apply_mod(
            mod_mod_id=self.mon_lib_srv.mod_mod_id_dict[apply_data_mod])
      self._set_incomplete_info()
      self.is_terminus = None
      self.monomer.set_classification()
      if (self.incomplete_info is None):
        self.resolve_unexpected()
    if (self.pdb_residue_id_str in apply_cif_links_mm_pdbres_dict):
      apply_cif_links_mm_pdbres_dict[self.pdb_residue_id_str].setdefault(
        self.i_conformer, []).append(self)

  def _collect_atom_names(self):
    self.ignored_atoms = {}
    self.active_atoms = []
    self.atom_names_given = []
    atom_id_str_pdbres_set = set()
    for atom in self.pdb_residue.atoms():
      atom_id_str_pdbres_set.add(atom.id_str(pdbres=True))
      if (atom.element.strip() == "Q"):
        self.ignored_atoms.setdefault(atom.name, []).append(atom)
      else:
        self.active_atoms.append(atom)
        self.atom_names_given.append(atom.name.replace(" ",""))
    return sorted(atom_id_str_pdbres_set)

  def _rna_sugar_pucker_analysis(self, params, next_pdb_residue):
    self.is_rna_dna = False
    self.is_rna2p = None
    if (self.monomer.is_peptide()): return
    from iotbx.pdb.rna_dna_detection import residue_analysis
    ra1 = residue_analysis(
      residue_atoms=self.pdb_residue.atoms(),
      distance_tolerance=params.bond_detection_distance_tolerance)
    if (ra1.problems is not None): return
    self.is_rna_dna = True
    if (not ra1.is_rna): return
    residue_2_p_atom = None
    if (next_pdb_residue is not None):
      residue_2_p_atom = next_pdb_residue.find_atom_by(name=" P  ")
    ana = rna_sugar_pucker_analysis.evaluate(
      params=params,
      residue_1_deoxy_ribo_atom_dict=ra1.deoxy_ribo_atom_dict,
      residue_1_c1p_outbound_atom=ra1.c1p_outbound_atom,
      residue_2_p_atom=residue_2_p_atom)
    self.is_rna2p = ana.is_2p
    if (self.is_rna2p): primary_mod_id = "rna2p"
    else:               primary_mod_id = "rna3p"
    self.monomer, chem_mod_ids = self.mon_lib_srv.get_comp_comp_id_mod(
      comp_comp_id=self.monomer,
      mod_ids=(primary_mod_id,))
    self._track_mods(chem_mod_ids=chem_mod_ids)

  def _get_mappings(self):
    self.monomer_atom_dict = atom_dict = self.monomer.atom_dict()
    deuterium_aliases = None
    processed_atom_names = {}
    self.expected_atoms = {}
    self.unexpected_atoms = {}
    self.duplicate_atoms = {}
    if (self.atom_name_interpretation is not None):
      replace_primes = False
    elif (self.is_rna_dna or self.monomer.is_rna_dna()):
      replace_primes = True
    else:
      n_primes = 0
      n_stars = 0
      for atom_name in self.atom_names_given:
        if (atom_name.find("'") >= 0): n_primes += 1
        if (atom_name.find("*") >= 0): n_stars += 1
      replace_primes = (n_primes != 0 and n_stars == 0)
    rna_dna_bb_cif_by_ref = None
    for i_atom,atom in enumerate(self.active_atoms):
      atom_name_given = self.atom_names_given[i_atom]
      if (self.mon_lib_names is None):
        atom_name = atom_name_given
      else:
        atom_name = self.mon_lib_names[i_atom]
        if (atom_name is None):
          atom_name = atom_name_given
      if (len(atom_name) != 0 and not atom_dict.has_key(atom_name)):
        auto_synomyms = []
        if (atom_name[0] in string.digits):
          auto_synomyms.append(atom_name[1:] + atom_name[0])
        elif (atom_name[-1] in string.digits):
          auto_synomyms.append(atom_name[-1] + atom_name[0:-1])
        if (replace_primes):
          atom_name = atom_name.replace("'", "*")
          if (atom_name != atom_name_given):
            auto_synomyms.append(atom_name)
            if (atom_name[0] in string.digits):
              auto_synomyms.append(atom_name[1:] + atom_name[0])
            elif (atom_name[-1] in string.digits):
              auto_synomyms.append(atom_name[-1] + atom_name[0:-1])
        for atom_name in auto_synomyms:
          if (atom_dict.has_key(atom_name)): break
        else:
          auto_synomyms.insert(0, atom_name_given)
          if (deuterium_aliases is None):
            deuterium_aliases = self.monomer.hydrogen_deuterium_aliases()
          for atom_name in auto_synomyms:
            atom_name = deuterium_aliases.get(atom_name)
            if (atom_name is not None): break
          else:
            for atom_name in auto_synomyms:
              atom_name = self.mon_lib_srv.comp_synonym_atom_list_dict.get(
                self.monomer.chem_comp.id, {}).get(atom_name, None)
              if (atom_name is not None): break
            else:
              atom_name = atom_name_given
      if (    len(atom_name) != 0
          and not atom_dict.has_key(atom_name)
          and self.is_rna_dna):
        aliases = pdb.rna_dna_atom_names_backbone_aliases
        if (rna_dna_bb_cif_by_ref is None):
          rna_dna_bb_cif_by_ref = {}
          for cif_name in atom_dict:
            ref_name = aliases.get(cif_name)
            if (ref_name is not None):
              rna_dna_bb_cif_by_ref[ref_name] = cif_name
        ref_name = aliases.get(atom_name)
        cif_name = rna_dna_bb_cif_by_ref.get(ref_name)
        if (cif_name is not None):
          atom_name = cif_name
      prev_atom = processed_atom_names.get(atom_name)
      if (prev_atom is None):
        processed_atom_names[atom_name] = atom
        if (atom_dict.has_key(atom_name)):
          self.expected_atoms[atom_name] = atom
        else:
          self.unexpected_atoms[atom_name] = atom
      else:
        self.duplicate_atoms.setdefault(atom_name, []).append(atom)
    if (    self.monomer.is_peptide()
        and self.atom_name_interpretation is None):
      self._rename_ot1_ot2("OXT" in atom_dict)
    self._set_missing_atoms()

  def _rename_ot1_ot2(self, oxt_in_atom_dict):
    if (not self.expected_atoms.has_key("O")):
      i_seq = self.unexpected_atoms.get("OT1", None)
      if (i_seq is not None):
        self.expected_atoms["O"] = i_seq
        del self.unexpected_atoms["OT1"]
    if (oxt_in_atom_dict):
      oxt_dict = self.expected_atoms
    else:
      oxt_dict = self.unexpected_atoms
    if (not oxt_dict.has_key("OXT")):
      i_seq = self.unexpected_atoms.get("OT2", None)
      if (i_seq is not None):
        oxt_dict["OXT"] = i_seq
        del self.unexpected_atoms["OT2"]

  def _set_missing_atoms(self):
    self.missing_non_hydrogen_atoms = {}
    self.missing_hydrogen_atoms = {}
    for atom in self.monomer.atom_list:
      if (not self.expected_atoms.has_key(atom.atom_id)):
        if (atom.type_symbol != "H"):
          self.missing_non_hydrogen_atoms[atom.atom_id] = atom
        else:
          self.missing_hydrogen_atoms[atom.atom_id] = atom

  def _get_incomplete_info(self):
    if (    len(self.unexpected_atoms) == 0
        and len(self.missing_non_hydrogen_atoms) > 0):
      if (self.monomer.is_peptide()):
        atom_ids = self.expected_atoms.keys()
        atom_ids.sort()
        atom_ids = " ".join(atom_ids)
        if (atom_ids == "CA"): return "c_alpha_only"
        if (atom_ids == "C CA N"): return "n_c_alpha_c_only"
        if (atom_ids == "C CA N O"): return "backbone_only"
        if (atom_ids == "C CA CB N O"): return "truncation_to_alanine"
      elif (self.is_rna_dna or self.monomer.is_rna_dna()):
        atom_ids = " ".join(self.expected_atoms.keys())
        if (atom_ids == "P"): return "p_only"
    return None

  def _set_incomplete_info(self):
    self.incomplete_info = self._get_incomplete_info()

  def resolve_unexpected(self):
    mod_mod_ids = []
    ani = self.atom_name_interpretation
    u = self.unexpected_atoms
    if (self.monomer.classification == "peptide"):
      if (ani is not None):
        u_mon_lib = {}
        for given_name,mon_lib_name in zip(ani.atom_names,
                                           self.mon_lib_names):
          i_seq = u.get(given_name)
          if (i_seq is None): continue
          if (mon_lib_name is None):
            u_mon_lib[given_name] = i_seq
          else:
            u_mon_lib[mon_lib_name] = i_seq
        u = u_mon_lib
      if ("HXT" in u):
        mod_mod_ids.append(self.mon_lib_srv.mod_mod_id_dict["COOH"])
      elif ("OXT" in u):
        mod_mod_ids.append(self.mon_lib_srv.mod_mod_id_dict["COO"])
      if (ani is not None):
        nitrogen_hydrogens = []
        for name in u.keys():
          # name is a mon_lib_name
          if (name in ["H1", "H2", "H3"]):
            nitrogen_hydrogens.append(name)
        nitrogen_hydrogen_translation = None
        if (len(nitrogen_hydrogens) == 3):
          if (self.monomer_atom_dict.get("H") is not None):
            mod_mod_ids.append(self.mon_lib_srv.mod_mod_id_dict["NH3"])
        elif (len(nitrogen_hydrogens) == 2):
          if (self.monomer.chem_comp.id == "PRO"):
            mod_mod_id = self.mon_lib_srv.mod_mod_id_dict["NH2"]
          else:
            mod_mod_id = self.mon_lib_srv.mod_mod_id_dict.get("NH2NOTPRO")
            if (mod_mod_id is None):
              raise RuntimeError("""\
A modified version of the monomer library is required to correctly
handle N-terminal hydrogens. The mod_NH2NOTPRO modification is missing.
This is a copy of mod_NH2, but without the HN2-N-CD angle.
Please contact cctbx@cci.lbl.gov for more information.""")
          mod_mod_ids.append(mod_mod_id)
          nitrogen_hydrogen_translation = ["HN1", "HN2"]
        elif (len(nitrogen_hydrogens) == 1):
          mod_mod_ids.append(self.mon_lib_srv.mod_mod_id_dict["NH1"])
          nitrogen_hydrogen_translation = ["HN"]
        if (nitrogen_hydrogen_translation is not None):
          j = 0
          for i,mon_lib_name in enumerate(self.mon_lib_names):
            if (not mon_lib_name in u): continue
            if (mon_lib_name in ["H1", "H2", "H3"]):
              self.mon_lib_names[i] = nitrogen_hydrogen_translation[j]
              j += 1
          assert j == len(nitrogen_hydrogen_translation)
      else:
        if (      ("H1" in u or "1H" in u)
              and ("H2" in u or "2H" in u)
              and ("H3" in u or "3H" in u)):
          if (self.monomer_atom_dict.get("H") is not None):
            mod_mod_ids.append(self.mon_lib_srv.mod_mod_id_dict["NH3"])
        elif (    ("HN1" in u or "1HN" in u)
              and ("HN2" in u or "2HN" in u)):
          mod_mod_ids.append(self.mon_lib_srv.mod_mod_id_dict["NH2"])
        elif (    "HN" in u):
          mod_mod_ids.append(self.mon_lib_srv.mod_mod_id_dict["NH1"])
        elif (    ("H1" in u or "1H" in u)
              or  ("H2" in u or "2H" in u)
              or  ("H3" in u or "3H" in u)):
          if (self.monomer_atom_dict.get("H") is not None):
            mod_mod_ids.append(self.mon_lib_srv.mod_mod_id_dict["NH3"])
        elif (    ("HN1" in u or "1HN" in u)
              or  ("HN2" in u or "2HN" in u)):
          mod_mod_ids.append(self.mon_lib_srv.mod_mod_id_dict["NH2"])
    elif (self.monomer.classification in ["RNA", "DNA"]):
      if (ani is not None):
        if (ani.have_op3_or_hop3):
          mod_mod_ids.append(self.mon_lib_srv.mod_mod_id_dict["p5*END"])
        elif (not ani.have_phosphate):
          mod_mod_ids.append(self.mon_lib_srv.mod_mod_id_dict["5*END"])
        if (ani.have_ho3prime):
          mod_mod_ids.append(self.mon_lib_srv.mod_mod_id_dict["3*END"])
      else:
        if ("O3T" in u):
          mod_mod_ids.append(self.mon_lib_srv.mod_mod_id_dict["p5*END"])
        else:
          e = self.expected_atoms
          if (    not "P"   in e
              and not "OP1" in e
              and not "OP2" in e):
            mod_mod_ids.append(self.mon_lib_srv.mod_mod_id_dict["5*END"])
        if ("HO3*" in u):
          mod_mod_ids.append(self.mon_lib_srv.mod_mod_id_dict["3*END"])
    for mod_mod_id in mod_mod_ids:
      self.apply_mod(mod_mod_id=mod_mod_id)

  def _track_mods(self, chem_mod_ids):
    for chem_mod_id in chem_mod_ids:
      self.chem_mod_ids.add(chem_mod_id)
      self.residue_name += "%" + chem_mod_id

  def apply_mod(self, mod_mod_id):
    if (mod_mod_id.chem_mod.id in self.chem_mod_ids):
      return
        # mod previously applied already, e.g. two links to same carbohydrate
    try:
      mod_mon = self.monomer.apply_mod(mod_mod_id)
    except Exception, e:
      import traceback
      msg = traceback.format_exc().splitlines()
      msg.extend([
        "apply_mod failure:",
        "  %s" % self.pdb_residue.id_str(),
        "  comp id: %s" % self.monomer.chem_comp.id,
        "  mod id: %s" % mod_mod_id.chem_mod.id])
      raise Sorry("\n".join(msg))
    self._track_mods(chem_mod_ids=[mod_mod_id.chem_mod.id])
    mod_mon.classification = self.monomer.classification
    self.monomer = mod_mon
    if (    mod_mod_id.chem_mod.name is not None
        and mod_mod_id.chem_mod.name.lower().find("terminus") >= 0):
      self.is_terminus = True # AD HOC manipulation
    self._get_mappings()

  def residue_altloc(self):
    result = self.residue_name
    if (self.conf_altloc != ""):
      result += ', conformer "%s"' % self.conf_altloc
    return result

  def is_unusual(self):
    m = self.monomer
    if (m is None): return True
    if (m.is_peptide()): return False
    if (m.is_rna_dna()): return False
    if (m.is_water()): return False
    return True

  def summary(self):
    if (self.monomer is None):
      classification = None
    else:
      classification = self.monomer.classification
    return monomer_mapping_summary(
      conf_altloc=self.conf_altloc,
      residue_name=self.residue_name,
      expected_atoms=self.expected_atoms.values(),
      unexpected_atoms=self.unexpected_atoms.values(),
      duplicate_atoms=flat_list(self.duplicate_atoms.values()),
      ignored_atoms=self.ignored_atoms.values(),
      classification=classification,
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

  def add_dihedral_proxies(self,
        dihedral_function_type,
        special_position_indices,
        dihedral_proxy_registry):
    self.dihedral_counters = add_dihedral_proxies(
      counters=counters(label="dihedral"),
      m_i=self,
      m_j=None,
      tor_list=self.monomer.tor_list,
      dihedral_function_type=dihedral_function_type,
      peptide_link_params=None,
      dihedral_proxy_registry=dihedral_proxy_registry,
      special_position_indices=special_position_indices).counters

  def add_chirality_proxies(self, special_position_indices,
                                  chirality_proxy_registry,
                                  chir_volume_esd):
    self.chirality_counters = add_chirality_proxies(
      counters=counters(label="chirality"),
      m_i=self,
      m_j=None,
      chir_list=self.monomer.chir_list,
      chirality_proxy_registry=chirality_proxy_registry,
      special_position_indices=special_position_indices,
      chir_volume_esd=chir_volume_esd).counters

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
    if (   chem_link_comp_id in [None, ""]
        or (chem_link_comp_id is not None
            and comp_id.lower() == chem_link_comp_id.lower())):
      self.is_comp_id_match = True
      if (chem_link_comp_id is None):
        self.len_comp_id_match = 0
      else:
        self.len_comp_id_match = len(chem_link_comp_id)
    else:
      self.is_comp_id_match = False
      self.len_comp_id_match = -1
    if (   chem_link_group_comp in [None, ""]
        or (comp_group is not None
            and comp_group.lower() == chem_link_group_comp.lower())):
      self.is_group_match = True
      if (chem_link_group_comp is None):
        self.len_group_match = 0
      else:
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

  def is_proper_match(self):
    return (
         self.len_comp_id_match_1 > 0
      or self.len_group_match_1 > 0
      or self.len_comp_id_match_2 > 0
      or self.len_group_match_2 > 0)

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
  if (m_j.expected_atoms.get("CN", None) is not None):
    link_id = "NM" + link_id
  elif (m_j.monomer.chem_comp.id == "PRO"):
    link_id = "P" + link_id
  return mon_lib_srv.link_link_id_dict[link_id]

def get_lib_link(mon_lib_srv, m_i, m_j):
  if (m_i.monomer.is_peptide() and m_j.monomer.is_peptide()):
    return get_lib_link_peptide(mon_lib_srv, m_i, m_j)
  elif (    (m_i.is_rna_dna or m_i.monomer.is_rna_dna())
        and (m_j.is_rna_dna or m_j.monomer.is_rna_dna())):
    if (m_i.is_rna2p):
      return mon_lib_srv.link_link_id_dict["rna2p"]
    return mon_lib_srv.link_link_id_dict["rna3p"]
  if (m_i.monomer.is_water() or m_j.monomer.is_water()): return None
  comp_id_1 = m_i.monomer.chem_comp.id
  comp_id_2 = m_j.monomer.chem_comp.id
  comp_1 = mon_lib_srv.get_comp_comp_id_direct(comp_id_1)
  comp_2 = mon_lib_srv.get_comp_comp_id_direct(comp_id_2)
  group_1 = comp_1.chem_comp.group
  group_2 = comp_2.chem_comp.group
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
      def get_atom(restr, comp_id, atom_id):
        ad = None
        if (comp_id == 1): ad = m_i.monomer_atom_dict
        if (comp_id == 2): ad = m_j.monomer_atom_dict
        if (ad is None):
          raise Sorry("""\
Corrupt CIF link definition:
  source info: %s
  link id: %s
  link name: %s
  %s atom atom_id: %s
  %s atom comp_id: %d (must be 1 or 2)""" % (
            str(match.link_link_id.source_info),
            chem_link.id,
            chem_link.name,
            restr, atom_id,
            restr, comp_id))
      match.n_unresolved_bonds = 0
      for bond in match.link_link_id.bond_list:
        atoms = [
          get_atom("bond", bond.atom_1_comp_id, bond.atom_id_1),
          get_atom("bond", bond.atom_2_comp_id, bond.atom_id_2)]
        if (None in atoms):
          match.n_unresolved_bonds += 1
      match.n_unresolved_angles = 0
      for angle in match.link_link_id.angle_list:
        atoms = [
          get_atom("angle", angle.atom_1_comp_id, angle.atom_id_1),
          get_atom("angle", angle.atom_2_comp_id, angle.atom_id_2),
          get_atom("angle", angle.atom_3_comp_id, angle.atom_id_3)]
        if (None in atoms):
          match.n_unresolved_angles += 1
      matches.append(match)
  if (len(matches) == 0): return None
  matches.sort()
  best_matches = []
  for m in matches:
    if (cmp(m, matches[0]) != 0): break
    best_matches.append(m)
  match = best_matches[0]
  if (not match.is_proper_match()):
    return None
  return match.link_link_id

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
  pdb_atoms = m_i.pdb_atoms
  atoms = [pdb_atoms[i_seq] for i_seq in i_seqs]
  if (not registry_process_result.is_new
      and not all_atoms_are_in_main_conf(atoms=atoms)):
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
    if (m_i.i_conformer != 0 and m_j.i_conformer != 0):
      assert m_i.i_conformer == m_j.i_conformer
    self.counters = counters
    self.broken_bond_i_seq_pairs = set()
    for bond in bond_list:
      if (   not m_i.monomer_atom_dict.has_key(bond.atom_id_1)
          or not m_j.monomer_atom_dict.has_key(bond.atom_id_2)):
        counters.corrupt_monomer_library_definitions += 1
        continue
      atoms = (m_i.expected_atoms.get(bond.atom_id_1, None),
               m_j.expected_atoms.get(bond.atom_id_2, None))
      if (None in atoms):
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
        i_seqs = [atom.i_seq for atom in atoms]
        proxy = geometry_restraints.bond_simple_proxy(
          i_seqs=i_seqs,
          distance_ideal=bond.value_dist,
          weight=1/bond.value_dist_esd**2)
        is_large_distance = False
        if (sites_cart is not None):
          r = geometry_restraints.bond(sites_cart=sites_cart, proxy=proxy)
          if (r.distance_model > distance_cutoff):
            is_large_distance = True
            self.broken_bond_i_seq_pairs.add(tuple(sorted(i_seqs)))
        if (not is_large_distance):
          if (    not m_i.is_first_conformer_in_chain
              and all_atoms_are_in_main_conf(atoms=atoms)):
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
        special_position_indices,
        broken_bond_i_seq_pairs=None):
    self.counters = counters
    if (m_j is None):
      m_1,m_2,m_3 = m_i,m_i,m_i
    elif (m_i.i_conformer != 0 and m_j.i_conformer != 0):
      assert m_i.i_conformer == m_j.i_conformer
    for angle in angle_list:
      if (m_j is not None):
        m_1,m_2,m_3 = [(m_i, m_j)[comp_id-1] for comp_id in (
          angle.atom_1_comp_id, angle.atom_2_comp_id, angle.atom_3_comp_id)]
      if (   not m_1.monomer_atom_dict.has_key(angle.atom_id_1)
          or not m_2.monomer_atom_dict.has_key(angle.atom_id_2)
          or not m_3.monomer_atom_dict.has_key(angle.atom_id_3)):
        counters.corrupt_monomer_library_definitions += 1
        continue
      atoms = (m_1.expected_atoms.get(angle.atom_id_1, None),
               m_2.expected_atoms.get(angle.atom_id_2, None),
               m_3.expected_atoms.get(angle.atom_id_3, None))
      if (None in atoms):
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
        i_seqs = [atom.i_seq for atom in atoms]
        if (involves_special_positions(special_position_indices, i_seqs)):
          counters.discarded_because_of_special_positions += 1
        elif (involves_broken_bonds(broken_bond_i_seq_pairs, i_seqs)):
          pass
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
        dihedral_function_type,
        peptide_link_params,
        dihedral_proxy_registry,
        special_position_indices,
        sites_cart=None,
        chem_link_id=None,
        broken_bond_i_seq_pairs=None):
    self.counters = counters
    self.chem_link_id = chem_link_id
    if (chem_link_id not in ["TRANS", "PTRANS", "NMTRANS"]):
      sites_cart = None
    if (m_j is None):
      m_1,m_2,m_3,m_4 = m_i,m_i,m_i,m_i
    elif (m_i.i_conformer != 0 and m_j.i_conformer != 0):
      assert m_i.i_conformer == m_j.i_conformer
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
      atoms = (m_1.expected_atoms.get(tor.atom_id_1, None),
               m_2.expected_atoms.get(tor.atom_id_2, None),
               m_3.expected_atoms.get(tor.atom_id_3, None),
               m_4.expected_atoms.get(tor.atom_id_4, None))
      if (None in atoms):
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
        i_seqs = [atom.i_seq for atom in atoms]
        if (involves_special_positions(special_position_indices, i_seqs)):
          counters.discarded_because_of_special_positions += 1
        elif (involves_broken_bonds(broken_bond_i_seq_pairs, i_seqs)):
          pass
        elif (    tor.id in ["psi", "phi"]
              and self.chem_link_id in ["TRANS", "PTRANS", "NMTRANS",
                                        "CIS",   "PCIS",   "NMCIS"]
              and peptide_link_params.discard_psi_phi):
          pass
        else:
          if (dihedral_function_type == "determined_by_sign_of_periodicity"):
            periodicity = tor.period
          elif (dihedral_function_type == "all_sinusoidal"):
            periodicity = max(1, tor.period)
          elif (dihedral_function_type == "all_harmonic"):
            periodicity = -abs(tor.period)
          else:
            raise RuntimeError(
              "Unknown dihedral_function_type: %s"
                % str(dihedral_function_type))
          try:
            if len(tor.alt_value_angle) == 0:
              alt_value_angle = None
            else:
              alt_value_angle = map(float,tor.alt_value_angle.split(","))
          except:
            alt_value_angle = None
          proxy = geometry_restraints.dihedral_proxy(
            i_seqs=i_seqs,
            angle_ideal=tor.value_angle,
            weight=1/tor.value_angle_esd**2,
            periodicity=periodicity, alt_angle_ideals=alt_value_angle)
          if (sites_cart is not None and tor.id == "omega"):
            assert abs(tor.value_angle - 180) < 1.e-6
            if (peptide_link_params.omega_esd_override_value is not None):
              assert peptide_link_params.omega_esd_override_value > 0
              proxy.weight = 1/peptide_link_params.omega_esd_override_value**2
            r = geometry_restraints.dihedral(
              sites_cart=sites_cart,
              proxy=proxy)
            if (abs(r.delta) > 180-peptide_link_params.cis_threshold):
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
        chir_volume_esd,
        lib_link=None,
        broken_bond_i_seq_pairs=None):
    self.counters = counters
    self.counters.unsupported_volume_sign = dicts.with_default_value(0)
    if (m_j is None):
      m_c,m_1,m_2,m_3 = m_i,m_i,m_i,m_i
    elif (m_i.i_conformer != 0 and m_j.i_conformer != 0):
      assert m_i.i_conformer == m_j.i_conformer
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
      atoms = (m_c.expected_atoms.get(chir.atom_id_centre, None),
               m_1.expected_atoms.get(chir.atom_id_1, None),
               m_2.expected_atoms.get(chir.atom_id_2, None),
               m_3.expected_atoms.get(chir.atom_id_3, None))
      if (None in atoms):
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
          i_seqs = [atom.i_seq for atom in atoms]
          if (involves_special_positions(special_position_indices, i_seqs)):
            counters.discarded_because_of_special_positions += 1
          elif (involves_broken_bonds(broken_bond_i_seq_pairs, i_seqs)):
            pass
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
        special_position_indices,
        broken_bond_i_seq_pairs=None):
    self.counters = counters
    self.counters.less_than_four_sites = dicts.with_default_value(0)
    if (    m_j is not None
        and m_i.i_conformer != 0 and m_j.i_conformer != 0):
      assert m_i.i_conformer == m_j.i_conformer
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
        atom = m_x.expected_atoms.get(plane_atom.atom_id, None)
        if (atom is None):
          if (m_x.monomer_atom_dict[plane_atom.atom_id].type_symbol == "H"):
            counters.unresolved_hydrogen += 1
          else:
            counters.unresolved_non_hydrogen += 1
            this_plane_has_unresolved_non_hydrogen = True
        elif (plane_atom.dist_esd in [None, 0]):
          counters.undefined += 1
        else:
          counters.resolved += 1
          i_seq = atom.i_seq
          if (special_position_indices is not None
              and i_seq in special_position_indices):
            counters.discarded_because_of_special_positions += 1
          else:
            i_seqs.append(i_seq)
            weights.append(1/plane_atom.dist_esd**2)
      if (len(i_seqs) < 4):
        if (this_plane_has_unresolved_non_hydrogen):
          counters.less_than_four_sites[plane.plane_id] += 1
      elif (involves_broken_bonds(broken_bond_i_seq_pairs, i_seqs)):
        pass
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
      assume_hydrogens_all_missing,
      factor_1_4_interactions,
      default_distance,
      minimum_distance):
  params = geometry_restraints.nonbonded_params(
    factor_1_4_interactions=factor_1_4_interactions,
    const_shrink_1_4_interactions=0,
    default_distance=default_distance,
    minimum_distance=minimum_distance)
  tables = {"": [], "h": []}
  for vdw in ener_lib.lib_vdw:
    assert vdw.H_flag in ["", "h"]
  if (assume_hydrogens_all_missing):
    reverse_prefs = ["", "h"]
  else:
    reverse_prefs = ["h", ""]
  for code in reverse_prefs:
    for vdw in tables[code]:
      atom_types = [vdw.atom_type_1, vdw.atom_type_2]
      atom_types.sort()
      params.distance_table.setdefault(
        atom_types[0])[atom_types[1]] = vdw.radius_min
  if (assume_hydrogens_all_missing):
    pref1, pref2 = ["vdwh_radius", "vdw_radius"]
  else:
    pref1, pref2 = ["vdw_radius", "vdwh_radius"]
  for atom_type,energy_lib_atom in ener_lib.lib_atom.items():
    if (len(atom_type) == 0): continue
    r = getattr(energy_lib_atom, pref1)
    if (r is None):
      r = getattr(energy_lib_atom, pref2)
    if (r is not None):
      params.radius_table[atom_type] = r
  return params

def is_same_model_as_before(model_type_indices, i_model, models):
  m_i = models[i_model]
  for j_model in xrange(0, i_model):
    if (model_type_indices[j_model] != j_model): continue
    if (m_i.is_identical_hierarchy(other=models[j_model])):
      model_type_indices[i_model] = j_model
      return True
  model_type_indices[i_model] = i_model
  return False

class build_chain_proxies(object):

  def __init__(self,
        mon_lib_srv,
        ener_lib,
        translate_cns_dna_rna_residue_names,
        rna_sugar_pucker_analysis_params,
        apply_cif_modifications,
        apply_cif_links_mm_pdbres_dict,
        link_distance_cutoff,
        not_linked_show_max,
        dihedral_function_type,
        chir_volume_esd,
        peptide_link_params,
        pdb_hierarchy,
        pdb_atoms,
        sites_cart,
        special_position_indices,
        keep_monomer_mappings,
        all_monomer_mappings,
        scattering_type_registry,
        nonbonded_energy_type_registry,
        geometry_proxy_registries,
        cystein_sulphur_i_seqs,
        cystein_monomer_mappings,
        is_unique_model,
        i_model,
        i_conformer,
        is_first_conformer_in_chain,
        conformer,
        conformation_dependent_restraints_list,
        log):
    self.conformation_dependent_restraints_list = \
      conformation_dependent_restraints_list
    unknown_residues = dicts.with_default_value(0)
    ad_hoc_single_atom_residues = dicts.with_default_value(0)
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
    mm_pairs_not_linked = []
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
    mm = None
    prev_mm = None
    prev_prev_mm = None
    pdb_residues = conformer.residues()
    for i_residue,residue in enumerate(pdb_residues):
      def _get_next_residue():
        j = i_residue + 1
        if (j == len(pdb_residues)): return None
        return pdb_residues[j]
      mm = monomer_mapping(
        pdb_atoms=pdb_atoms,
        mon_lib_srv=mon_lib_srv,
        translate_cns_dna_rna_residue_names
          =translate_cns_dna_rna_residue_names,
        rna_sugar_pucker_analysis_params=rna_sugar_pucker_analysis_params,
        apply_cif_modifications=apply_cif_modifications,
        apply_cif_links_mm_pdbres_dict=apply_cif_links_mm_pdbres_dict,
        i_model=i_model,
        i_conformer=i_conformer,
        is_first_conformer_in_chain=is_first_conformer_in_chain,
        conf_altloc=conformer.altloc,
        pdb_residue=residue,
        next_pdb_residue=_get_next_residue())
      if (mm.monomer is None):
        def use_scattering_type_if_available_to_define_nonbonded_type():
          if (   residue.atoms_size() != 1
              or len(mm.active_atoms) != 1): return False
          atom = mm.active_atoms[0]
          ad_hoc = ad_hoc_single_atom_residue(
            residue_name=residue.resname,
            atom_name=atom.name,
            atom_element=atom.element)
          if (ad_hoc.scattering_type is None): return False
          entry = ener_lib.lib_atom.get(ad_hoc.energy_type, None)
          if (entry is None): return False
          i_seq = atom.i_seq
          scattering_type_registry.assign_directly(
            i_seq=i_seq, symbol=ad_hoc.scattering_type)
          nonbonded_energy_type_registry.assign_directly(
            i_seq=i_seq, symbol=ad_hoc.energy_type)
          ad_hoc_single_atom_residues[mm.residue_name] += 1
          return True
        if (not use_scattering_type_if_available_to_define_nonbonded_type()):
          unknown_residues[mm.residue_name] += 1
        n_chain_breaks += 1
      elif (prev_mm is not None and not residue.link_to_previous):
        n_chain_breaks += 1
      else:
        if (prev_mm is not None and prev_mm.monomer is not None):
          prev_mm.lib_link = get_lib_link(
            mon_lib_srv=mon_lib_srv,
            m_i=prev_mm,
            m_j=mm)
          if (prev_mm.lib_link is None):
            link_ids[None] += 1
            mm_pairs_not_linked.append((prev_mm, mm))
          else:
            mod_id = prev_mm.lib_link.chem_link.mod_id_1
            if (mod_id not in [None, ""]):
              mod_mod_id = mon_lib_srv.mod_mod_id_dict[mod_id]
              prev_mm.apply_mod(mod_mod_id=mod_mod_id)
              prev_mm.resolve_unexpected()
            mod_id = prev_mm.lib_link.chem_link.mod_id_2
            if (mod_id not in [None, ""]):
              mod_mod_id = mon_lib_srv.mod_mod_id_dict[mod_id]
              mm.apply_mod(mod_mod_id=mod_mod_id)
              mm.resolve_unexpected()
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
            broken_bond_i_seq_pairs = link_resolution.broken_bond_i_seq_pairs
            n_chain_breaks += len(broken_bond_i_seq_pairs)
            link_resolution = add_angle_proxies(
              counters=counters(label="link_angle"),
              m_i=prev_mm,
              m_j=mm,
              angle_list=prev_mm.lib_link.angle_list,
              angle_proxy_registry=geometry_proxy_registries.angle,
              special_position_indices=special_position_indices,
              broken_bond_i_seq_pairs=broken_bond_i_seq_pairs)
            n_unresolved_chain_link_angles \
              += link_resolution.counters.unresolved_non_hydrogen
            link_resolution = add_dihedral_proxies(
              counters=counters(label="link_dihedral"),
              m_i=prev_mm,
              m_j=mm,
              tor_list=prev_mm.lib_link.tor_list,
              dihedral_function_type=dihedral_function_type,
              peptide_link_params=peptide_link_params,
              dihedral_proxy_registry=geometry_proxy_registries.dihedral,
              special_position_indices=special_position_indices,
              sites_cart=sites_cart,
              chem_link_id=prev_mm.lib_link.chem_link.id,
              broken_bond_i_seq_pairs=broken_bond_i_seq_pairs)
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
              chir_volume_esd=chir_volume_esd,
              lib_link=prev_mm.lib_link,
              broken_bond_i_seq_pairs=broken_bond_i_seq_pairs)
            n_unresolved_chain_link_chiralities \
              += link_resolution.counters.unresolved_non_hydrogen
            link_resolution = add_planarity_proxies(
              counters=counters(label="link_planarity"),
              m_i=prev_mm,
              m_j=mm,
              plane_list=prev_mm.lib_link.get_planes(),
              planarity_proxy_registry=geometry_proxy_registries.planarity,
              special_position_indices=special_position_indices,
              broken_bond_i_seq_pairs=broken_bond_i_seq_pairs)
            n_unresolved_chain_link_planarities \
              += link_resolution.counters.unresolved_non_hydrogen
      if (mm.monomer is not None):
        if (mm.is_unusual()):
          unusual_residues[mm.residue_name] += 1
        if (    mm.is_terminus == True
            and i_residue > 0
            and i_residue < conformer.residues_size()-1):
          inner_chain_residues_flagged_as_termini.append(residue.id_str())
        n_expected_atoms += len(mm.expected_atoms)
        for atom_name in mm.unexpected_atoms.keys():
          unexpected_atoms[mm.residue_name+","+atom_name] += 1
        for atom_name,i_seqs in mm.ignored_atoms.items():
          ignored_atoms[mm.residue_name+","+atom_name] += len(i_seqs)
        for atom_name,i_seqs in mm.duplicate_atoms.items():
          duplicate_atoms[mm.residue_name+","+atom_name] += len(i_seqs)
        if (mm.incomplete_info is not None):
          incomplete_infos[mm.incomplete_info] += 1
        if (mm.monomer.classification is not None):
          classifications[mm.monomer.classification] += 1
        for chem_mod_id in mm.chem_mod_ids:
          modifications_used[chem_mod_id] += 1
        scattering_type_registry.assign_from_monomer_mapping(
          conf_altloc=conformer.altloc, mm=mm)
        nonbonded_energy_type_registry.assign_from_monomer_mapping(
          conf_altloc=conformer.altloc, mm=mm)
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
          dihedral_function_type=dihedral_function_type,
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
          chirality_proxy_registry=geometry_proxy_registries.chirality,
          chir_volume_esd=chir_volume_esd)
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
          sulphur_atom = mm.expected_atoms.get("SG", None)
          # XXX keep track of weights
          if (sulphur_atom is not None
              and sulphur_atom.i_seq not in cystein_sulphur_i_seqs):
            cystein_sulphur_i_seqs.append(sulphur_atom.i_seq)
            cystein_monomer_mappings.append(mm)
      if (conformation_dependent_restraints.is_available):
        cdr = conformation_dependent_restraints \
                .build_conformation_dependent_angle_proxies(
          angle_proxy_registry=geometry_proxy_registries.angle,
          dihedral_proxy_registry=geometry_proxy_registries.dihedral,
          monomer_mappings=(prev_prev_mm, prev_mm, mm),
          connectivity_i_j=True,
          connectivity_j_k=True,
          sites_cart=sites_cart)
        self.conformation_dependent_restraints_list.append(cdr)
      if (keep_monomer_mappings):
        all_monomer_mappings.append(mm)
      else:
        all_monomer_mappings.append(mm.summary())
      prev_prev_mm = prev_mm
      prev_mm = mm
      prev_mm.lib_link = None
    #
    if (is_unique_model and log is not None):
      print >> log, "        Number of residues, atoms: %d, %d" % (
        conformer.residues_size(),
        n_expected_atoms + flex.sum(flex.long(unexpected_atoms.values())))
      if (len(unknown_residues) > 0):
        print >> log, "          Unknown residues:", unknown_residues
      if (len(ad_hoc_single_atom_residues) > 0):
        print >> log, "          Ad-hoc single atom residues:", \
          ad_hoc_single_atom_residues
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
        if (len(link_ids) != 1):
          if (not_linked_show_max is None):
            show_max = len(mm_pairs_not_linked)
          else:
            show_max = not_linked_show_max
          n_not_shown = max(0, len(mm_pairs_not_linked) - show_max)
          if (n_not_shown == 1):
            show_max += 1
            n_not_shown = 0
          for pair in mm_pairs_not_linked[:show_max]:
            print >> log, "            Not linked:"
            for mm in pair:
              print >> log, "              %s" % mm.pdb_residue.id_str()
          if (n_not_shown != 0):
            print >> log, \
              "            ... (remaining %d not shown)" % n_not_shown
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
        params=None,
        file_name=None,
        raw_records=None,
        pdb_inp=None,
        special_position_settings=None,
        crystal_symmetry=None,
        force_symmetry=False,
        substitute_non_crystallographic_unit_cell_if_necessary=False,
        strict_conflict_handling=True,
        keep_monomer_mappings=False,
        max_atoms=None,
        log=None,
        for_dihedral_reference=False):
    assert special_position_settings is None or crystal_symmetry is None
    if (params is None): params = master_params.extract()
    self.params = params
    if for_dihedral_reference:
      self.params.peptide_link.discard_psi_phi=False
    timer = user_plus_sys_time()
    self.time_building_chain_proxies = None
    if (log is not None and file_name is not None):
      print >> log, file_name
    if (pdb_inp is not None):
      assert raw_records is None
      self.pdb_inp = pdb_inp
    elif (file_name is not None):
      assert raw_records is None
      self.pdb_inp = pdb.input(file_name=file_name)
    else:
      if (isinstance(raw_records, str)):
        raw_records = flex.split_lines(raw_records)
      elif (not isinstance(raw_records, flex.std_string)):
        raw_records = flex.std_string(raw_records)
      self.pdb_inp = pdb.input(source_info=None, lines=raw_records)
    self.pdb_hierarchy = self.pdb_inp.construct_hierarchy()
    self.pdb_atoms = self.pdb_hierarchy.atoms()
    self.pdb_atoms.reset_i_seq()
    if (log is not None):
      print >> log, "  Monomer Library directory:"
      print >> log, "   ", show_string(mon_lib_srv.root_path)
      print >> log, "  Total number of atoms:", self.pdb_atoms.size()
    self.special_position_settings = None
    self._site_symmetry_table = None
    self.sites_cart = None
    self._sites_cart_exact = None
    if (max_atoms is not None
        and self.pdb_atoms.size() > max_atoms):
      if (log is not None):
        print >> log, "  More than %d atoms: no processing." % max_atoms
        return
    self.sites_cart = self.pdb_atoms.extract_xyz()
    models = self.pdb_hierarchy.models()
    if (log is not None):
      print >> log, "  Number of models:", len(models)
    n_seq = self.pdb_atoms.size()
    def set_model_indices():
      self.model_indices = flex.size_t(n_seq, n_seq)
      for i_model,model in enumerate(models):
        self.model_indices.set_selected(model.atoms().extract_i_seq(), i_model)
      assert self.model_indices.count(n_seq) == 0
    set_model_indices()
    altloc_i_conformer = {}
    def set_conformer_indices():
      self.conformer_indices = flex.size_t(n_seq, 0)
      altloc_indices = self.pdb_hierarchy.altloc_indices()
      if ("" in altloc_indices): p = 0
      else:                      p = 1
      altlocs = sorted(altloc_indices.keys())
      for i,altloc in enumerate(altlocs):
        if (altloc == ""): continue
        self.conformer_indices.set_selected(altloc_indices[altloc], i+p)
        altloc_i_conformer[altloc] = i+p
      altloc_i_conformer[""] = 0
    set_conformer_indices()
    sym_excl_residue_groups = []
    def set_sym_excl_indices():
      self.sym_excl_indices = flex.size_t(n_seq, 0)
      for rg in self.pdb_hierarchy.residue_groups():
        rg_atoms = rg.atoms()
        if (rg_atoms.extract_occ().all_lt(1.0)):
          sym_excl_residue_groups.append(rg)
          self.sym_excl_indices.set_selected(
            rg_atoms.extract_i_seq(), len(sym_excl_residue_groups))
    set_sym_excl_indices()
    #
    if (    special_position_settings is None
        and crystal_symmetry is not None):
      special_position_settings = crystal_symmetry.special_position_settings()
    self.special_position_settings = self.pdb_inp.special_position_settings(
      special_position_settings=special_position_settings,
      weak_symmetry=not force_symmetry)
    if (self.special_position_settings is not None
        and (   self.special_position_settings.unit_cell() is None
             or self.special_position_settings.space_group_info() is None)):
      self.special_position_settings = None
    if (    self.special_position_settings is None
        and substitute_non_crystallographic_unit_cell_if_necessary):
      self.special_position_settings = crystal.non_crystallographic_symmetry(
        sites_cart=self.sites_cart).special_position_settings()
    if (self.special_position_settings is None):
      self.special_position_indices = None
    else:
      self.special_position_indices = \
        self.site_symmetry_table().special_position_indices()
    self.process_apply_cif_modification(mon_lib_srv=mon_lib_srv, log=log)
    self.process_apply_cif_link(mon_lib_srv=mon_lib_srv, log=log)
    self.conformation_dependent_restraints_list = []
    model_type_indices = [-1] * len(models)
    self.all_monomer_mappings = []
    self.scattering_type_registry = scattering_type_registry(
      # XXX should be same as in pdb_inp.xray_structure_simple
      scattering_types=self.pdb_atoms.extract_element(strip=True),
      strict_conflict_handling=strict_conflict_handling)
    self.nonbonded_energy_type_registry = nonbonded_energy_type_registry(
      n_seq=self.pdb_atoms.size(),
      strict_conflict_handling=strict_conflict_handling)
    self.geometry_proxy_registries = geometry_restraints_proxy_registries(
      n_seq=self.pdb_atoms.size(),
      strict_conflict_handling=strict_conflict_handling)
    self.cystein_sulphur_i_seqs = flex.size_t()
    self.cystein_monomer_mappings = []
    n_unique_models = 0
    for i_model,model in enumerate(models):
      if (log is not None):
        print >> log, '  Model: "%s"' % model.id
      is_unique_model = not is_same_model_as_before(
        model_type_indices, i_model, models)
      if (is_unique_model):
        n_unique_models += 1
      elif (log is not None):
        print >> log, "    Same as model", \
          models[model_type_indices[i_model]].id
      if (is_unique_model and log is not None):
        print >> log, "    Number of chains:", model.chains_size()
      flush_log(log)
      self.geometry_proxy_registries.initialize_tables()
      apply_cif_links_mm_pdbres_dict = dict(
        self.empty_apply_cif_links_mm_pdbres_dict)
      for chain in model.chains():
        conformers = chain.conformers()
        if (is_unique_model and log is not None):
          print >> log, '    Chain: "%s"' % chain.id
          print >> log, "      Number of atoms:", chain.atoms_size()
          print >> log, "      Number of conformers:", len(conformers)
          flush_log(log)
        for j_conformer,conformer in enumerate(conformers):
          if (is_unique_model and log is not None):
            print >> log, '      Conformer: "%s"' % conformer.altloc
            flush_log(log)
          i_conformer = altloc_i_conformer[conformer.altloc]
          chain_proxies = build_chain_proxies(
            mon_lib_srv=mon_lib_srv,
            ener_lib=ener_lib,
            translate_cns_dna_rna_residue_names
              =self.params.translate_cns_dna_rna_residue_names,
            rna_sugar_pucker_analysis_params
              =self.params.rna_sugar_pucker_analysis,
            apply_cif_modifications=self.apply_cif_modifications,
            apply_cif_links_mm_pdbres_dict=apply_cif_links_mm_pdbres_dict,
            link_distance_cutoff=self.params.link_distance_cutoff,
            not_linked_show_max=self.params.show_max_items.not_linked,
            dihedral_function_type=self.params.dihedral_function_type,
            chir_volume_esd=self.params.chir_volume_esd,
            peptide_link_params=self.params.peptide_link,
            pdb_hierarchy=self.pdb_hierarchy,
            pdb_atoms=self.pdb_atoms,
            sites_cart=self.sites_cart,
            special_position_indices=self.special_position_indices,
            keep_monomer_mappings=keep_monomer_mappings,
            all_monomer_mappings=self.all_monomer_mappings,
            scattering_type_registry=self.scattering_type_registry,
            nonbonded_energy_type_registry=self.nonbonded_energy_type_registry,
            geometry_proxy_registries=self.geometry_proxy_registries,
            cystein_sulphur_i_seqs=self.cystein_sulphur_i_seqs,
            cystein_monomer_mappings=self.cystein_monomer_mappings,
            is_unique_model=is_unique_model,
            i_model=i_model,
            i_conformer=i_conformer,
            is_first_conformer_in_chain=(j_conformer == 0),
            conformer=conformer,
            conformation_dependent_restraints_list=
              self.conformation_dependent_restraints_list,
            log=log)
          self.conformation_dependent_restraints_list = \
            chain_proxies.conformation_dependent_restraints_list
          del chain_proxies
          flush_log(log)
      #
      n_unresolved_apply_cif_link_bonds = 0
      n_unresolved_apply_cif_link_angles = 0
      n_unresolved_apply_cif_link_dihedrals = 0
      n_unresolved_apply_cif_link_chiralities = 0
      n_unresolved_apply_cif_link_planarities = 0
      for apply in self.apply_cif_links:
        if (apply.was_used): continue
        mms = []
        for pdbres in apply.pdbres_pair:
          mms.append(apply_cif_links_mm_pdbres_dict[pdbres])
        for m_i_list in mms[0].values():
          assert len(m_i_list) == 1
          m_i = m_i_list[0]
          def raise_missing_cif(i_pair):
            raise Sorry(
              "Error processing apply_cif_link:\n"
              + "  data_link: %s\n" % show_string(apply.data_link)
              + "  Missing CIF file for residue: %s"
                  % apply.pdbres_pair[i_pair])
          if (m_i.monomer is None):
            raise_missing_cif(i_pair=0)
          for m_j_list in mms[1].values():
            assert len(m_j_list) == 1
            m_j = m_j_list[0]
            if (m_j.monomer is None):
              raise_missing_cif(i_pair=1)
            if (    m_i.i_conformer != 0
                and m_j.i_conformer != 0
                and m_i.i_conformer != m_j.i_conformer):
              continue
            apply.was_used = True
            def raise_if_corrupt(link_resolution):
              counters = link_resolution.counters
              if (counters.corrupt_monomer_library_definitions != 0):
                raise Sorry(
                  "Error processing %s:\n" % counters.label
                  + "  data_link: %s\n" % show_string(apply.data_link)
                  + "  residue 1: %s\n" % apply.pdbres_pair[0]
                  + "  residue 2: %s\n" % apply.pdbres_pair[1]
                  + "  Possible problems giving rise to this error include:\n"
                  + "    - residue selections in apply_cif_link block swapped\n"
                  + "    - corrupt CIF link definition\n"
                  + "    - corrupt or missing CIF modifications associated with"
                      " this link\n"
                  + "  If none of this applies, send email to:\n"
                  + "    bugs@phenix-online.org")
            link = mon_lib_srv.link_link_id_dict[apply.data_link]
            link_resolution = add_bond_proxies(
              counters=counters(label="apply_cif_link_bond"),
              m_i=m_i,
              m_j=m_j,
              bond_list=link.bond_list,
              bond_simple_proxy_registry=self.geometry_proxy_registries
                .bond_simple,
              sites_cart=self.sites_cart,
              distance_cutoff=self.params.link_distance_cutoff)
            raise_if_corrupt(link_resolution)
            n_unresolved_apply_cif_link_bonds \
              += link_resolution.counters.unresolved_non_hydrogen
            link_resolution = add_angle_proxies(
              counters=counters(label="apply_cif_link_angle"),
              m_i=m_i,
              m_j=m_j,
              angle_list=link.angle_list,
              angle_proxy_registry=self.geometry_proxy_registries.angle,
              special_position_indices=self.special_position_indices)
            raise_if_corrupt(link_resolution)
            n_unresolved_apply_cif_link_angles \
              += link_resolution.counters.unresolved_non_hydrogen
            link_resolution = add_dihedral_proxies(
              counters=counters(label="apply_cif_link_dihedral"),
              m_i=m_i,
              m_j=m_j,
              tor_list=link.tor_list,
              dihedral_function_type=self.params.dihedral_function_type,
              peptide_link_params=self.params.peptide_link,
              dihedral_proxy_registry=self.geometry_proxy_registries.dihedral,
              special_position_indices=self.special_position_indices,
              sites_cart=self.sites_cart,
              chem_link_id=link.chem_link.id)
            raise_if_corrupt(link_resolution)
            n_unresolved_apply_cif_link_dihedrals \
              += link_resolution.counters.unresolved_non_hydrogen
            link_resolution = add_chirality_proxies(
              counters=counters(label="apply_cif_link_chirality"),
              m_i=m_i,
              m_j=m_j,
              chir_list=link.chir_list,
              chirality_proxy_registry=self.geometry_proxy_registries.chirality,
              special_position_indices=self.special_position_indices,
              chir_volume_esd=self.params.chir_volume_esd,
              lib_link=link)
            raise_if_corrupt(link_resolution)
            n_unresolved_apply_cif_link_chiralities \
              += link_resolution.counters.unresolved_non_hydrogen
            link_resolution = add_planarity_proxies(
              counters=counters(label="apply_cif_link_planarity"),
              m_i=m_i,
              m_j=m_j,
              plane_list=link.get_planes(),
              planarity_proxy_registry=self.geometry_proxy_registries.planarity,
              special_position_indices=self.special_position_indices)
            raise_if_corrupt(link_resolution)
            n_unresolved_apply_cif_link_planarities \
              += link_resolution.counters.unresolved_non_hydrogen
      if (log is not None):
        if (n_unresolved_apply_cif_link_bonds > 0):
          print >> log, \
            "          Unresolved apply_cif_link bonds:", \
            n_unresolved_apply_cif_link_bonds
        if (n_unresolved_apply_cif_link_angles > 0):
          print >> log, \
            "          Unresolved apply_cif_link angles:", \
            n_unresolved_apply_cif_link_angles
        if (n_unresolved_apply_cif_link_dihedrals > 0):
          print >> log, \
            "          Unresolved apply_cif_link dihedrals:", \
            n_unresolved_apply_cif_link_dihedrals
        if (n_unresolved_apply_cif_link_chiralities > 0):
          print >> log, \
            "          Unresolved apply_cif_link chiralities:", \
            n_unresolved_apply_cif_link_chiralities
        if (n_unresolved_apply_cif_link_planarities > 0):
          print >> log, \
            "          Unresolved apply_cif_link planarities:", \
            n_unresolved_apply_cif_link_planarities
        flush_log(log)
    for apply in self.apply_cif_links:
      if (not apply.was_used):
        raise RuntimeError(
          "Unused apply_cif_link: %s %s" % (
            apply.data_link, str(apply.pdbres_pair)))
    self.geometry_proxy_registries.discard_tables()
    self.scattering_type_registry.discard_tables()
    self.nonbonded_energy_type_registry.discard_tables()
    if (log is not None):
      if (n_unique_models != 1):
        print >> log, "  Number of unique models:", n_unique_models
      if (len(sym_excl_residue_groups) != 0):
        print >> log, \
          "  Residues with excluded nonbonded symmetry interactions:", \
          len(sym_excl_residue_groups)
        show_residue_groups(
          residue_groups=sym_excl_residue_groups,
          log=log,
          prefix="    ",
          max_items=self.params.show_max_items
            .residues_with_excluded_nonbonded_symmetry_interactions)
      self.geometry_proxy_registries.report(log=log, prefix="  ")
      self.fatal_problems_report(prefix="  ", log=log)
    self.time_building_chain_proxies = timer.elapsed()

  def fatal_problems_report(self, prefix="", log=None):
    self.scattering_type_registry.report(
      pdb_atoms=self.pdb_atoms, log=log, prefix=prefix)
    self.nonbonded_energy_type_registry.report(
      pdb_atoms=self.pdb_atoms, log=log, prefix=prefix)

  def fatal_problems_message(self,
        ignore_unknown_scattering_types=False,
        ignore_unknown_nonbonded_energy_types=False):
    result = ["Fatal problems interpreting PDB file:"]
    for reg,ignore in [
          (self.scattering_type_registry,
            ignore_unknown_scattering_types),
          (self.nonbonded_energy_type_registry,
             ignore_unknown_nonbonded_energy_types)]:
      if (not ignore):
        n_unknown = reg.n_unknown_type_symbols()
        if (n_unknown != 0):
          result.append("%s: %d" % (reg.report_unknown_message(), n_unknown))
    if (len(result) == 1): return None
    result.extend([
      "  Please edit the PDB file to resolve the problems and/or supply a",
      "  CIF file with matching restraint definitions, along with",
      "  apply_cif_modification and apply_cif_link parameter definitions",
      "  if necessary."])
    if (self.scattering_type_registry.n_unknown_type_symbols() != 0):
      result.extend([
        "  It is best practice to define the element names in",
        "  columns 77-78 of the PDB file."])
    result.extend([
      "  Also note that phenix.ready_set and phenix.elbow are available",
      "  for creating restraint definitions (CIF files)."])
    return "\n  ".join(result)

  def extract_secondary_structure (self) :
    return self.pdb_inp.extract_secondary_structure()

  def site_symmetry_table(self):
    if (self._site_symmetry_table is None):
      assert self.special_position_settings is not None
      assert self.sites_cart is not None
      self._site_symmetry_table = \
        self.special_position_settings.site_symmetry_table(
          sites_cart=self.sites_cart,
          unconditional_general_position_flags=(
            self.pdb_atoms.extract_occ() != 1))
    return self._site_symmetry_table

  def sites_cart_exact(self):
    if (self._sites_cart_exact is None):
      self._sites_cart_exact = self.site_symmetry_table().apply_symmetry_sites(
        unit_cell=self.special_position_settings.unit_cell(),
        sites_cart=self.sites_cart)
    return self._sites_cart_exact

  def sel_classification(self, classification):
    result = flex.bool(self.pdb_atoms.size(), False)
    for summary in self.all_monomer_mappings:
      if (summary.classification == classification):
        result.set_selected(summary.all_associated_i_seqs(), True)
    return result

  def sel_backbone_or_sidechain(self, backbone_flag, sidechain_flag):
    result = flex.bool(self.pdb_atoms.size(), False)
    for summary in self.all_monomer_mappings:
      if (summary.classification == "peptide"):
        for atom in summary.expected_atoms:
          # XXX hydrogens not included
          if (atom.name.strip() in ["N", "CA", "C", "O"]):
            result[atom.i_seq] = backbone_flag
          else:
            result[atom.i_seq] = sidechain_flag
      elif (summary.classification in ["RNA", "DNA"]):
        for atom in summary.expected_atoms:
          # XXX hydrogens not included
          if (atom.name.strip()
                in ["P", "O1P", "O2P", "O3'", "O5'",
                         "OP1", "OP2", "O3*", "O5*",
                                       "O2*", "O2'",
                    "O4'", "C1'", "C2'", "C3'", "C4'", "C5'",
                    "O4*", "C1*", "C2*", "C3*", "C4*", "C5*"]):
            result[atom.i_seq] = backbone_flag
          else:
            result[atom.i_seq] = sidechain_flag
    return result

  def sel_backbone(self):
    return self.sel_backbone_or_sidechain(True, False)

  def sel_sidechain(self):
    return self.sel_backbone_or_sidechain(False, True)

  def sel_phosphate(self):
    result = flex.bool(self.pdb_atoms.size(), False)
    for summary in self.all_monomer_mappings:
      if (summary.classification in ["RNA", "DNA"]):
        for atom in summary.expected_atoms:
          if (atom.name.strip()
                in ["P", "O1P", "O2P", "O3'", "O5'",
                         "OP1", "OP2", "O3*", "O5*"]):
            result[atom.i_seq] = True
    return result

  def sel_ribose(self):
    result = flex.bool(self.pdb_atoms.size(), False)
    for summary in self.all_monomer_mappings:
      if (summary.classification in ["RNA", "DNA"]):
        for atom in summary.expected_atoms:
          if (atom.name.strip()
                in ["O4'", "C1'", "C2'", "C3'", "C4'", "C5'", "O2'",
                    "O4*", "C1*", "C2*", "C3*", "C4*", "C5*", "O2*"]):
            result[atom.i_seq] = True
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
      sel = self.pdb_hierarchy.atom_selection_cache().selection_parser(
        word_iterator=word_iterator,
        callback=self._selection_callback,
        expect_nonmatching_closing_parenthesis=True)
      result_stack.append(self.sel_within(radius=radius,primary_selection=sel))
    else:
      return False
    return True

  def selection(self, string, cache=None):
    if (cache is None): cache = self.pdb_hierarchy.atom_selection_cache()
    return cache.selection(
      string=string,
      callback=self._selection_callback)

  def iselection(self, string, cache=None):
    return self.selection(string=string, cache=cache).iselection()

  def process_apply_cif_modification(self, mon_lib_srv, log):
    self.apply_cif_modifications = {}
    atoms = self.pdb_atoms
    sel_cache = None
    for apply in self.params.apply_cif_modification:
      if (apply.data_mod is None): continue
      print >> log, "  apply_cif_modification:"
      print >> log, "    data_mod:", apply.data_mod
      mod = mon_lib_srv.mod_mod_id_dict.get(apply.data_mod)
      if (mod is None):
        print >> log
        raise Sorry(
          "Missing CIF modification: data_mod_%s\n" % apply.data_mod
          + "  Please check for spelling errors or specify the file name\n"
          + "  with the modification as an additional argument.")
      print >> log, "    residue_selection:", apply.residue_selection
      if (sel_cache is None):
        sel_cache = self.pdb_hierarchy.atom_selection_cache()
      iselection = self.phil_atom_selection(
        cache=sel_cache,
        scope_extract=apply,
        attr="residue_selection").iselection()
      pdbres_set = set()
      for i_seq in iselection:
        pdbres_set.add(atoms[i_seq].id_str(pdbres=True))
      for pdbres in pdbres_set:
        self.apply_cif_modifications.setdefault(pdbres, []).append(
          apply.data_mod)

  def process_apply_cif_link(self, mon_lib_srv, log):
    self.apply_cif_links = []
    self.empty_apply_cif_links_mm_pdbres_dict = {}
    atoms = self.pdb_atoms
    sel_cache = None
    for apply in self.params.apply_cif_link:
      if (apply.data_link is None): continue
      print >> log, "  apply_cif_link:"
      print >> log, "    data_link:", apply.data_link
      link = mon_lib_srv.link_link_id_dict.get(apply.data_link)
      if (link is None):
        print >> log
        raise Sorry(
          "Missing CIF link: data_link_%s\n" % apply.data_link
          + "  Please check for spelling errors or specify the file name\n"
          + "  with the link as an additional argument.")
      mod_ids = []
      for mod_attr in ["mod_id_1", "mod_id_2"]:
        mod_id = getattr(link.chem_link, mod_attr)
        if (mod_id == ""): mod_id = None
        mod_ids.append(mod_id)
        if (mod_id is not None):
          print >> log, "      %s:" % mod_attr, mod_id
          mod = mon_lib_srv.mod_mod_id_dict.get(mod_id)
          if (mod is None):
            print >> log
            raise Sorry(
              "Missing CIF modification: data_mod_%s\n" % mod_id
              + "  Please check for spelling errors or specify the file name\n"
              + "  with the modification as an additional argument.")
      sel_attrs = ["residue_selection_"+n for n in ["1", "2"]]
      pdbres_pair = []
      for attr in sel_attrs:
        print >> log, "    %s:" % attr, getattr(apply, attr)
        if (sel_cache is None):
          sel_cache = self.pdb_hierarchy.atom_selection_cache()
        iselection = self.phil_atom_selection(
          cache=sel_cache,
          scope_extract=apply,
          attr=attr).iselection()
        pdbres_set = set()
        pdbres_set_no_segid = set()
        for i_seq in iselection:
          atom = atoms[i_seq]
          pdbres_set.add(atom.id_str(pdbres=True))
          pdbres_set_no_segid.add(atom.id_str(pdbres=True,suppress_segid=True))
        if (len(pdbres_set) != 1):
          # XXX models?
          if (len(pdbres_set_no_segid) != 1):
            raise Sorry("Not exactly one residue selected.")
          raise Sorry(
            "Selected residue has multiple segid. This is not supported.")
        pdbres_pair.append(list(pdbres_set)[0])
      for pdbres,mod_id in zip(pdbres_pair, mod_ids):
        if (mod_id is not None):
          self.apply_cif_modifications.setdefault(pdbres, []).append(mod_id)
      self.apply_cif_links.append(group_args(
        pdbres_pair=pdbres_pair,
        data_link=apply.data_link,
        was_used=False))
      for pdbres in pdbres_pair:
        self.empty_apply_cif_links_mm_pdbres_dict[pdbres] = {}

  def create_disulfides(self, disulfide_distance_cutoff, log=None):
    if (self.model_indices is not None):
      model_indices = self.model_indices.select(self.cystein_sulphur_i_seqs)
    conformer_indices = self.conformer_indices.select(
      self.cystein_sulphur_i_seqs)
    sym_excl_indices = self.sym_excl_indices.select(
      self.cystein_sulphur_i_seqs)
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
      sym_excl_indices=sym_excl_indices,
      nonbonded_params=geometry_restraints.nonbonded_params(
        default_distance=1),
      nonbonded_types=flex.std_string(conformer_indices.size()),
      nonbonded_distance_cutoff_plus_buffer=disulfide_distance_cutoff,
      min_cubicle_edge=5,
      shell_asu_tables=[pair_asu_table])
    labels = [self.pdb_atoms[i_seq].id_str()
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

  def atom_selection(self, parameter_name, string, cache=None):
    try:
      return self.selection(string=string, cache=cache)
    except KeyboardInterrupt: raise
    except Exception, e: # keep e alive to avoid traceback
      fe = format_exception()
      raise Sorry('Invalid atom selection:\n  %s="%s"\n  (%s)' % (
        parameter_name, string, fe))

  def phil_atom_selection(self,
        cache,
        scope_extract,
        attr,
        allow_none=False,
        allow_auto=False,
        raise_if_empty_selection=True):
    def parameter_name():
      return scope_extract.__phil_path__(object_name=attr)
    string = getattr(scope_extract, attr)
    if (string is None):
      if (allow_none): return None
      raise Sorry('Atom selection cannot be None:\n  %s=None' % (
        parameter_name()))
    elif (string is Auto):
      if (allow_auto): return Auto
      raise Sorry('Atom selection cannot be Auto:\n  %s=Auto' % (
        parameter_name()))
    try: result = self.selection(string=string, cache=cache)
    except KeyboardInterrupt: raise
    except Exception, e: # keep e alive to avoid traceback
      fe = format_exception()
      raise Sorry('Invalid atom selection:\n  %s=%s\n  (%s)' % (
        parameter_name(), show_string(string), fe))
    if (raise_if_empty_selection and result.count(True) == 0):
      raise Sorry('Empty atom selection:\n  %s=%s' % (
        parameter_name(), show_string(string)))
    return result

  def phil_atom_selection_multiple(self,
        cache,
        scope_extract,
        attr,
        allow_none=False,
        allow_auto=False,
        raise_if_empty_selection=True):
    result = []
    def parameter_name():
      return scope_extract.__phil_path__(object_name=attr)
    string_list = getattr(scope_extract, attr)
    for string in string_list:
      if (string is None):
        if (allow_none): return None
        raise Sorry('Atom selection cannot be None:\n  %s=None' % (
          parameter_name()))
      elif (string is Auto):
        if (allow_auto): return Auto
        raise Sorry('Atom selection cannot be Auto:\n  %s=Auto' % (
          parameter_name()))
      try:
          result.append(self.selection(string=string, cache=cache).iselection())
      except KeyboardInterrupt: raise
      except Exception, e: # keep e alive to avoid traceback
        fe = format_exception()
        raise Sorry('Invalid atom selection:\n  %s=%s\n  (%s)' % (
          parameter_name(), show_string(string), fe))
      if (raise_if_empty_selection and result.count(True) == 0):
        raise Sorry('Empty atom selection:\n  %s=%s' % (
          parameter_name(), show_string(string)))
    return result

  def phil_atom_selections_as_i_seqs(self, cache, scope_extract, sel_attrs):
    result = []
    for attr in sel_attrs:
      iselection = self.phil_atom_selection(
        cache=cache,
        scope_extract=scope_extract,
        attr=attr,
        raise_if_empty_selection=False).iselection()
      if (iselection.size() != 1):
        atom_sel = getattr(scope_extract, attr)
        if (iselection.size() == 0):
          raise Sorry("No atom selected: %s" % show_string(atom_sel))
        else:
          raise Sorry(
            "More than one atom selected: %s\n"
            "  Number of selected atoms: %d" % (
              show_string(atom_sel), iselection.size()))
      result.append(iselection[0])
    return result

  def phil_atom_selections_as_i_seqs_multiple(self,
                                              cache,
                                              scope_extract,
                                              sel_attrs):
    result = []
    for attr in sel_attrs:
        iselection = self.phil_atom_selection_multiple(
          cache=cache,
          scope_extract=scope_extract,
          attr=attr,
          raise_if_empty_selection=False)
        atom_sel = getattr(scope_extract, attr)
        for i in iselection:
          if (i.size() == 0):
            raise Sorry("No atom selected: %s" % show_string(atom_sel))
          for atom in i:
            result.append(atom)
    return result

  def process_geometry_restraints_remove(self,
        params,
        geometry_restraints_manager):
    path = params.__phil_path__() + "."
    grm = geometry_restraints_manager
    cache = self.pdb_hierarchy.atom_selection_cache()
    for atom_selection in params.angles:
      grm.remove_angles_in_place(selection=self.atom_selection(
        parameter_name=path+"angles", string=atom_selection, cache=cache))
    for atom_selection in params.dihedrals:
      grm.remove_dihedrals_in_place(selection=self.atom_selection(
        parameter_name=path+"dihedrals", string=atom_selection, cache=cache))
    for atom_selection in params.chiralities:
      grm.remove_chiralities_in_place(selection=self.atom_selection(
        parameter_name=path+"chiralities", string=atom_selection, cache=cache))
    for atom_selection in params.planarities:
      grm.remove_planarities_in_place(selection=self.atom_selection(
        parameter_name=path+"planarities", string=atom_selection, cache=cache))

  def process_geometry_restraints_edits_bond(self, sel_cache, params, log):
    bond_sym_proxies = []
    bond_distance_model_max = 0
    if (len(params.bond) == 0):
      return group_args(
        bond_sym_proxies=bond_sym_proxies,
        bond_distance_model_max=bond_distance_model_max)
    print >> log, "  Custom bonds:"
    atoms = self.pdb_atoms
    unit_cell = self.special_position_settings.unit_cell()
    space_group = self.special_position_settings.space_group()
    uc_shortest_vector = unit_cell.shortest_vector_sq()**0.5
    max_bond_length = uc_shortest_vector
    ebdl = params.excessive_bond_distance_limit
    if (ebdl not in [None, Auto] and ebdl > 0 and ebdl < max_bond_length):
      max_bond_length = ebdl
    n_excessive = 0
    sel_attrs = ["atom_selection_"+n for n in ["1", "2"]]
    for bond in params.bond:
      def show_atom_selections():
        for attr in sel_attrs:
          print >> log, "      %s = %s" % (
            attr, show_string(getattr(bond, attr, None)))
      slack = bond.slack
      if (slack is None or slack < 0):
        slack = 0
      if (bond.distance_ideal is None):
        print >> log, "    Warning: Ignoring bond with distance_ideal = None:"
        show_atom_selections()
      elif (bond.distance_ideal < 0):
        print >> log, "    Warning: Ignoring bond with distance_ideal < 0:"
        show_atom_selections()
        print >> log, "      distance_ideal = %.6g" % bond.distance_ideal
      elif (bond.sigma is None):
        print >> log, "    Warning: Ignoring bond with sigma = None:"
        show_atom_selections()
        print >> log, "      distance_ideal = %.6g" % bond.distance_ideal
      elif (bond.sigma <= 0):
        print >> log, "    Warning: Ignoring bond with sigma <= 0:"
        show_atom_selections()
        print >> log, "      distance_ideal = %.6g" % bond.distance_ideal
        print >> log, "      sigma = %.6g" % bond.sigma
        print >> log, "      slack = %.6g" % slack
      elif (bond.action != "add"):
        raise Sorry("%s = %s not implemented." %
          bond.__phil_path_and_value__(object_name="action"))
      else:
        i_seqs = self.phil_atom_selections_as_i_seqs(
          cache=sel_cache, scope_extract=bond, sel_attrs=sel_attrs)
        if (bond.symmetry_operation is None):
          s = "x,y,z"
        else:
          s = bond.symmetry_operation
        rt_mx_ji = sgtbx.rt_mx(symbol=s, t_den=space_group.t_den())
        p = geometry_restraints.bond_sym_proxy(
          i_seqs=i_seqs,
          distance_ideal=bond.distance_ideal,
          weight=geometry_restraints.sigma_as_weight(sigma=bond.sigma),
          slack=slack,
          rt_mx_ji=rt_mx_ji)
        bond_sym_proxies.append(p)
        b = geometry_restraints.bond(
          unit_cell=unit_cell,
          sites_cart=self.sites_cart,
          proxy=p)
        print >> log, "    bond:"
        for i in [0,1]:
          print >> log, "      atom %d:" % (i+1), atoms[p.i_seqs[i]].quote()
        print >> log, "      symmetry operation:", str(p.rt_mx_ji)
        if (not space_group.contains(smx=p.rt_mx_ji)):
          raise Sorry(
            'The bond symmetry operation "%s" is not compatible'
            ' with space group %s.' % (
              str(p.rt_mx_ji),
              self.special_position_settings.space_group_info()
                .symbol_and_number()))
        print >> log, "      distance_model: %7.3f" % b.distance_model
        print >> log, "      distance_ideal: %7.3f" % b.distance_ideal
        print >> log, "      ideal - model:  %7.3f" % b.delta
        print >> log, "      slack:          %7.3f" % b.slack
        print >> log, "      delta_slack:    %7.3f" % b.delta_slack
        print >> log, "      sigma:          %8.4f" % \
          geometry_restraints.weight_as_sigma(weight=b.weight)
        if (bond_distance_model_max < b.distance_model):
          bond_distance_model_max = b.distance_model
        if (b.distance_model > max_bond_length):
          print >> log, "      *** WARNING: EXCESSIVE BOND LENGTH. ***"
          n_excessive += 1
    if (n_excessive != 0):
      if (max_bond_length == uc_shortest_vector):
        print >> log, "  Excessive bond length limit at hard upper bound:" \
          " length of shortest vector between unit cell lattice points: %.6g" \
            % uc_shortest_vector
      else:
        print >> log, "  %s = %.6g" % \
          params.__phil_path_and_value__("excessive_bond_distance_limit")
        print >> log, \
          "    Please assign a larger value to this parameter if necessary."
      raise Sorry(
        "Custom bonds with excessive length: %d\n"
        "  Please check the log file for details." % n_excessive)
    print >> log, "    Total number of custom bonds:", len(bond_sym_proxies)
    return group_args(
      bond_sym_proxies=bond_sym_proxies,
      bond_distance_model_max=bond_distance_model_max)

  def process_geometry_restraints_edits_angle(self, sel_cache, params, log):
    result = []
    if (len(params.angle) == 0): return result
    if (self.special_position_indices is None):
      special_position_indices = []
    else:
      special_position_indices = self.special_position_indices
    print >> log, "  Custom angles:"
    atoms = self.pdb_atoms
    sel_attrs = ["atom_selection_"+n for n in ["1", "2", "3"]]
    for angle in params.angle:
      def show_atom_selections():
        for attr in sel_attrs:
          print >> log, "      %s = %s" % (
            attr, show_string(getattr(angle, attr, None)))
      if (angle.angle_ideal is None):
        print >> log, "    Warning: Ignoring angle with angle_ideal = None:"
        show_atom_selections()
      elif (angle.sigma is None):
        print >> log, "    Warning: Ignoring angle with sigma = None:"
        show_atom_selections()
        print >> log, "      angle_ideal = %.6g" % angle.angle_ideal
      elif (angle.sigma is None or angle.sigma <= 0):
        print >> log, "    Warning: Ignoring angle with sigma <= 0:"
        show_atom_selections()
        print >> log, "      angle_ideal = %.6g" % angle.angle_ideal
        print >> log, "      sigma = %.6g" % angle.sigma
      elif (angle.action != "add"):
        raise Sorry("%s = %s not implemented." %
          angle.__phil_path_and_value__("action"))
      else:
        i_seqs = self.phil_atom_selections_as_i_seqs(
          cache=sel_cache, scope_extract=angle, sel_attrs=sel_attrs)
        p = geometry_restraints.angle_proxy(
          i_seqs=i_seqs,
          angle_ideal=angle.angle_ideal,
          weight=geometry_restraints.sigma_as_weight(sigma=angle.sigma))
        a = geometry_restraints.angle(
          sites_cart=self.sites_cart,
          proxy=p)
        print >> log, "    angle:"
        n_special = 0
        for i,i_seq in enumerate(p.i_seqs):
          print >> log, "      atom %d:" % (i+1), atoms[i_seq].quote(),
          if (i_seq in special_position_indices):
            n_special += 1
            print >> log, "# SPECIAL POSITION",
          print >> log
        print >> log, "      angle_model: %7.2f" % a.angle_model
        print >> log, "      angle_ideal: %7.2f" % a.angle_ideal
        print >> log, "      ideal - model:  %7.2f" % a.delta
        print >> log, "      sigma: %.6g" % \
          geometry_restraints.weight_as_sigma(weight=a.weight)
        if (n_special != 0):
          raise Sorry(
            "Custom angle involves %d special position%s:\n"
            "  Please inspect the output for details."
              % plural_s(n_special))
        result.append(p)
    print >> log, "    Total number of custom angles:", len(result)
    return result

  def process_geometry_restraints_edits_planarity(self,
                                                  sel_cache,
                                                  params,
                                                  log):
    result = []
    if (len(params.planarity) == 0): return result
    if (self.special_position_indices is None):
      special_position_indices = []
    else:
      special_position_indices = self.special_position_indices
    sel_attrs = ["atom_selection"]
    print >> log, "  Custom planarities:"
    for planarity in params.planarity:
       i_seqs = self.phil_atom_selections_as_i_seqs_multiple(
         cache=sel_cache, scope_extract=planarity, sel_attrs=sel_attrs)
       weights = []
       for i_seq in i_seqs:
         weights.append(geometry_restraints.sigma_as_weight(sigma=planarity.sigma))
       proxy = geometry_restraints.planarity_proxy(
         i_seqs=i_seqs,
         weights=flex.double(weights))
       plane = geometry_restraints.planarity(
         sites_cart=self.sites_cart,
         proxy=proxy)
       result.append(proxy)
    print >> log, "    Total number of custom planarities:", len(result)
    return result

  def process_geometry_restraints_edits(self, params, log):
    sel_cache = self.pdb_hierarchy.atom_selection_cache()
    result = self.process_geometry_restraints_edits_bond(
        sel_cache=sel_cache, params=params, log=log)
    result.angle_proxies=self.process_geometry_restraints_edits_angle(
        sel_cache=sel_cache, params=params, log=log)
    result.planarity_proxies=self.process_geometry_restraints_edits_planarity(
        sel_cache=sel_cache, params=params, log=log)
    return result

  def process_hydrogen_bonds (self, bonds_table, log, verbose=False) :
    atoms = self.pdb_atoms
    def show_atoms (i_seqs, log) :
      for i_seq in i_seqs :
        print >> log, "     %s" % atoms[i_seq].fetch_labels().quote()
    unit_cell = self.special_position_settings.unit_cell()
    space_group = self.special_position_settings.space_group()
    bond_sym_proxies = []
    for bond in bonds_table.get_bond_restraint_data() :
      i_seqs = [bond.donor_i_seq, bond.acceptor_i_seq]
      slack = bond.slack
      if (slack is None or slack < 0):
        slack = 0
      if (bond.distance_ideal is None):
        print >> log, "    Warning: Ignoring bond with distance_ideal = None:"
        show_atoms(i_seqs, log)
      elif (bond.distance_ideal <= 0):
        print >> log, "    Warning: Ignoring bond with distance_ideal <= 0:"
        show_atoms(i_seqs, log)
        print >> log, "      distance_ideal = %.6g" % bond.distance_ideal
      elif (bond.sigma is None):
        print >> log, "    Warning: Ignoring bond with sigma = None:"
        show_atoms(i_seqs, log)
        print >> log, "      distance_ideal = %.6g" % bond.distance_ideal
      elif (bond.sigma <= 0):
        print >> log, "    Warning: Ignoring bond with sigma <= 0:"
        show_atoms(i_seqs, log)
        print >> log, "      distance_ideal = %.6g" % bond.distance_ideal
        print >> log, "      sigma = %.6g" % bond.sigma
        print >> log, "      slack = %.6g" % slack
      else:
        rt_mx_ji = sgtbx.rt_mx(symbol="x,y,z", t_den=space_group.t_den())
        p = geometry_restraints.bond_sym_proxy(
          i_seqs=i_seqs,
          distance_ideal=bond.distance_ideal,
          weight=geometry_restraints.sigma_as_weight(sigma=bond.sigma),
          slack=slack,
          rt_mx_ji=rt_mx_ji)
        bond_sym_proxies.append(p)
        b = geometry_restraints.bond(
          unit_cell=unit_cell,
          sites_cart=self.sites_cart,
          proxy=p)
        if verbose :
          print >> log, "    hydrogen bond:"
          for i in [0,1]:
            print >> log, "      atom %d:" % (i+1), atoms[p.i_seqs[i]].quote()
          print >> log, "      distance_model: %7.3f" % b.distance_model
          print >> log, "      distance_ideal: %7.3f" % b.distance_ideal
          print >> log, "      ideal - model:  %7.3f" % b.delta
          print >> log, "      slack:          %7.3f" % b.slack
          print >> log, "      delta_slack:    %7.3f" % b.delta_slack
          print >> log, "      sigma:          %8.4f" % \
            geometry_restraints.weight_as_sigma(weight=b.weight)
    print >> log, "  Total number of hydrogen bonds:", len(bond_sym_proxies)
    return bond_sym_proxies

  def construct_geometry_restraints_manager(self,
        ener_lib,
        disulfide_link,
        plain_pairs_radius=None,
        params_edits=None,
        params_remove=None,
        h_bond_table=None,
        assume_hydrogens_all_missing=True,
        log=None):
    assert self.special_position_settings is not None
    timer = user_plus_sys_time()
    bond_params_table = geometry_restraints.extract_bond_params(
      n_seq=self.sites_cart.size(),
      bond_simple_proxies=self.geometry_proxy_registries.bond_simple.proxies)
    bond_distances_model = geometry_restraints.bond_distances_model(
      sites_cart=self.sites_cart,
      proxies=self.geometry_proxy_registries.bond_simple.proxies)
    if (bond_distances_model.size() > 0):
      excessive_bonds = (
        bond_distances_model > self.special_position_settings.unit_cell()
          .shortest_vector_sq()**.5).iselection()
      if (excessive_bonds.size() > 0 and not
          self.params.proceed_with_excessive_length_bonds):
        atoms = self.pdb_atoms
        proxies = self.geometry_proxy_registries.bond_simple.proxies
        print >> log, "  Bonds with excessive lengths:"
        for i_proxy in excessive_bonds:
          proxy = proxies[i_proxy]
          bond = geometry_restraints.bond(
            sites_cart=self.sites_cart, proxy=proxy)
          print >> log, "    Distance model: %.6g (ideal: %.6g)" % (
            bond.distance_model, bond.distance_ideal)
          for i_seq in proxy.i_seqs:
            print >> log, "      %s" % atoms[i_seq].format_atom_record()
        raise Sorry("Number of bonds with excessive lengths: %d" %
          excessive_bonds.size())
    disulfide_sym_table, max_disulfide_bond_distance = \
      self.create_disulfides(
        disulfide_distance_cutoff=self.params.disulfide_distance_cutoff,
        log=log)
    max_bond_distance = max_disulfide_bond_distance
    if (bond_distances_model.size() > 0):
      max_bond_distance = max(max_bond_distance,
        flex.max(bond_distances_model))
    if (params_edits is None):
      processed_edits = None
    else:
      processed_edits = self.process_geometry_restraints_edits(
        params=params_edits, log=log)
      max_bond_distance = max(max_bond_distance,
        processed_edits.bond_distance_model_max)
    if (h_bond_table is None) :
      hydrogen_bonds = None
    else :
      hydrogen_bonds = self.process_hydrogen_bonds(h_bond_table,
        log=log)
    asu_mappings = self.special_position_settings.asu_mappings(
      buffer_thickness=max_bond_distance*3)
        # factor 3 is to reach 1-4 interactions
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
      bond_params_table.update(
        i_seq=i_seq,
        j_seq=j_seq,
        params=geometry_restraints.bond_params(
          distance_ideal=disulfide_bond.value_dist,
          weight=1/disulfide_bond.value_dist_esd**2))
      bond_asu_table.add_pair(
        i_seq=i_seq,
        j_seq=j_seq,
        rt_mx_ji=sym_pair.rt_mx_ji)
    #
    if (processed_edits is not None):
      for proxy in processed_edits.bond_sym_proxies:
        if (proxy.weight <= 0): continue
        i_seq, j_seq = proxy.i_seqs
        bond_params_table.update(i_seq=i_seq, j_seq=j_seq, params=proxy)
        bond_asu_table.add_pair(
          i_seq=i_seq,
          j_seq=j_seq,
          rt_mx_ji=proxy.rt_mx_ji)
      for proxy in processed_edits.angle_proxies:
        self.geometry_proxy_registries.angle.append_custom_proxy(proxy=proxy)
      for proxy in processed_edits.planarity_proxies:
        self.geometry_proxy_registries.planarity.append_custom_proxy(proxy=proxy)
    #
    if (hydrogen_bonds is not None) :
      for proxy in hydrogen_bonds :
        if (proxy.weight <= 0): continue
        i_seq, j_seq = proxy.i_seqs
        bond_params_table.update(i_seq=i_seq, j_seq=j_seq, params=proxy)
        bond_asu_table.add_pair(
          i_seq=i_seq,
          j_seq=j_seq,
          rt_mx_ji=proxy.rt_mx_ji)
    #
    shell_asu_tables = crystal.coordination_sequences.shell_asu_tables(
      pair_asu_table=bond_asu_table,
      max_shell=3)
    shell_sym_tables = [shell_asu_table.extract_pair_sym_table()
      for shell_asu_table in shell_asu_tables]
    #
    nonbonded_params = ener_lib_as_nonbonded_params(
      ener_lib=ener_lib,
      assume_hydrogens_all_missing=assume_hydrogens_all_missing,
      factor_1_4_interactions=self.params.vdw_1_4_factor,
      default_distance=self.params.default_vdw_distance,
      minimum_distance=self.params.min_vdw_distance)
    result = geometry_restraints.manager.manager(
      crystal_symmetry=self.special_position_settings,
      model_indices=self.model_indices,
      conformer_indices=self.conformer_indices,
      sym_excl_indices=self.sym_excl_indices,
      site_symmetry_table=self.site_symmetry_table(),
      bond_params_table=bond_params_table,
      shell_sym_tables=shell_sym_tables,
      nonbonded_params=nonbonded_params,
      nonbonded_types=self.nonbonded_energy_type_registry.symbols,
      nonbonded_function=geometry_restraints.prolsq_repulsion_function(),
      nonbonded_distance_cutoff=self.params.nonbonded_distance_cutoff,
      nonbonded_buffer=self.params.nonbonded_buffer,
      angle_proxies=self.geometry_proxy_registries.angle.proxies,
      dihedral_proxies=self.geometry_proxy_registries.dihedral.proxies,
      chirality_proxies=self.geometry_proxy_registries.chirality.proxies,
      planarity_proxies=self.geometry_proxy_registries.planarity.proxies,
      max_reasonable_bond_distance=self.params.max_reasonable_bond_distance,
      plain_pairs_radius=plain_pairs_radius)
    if (params_remove is not None):
      self.process_geometry_restraints_remove(
        params=params_remove, geometry_restraints_manager=result)
    self.time_building_geometry_restraints_manager = timer.elapsed()
    return result

  def extract_xray_structure(self, unknown_scattering_type_substitute = "?"):
    from cctbx import xray
    from cctbx import adptbx
    from cctbx import eltbx
    import cctbx.eltbx.xray_scattering
    from itertools import count
    assert self.special_position_settings is not None
    result = xray.structure(
      special_position_settings=self.special_position_settings)
    sites_frac = result.unit_cell().fractionalize(
      sites_cart=self.sites_cart_exact())
    site_symmetry_table = self.site_symmetry_table()
    if (site_symmetry_table is not None):
      assert site_symmetry_table.indices().size() == sites_frac.size()
    scattering_types = self.scattering_type_registry.symbols
    if (scattering_types is None):
      scattering_types = self.get_element_symbols(strip_symbols=True)
    site_symmetry_ops = None
    for i_seq,atom,site_frac,scattering_type in zip(
          count(),
          self.pdb_atoms,
          sites_frac,
          scattering_types):
      assert atom.i_seq == i_seq
      from cctbx.eltbx.xray_scattering import get_standard_label
      scattering_type = get_standard_label(
        label=scattering_type, exact=True, optional=True)
      if (scattering_type is not None):
        charge = atom.charge_tidy(strip=True)
        if (charge is not None and len(charge) != 0):
          scattering_type_with_charge = get_standard_label(
            label=scattering_type+charge, exact=False, optional=True)
          if (scattering_type_with_charge is not None):
            scattering_type = scattering_type_with_charge
      if (scattering_type is None):
        if (unknown_scattering_type_substitute is None):
          raise RuntimeError("Unknown scattering type: %s" %
            atom.format_atom_record(cut_after_label_columns=True))
        scattering_type = unknown_scattering_type_substitute
      if (not atom.uij_is_defined()):
        u = adptbx.b_as_u(atom.b)
      else:
        u = adptbx.u_cart_as_u_star(result.unit_cell(), atom.uij)
      if (site_symmetry_table is not None):
        site_symmetry_ops = site_symmetry_table.get(i_seq=i_seq)
      result.add_scatterer(
        scatterer=xray.scatterer(
          label=atom.id_str(),
          site=site_frac,
          u=u,
          occupancy=atom.occ,
          scattering_type=scattering_type),
        site_symmetry_ops=site_symmetry_ops)
    return result

def show_residue_groups(residue_groups, log, prefix, max_items):
  if (len(residue_groups) == max_items + 1):
    max_items += 1
  for rg in residue_groups[:max_items]:
    if (rg.unique_resnames().size() == 1):
      print >> log, prefix+"residue:"
    else:
      print >> log, prefix+"residue group:"
    rg_atoms = rg.atoms()
    def show_atom(i):
      a = rg_atoms[i]
      print >> log, prefix+"  %s occ=%.2f" % (a.id_str(), a.occ)
    show_atom(0)
    n = rg_atoms.size()
    if (n > 3):
      print >> log, prefix+"  ... (%d atoms not shown)" % (n-2)
    elif (n == 3):
      show_atom(1)
    if (n > 1):
      show_atom(-1)
  n = len(residue_groups) - max_items
  if (n > 0):
    print >> log, prefix+"... (remaining %d not shown)" % n

class process(object):

  def __init__(self,
        mon_lib_srv,
        ener_lib,
        params=None,
        file_name=None,
        raw_records=None,
        pdb_inp=None,
        strict_conflict_handling=False,
        special_position_settings=None,
        crystal_symmetry=None,
        force_symmetry=False,
        substitute_non_crystallographic_unit_cell_if_necessary=False,
        keep_monomer_mappings=False,
        max_atoms=None,
        log=None,
        for_dihedral_reference=False):
    self.mon_lib_srv = mon_lib_srv
    self.ener_lib = ener_lib
    self.log = log
    self.for_dihedral_reference=for_dihedral_reference
    self.all_chain_proxies = build_all_chain_proxies(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      params=params,
      file_name=file_name,
      raw_records=raw_records,
      pdb_inp=pdb_inp,
      special_position_settings=special_position_settings,
      crystal_symmetry=crystal_symmetry,
      force_symmetry=force_symmetry,
      substitute_non_crystallographic_unit_cell_if_necessary
        =substitute_non_crystallographic_unit_cell_if_necessary,
      strict_conflict_handling=strict_conflict_handling,
      keep_monomer_mappings=keep_monomer_mappings,
      max_atoms=max_atoms,
      log=log,
      for_dihedral_reference=for_dihedral_reference)
    if (log is not None
        and self.all_chain_proxies.time_building_chain_proxies is not None):
      print >> log, \
        "  Time building chain proxies: %.2f, per 1000 atoms: %.2f" % (
          self.all_chain_proxies.time_building_chain_proxies,
          self.all_chain_proxies.time_building_chain_proxies * 1000
            / max(1,self.all_chain_proxies.pdb_atoms.size()))
    self._geometry_restraints_manager = None
    self._xray_structure = None

  def geometry_restraints_manager(self,
        plain_pairs_radius=None,
        params_edits=None,
        params_remove=None,
        h_bond_table=None,
        assume_hydrogens_all_missing=True,
        show_energies=True,
        hard_minimum_bond_distance_model=0.001):
    if (    self.all_chain_proxies.sites_cart is not None
        and self.all_chain_proxies.special_position_settings is not None
        and self._geometry_restraints_manager is None):
      self._geometry_restraints_manager \
        = self.all_chain_proxies.construct_geometry_restraints_manager(
            ener_lib=self.ener_lib,
            disulfide_link=self.mon_lib_srv.link_link_id_dict["SS"],
            plain_pairs_radius=plain_pairs_radius,
            params_edits=params_edits,
            params_remove=params_remove,
            h_bond_table=h_bond_table,
            assume_hydrogens_all_missing=assume_hydrogens_all_missing,
            log=self.log)
      if (self.log is not None):
        print >> self.log, \
          "  Time building geometry restraints manager: %.2f seconds" % (
            self.all_chain_proxies.time_building_geometry_restraints_manager)
        print >> self.log
        def note_geo():
          print >> self.log, """\
  NOTE: a complete listing of the restaints can be found in the
        .geo file."""
        note_geo()
        print >> self.log
        flush_log(self.log)
        site_labels = [atom.id_str()
          for atom in self.all_chain_proxies.pdb_atoms]
        pair_proxies = self._geometry_restraints_manager.pair_proxies(
          sites_cart=self.all_chain_proxies.sites_cart_exact())
        params = self.all_chain_proxies.params
        pair_proxies.bond_proxies.show_histogram_of_model_distances(
          sites_cart=self.all_chain_proxies.sites_cart_exact(),
          n_slots=params.show_histogram_slots.bond_lengths,
          f=self.log,
          prefix="  ")
        smallest_distance_model = \
          pair_proxies.bond_proxies.show_sorted(
            by_value="residual",
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            site_labels=site_labels,
            f=self.log,
            prefix="  ",
            max_items=params.show_max_items.bond_restraints_sorted_by_residual)
        if (    smallest_distance_model is not None
            and hard_minimum_bond_distance_model is not None
            and smallest_distance_model < hard_minimum_bond_distance_model):
          raise Sorry("""Bond restraint model distance < %.6g:
  Please inspect the output above and correct the input PDB file.""" % (
            hard_minimum_bond_distance_model))
        print >> self.log
        self._geometry_restraints_manager.angle_proxies \
          .show_histogram_of_deltas(
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            n_slots=params.show_histogram_slots
              .bond_angle_deviations_from_ideal,
            f=self.log,
            prefix="  ")
        self._geometry_restraints_manager.angle_proxies \
          .show_sorted(
            by_value="residual",
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            site_labels=site_labels,
            f=self.log,
            prefix="  ",
            max_items=params.show_max_items
              .bond_angle_restraints_sorted_by_residual)
        print >> self.log
        self._geometry_restraints_manager.dihedral_proxies \
          .show_histogram_of_deltas(
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            n_slots=params.show_histogram_slots
              .dihedral_angle_deviations_from_ideal,
            f=self.log,
            prefix="  ")
        self._geometry_restraints_manager.dihedral_proxies \
          .show_sorted(
            by_value="residual",
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            site_labels=site_labels,
            f=self.log,
            prefix="  ",
            max_items=params.show_max_items
              .dihedral_angle_restraints_sorted_by_residual)
        print >> self.log
        self._geometry_restraints_manager.chirality_proxies \
          .show_histogram_of_deltas(
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            n_slots=params.show_histogram_slots
              .chiral_volume_deviations_from_ideal,
            f=self.log,
            prefix="  ")
        self._geometry_restraints_manager.chirality_proxies \
          .show_sorted(
            by_value="residual",
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            site_labels=site_labels,
            f=self.log,
            prefix="  ",
            max_items=params.show_max_items
              .chirality_restraints_sorted_by_residual)
        print >> self.log
        self._geometry_restraints_manager.planarity_proxies \
          .show_sorted(
            by_value="residual",
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            site_labels=site_labels,
            f=self.log,
            prefix="  ",
            max_items=params.show_max_items
              .planarity_restraints_sorted_by_residual)
        print >> self.log
        pair_proxies.nonbonded_proxies.show_histogram_of_model_distances(
          sites_cart=self.all_chain_proxies.sites_cart_exact(),
          n_slots=params.show_histogram_slots.nonbonded_interaction_distances,
          f=self.log,
          prefix="  ")
        pair_proxies.nonbonded_proxies.show_sorted(
          by_value="delta",
          sites_cart=self.all_chain_proxies.sites_cart_exact(),
          site_labels=site_labels,
          f=self.log,
          prefix="  ",
          max_items=params.show_max_items
            .nonbonded_interactions_sorted_by_model_distance)
        print >> self.log
        note_geo()
        flush_log(self.log)
        if (show_energies):
          print >> self.log
          timer = user_plus_sys_time()
          energies = self._geometry_restraints_manager.energies_sites(
            sites_cart=self.all_chain_proxies.sites_cart_exact())
          energies.show(f=self.log, prefix="  ")
          print >> self.log, "  Time first energy calculation" \
                             " (mainly nonbonded setup): %.2f" % (
            timer.elapsed())
          flush_log(self.log)
        if not self.for_dihedral_reference:
          self.clash_guard()
    return self._geometry_restraints_manager

  def clash_guard(self, hard_minimum_nonbonded_distance=0.001):
    params = self.all_chain_proxies.params.clash_guard
    if (params.nonbonded_distance_threshold is None): return
    geo = self._geometry_restraints_manager
    n_below_threshold = (
      geo.nonbonded_model_distances() < params.nonbonded_distance_threshold) \
        .count(True)
    if ((params.max_number_of_distances_below_threshold is not None
         and n_below_threshold
               > params.max_number_of_distances_below_threshold)
        or
        (params.max_fraction_of_distances_below_threshold is not None
         and n_below_threshold
               / max(1,geo.sites_cart_used_for_pair_proxies().size())
                 > params.max_fraction_of_distances_below_threshold)):
      phil_path = params.__phil_path__()
      raise Sorry("""%s failure:
  Number of nonbonded interaction distances < %.6g: %d
    Please inspect the histogram of nonbonded interaction distances above.
    To disable this error, run the same command again with the following
    additional argument:
      %s.nonbonded_distance_threshold=None""" %
        (phil_path, params.nonbonded_distance_threshold,
         n_below_threshold, phil_path))
    #
    n_below_hard_minimum_nonbonded_distance = (
      geo.nonbonded_model_distances() < hard_minimum_nonbonded_distance) \
        .count(True)
    if (n_below_hard_minimum_nonbonded_distance != 0 and
      params.nonbonded_distance_threshold >=0):
      raise Sorry("""Number of nonbonded interaction distances < %.6g: %d
  Please inspect the output above and correct the input PDB file.""" % (
        hard_minimum_nonbonded_distance,
        n_below_hard_minimum_nonbonded_distance))

  def xray_structure(self, show_summary = True):
    log = self.log
    if (    self.all_chain_proxies.sites_cart is not None
        and self.all_chain_proxies.special_position_settings is not None
        and self._xray_structure is None):
      self._xray_structure = self.all_chain_proxies.extract_xray_structure()
      self._xray_structure.scattering_type_registry(
        types_without_a_scattering_contribution=["?"])
      if (log is not None and show_summary):
        self._xray_structure.show_summary(f = log, prefix="  ")
        self._xray_structure.show_special_position_shifts(
          sites_cart_original=self.all_chain_proxies.sites_cart,
          out= log, prefix="  ")
        self._xray_structure.scattering_type_registry().show(
          show_gaussians=False, out = log, prefix="  ")
        flush_log(log)
    return self._xray_structure

  def show_selected_atoms(self,
        selection,
        header_lines=None,
        out=None,
        prefix=""):
    if (out is None): out = sys.stdout
    if (header_lines is not None):
      for line in header_lines:
        print >> out, prefix+line
    sub_hierarchy = self.all_chain_proxies.pdb_hierarchy.select(
      atom_selection=selection)
    s = sub_hierarchy.as_pdb_string()
    if (len(s) == 0 and header_lines is not None):
      s = "  None\n"
    if (prefix == ""):
      out.write(s)
    else:
      for line in s.splitlines():
        print >> out, prefix+line

  def show_atoms_without_ncs_restraints(self,
        ncs_restraints_groups,
        out=None,
        prefix=""):
    self.show_selected_atoms(
      selection=~ncs_restraints_groups.selection_restrained(
        n_seq=self.all_chain_proxies.pdb_atoms.size()),
      header_lines=["Atoms without NCS restraints:"],
      out=out,
      prefix=prefix)

def run(
      args,
      params=None,
      strict_conflict_handling=True,
      substitute_non_crystallographic_unit_cell_if_necessary=False,
      return_all_processed_pdb_files=False,
      max_atoms=None,
      log=None):
  if (log is None): log = sys.stdout
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  pdb_file_names = []
  cif_file_names = []
  for arg in args:
    if (not os.path.isfile(arg)):
      raise Sorry("No such file: %s" % show_string(arg))
    if (iotbx.pdb.is_pdb_file(file_name=arg)):
      pdb_file_names.append(arg)
    else:
      try:
        cif_object = server.read_cif(file_name=arg)
      except KeyboardInterrupt: raise
      except:
        raise Sorry("Unknown file format: %s" % show_string(arg))
      else:
        print >> log, "Processing CIF file: %s" % show_string(arg)
        for srv in [mon_lib_srv, ener_lib]:
          srv.process_cif_object(cif_object=cif_object, file_name=arg)
  all_processed_pdb_files = []
  for file_name in pdb_file_names:
    processed_pdb_file = process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      params=params,
      file_name=file_name,
      strict_conflict_handling=strict_conflict_handling,
      substitute_non_crystallographic_unit_cell_if_necessary
        =substitute_non_crystallographic_unit_cell_if_necessary,
      max_atoms=max_atoms,
      log=log)
    processed_pdb_file.geometry_restraints_manager()
    processed_pdb_file.xray_structure()
    if (return_all_processed_pdb_files):
      all_processed_pdb_files.append(processed_pdb_file)
  if (return_all_processed_pdb_files):
    return all_processed_pdb_files
  if (len(pdb_file_names) > 0):
    return processed_pdb_file
  return None

if (__name__ == "__main__"):
  run(sys.argv[1:])
