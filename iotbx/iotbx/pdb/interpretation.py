from iotbx import pdb
import iotbx.pdb.parser
import iotbx.pdb.cryst1_interpretation
import iotbx.pdb.atom
from cctbx.array_family import flex
from libtbx.itertbx import count
import string

class empty: pass

class model:

  def __init__(self, stage_1, serial):
    self.stage_1 = stage_1
    self.serial = serial
    self.altLocs = {}
    self.conformers = []

  def add_conformer(self, conformer):
    assert not self.altLocs.has_key(conformer.altLoc)
    self.altLocs[conformer.altLoc] = len(self.conformers)
    self.conformers.append(conformer)

class conformer:

  def __init__(self, model, altLoc, iselection):
    self.model = model
    self.altLoc = altLoc
    self.iselection = iselection

  def iselection_common_atoms(self, other):
    assert other.model.stage_1 is self.model.stage_1
    return self.model.stage_1.selection_cache().intersection(
      iselections=[self.iselection, other.iselection]).iselection()

  def get_chains(self):
    stage_1 = self.model.stage_1
    ter_block_identifiers = stage_1.get_ter_block_identifiers()
    chains = []
    residue_iselections = []
    isel_residue = flex.size_t()
    prev_atom = None
    prev_residue_labels = ""
    prev_ter_block_identifier = -1
    for i_seq in self.iselection:
      atom = stage_1.atom_attributes_list[i_seq]
      residue_labels = atom.residue_labels()
      ter_block_identifier = ter_block_identifiers[i_seq]
      if (prev_atom is not None):
        if (   ter_block_identifier != prev_ter_block_identifier
            or not atom.is_in_same_chain(prev_atom)):
          if (isel_residue.size() > 0):
            residue_iselections.append(isel_residue)
            isel_residue = flex.size_t()
          if (len(residue_iselections) > 0):
            chains.append(pdb.interpretation.chain(
              conformer=self, residue_iselections=residue_iselections))
            residue_iselections = []
        elif (residue_labels != prev_residue_labels):
          if (isel_residue.size() > 0):
            residue_iselections.append(isel_residue)
            isel_residue = flex.size_t()
      isel_residue.append(i_seq)
      prev_atom = atom
      prev_residue_labels = residue_labels
      prev_ter_block_identifier = ter_block_identifier
    if (isel_residue.size() > 0):
      residue_iselections.append(isel_residue)
    if (len(residue_iselections) > 0):
      chains.append(pdb.interpretation.chain(
        conformer=self, residue_iselections=residue_iselections))
    return chains

class chain:

  def __init__(self, conformer, residue_iselections):
    self.conformer = conformer
    self.residues = []
    for iselection in residue_iselections:
      self.residues.append(
        residue(chain=self, iselection=iselection))

class residue:

  def __init__(self, chain, iselection):
    self.chain = chain
    self.iselection = iselection

class stage_1:

  def __init__(self, file_name=None, raw_records=None):
    assert [file_name, raw_records].count(None) == 1
    if (raw_records is None):
      raw_records = open(file_name)
    columns_73_76_eval = pdb.parser.columns_73_76_evaluator(
      raw_records=raw_records)
    raw_records = columns_73_76_eval.raw_records
    self.ignore_columns_73_and_following = columns_73_76_eval.is_old_style
    self.crystal_symmetry = None
    self.remark_290_records = []
    self.ter_indices = []
    self.break_indices = []
    self.atom_attributes_list = []
    self.conect_records = []
    self.link_records = []
    self.ssbond_records = []
    self.sltbrg_records = []
    self.model_serial_list = []
    model_serial = None
    altLoc_dict = {}
    state = empty()
    self.state = state
    for state.line_number,state.raw_record in zip(count(1), raw_records):
      record_name = state.raw_record[:6]
      if (record_name == "CRYST1"):
        self.crystal_symmetry = pdb.cryst1_interpretation.crystal_symmetry(
          cryst1_record=state.raw_record,
          line_number=state.line_number)
      elif (state.raw_record.startswith("REMARK 290 ")):
        self.remark_290_records.append(state.raw_record.rstrip())
      elif (record_name == "MODEL "):
        model_serial = self.parse_record().serial
        self.model_serial_list.append(model_serial)
      elif (record_name == "ENDMDL"):
        model_serial = None
      elif (record_name.rstrip() == "TER"):
        self.ter_indices.append(len(self.atom_attributes_list))
      elif (record_name.rstrip() == "BREAK"):
        self.break_indices.append(len(self.atom_attributes_list))
      elif (record_name in ("ATOM  ", "HETATM")):
        atom_attributes = pdb.atom.attributes(line_number=state.line_number)
        atom_attributes.set_from_ATOM_record(self.parse_record())
        if (model_serial is None):
          atom_attributes.MODELserial = -1
        else:
          atom_attributes.MODELserial = len(self.model_serial_list)-1
        self.atom_attributes_list.append(atom_attributes)
        altLoc_dict[atom_attributes.altLoc] = 0
      elif (record_name == "SIGATM"):
        if (len(self.atom_attributes_list) > 0):
          self.atom_attributes_list[-1].set_from_SIGATM_record(
            self.parse_record())
      elif (record_name == "ANISOU"):
        if (len(self.atom_attributes_list) > 0):
          self.atom_attributes_list[-1].set_from_ANISOU_record(
            self.parse_record())
      elif (record_name == "SIGUIJ"):
        if (len(self.atom_attributes_list) > 0):
          self.atom_attributes_list[-1].set_from_SIGUIJ_record(
            self.parse_record())
      elif (record_name == "CONECT"):
        self.conect_records.append(state.raw_record)
      elif (record_name == "LINK  "):
        self.link_records.append(state.raw_record)
      elif (record_name == "SSBOND"):
        self.ssbond_records.append(state.raw_record)
      elif (record_name == "SLTBRG"):
        self.sltbrg_records.append(state.raw_record)
    del self.state
    self._sites_cart = None
    self._ter_block_identifiers = None
    self._break_block_identifiers = None
    self._selection_cache = None
    self._clean_model_serial_list()
    self._fix_false_blank_altLoc_identifiers(altLoc_dict)

  def parse_record(self):
    return pdb.parser.pdb_record(
      raw_record=self.state.raw_record,
      line_number=self.state.line_number,
      ignore_columns_73_and_following=self.ignore_columns_73_and_following)

  def get_sites_cart(self):
    if (self._sites_cart is None):
      self._sites_cart = flex.vec3_double()
      for atom in self.atom_attributes_list:
        self._sites_cart.append(atom.coordinates)
    return self._sites_cart

  def get_block_identifiers(self, block_indices):
    if (len(block_indices) == 0):
      return flex.size_t(len(self.atom_attributes_list), 0)
    result = flex.size_t()
    for i,j in enumerate(block_indices):
      result.resize(j, i)
    result.resize(len(self.atom_attributes_list), result.back()+1)
    return result

  def get_ter_block_identifiers(self):
    if (self._ter_block_identifiers is None):
      self._ter_block_identifiers = self.get_block_identifiers(
        self.ter_indices)
    return self._ter_block_identifiers

  def get_break_block_identifiers(self):
    if (self._break_block_identifiers is None):
      self._break_block_identifiers = self.get_block_identifiers(
        self.break_indices)
    return self._break_block_identifiers

  def selection_cache(self):
    if (self._selection_cache is None):
      self._selection_cache = pdb.atom.selection_cache(
        atom_attributes_list=self.atom_attributes_list)
    return self._selection_cache

  def _clean_model_serial_list(self):
    self.n_model_serial_numbers_changed = 0
    clean_model_serial_list = []
    for model_serial in self.model_serial_list:
      if (model_serial in clean_model_serial_list):
        model_serial = max(clean_model_serial_list) + 1
        if (model_serial in self.model_serial_list):
          model_serial = max(self.model_serial_list) + 1
        self.n_model_serial_numbers_changed += 1
      clean_model_serial_list.append(model_serial)
    self.model_serial_list = clean_model_serial_list
    if (len(self.model_serial_list) > 0):
      for atom_attributes in self.atom_attributes_list:
        if (atom_attributes.MODELserial < 0): continue
        atom_attributes.MODELserial \
          = self.model_serial_list[atom_attributes.MODELserial]
    if (self.n_model_serial_numbers_changed > 0):
      self._selection_cache = None

  def _fix_false_blank_altLoc_identifiers(self, altLoc_dict):
    self.n_patched_altLocs = 0
    if (len(altLoc_dict) == 1 or " " not in altLoc_dict): return
    sel_cache = self.selection_cache()
    is_processed = flex.bool(sel_cache.n_seq, 00000)
    for MODELserial,isel_model in sel_cache.MODELserial.items():
      false_blank_related = {}
      altLoc_groups = altLoc_grouping()
      for altLoc_pivot,isel_altLoc_pivot in sel_cache.altLoc.items():
        if (altLoc_pivot == " "): continue
        isel_pivot = sel_cache.intersection(
          iselections=[isel_model, isel_altLoc_pivot]).iselection()
        altLoc_group = {}
        for i_seq_pivot in isel_pivot:
          if (is_processed[i_seq_pivot]): continue
          atom = self.atom_attributes_list[i_seq_pivot]
          isel_related_atoms = sel_cache.intersection(sel_cache.get_labels(
            name=atom.name,
            altLoc=None,
            resName=atom.resName,
            chainID=atom.chainID,
            resSeq=atom.resSeq,
            iCode=atom.iCode,
            segID=atom.segID,
            MODELserial=atom.MODELserial)).iselection()
          for i_seq in isel_related_atoms:
            atom = self.atom_attributes_list[i_seq]
            is_processed[i_seq] = 0001
            altLoc_group[atom.altLoc] = 0
            if (atom.altLoc == " "):
              false_blank_related[i_seq] = altLoc_pivot
        altLoc_groups.add_group(altLoc_group.keys())
      membership_indices = altLoc_groups.get_membership_indices()
      altLoc_replacements = altLoc_groups.get_false_blank_altLoc_replacements()
      for i_seq,altLoc in false_blank_related.items():
        atom = self.atom_attributes_list[i_seq]
        assert atom.altLoc == " "
        atom.altLoc = altLoc_replacements[membership_indices[altLoc]]
        self.n_patched_altLocs += 1
    if (self.n_patched_altLocs > 0):
      self._selection_cache = None

  def get_models_and_conformers(self):
    sel_cache = self.selection_cache()
    altLoc_unions = {}
    isel_blank_altLoc = sel_cache.altLoc.get(" ", None)
    if (isel_blank_altLoc is None):
      altLoc_unions = sel_cache.altLoc
    else:
      for altLoc,isel_altLoc in sel_cache.altLoc.items():
        altLoc_unions[altLoc] = sel_cache.union(
          iselections=[isel_blank_altLoc, isel_altLoc]).iselection()
    models = []
    for MODELserial,isel_model in sel_cache.MODELserial.items():
      conformer_dict = {}
      for altLoc,isel_altLoc_union in altLoc_unions.items():
        isel_conformer = sel_cache.intersection(
          iselections=[isel_model, isel_altLoc_union]).iselection()
        if (isel_conformer.size() > 0):
          conformer_dict[altLoc] = isel_conformer
      if (len(conformer_dict) > 1 and " " in conformer_dict):
        del conformer_dict[" "]
      altLocs = conformer_dict.keys()
      altLocs.sort()
      model_ = model(stage_1=self, serial=MODELserial)
      for altLoc in altLocs:
        conformer_ = conformer(
          model=model_,
          altLoc=altLoc,
          iselection=conformer_dict[altLoc])
        model_.add_conformer(conformer_)
      models.append(model_)
    return models

class altLoc_grouping:

  def __init__(self):
    self.group_list = []

  def get_membership_indices(self):
    result = {}
    for i,group in enumerate(self.group_list):
      for altLoc in group:
        if (altLoc == " "): continue
        assert not result.has_key(altLoc)
        result[altLoc] = i
    return result

  def add_group(self, new_group):
    if (len(new_group) == 0): return
    membership_indices = self.get_membership_indices()
    indices_existing = {}
    for altLoc in new_group:
      i = membership_indices.get(altLoc, None)
      if (i is not None):
        indices_existing[i] = 0
    if (len(indices_existing) == 0):
      self.group_list.append(new_group[:])
      self.group_list[-1].sort()
    else:
      indices_existing = indices_existing.keys()
      indices_existing.sort()
      if (len(indices_existing) > 1):
        sel = flex.bool(len(self.group_list), 0001)
        merged_group = self.group_list[indices_existing[0]]
        for i in indices_existing[1:]:
          merged_group.extend(self.group_list[i])
          sel[i] = 00000
        self.group_list = flex.select(sequence=self.group_list, flags=sel)
      group = self.group_list[indices_existing[0]]
      group.extend(new_group)
      group = dict([(altLoc,0) for altLoc in group]).keys()
      group.sort()
      self.group_list[indices_existing[0]] = group

  def get_false_blank_altLoc_replacements(self):
    result = {}
    all_altLocs = []
    for group in self.group_list:
      all_altLocs.extend(group)
    for i,group in enumerate(self.group_list):
      if (" " in group):
        replacement = get_false_blank_altLoc_replacement(
          altLoc_group=group,
          all_altLocs=all_altLocs)
        result[i] = replacement
        all_altLocs.append(replacement)
    return result

def get_false_blank_altLoc_replacement(
      altLoc_group,
      all_altLocs,
      letter_priorities="XYZUVWQRSTIJKLMNOPHGFEDCBA",
      digit_priorities="0987654321",
      other_priorities="#.:=_^~@$%&*+-!|/\'\"`?<>[]{}(),;"):
  altLoc_group = altLoc_group[:]
  altLoc_group.sort()
  n_letters = 0
  n_digits = 0
  for altLoc in altLoc_group:
    if (altLoc in string.letters):
      n_letters += 1
    elif (altLoc in string.digits):
      n_digits += 1
  if (n_letters >= n_digits):
    priorities = letter_priorities + digit_priorities
  else:
    priorities = digit_priorities + letter_priorities
  priorities += other_priorities + letter_priorities.lower()
  for result in priorities:
    if (not result in all_altLocs):
      return result
  raise RuntimeError(
    "Cannot find a replacement for false blank altLoc identifiers.")
