
# Implemented based on PDB v3.2 specification at:
#   http://www.wwpdb.org/documentation/format32/sect5.html

# NOTE: all hydrogen bond information is returned as atom pairs, donor first.

from libtbx.utils import Sorry, Usage
from libtbx import smart_open
import libtbx.phil
from libtbx import group_args
import string, sys, os

ss_input_params_str = """
  file_name = None
    .type = path
    .multiple = True
    .optional = True
  use_hydrogens = True
    .type = bool
  include_helices = True
    .type = bool
  include_sheets = True
    .type = bool
"""
ss_input_params = libtbx.phil.parse(ss_input_params_str)

ss_group_params_str = """
helix
  .multiple = True
  .optional = True
{
  selection = None
    .type = str
    .style = bold selection
  helix_class = *1 2 3 4 5 6 7 8 9 10
    .type = choice
    .caption = alpha other pi other 3_10 other other other other other
    .help = Type of helix, defaults to alpha (1).  Only alpha, pi, and 3_10 \
      helices are used for hydrogen-bond restraints.
    .style = bold
  restraint_sigma = None
    .type = float
  backbone_only = False
    .type = bool
    .help = Only applies to rigid-body groupings, and not H-bond restraints \
      which are already backbone-only.
}
sheet
  .multiple = True
  .optional = True
{
  first_strand = None
    .type = str
    .style = bold selection
  strand
    .multiple = True
    .optional = False
  {
    selection = None
      .type = str
      .style = bold selection
    sense = parallel antiparallel *unknown
      .type = choice
      .style = bold
    bond_start_current = None
      .type = str
      .style = bold selection
    bond_start_previous = None
      .type = str
      .style = bold selection
  }
  restraint_sigma = None
    .type = float
  backbone_only = False
    .type = bool
    .help = Only applies to rigid-body groupings, and not H-bond restraints \
      which are already backbone-only.
}
"""

class structure_base (object) :
  def extract_h_bonds (self, params) :
    return []

  def as_pdb_str (self) :
    return None

  def __str__ (self) :
    return self.as_pdb_str()

  def as_pymol_dashes (self, params, object_name=None) :
    cmds = []
    prefix = ""
    if object_name is not None :
      prefix = "%s and " % object_name
    bonded_atoms = self.extract_h_bonds(params)
    for (atom1, atom2) in bonded_atoms :
      sele1 = "(%schain '%s' and resi %d and name %s)" % (prefix,
        atom1.chain_id, atom1.resseq, atom1.name)
      sele2 = "(%schain '%s' and resi %d and name %s)" % (prefix,
        atom2.chain_id, atom2.resseq, atom2.name)
      cmd = "dist %s, %s" % (sele1, sele2)
      cmds.append(cmd)
    return "\n".join(cmds)

class annotation (structure_base) :
  def __init__ (self, records) :
    self.helices = parse_helix_records(records)
    self.sheets = parse_sheet_records(records)

  def as_pdb_str (self) :
    records = []
    for helix in self.helices :
      records.append(helix.as_pdb_str())
    for sheet in self.sheets :
      records.append(sheet.as_pdb_str())
    return "\n".join(records)

  def extract_h_bonds (self, params) :
    bonded_atoms = []
    if params.include_helices :
      for helix in self.helices :
        helix_bonds = helix.extract_h_bonds(params)
        bonded_atoms.extend(helix_bonds)
    if params.include_sheets :
      for sheet in self.sheets :
        sheet_bonds = sheet.extract_h_bonds(params)
        bonded_atoms.extend(sheet_bonds)
    return bonded_atoms

  def as_atom_selections (self, params) :
    selections = []
    if params.include_helices :
      for helix in self.helices :
        try :
          selections.extend(helix.as_atom_selections(params))
        except RuntimeError, e :
          pass
    if params.include_sheets :
      for sheet in self.sheets :
        selections.extend(sheet.as_atom_selections(params))
    return selections

  def overall_helix_selection (self, params=ss_input_params) :
    selections = []
    for helix in self.helices :
      try :
        selections.extend(helix.as_atom_selections(params))
      except RuntimeError, e :
        pass
    return "(" + ") or (".join(selections) + ")"

  def overall_sheet_selection (self, params=ss_input_params) :
    for sheet in self.sheets :
      selections.extend(sheet.as_atom_selections(params))
    return "(" + ") or (".join(selections) + ")"

  def as_bond_selections (self, params) :
    bonded_atoms = self.extract_h_bonds(params)
    selections = []
    for (atom1, atom2) in bonded_atoms :
      selection_1 = "name %s and chain '%s' and resseq %d and icode '%s'" % (
        atom1.name, atom1.chain_id, atom1.resseq, atom1.icode)
      selection_2 = "name %s and chain '%s' and resseq %d and icode '%s'" % (
        atom2.name, atom2.chain_id, atom2.resseq, atom2.icode)
      selections.append((selection_1, selection_2))
    return selections

  def as_restraint_groups (self, log=sys.stderr, prefix_scope="") :
    phil_strs = []
    for helix in self.helices :
      helix_phil = helix.as_restraint_group(log, prefix_scope)
      if helix_phil is not None :
        phil_strs.append(helix_phil)
    for sheet in self.sheets :
      sheet_phil = sheet.as_restraint_group(log, prefix_scope)
      if sheet_phil is not None :
        phil_strs.append(sheet_phil)
    return "\n".join(phil_strs)

#-----------------------------------------------------------------------
class pdb_helix (structure_base, group_args) :
  def as_pdb_str (self) :
    format = "HELIX  %3d %3s %3s %1s %4d%1s %3s %1s %4d%1s%2d%30s %5d"
    out = format % (self.serial, self.helix_id, self.start_resname,
      self.start_chain_id, self.start_resseq, self.start_icode,
      self.end_resname, self.end_chain_id, self.end_resseq, self.end_icode,
      self.helix_class, self.comment, self.length)
    return out.strip()

  def continuity_check (self) :
    if self.start_icode != self.end_icode :
      raise RuntimeError("Don't know how to deal with helices with multiple "+
        "insertion codes ('%s' vs. '%s')." % (self.start_icode, self.end_icode))
    if self.start_chain_id != self.end_chain_id :
      raise RuntimeError("Don't know how to deal with helices with multiple "+
        "chain IDs ('%s' vs. '%s')." % (self.start_chain_id, self.end_chain_id))

  def as_atom_selections (self, params) :
    self.continuity_check()
    sele = "chain '%s' and resseq %d:%d and icode '%s'" % (self.start_chain_id,
      self.start_resseq, self.end_resseq, self.end_icode)
    return [sele]

  def as_restraint_group (self, log=sys.stderr, prefix_scope="") :
    if self.start_chain_id != self.end_chain_id :
      print >> log, "Helix chain ID mismatch: starts in %s, ends in %s" % (
        self.start_chain_id, self.end_chain_id)
      return None
    resid_start = "%d%s" % (self.start_resseq, self.start_icode)
    resid_end = "%d%s" % (self.end_resseq, self.end_icode)
    sele = "chain '%s' and resid %s through %s" % (self.start_chain_id,
      resid_start, resid_end)
    if prefix_scope != "" and not prefix_scope.endswith(".") :
      prefix_scope += "."
    rg = """\
%shelix {
  selection = %s
  helix_class = %d
}""" % (prefix_scope, sele, self.helix_class)
    return rg

  def extract_h_bonds (self, params) :
    self.continuity_check()
    bonded_atoms = []
    i = 0
    if self.helix_class == 1 : # alpha
      j = 4
    else :
      if self.helix_class == 5 : # 3_10
        j = 3
      elif self.helix_class == 3 : # pi
        j = 5
      else :
        raise RuntimeError("Don't know how to deal with helix class %d." %
          self.helix_class)
    acceptor_name = "O"
    donor_name = "N"
    if params.use_hydrogens :
      donor_name = "H"
    while j <= self.length :
      resseq1 = self.start_resseq + i
      resseq2 = self.start_resseq + j
      i += 1
      j += 1
      #print resseq1, resseq2, self.end_resseq, self.helix_class
      if not resseq2 <= self.end_resseq :
        break
      acceptor = group_args(
        chain_id=self.start_chain_id,
        resseq=resseq1,
        name=acceptor_name,
        icode=self.start_icode)
      donor = group_args(
        chain_id=self.start_chain_id,
        resseq=resseq2,
        name=donor_name,
        icode=self.start_icode)
      bonded_atoms.append((donor, acceptor))
    return bonded_atoms

def parse_helix_records (records) :
  helices = []
  for line in records :
    if not line.startswith("HELIX") :
      continue
    current_helix = pdb_helix(
      serial=string.atoi(line[7:10]),
      helix_id=line[11:14].strip(),
      start_resname=line[15:18],
      start_chain_id=line[19],
      start_resseq=string.atoi(line[21:25]),
      start_icode=line[25],
      end_resname=line[27:30],
      end_chain_id=line[31],
      end_resseq=string.atoi(line[33:37]),
      end_icode=line[37],
      helix_class=string.atoi(line[38:40]),
      comment=line[40:70],
      length=string.atoi(line[71:76]))
    helices.append(current_helix)
  return helices

def restraint_groups_as_pdb_helices (pdb_hierarchy, helices, log=sys.stderr) :
  isel = pdb_hierarchy.atom_selection_cache().iselection
  atoms = [ a for a in pdb_hierarchy.atoms_with_labels() ]
  pdb_helices = []
  for i, helix_params in enumerate(helices) :
    if helix_params.selection is None :
      print >> log, "Empty helix at serial %d." % (i+1)
      continue
    sele_str = ("(%s) and name N and (altloc 'A' or altloc ' ')" %
                helix_params.selection)
    amide_isel = isel(sele_str)
    start_atom = atoms[amide_isel[0]]
    end_atom = atoms[amide_isel[-1]]
    current_helix = pdb_helix(
      serial=i+1,
      helix_id=i+1,
      start_resname=start_atom.resname,
      start_chain_id=start_atom.chain_id,
      start_resseq=start_atom.resseq,
      start_icode=start_atom.icode,
      end_resname=end_atom.resname,
      end_chain_id=end_atom.chain_id,
      end_resseq=end_atom.resseq,
      end_icode=end_atom.icode,
      helix_class=int(helix_params.helix_class),
      comment="",
      length=amide_isel.size())
    pdb_helices.append(current_helix)
  return pdb_helices

#-----------------------------------------------------------------------
class pdb_sheet (structure_base, group_args) :
  def add_strand (self, strand) :
    self.strands.append(strand)

  def add_registration (self, registration) :
    self.registrations.append(registration)

  def as_atom_selections (self, params) :
    strand_selections = []
    for strand in self.strands :
      if strand.start_icode != strand.end_icode :
        continue
      sele = "chain '%s' and resseq %d:%d and icode '%s'" % (
        strand.start_chain_id, strand.start_resseq, strand.end_resseq,
        strand.start_icode)
      strand_selections.append(sele)
    return strand_selections

  def as_restraint_group (self, log=sys.stderr, prefix_scope="") :
    if len(self.strands) == 0 :
      return None
    selections = []
    senses = []
    reg_curr = []
    reg_prev = []
    for (strand,registration) in zip(self.strands, self.registrations) :
      resid_start = "%d%s" % (strand.start_resseq, strand.start_icode)
      resid_end = "%d%s" % (strand.end_resseq, strand.end_icode)
      sele = "chain '%s' and resid %s through %s" % (strand.start_chain_id,
        resid_start, resid_end)
      selections.append(sele)
      if strand.sense == 0 :
        senses.append("unknown")
      elif strand.sense == -1 :
        senses.append("antiparallel")
      elif strand.sense == 1 :
        senses.append("parallel")
      else :
        raise Sorry("Sense must be 0, 1, or -1.")
      if registration is not None :
        sele_base = "chain '%s' and resid %s and name %s"
        resid_curr = "%d%s" % (registration.cur_resseq,registration.cur_icode)
        resid_prev = "%d%s" % (registration.prev_resseq,registration.prev_icode)
        reg_curr.append(sele_base % (registration.cur_chain_id,resid_curr,"N"))
        reg_prev.append(sele_base % (registration.prev_chain_id,resid_prev,"O"))
      else :
        reg_curr.append(None)
        reg_prev.append(None)
    n = 0
    first_strand = None
    strands = []
    for (sele, sense, curr, prev) in zip(selections,senses,reg_curr,reg_prev) :
      if n == 0 :
        first_strand = sele
      else :
        strands.append("""\
  strand {
    selection = %s
    sense = %s
    bond_start_current = %s
    bond_start_previous = %s
  }""" % (sele, sense, curr, prev))
      n += 1
    assert first_strand is not None
    if prefix_scope != "" and not prefix_scope.endswith(".") :
      prefix_scope += "."
    phil_str = """
%ssheet {
  first_strand = %s
%s
}""" % (prefix_scope, first_strand, "\n".join(strands))
    return phil_str

  def extract_h_bonds (self, params) :
    assert len(self.strands) == len(self.registrations)
    bonded_atoms = []
    errors = 0
    donor_name = "N"
    acceptor_name = "O"
    if params.use_hydrogens :
      donor_name = "H"
    for i, strand in enumerate(self.strands) :
      registration = self.registrations[i]
      prev_strand = self.strands[i - 1]
      if registration is None : # usually the first strand, but not always!
        continue
      if ((strand.start_icode != strand.end_icode) or
          (prev_strand.start_icode != prev_strand.end_icode)) :
        errors += 1
        continue # don't raise exception - other strands may be okay
      cur_resseq = registration.cur_resseq
      prev_resseq = registration.prev_resseq
      if strand.sense == -1 :
        while ((prev_resseq <= prev_strand.end_resseq) and
               (cur_resseq >= strand.start_resseq)) :
          # O (current) --> H/N (previous)
          acceptor1 = group_args(
            chain_id=strand.start_chain_id,
            resseq=cur_resseq,
            name=acceptor_name,
            icode=strand.start_icode)
          donor1 = group_args(
            chain_id=prev_strand.start_chain_id,
            resseq=prev_resseq,
            name=donor_name,
            icode=prev_strand.start_icode)
          bonded_atoms.append((donor1, acceptor1))
          # H/N (current) --> O (previous)
          donor2 = group_args(
            chain_id=strand.start_chain_id,
            resseq=cur_resseq,
            name=donor_name,
            icode=strand.start_icode)
          acceptor2 = group_args(
            chain_id=prev_strand.start_chain_id,
            resseq=prev_resseq,
            name=acceptor_name,
            icode=prev_strand.start_icode)
          bonded_atoms.append((donor2, acceptor2))
          prev_resseq += 2
          cur_resseq -= 2
      elif strand.sense == 1 :
        while ((prev_resseq <= prev_strand.end_resseq) and
               (cur_resseq <= strand.end_resseq)) :
          # O (current) --> H/N (previous)
          acceptor1 = group_args(
            chain_id=strand.start_chain_id,
            resseq=cur_resseq,
            name=acceptor_name,
            icode=strand.start_icode)
          donor1 = group_args(
            chain_id=prev_strand.start_chain_id,
            resseq=prev_resseq,
            name=donor_name,
            icode=prev_strand.start_icode)
          bonded_atoms.append((donor1, acceptor1))
          if (cur_resseq + 2) > strand.end_resseq :
            break
          # H/N (current + 2) --> O (previous)
          donor2 = group_args(
            chain_id=strand.start_chain_id,
            resseq=cur_resseq + 2,
            name=donor_name,
            icode=strand.start_icode)
          acceptor2 = group_args(
            chain_id=prev_strand.start_chain_id,
            resseq=prev_resseq,
            name=acceptor_name,
            icode=prev_strand.start_icode)
          bonded_atoms.append((donor2, acceptor2))
          cur_resseq += 2
          prev_resseq += 2
      else :
        raise RuntimeError("Strand sense must be either -1 or 1, except for "+
          "the first strand in a sheet.")
    return bonded_atoms

  def as_pdb_str (self) :
    assert len(self.strands) == len(self.registrations)
    lines = []
    for strand, reg in zip(self.strands, self.registrations) :
      format1 = "SHEET  %3d %3s%2d %3s %1s%4d%1s %3s %1s%4d%1s%2d"
      format2 = "%4s%3s %1s%4d%1s %4s%3s %1s%4d%1s"
      line = format1 % (strand.strand_id, self.sheet_id, self.n_strands,
        strand.start_resname, strand.start_chain_id, strand.start_resseq,
        strand.start_icode, strand.end_resname, strand.end_chain_id,
        strand.end_resseq, strand.end_icode, strand.sense)
      if reg is not None :
        line += " "
        line += format2 % (reg.cur_atom, reg.cur_resname, reg.cur_chain_id,
          reg.cur_resseq, reg.cur_icode, reg.prev_atom, reg.prev_resname,
          reg.prev_chain_id, reg.prev_resseq, reg.prev_icode)
      else :
        pass
        #assert strand.sense == 0
      lines.append(line.strip())
    return "\n".join(lines)

def parse_sheet_records (records) :
  sheets = []
  current_sheet = None
  current_sheet_id = None
  for line in records :
    if not line.startswith("SHEET") :
      continue
    line = "%-80s" % line # XXX: flex.split_lines strips each line
    sheet_id = line[11:14]
    n_strands = string.atoi(line[14:16])
    if sheet_id != current_sheet_id :
      if current_sheet is not None :
        # XXX: n_strands is frequently incorrect!
        assert (len(current_sheet.strands) == len(current_sheet.registrations))
        sheets.append(current_sheet)
      current_sheet = pdb_sheet(
        sheet_id=sheet_id,
        n_strands=n_strands,
        strands=[],
        registrations=[])
      current_sheet_id = sheet_id
    sense = string.atoi(line[38:40])
    current_strand = group_args(
      sheet_id=sheet_id,
      strand_id=string.atoi(line[7:10]),
      start_resname=line[17:20],
      start_chain_id=line[21],
      start_resseq=string.atoi(line[22:26]),
      start_icode=line[26],
      end_resname=line[28:31],
      end_chain_id=line[32],
      end_resseq=string.atoi(line[33:37]),
      end_icode=line[37],
      sense=sense)
    current_sheet.add_strand(current_strand)
    if sense == 0 :
      current_sheet.add_registration(None)
    else :
      if line[41:].strip() == "" :
        registration = None
      else :
        try :
          registration = group_args(
            cur_atom=line[41:45],
            cur_resname=line[45:48],
            cur_chain_id=line[49],
            cur_resseq=string.atoi(line[50:54]),
            cur_icode=line[54],
            prev_atom=line[56:60],
            prev_resname=line[60:63],
            prev_chain_id=line[64],
            prev_resseq=string.atoi(line[65:69]),
            prev_icode=line[69])
        except ValueError :
          registration = None
      current_sheet.add_registration(registration)
  if current_sheet is not None :
    sheets.append(current_sheet)
  return sheets

def restraint_groups_as_pdb_sheets (pdb_hierarchy, sheets, log=sys.stderr) :
  isel = pdb_hierarchy.atom_selection_cache().iselection
  atoms = [ a for a in pdb_hierarchy.atoms_with_labels() ]
  pdb_sheets = []
  for i, sheet in enumerate(sheets) :
    sheet_id = string.uppercase[i]
    if sheet.first_strand is None :
      print >> log, "Missing first strand in sheet %s" % sheet_id
    current_sheet = pdb_sheet(
      sheet_id=sheet_id,
      n_strands=1+len(sheet.strand),
      strands=[],
      registrations=[])
    first_strand = __strand_group_as_pdb_strand(isel=isel,
      selection=sheet.first_strand,
      atoms=atoms,
      log=log,
      sense=None)
    current_sheet.add_strand(first_strand)
    current_sheet.add_registration(None)
    base_sele = "(%s) and name N and (altloc 'A' or altloc ' ')"
    for strand in sheet.strand :
      pdb_strand = __strand_group_as_pdb_strand(isel=isel,
        selection=strand.selection,
        atoms=atoms,
        log=log,
        sense=strand.sense)
      current_sheet.add_strand(pdb_strand)
      s1 = base_sele % strand.bond_start_current
      s2 = base_sele % strand.bond_start_previous
      reg_curr_isel = isel(s1)
      reg_prev_isel = isel(s2)
      if reg_curr_isel.size() == 0 or reg_prev_isel.size() == 0 :
        current_sheet.add_registration(None)
        continue
      reg_curr_atom = atoms[reg_curr_isel[0]]
      reg_prev_atom = atoms[reg_prev_isel[0]]
      registration = group_args(
        cur_atom="N", #reg_curr_atom.name,
        cur_resname=reg_curr_atom.resname,
        cur_chain_id=reg_curr_atom.chain_id,
        cur_resseq=reg_curr_atom.resseq,
        cur_icode=reg_curr_atom.icode,
        prev_atom=reg_prev_atom.name,
        prev_resname="O", #reg_prev_atom.resname,
        prev_chain_id=reg_prev_atom.chain_id,
        prev_resseq=reg_prev_atom.resseq,
        prev_icode=reg_prev_atom.icode)
      current_sheet.add_registration(registration)
    pdb_sheets.append(current_sheet)
  return pdb_sheets

def __strand_group_as_pdb_strand (isel, selection, atoms, log, sense) :
  if sense is None or sense == "unknown" :
    int_sense = 0
  elif sense == "parallel" :
    int_sense = 1
  elif sense == "antiparallel" :
    int_sense = -1
  strand_isel = isel("(%s) and name N and (altloc 'A' or altloc ' ')" % (
    selection))
  start_atom = atoms[strand_isel[0]]
  end_atom = atoms[strand_isel[-1]]
  pdb_strand = group_args(
    sheet_id=sheet_id,
    strand_id=i+1,
    start_resname=start_atom.resname,
    start_chain_id=start_atom.chain_id,
    start_resseq=start_atom.resseq,
    start_icode=start_atom.icode,
    end_resname=end_atom.resname,
    end_chain_id=end_atom.chain_id,
    end_resseq=end_atom.resseq,
    end_icode=end_atom.icode,
    sense=int_sense)
  return pdb_strand

#-----------------------------------------------------------------------
def process_files (params, records=None) :
  assert len(params.file_name) > 0 or records is not None
  if records is None :
    records = []
    for file_name in params.file_name :
      lines = smart_open.for_reading(file_name).readlines()
      for line in lines :
        if line.startswith("HELIX") or line.startswith("SHEET") :
          records.append(line)
  secondary_structure = annotation(records=records)
  return secondary_structure

def run (args, out=sys.stdout, log=sys.stderr) :
  import iotbx.pdb
  master_phil = libtbx.phil.parse("""
  show_restraints = True
    .type = bool
  echo_pdb_records = False
    .type = bool
  echo_pymol_cmds = False
    .type = bool
  %s
""" % ss_input_params_str)
  user_phil = []
  for arg in args :
    if os.path.isfile(arg) :
      file_name = os.path.abspath(arg)
      base, ext = os.path.splitext(file_name)
      if ext in [".pdb", ".ent", ".gz"] :
        user_phil.append(libtbx.phil.parse("file_name=\"%s\"" % file_name))
      else :
        user_phil.append(libtbx.phil.parse(file_name=file_name))
    else :
      if arg.startswith("--") :
        arg = arg[2:] + "=True"
      try :
        cmdline_phil = libtbx.phil.parse(arg)
      except RuntimeError, e :
        print >> log, str(e)
      else :
        user_phil.append(cmdline_phil)
  working_phil = master_phil.fetch(sources=user_phil)
  params = working_phil.extract()
  if len(params.file_name) == 0 :
    raise Usage("Please supply at least one PDB file.")
  secondary_structure = process_files(params)
  if params.echo_pdb_records :
    print >> out, secondary_structure.as_pdb_str()
  elif params.echo_pymol_cmds :
    print >> out, secondary_structure.as_pymol_dashes(params)
  elif params.show_restraints :
    print >> out, secondary_structure.as_restraint_groups(log=log)

#-----------------------------------------------------------------------
def exercise_single () :
  from iotbx import file_reader
  from scitbx.array_family import flex
  from libtbx import test_utils
  # XXX: the PDB's annotation for 1ywf is simply wrong - the registers for the
  # last two strands are offset, leading to "bonds" with a distance of > 6 A.
  # The records below are correct.  However, the 3_10 helices will not result
  # in any H-bonds, since the only bonds present in the structure involve
  # adjacent residues not included in those specific HELIX records.  (ksdssp
  # simply ignores the 3_10 helix at 192-194 and combines it with the adjacent
  # alpha helices, which is also wrong.)
  ptpb_1ywf_records = """\
HELIX    1   1 ALA A   16  THR A   18  5                                   3
HELIX    2   2 ASP A   37  GLY A   48  1                                  12
HELIX    3   3 SER A   57  GLY A   65  1                                   9
HELIX    4   4 ASN A  119  PHE A  133  1                                  15
HELIX    5   5 PRO A  134  ARG A  136  5                                   3
HELIX    6   6 GLY A  138  ALA A  152  1                                  15
HELIX    7   7 ASP A  165  VAL A  178  1                                  14
HELIX    8   8 ASP A  181  ARG A  191  1                                  11
HELIX    9   9 SER A  192  ASP A  194  5                                   3
HELIX   10  10 SER A  195  GLN A  209  1                                  15
HELIX   11  11 ALA A  216  ALA A  225  1                                  10
HELIX   12  12 SER A  228  GLY A  233  1                                   6
HELIX   13  13 ARG A  235  GLY A  251  1                                  17
HELIX   14  14 SER A  252  ALA A  260  1                                   9
HELIX   15  15 SER A  263  LEU A  275  1                                  13
SHEET    1   A 5 ARG A  13  ASP A  14  0
SHEET    2   A 5 LEU A  27  SER A  30 -1  O  ARG A  29   N  ARG A  13
SHEET    3   A 5 VAL A 156  HIS A 159  1  O  VAL A 156   N  PHE A  28
SHEET    4   A 5 ASP A  51  ASP A  54  1  N  ALA A  51   O  LEU A 157
SHEET    5   A 5 ASP A  74  LEU A  77  1  O  HIS A  74   N  VAL A  52"""

  lines = flex.std_string()
  lines.extend(flex.split_lines(ptpb_1ywf_records))
  params = ss_input_params.extract()
  ss = annotation(records=lines)
  ss_out = ss.as_pdb_str()
  assert not test_utils.show_diff(ss_out, ptpb_1ywf_records)
  pml_out = ss.as_pymol_dashes(params=params)
  assert len(pml_out.splitlines()) == 109
  params.include_helices = False
  pml_out = ss.as_pymol_dashes(params=params)
  assert len(pml_out.splitlines()) == 11
  params.include_helices = True
  params.include_sheets = False
  pml_out = ss.as_pymol_dashes(params=params)
  assert len(pml_out.splitlines()) == 98
  assert (pml_out.splitlines()[0] ==
    """dist (chain 'A' and resi 41 and name H), (chain 'A' and resi 37 and name O)""")
  params.include_sheets = True
  assert len(ss.as_atom_selections(params=params)) == 20
  assert (ss.as_atom_selections(params=params)[0] ==
    """chain 'A' and resseq 16:18 and icode ' '""")
  print "OK"

def tst_pdb_file () :
  from iotbx import file_reader
  pdb_in = file_reader.any_file(file_name, force_type="pdb")
  old_ss = pdb_in.file_object.secondary_structure_section()
  structure = pdb_in.file_object.extract_secondary_structure()
  new_ss = structure.as_pdb_str()
  old_ss = "\n".join(old_ss)
  assert not test_utils.show_diff(new_ss, old_ss)

if __name__ == "__main__" :
  exercise_single()
