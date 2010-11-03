#
# XXX: see also mmtbx.monomer_library.secondary_structure
#
# Implemented based on PDB v3.2 specification at:
#   http://www.wwpdb.org/documentation/format32/sect5.html
#
# NOTE: the routines for extracting hydrogen bonds are not used in phenix.
# I have left them here because they're simpler (and probably still accurate
# in most cases), and useful for pymol-related I/O.

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
    .style = hidden
  use_hydrogens = True
    .type = bool
    .style = hidden
  include_helices = True
    .type = bool
    .style = hidden
  include_sheets = True
    .type = bool
    .style = hidden
"""
ss_input_params = libtbx.phil.parse(ss_input_params_str)

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
    selections = []
    for sheet in self.sheets :
      try:
        selections.extend(sheet.as_atom_selections(params))
      except RuntimeError, e :
        pass
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

#-----------------------------------------------------------------------
class pdb_helix (structure_base, group_args) :
  def as_pdb_str (self) :
    format = "HELIX  %3d %3s %3s%2s %4d%1s %3s%2s %4d%1s%2d%30s %5d"
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

def parse_chain_id (chars) :
  assert len(chars) == 2
  if chars == "  " :
    return " "
  else :
    return chars.strip()

def parse_helix_records (records) :
  helices = []
  for line in records :
    if not line.startswith("HELIX") :
      continue
    try :
      length = string.atoi(line[71:76])
    except ValueError :
      length = 0
    current_helix = pdb_helix(
      serial=string.atoi(line[7:10]),
      helix_id=line[11:14].strip(),
      start_resname=line[15:18],
      start_chain_id=parse_chain_id(line[18:20]),
      start_resseq=string.atoi(line[21:25]),
      start_icode=line[25],
      end_resname=line[27:30],
      end_chain_id=parse_chain_id(line[30:32]),
      end_resseq=string.atoi(line[33:37]),
      end_icode=line[37],
      helix_class=string.atoi(line[38:40]),
      comment=line[40:70],
      length=length) #string.atoi(line[71:76]))
    helices.append(current_helix)
  return helices

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

  def extract_h_bonds_as_i_seqs (self) :
    pass

  def as_pdb_str (self) :
    assert len(self.strands) == len(self.registrations)
    lines = []
    for strand, reg in zip(self.strands, self.registrations) :
      format1 = "SHEET  %3d %3s%2d %3s%2s%4d%1s %3s%2s%4d%1s%2d"
      format2 = "%4s%3s%2s%4d%1s %4s%3s%2s%4d%1s"
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
      start_chain_id=parse_chain_id(line[20:22]),
      start_resseq=string.atoi(line[22:26]),
      start_icode=line[26],
      end_resname=line[28:31],
      end_chain_id=parse_chain_id(line[31:33]),
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
            cur_chain_id=parse_chain_id(line[48:50]),
            cur_resseq=string.atoi(line[50:54]),
            cur_icode=line[54],
            prev_atom=line[56:60],
            prev_resname=line[60:63],
            prev_chain_id=parse_chain_id(line[63:65]),
            prev_resseq=string.atoi(line[65:69]),
            prev_icode=line[69])
        except ValueError :
          registration = None
      current_sheet.add_registration(registration)
  if current_sheet is not None :
    sheets.append(current_sheet)
  return sheets

#-----------------------------------------------------------------------
def process_records (records=None, pdb_files=None, allow_none=True) :
  assert records is not None or pdb_files is not None
  if records is None :
    records = []
    for file_name in pdb_files :
      lines = smart_open.for_reading(file_name).readlines()
      for line in lines :
        if line.startswith("HELIX") or line.startswith("SHEET") :
          records.append(line)
  if len(records) == 0 :
    if allow_none :
      return None
    else :
      raise Sorry("No secondary structure found.")
  secondary_structure = annotation(records=records)
  return secondary_structure

def run (args, out=sys.stdout, log=sys.stderr, cmd_params_str="") :
  import iotbx.pdb
  master_phil = libtbx.phil.parse("""%s
  echo_pdb_records = False
    .type = bool
  echo_pymol_cmds = False
    .type = bool
  %s
""" % (cmd_params_str, ss_input_params_str))
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
  secondary_structure = process_records(pdb_files=params.file_name)
  if secondary_structure is not None :
    if params.echo_pdb_records :
      print >> out, secondary_structure.as_pdb_str()
    elif params.echo_pymol_cmds :
      print >> out, secondary_structure.as_pymol_dashes(params)
  return secondary_structure

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
  two_char_chain_records = """\
HELIX    1   1 THRA1   11  ASPA1   39  1                                  29
HELIX    2   2 GLUA1   46  ARGA1   73  1                                  28
HELIX    3   3 THRA1   93  ALAA1  120  1                                  28
HELIX    4   4 PROA1  124  HISA1  133  1                                  10
HELIX    5   5 LEUA1  135  ARGA1  154  1                                  20
HELIX    6   6 THRa1   11  TYRa1   37  1                                  27
HELIX    7   7 GLUa1   46  GLNa1   72  1                                  27
HELIX    8   8 THRa1   93  ALAa1  120  1                                  28
HELIX    9   9 PROa1  124  HISa1  133  1                                  10"""
  lines = flex.std_string()
  lines.extend(flex.split_lines(two_char_chain_records))
  ss = annotation(records=lines)
  ss_out = ss.as_pdb_str()
  assert not test_utils.show_diff(ss_out, two_char_chain_records)
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
