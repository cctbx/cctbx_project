from __future__ import division
#
# Implemented based on PDB v3.2 specification at:
#   http://www.wwpdb.org/documentation/format32/sect5.html
#
# NOTE: the routines for extracting hydrogen bonds are not used in phenix.
# I have left them here because they're simpler (and probably still accurate
# in most cases), and useful for pymol-related I/O.
#
# Oleg on 3-26-2014:
# I'm commenting them out because they contain some bugs and duplicate the
# existing functionality (tested and used) from proteins.py
# pymol-related I/O used in phenix.secondary_structure_restraints is
# implemented and used from
#   mmtbx.geometry_restraints.hbond.as_pymol_dashes(...)
#
# The great part of relevant code which is actually in use
# (generation of hydrogen bonds, phil records) written as additional methods
# for annotation class (phil generation) and as stand-alone functions for
# h-bonds generation are in:
#   mmtbx/secondary_structure/proteins.py
#
# It might be useful to move all the relevant secondary-structure code
# somewhere in one place.
#
# Oleg on 11-25-2014
# Carried out massive refactoring of secondary_structure objects.
# proteins.py now contains only h-bond generation. Probably this is also
# should be methods of secondary structure classes. Then proteins.py could be
# eliminated completely.
#
#
# Oleg on 5-27-2015
# The purpose of these classes - store as precisely as possible the
# knowledge about secondary structure itself, without connection to a particular
# model file. Therefore these classes are just containers of information with
# methods to convert it to various formats, e.g. pdb HELIX/SHEET records,
# phil files, cif format (in future). No connection to model means no i_seqs,
# no integer or bool selections etc. Object of annotation should be passed
# whenever the SS information is needed to the places where particular pdb
# model is available.
# The outcome from this rationale is:
# 1. No phil parameters in the module.
# 2. Essentially only constructors from various formats and exports to
#    various formats here
# 3. Convenient access to the various data to work with in different parts
#    of code.
#

from libtbx.utils import Sorry
import libtbx.phil
from libtbx import adopt_init_args
import sys

# This switch representation of seletions used for phil output from
# "resid 55 through 66" to "resseq 55:66"
use_resids = False # XXX: for debugging purposes only

class structure_base (object) :

  def as_pdb_str (self) :
    return None

  def __str__ (self) :
    return self.as_pdb_str()

  @staticmethod
  def parse_chain_id (chars) :
    assert len(chars) == 2
    if chars == "  " :
      return " "
    else :
      return chars.strip()

  @staticmethod
  def filter_helix_records(lines):
    result = []
    for line in lines:
      if line.startswith("HELIX"):
        result.append(line)
    return result

  @staticmethod
  def filter_and_split_sheet_records(lines):
    """Filter out only SHEET records and split them by sheet identifier.
    returns [[lines with equal sheetID], ... ,[lines with equal sheetID]]
    """
    result = []
    current_sh_lines = []
    current_sh_id = ""
    for line in lines:
      if line.startswith("SHEET"):
        line = "%-80s" % line # XXX: flex.split_lines strips each line
        sheet_id = line[11:14]
        if sheet_id == current_sh_id:
          current_sh_lines.append(line)
        else:
          if current_sh_lines != []:
            result.append(current_sh_lines)
            current_sh_lines = []
          current_sh_lines = [line]
          current_sh_id = sheet_id
    if current_sh_lines != []:
      result.append(current_sh_lines)
    return result

class annotation(structure_base):
  def __init__(self, helices=None, sheets=None):
    assert (not None in [helices, sheets])
    self.helices = helices
    self.sheets = sheets

  @classmethod
  def from_records(cls, records=None, log=None):
    "Initialize annotation from pdb HELIX/SHEET records"
    helices = []
    sheets = []
    for line in cls.filter_helix_records(records):
      try:
        h = pdb_helix.from_pdb_record(line)
      except ValueError:
        print >> log, "Bad HELIX record, was skipped:\n%s" % line
      else:
        helices.append(h)
    for sh_lines in cls.filter_and_split_sheet_records(records):
      try:
        sh = pdb_sheet.from_pdb_records(sh_lines)
      except ValueError:
        print >> log, "Bad SHEET records, was skipped:\n"
        for l in sh_lines:
          print >> log, "  %s" % l
      else:
        sheets.append(sh)

    return cls(helices=helices, sheets=sheets)

  @classmethod
  def from_phil(cls, phil_helices, phil_sheets, pdb_hierarchy, log=None):
    helices = []
    sheets = []
    for i, helix_param in enumerate(phil_helices):
      if helix_param.selection is not None:
        h = pdb_helix.from_phil_params(helix_param, pdb_hierarchy, i, log)
        helices.append(h)
    for i, sheet_param in enumerate(phil_sheets):
      if sheet_param.first_strand is not None:
        sh = pdb_sheet.from_phil_params(sheet_param, pdb_hierarchy, log)
        sheets.append(sh)
    return cls(helices=helices, sheets=sheets)

  def as_pdb_str (self) :
    records = []
    for helix in self.helices :
      records.append(helix.as_pdb_str())
    for sheet in self.sheets :
      records.append(sheet.as_pdb_str())
    return "\n".join(records)

  def as_restraint_groups (self, log=sys.stdout, prefix_scope="",
      add_segid=None) :
    phil_strs = []
    for helix in self.helices :
      helix_phil = helix.as_restraint_group(log, prefix_scope, add_segid)
      if helix_phil is not None :
        phil_strs.append(helix_phil)
    for sheet in self.sheets :
      sheet_phil = sheet.as_restraint_group(log, prefix_scope, add_segid)
      if sheet_phil is not None :
        phil_strs.append(sheet_phil)
    return "\n".join(phil_strs)

  def as_atom_selections (self):
    selections = []
    for helix in self.helices :
      try :
        selections.extend(helix.as_atom_selections())
      except RuntimeError, e :
        pass
    for sheet in self.sheets :
      selections.extend(sheet.as_atom_selections())
    return selections

  def overall_helix_selection (self) :
    selections = []
    for helix in self.helices :
      try :
        selections.extend(helix.as_atom_selections())
      except RuntimeError, e :
        pass
    return "(" + ") or (".join(selections) + ")"

  def overall_sheet_selection (self) :
    selections = []
    for sheet in self.sheets :
      try:
        selections.extend(sheet.as_atom_selections())
      except RuntimeError, e :
        pass
    return "(" + ") or (".join(selections) + ")"


  def as_bond_selections (self) :
    assert 0, "Probably is not used anywhere"
    bonded_atoms = self.extract_h_bonds(params)
    selections = []
    for (atom1, atom2) in bonded_atoms :
      selection_1 = "name %s and chain '%s' and resseq %d and icode '%s'" % (
        atom1.name, atom1.chain_id, atom1.resseq, atom1.icode)
      selection_2 = "name %s and chain '%s' and resseq %d and icode '%s'" % (
        atom2.name, atom2.chain_id, atom2.resseq, atom2.icode)
      selections.append((selection_1, selection_2))
    return selections

  def get_n_helices(self):
    return len(self.helices)

  def get_n_sheets(self):
    return len(self.sheets)



  def get_n_defined_hbonds(self):
    n_hb = 0
    if self.get_n_helices() > 0:
      for h in self.helices:
        n_hb += h.get_n_defined_hbonds()
    if self.get_n_sheets() > 0:
      for sh in self.sheets:
        n_hb += sh.get_n_defined_hbonds()
    return n_hb

#=============================================================================
#        88        88 88888888888 88          88 8b        d8
#        88        88 88          88          88  Y8,    ,8P
#        88        88 88          88          88   `8b  d8'
#        88aaaaaaaa88 88aaaaa     88          88     Y88P
#        88""""""""88 88"""""     88          88     d88b
#        88        88 88          88          88   ,8P  Y8,
#        88        88 88          88          88  d8'    `8b
#        88        88 88888888888 88888888888 88 8P        Y8
#=============================================================================

class pdb_helix (structure_base) :
  _helix_class_array = ['unknown','alpha', 'unknown', 'pi', 'unknown',
        '3_10', 'unknown', 'unknown', 'unknown', 'unknown', 'unknown']

  def __init__ (self,
        serial,
        helix_id,
        start_resname,
        start_chain_id,
        start_resseq,
        start_icode,
        end_resname,
        end_chain_id,
        end_resseq,
        end_icode,
        helix_class,
        comment,
        length,
        hbond_list=[], # list of (donor, acceptor) selecitons
        helix_selection=None,
        enabled=True,
        sigma=0.05,
        slack=0,
        top_out=False,
        ):
    adopt_init_args(self, locals())
    assert (length > 0), "Bad helix length"
    if isinstance(self.helix_class, int):
      self.helix_class = self._helix_class_array[helix_class]
    assert (self.helix_class in self._helix_class_array), \
        "Bad helix class: %s" % helix_class

    if isinstance(self.helix_id, int):
      self.helix_id = "%s" % self.helix_id
      self.helix_id = self.helix_id[:3]
    elif self.helix_id is None:
      self.helix_id = "%s" % self.serial
      self.helix_id = self.helix_id[:3]
    else:
      assert isinstance(self.helix_id, str)

  @classmethod
  def helix_class_to_int(cls, h_class):
    return cls._helix_class_array.index(h_class)

  @classmethod
  def helix_class_to_str(cls, h_class):
    return cls._helix_class_array[h_class]

  @classmethod
  def from_pdb_record(cls, line):
    "Raises ValueError in case of corrupted HELIX record!!!"
    return cls(
      serial=int(line[7:10]),
      helix_id=line[11:14].strip(),
      start_resname=line[15:18],
      start_chain_id=cls.parse_chain_id(line[18:20]),
      start_resseq=int(line[21:25]),
      start_icode=line[25],
      end_resname=line[27:30],
      end_chain_id=cls.parse_chain_id(line[30:32]),
      end_resseq=int(line[33:37]),
      end_icode=line[37],
      helix_class=cls.helix_class_to_str(int(line[38:40])),
      helix_selection=None,
      comment=line[40:70],
      length=int(line[71:76])) #string.atoi(line[71:76]))

  @classmethod
  def from_phil_params(cls, helix_params, pdb_hierarchy, serial=0, log=None):
    isel = pdb_hierarchy.atom_selection_cache().iselection
    atoms = [ a for a in pdb_hierarchy.atoms_with_labels() ]
    if helix_params.selection is None :
      print >> log, "Empty helix at serial %d." % (serial)
      # continue
    sele_str = ("(%s) and (name N) and (altloc 'A' or altloc ' ')" %
                helix_params.selection)
    if helix_params.serial_number is not None:
      serial = helix_params.serial_number
    amide_isel = isel(sele_str)
    start_atom = atoms[amide_isel[0]]
    end_atom = atoms[amide_isel[-1]]
    hbonds = []
    for hb in helix_params.hbond:
      if hb.donor is None:
        print >> log, "Donor selection in hbond cannot be None"
        continue
      if hb.acceptor is None:
        print >> log, "Acceptor selection in hbond cannot be None"
        continue
      hbonds.append((hb.donor, hb.acceptor))
    return cls(
      serial=serial,
      helix_id=helix_params.helix_identifier,
      start_resname=start_atom.resname,
      start_chain_id=start_atom.chain_id,
      start_resseq=int(start_atom.resseq),
      start_icode=start_atom.icode,
      end_resname=end_atom.resname,
      end_chain_id=end_atom.chain_id,
      end_resseq=int(end_atom.resseq),
      end_icode=end_atom.icode,
      helix_class=helix_params.helix_type,
      comment="",
      length=amide_isel.size(),
      hbond_list=hbonds,
      helix_selection=helix_params.selection,
      enabled=helix_params.enabled,
      sigma=helix_params.sigma,
      slack=helix_params.slack,
      top_out=helix_params.top_out,
      )

  def as_restraint_group(self, log=sys.stdout, prefix_scope="",
      add_segid=None, show_hbonds=False):
    if self.start_chain_id != self.end_chain_id :
      print >> log, "Helix chain ID mismatch: starts in %s, ends in %s" % (
        self.start_chain_id, self.end_chain_id)
      return None
    segid_extra = ""
    if add_segid is not None :
      segid_extra = "and segid '%s' " % add_segid
    if use_resids :
      resid_start = "%d%s" % (self.start_resseq, self.start_icode)
      resid_end = "%d%s" % (self.end_resseq, self.end_icode)
      sele = "chain '%s' %sand resid %s through %s" % (self.start_chain_id,
        segid_extra, resid_start, resid_end)
    else :
      sele = "chain '%s' %sand resseq %d:%d" % (self.start_chain_id,
        segid_extra, self.start_resseq, self.end_resseq)
    if prefix_scope != "" and not prefix_scope.endswith(".") :
      prefix_scope += "."
    serial_and_id = ""
    if self.serial is not None and self.serial > 0:
      serial_and_id += "\n  serial_number = %d" % self.serial
    if self.helix_id is not None:
      serial_and_id += "\n  helix_identifier = %s" % self.helix_id
    # if serial_and_id != "":
    #   serial_and_id += "\n"
    hbond_restr = ""
    if show_hbonds:
      if self.get_n_defined_hbonds() > 0:
        for hb in self.hbond_list:
          hb_str = "\n  hbond {\n    donor = \"%s\"\n    acceptor = \"%s\"\n  }" % (
            hb[0], hb[1])
          hbond_restr += hb_str
    rg = """\
%sprotein.helix {%s
  selection = "%s"
  helix_type = %s%s
}""" % (prefix_scope, serial_and_id, sele,
        self.helix_class, hbond_restr)
    return rg

  def as_pdb_str (self):
    def h_class_to_pdb_int(h_class):
      h_class_int = self.helix_class_to_int(h_class)
      if h_class_int == 0:
        return 1
      return h_class_int
    format = "HELIX  %3d %3s %3s%2s %4d%1s %3s%2s %4d%1s%2d%30s %5d"
    out = format % (self.serial,
      self.serial if self.helix_id is None else self.helix_id[:3],
      self.start_resname,
      self.start_chain_id, self.start_resseq, self.start_icode,
      self.end_resname, self.end_chain_id, self.end_resseq, self.end_icode,
      h_class_to_pdb_int(self.helix_class), self.comment, self.length)
    return out.strip()

  def continuity_check (self) :
    if self.start_icode != self.end_icode :
      raise RuntimeError("Don't know how to deal with helices with multiple "+
        "insertion codes ('%s' vs. '%s')." % (self.start_icode, self.end_icode))
    if self.start_chain_id != self.end_chain_id :
      raise RuntimeError("Don't know how to deal with helices with multiple "+
        "chain IDs ('%s' vs. '%s')." % (self.start_chain_id, self.end_chain_id))

  def as_atom_selections (self) :
    self.continuity_check()
    sele = "chain '%s' and resseq %d:%d and icode '%s'" % (self.start_chain_id,
      self.start_resseq, self.end_resseq, self.end_icode)
    return [sele]

  def get_n_defined_hbonds(self):
    if self.hbond_list is not None:
      return len(self.hbond_list)
    return 0

  def get_n_maximum_hbonds(self):
    if self.helix_class == 'alpha':
      return self.length-4
    elif self.helix_class=='3_10':
      return self.length-3
    elif helix.helix_class=='pi':
      return self.length-5
    elif helix.helix_class=='unknown':
      return 0
    else:
      # Should never happen
      assert 0, "Wrong helix_class creeped in object fields:" % self.helix_class

#=============================================================================
#       ad88888ba  88        88 88888888888 88888888888 888888888888
#      d8"     "8b 88        88 88          88               88
#      Y8,         88        88 88          88               88
#      `Y8aaaaa,   88aaaaaaaa88 88aaaaa     88aaaaa          88
#        `"""""8b, 88""""""""88 88"""""     88"""""          88
#              `8b 88        88 88          88               88
#      Y8a     a8P 88        88 88          88               88
#       "Y88888P"  88        88 88888888888 88888888888      88
#=============================================================================

class pdb_strand(structure_base):
  def __init__ (self,
      sheet_id,
      strand_id,
      start_resname,
      start_chain_id,
      start_resseq,
      start_icode,
      end_resname,
      end_chain_id,
      end_resseq,
      end_icode,
      sense) :
    adopt_init_args(self, locals())
    assert (sheet_id > 0) and (strand_id > 0)
    assert (sense in [-1, 0, 1]), "Bad sense."
    if start_icode != end_icode:
      raise Sorry("Different insertion codes for the beginning and the end of \
beta strand are not supported.")

  @classmethod
  def from_pdb_record(cls, line):
    return cls(sheet_id=line[11:14],
        strand_id=int(line[7:10]),
        start_resname=line[17:20],
        start_chain_id=cls.parse_chain_id(line[20:22]),
        start_resseq=int(line[22:26]),
        start_icode=line[26],
        end_resname=line[28:31],
        end_chain_id=cls.parse_chain_id(line[31:33]),
        end_resseq=int(line[33:37]),
        end_icode=line[37],
        sense=int(line[38:40]))


  def as_atom_selections(self):
    return "chain '%s' and resseq %d:%d and icode '%s'" % (
        self.start_chain_id, self.start_resseq, self.end_resseq,
        self.start_icode)

class pdb_strand_register(structure_base):
  def __init__ (self,
      cur_atom,
      cur_resname,
      cur_chain_id,
      cur_resseq,
      cur_icode,
      prev_atom,
      prev_resname,
      prev_chain_id,
      prev_resseq,
      prev_icode) :
    adopt_init_args(self, locals())

  @classmethod
  def from_pdb_record(cls, line):
    return cls(cur_atom=line[41:45],
        cur_resname=line[45:48],
        cur_chain_id=cls.parse_chain_id(line[48:50]),
        cur_resseq=int(line[50:54]),
        cur_icode=line[54],
        prev_atom=line[56:60],
        prev_resname=line[60:63],
        prev_chain_id=cls.parse_chain_id(line[63:65]),
        prev_resseq=int(line[65:69]),
        prev_icode=line[69])

class pdb_sheet(structure_base):
  def __init__ (self,
      sheet_id,
      n_strands,
      strands,
      registrations,
      hbond_list=[], # list of (donor, acceptor) selecitons
      ):
    adopt_init_args(self, locals())
    if isinstance(self.sheet_id, int):
      self.sheet_id = "%s" % self.sheet_id
      self.sheet_id = self.sheet_id[:3]
    else:
      assert isinstance(self.sheet_id, str)

  @classmethod
  def from_pdb_records(cls,records):
    "Raises ValueError in case of bad SHEET record"
    assert len(records) > 0
    sheet_id = records[0][11:14]
    n_strands = int(records[0][14:16])
    if n_strands != len(records):
      # correcting wrong n_strands
      n_strands = len(records)
    strands = []
    registrations = []
    for r in records:
      s = pdb_strand.from_pdb_record(r)
      reg = None
      sense = int(r[38:40])
      if sense != 0:
        reg = pdb_strand_register.from_pdb_record(r)
      strands.append(s)
      registrations.append(reg)
    return cls(sheet_id=sheet_id,
               n_strands=n_strands,
               strands=strands,
               registrations=registrations)

  @classmethod
  def from_phil_params(cls, sheet_params, pdb_hierarchy, log=None):
    isel = pdb_hierarchy.atom_selection_cache().iselection
    atoms = [a for a in pdb_hierarchy.atoms_with_labels()]
    if sheet_params.first_strand is None:
      raise Sorry("Empty first strand selection")
    sheet_id="1"
    if sheet_params.sheet_id is not None:
      sheet_id =  "%3s" % sheet_params.sheet_id[:3]
    n_strands = len(sheet_params.strand) + 1
    sele_str = ("(%s) and (name N) and (altloc 'A' or altloc ' ')" %
                    sheet_params.first_strand)
    amide_isel = isel(sele_str)
    start_atom = atoms[amide_isel[0]]
    end_atom = atoms[amide_isel[-1]]
    first_strand = pdb_strand(
        sheet_id=sheet_id,
        strand_id=1,
        start_resname=start_atom.resname,
        start_chain_id=start_atom.chain_id,
        start_resseq=int(start_atom.resseq),
        start_icode=start_atom.icode,
        end_resname=end_atom.resname,
        end_chain_id=end_atom.chain_id,
        end_resseq=int(end_atom.resseq),
        end_icode=end_atom.icode,
        sense=0)
    strands = [first_strand]
    registrations = [None]
    for i, strand_param in enumerate(sheet_params.strand):
      sele_str = ("(%s) and (name N) and (altloc 'A' or altloc ' ')" %
                      strand_param.selection)
      amide_isel = isel(sele_str)
      start_atom = atoms[amide_isel[0]]
      end_atom = atoms[amide_isel[-1]]
      sense = cls.sense_to_int(strand_param.sense)
      strand = pdb_strand(
          sheet_id=sheet_id,
          strand_id=i+2,
          start_resname=start_atom.resname,
          start_chain_id=start_atom.chain_id,
          start_resseq=int(start_atom.resseq),
          start_icode=start_atom.icode,
          end_resname=end_atom.resname,
          end_chain_id=end_atom.chain_id,
          end_resseq=int(end_atom.resseq),
          end_icode=end_atom.icode,
          sense=sense)
      reg_cur_sel = isel(strand_param.bond_start_current)
      if len(reg_cur_sel) == 0:
        raise Sorry("This bond_start_current yields 0 atoms:\n %s" % strand_param.bond_start_current)
      reg_cur_atom = atoms[reg_cur_sel[0]]
      reg_prev_sel = isel(strand_param.bond_start_previous)
      if len(reg_prev_sel) == 0:
        raise Sorry("This bond_start_current yields 0 atoms:\n %s" % strand_param.bond_start_previous)
      reg_prev_atom = atoms[reg_prev_sel[0]]
      reg = pdb_strand_register(
          cur_atom=reg_cur_atom.name,
          cur_resname=reg_cur_atom.resname,
          cur_chain_id=reg_cur_atom.chain_id,
          cur_resseq=int(reg_cur_atom.resseq),
          cur_icode=reg_cur_atom.icode,
          prev_atom=reg_prev_atom.name,
          prev_resname=reg_prev_atom.resname,
          prev_chain_id=reg_prev_atom.chain_id,
          prev_resseq=int(reg_prev_atom.resseq),
          prev_icode=reg_prev_atom.icode)
      strands.append(strand)
      registrations.append(reg)
    hbonds = []
    for hb in sheet_params.hbond:
      if hb.donor is None:
        print >> log, "Donor selection in hbond cannot be None"
        continue
      if hb.acceptor is None:
        print >> log, "Acceptor selection in hbond cannot be None"
        continue
      hbonds.append((hb.donor, hb.acceptor))

    return cls(sheet_id=sheet_id,
               n_strands=n_strands,
               strands=strands,
               registrations=registrations,
               hbond_list=hbonds)
  @staticmethod
  def sense_to_int(str_sense):
    sense = 0
    if str_sense == "parallel":
      sense = 1
    elif str_sense == "antiparallel":
      sense = -1
    return sense

  def add_strand (self, strand) :
    self.strands.append(strand)

  def add_registration (self, registration) :
    self.registrations.append(registration)

  def as_atom_selections (self) :
    strand_selections = []
    for strand in self.strands :
      strand_selections.append(strand.as_atom_selections())
    return strand_selections

  def get_n_defined_hbonds(self):
    if self.hbond_list is not None:
      return len(self.hbond_list)
    return 0

  def as_pdb_str (self) :
    assert len(self.strands) == len(self.registrations)
    lines = []
    # print self.strands
    # print self.registrations
    # STOP()
    for strand, reg in zip(self.strands, self.registrations) :
      format1 = "SHEET  %3d %3s%2d %3s%2s%4d%1s %3s%2s%4d%1s%2d"
      format2 = "%4s%3s%2s%4d%1s %4s%3s%2s%4d%1s"
      # print "STRAND, REG", strand, reg
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

  def particular_strand_as_pdb_str(self, strand_id=None):
    if strand_id is None:
      return ""
    line = ""
    for strand, reg in zip(self.strands, self.registrations):
      if strand.strand_id == strand_id:
        format1 = "SHEET  %3d %3s%2d %3s%2s%4d%1s %3s%2s%4d%1s%2d"
        format2 = "%4s%3s%2s%4d%1s %4s%3s%2s%4d%1s"
        # print "STRAND, REG", strand, reg
        line = format1 % (strand.strand_id, self.sheet_id, self.n_strands,
          strand.start_resname, strand.start_chain_id, strand.start_resseq,
          strand.start_icode, strand.end_resname, strand.end_chain_id,
          strand.end_resseq, strand.end_icode, strand.sense)
        if reg is not None :
          line += " "
          line += format2 % (reg.cur_atom, reg.cur_resname, reg.cur_chain_id,
            reg.cur_resseq, reg.cur_icode, reg.prev_atom, reg.prev_resname,
            reg.prev_chain_id, reg.prev_resseq, reg.prev_icode)
    return line

  def as_restraint_group(self, log=sys.stdout, prefix_scope="",
      add_segid=None, show_hbonds=False):
    if len(self.strands) == 0 :
      return None
    selections = []
    senses = []
    reg_curr = []
    reg_prev = []
    segid_extra = ""
    if add_segid is not None :
      segid_extra = "and segid '%s' " % add_segid
    for (strand,registration) in zip(self.strands, self.registrations) :
      if use_resids :
        resid_start = "%d%s" % (strand.start_resseq, strand.start_icode)
        resid_end = "%d%s" % (strand.end_resseq, strand.end_icode)
        sele = "chain '%s' %sand resid %s through %s" % (strand.start_chain_id,
          segid_extra, resid_start, resid_end)
      else :
        sele = "chain '%s' %sand resseq %d:%d" % (strand.start_chain_id,
          segid_extra, strand.start_resseq, strand.end_resseq)
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
        if use_resids :
          sele_base = "chain '%s' %sand resid %s and name %s"
          resid_curr = "%d%s" % (registration.cur_resseq,registration.cur_icode)
          resid_prev = "%d%s" % (registration.prev_resseq,
              registration.prev_icode)
          reg_curr.append(sele_base % (registration.cur_chain_id,segid_extra,
              resid_curr, registration.cur_atom.strip()))
          reg_prev.append(sele_base % (registration.prev_chain_id,segid_extra,
              resid_prev, registration.prev_atom.strip()))
        else :
          reg_curr.append("chain '%s' %sand resseq %d and name %s" % (
              registration.cur_chain_id, segid_extra, registration.cur_resseq,
              registration.cur_atom.strip()))
          reg_prev.append("chain '%s' %sand resseq %d and name %s" % (
              registration.prev_chain_id, segid_extra, registration.prev_resseq,
              registration.prev_atom.strip()))
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
    selection = "%s"
    sense = %s
    bond_start_current = "%s"
    bond_start_previous = "%s"
  }""" % (sele, sense, curr, prev))
      n += 1
    assert first_strand is not None
    if prefix_scope != "" and not prefix_scope.endswith(".") :
      prefix_scope += "."
    hbond_restr = ""
    if show_hbonds:
      if self.get_n_defined_hbonds() > 0:
        for hb in self.hbond_list:
          hb_str = "\n  hbond {\n    donor = \"%s\"\n    acceptor = \"%s\"\n  }" % (
            hb[0], hb[1])
          hbond_restr += hb_str
    phil_str = """
%sprotein.sheet {
  sheet_id = "%s"
  first_strand = "%s"
%s%s
}""" % (prefix_scope, self.sheet_id, first_strand, "\n".join(strands), hbond_restr)
    return phil_str

if __name__ == "__main__" :
  pass
