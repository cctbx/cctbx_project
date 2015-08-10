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


def segments_are_similar(atom_selection_1=None,
   atom_selection_2=None,
   hierarchy=None,
   maximum_length_difference=None,
   minimum_overlap=None):

    "Return True if two segments selected are similar"
    "Different lengths optional."
    if minimum_overlap is None:
      minimum_overlap=4 # set default

    # Check to see if they overlap sufficiently and have similar lengths

    asc=hierarchy.atom_selection_cache()
    sel=asc.selection(string = atom_selection_1)
    try:
      h1=hierarchy.deep_copy().select(sel)  # keep original hierarchy too
      number_self=h1.overall_counts().n_residues
    except Exception, e:
      return False

    asc=hierarchy.atom_selection_cache()
    sel=asc.selection(string = atom_selection_2)
    try:
      h2=hierarchy.deep_copy().select(sel)
      number_other=h2.overall_counts().n_residues
    except Exception, e:
      return False

    if maximum_length_difference is not None and \
        abs(number_self-number_other)>maximum_length_difference:
      return False

    asc=hierarchy.atom_selection_cache()
    atom_selection="(%s) and (%s)" %(atom_selection_1,atom_selection_2)
    sel=asc.selection(string = atom_selection)
    try:
      h12=hierarchy.deep_copy().select(sel)
      number_both=h12.overall_counts().n_residues
    except Exception, e:
      return False

    if number_both<minimum_overlap:
      return False

    return True


class one_strand_pair_registration_atoms:
  def __init__(self,
      strand_a_atom=None,
      strand_b_atom=None,
      strand_a_position=None,
      strand_b_position=None,
      sense=None,
      log=sys.stdout):

    self.ok=True
    self.sense=sense
    self.std_position_a=None
    self.std_position_b=None

    #  If CA residue i of strand n matches with residue i' of strand n+1:
    #  For both parallel and antiparallel: N,O alternately H-bond to strands
    #    to right and to left

    #  For parallel strands:
    #  O of residue i in strand n H-bonds to N of residue i'+1 in strand n+1.
    #  N of residue i in strand n H-bonds to O of residue i'-1 in strand n+1.
    #  So... O 31 -- N 72  is equivalent to: (i=31  i'=71)
    #        O 33 -- N 74  (just naming residues 2 away)
    #        N 31 -- O 70  (naming N instead of O)
    #  and not equivalent to:
    #        O 32 -- N 73  (points the other way)

    #  For antiparallel strands:
    #  O of residue i in strand n H-bonds to N of residue i' in strand n+1
    #  N of residue i in strand n H-bonds to O of residue i' in strand n+1

    if not strand_a_atom or not strand_b_atom:
      return

    strand_a_atom=strand_a_atom.replace(" ","")
    strand_b_atom=strand_b_atom.replace(" ","")

    if (strand_a_atom=="O" and strand_b_atom != "N") or \
       (strand_a_atom=="N" and strand_b_atom != "O") or \
       (strand_a_atom!="N" and strand_b_atom != "N") or \
       (strand_a_atom!="O" and strand_b_atom != "O"):
          print >>log,"Cannot interpret H-bonding of %s to %s " %(
            strand_a_atom,strand_b_atom)
          ok=False
          return

    if not sense in [-1,1]:
      print >>log,"Cannot interpret bonding with sense not -1 or 1:"
      self.ok=False
      return

    if strand_a_position is None or strand_b_position is None:
      self.ok=False # can't do anything
      return

    if sense==1:  # parallel
      if strand_a_atom=="O":
        self.std_position_a=strand_a_position
        self.std_position_b=strand_b_position
      else:
        self.std_position_a=strand_a_position
        self.std_position_b=strand_b_position+2  # equiv to O bonding here
    else:  # antiparallel
      self.std_position_a=strand_a_position
      self.std_position_b=strand_b_position

  def is_equivalent_to(self,other=None):
    if not self.ok: return None  # can't tell

    #  Equivalent if std_position_a and std_position_b are the same for both
    #  or both offset by the same multiple of 2 (opposite multiple of 2 if
    #  antiparallel)

    if self.std_position_a is None and other.std_position_a is None:
      return None  # both are None
    elif self.std_position_a is None or other.std_position_a is None:
      return False # one is None
    if self.std_position_b is None and other.std_position_b is None:
      return None  # both are None
    elif self.std_position_b is None or other.std_position_b is None:
      return False # one is None

    delta_a=self.std_position_a-other.std_position_a
    if delta_a-2*(delta_a//2)!=0:  # must be multiple of 2
      return False
    delta_b=self.std_position_b-other.std_position_b
    if delta_b-2*(delta_b//2)!=0:  # must be multiple of 2
      return False
    if self.sense==1 and delta_a != delta_b:
      return False
    elif self.sense==-1 and delta_a != -delta_b:
      return False

    return True # all ok

class registration_atoms:
  def __init__(self,hierarchy=None,
      strand_1a=None,strand_1b=None,
      strand_2a=None,strand_2b=None,
      registration_1=None,registration_2=None,
      sense=None,
      log=sys.stdout):
    adopt_init_args(self, locals())

    self.ok=True
    # 1a,1b are the prev and cur strands for registration 1; 2a,2b for reg 2

    ok_a,strand_1a_position,strand_2a_position,\
      strand_1a_atom,strand_2a_atom=\
      self.get_positions_and_atoms_for_one_strand(prev=True)
    if not ok_a: ok=False

    ok_b,strand_1b_position,strand_2b_position, \
      strand_1b_atom,strand_2b_atom=\
      self.get_positions_and_atoms_for_one_strand(prev=False)
    if not ok_b: ok=False

    # pair_1 are the atoms and positions for one pair of strands
    registration_pair_1=one_strand_pair_registration_atoms(
      strand_a_atom=strand_1a_atom,
      strand_b_atom=strand_1b_atom,
      strand_a_position=strand_1a_position,
      strand_b_position=strand_1b_position,
      sense=sense,log=log)
    if not registration_pair_1.ok: self.ok=False


    registration_pair_2=one_strand_pair_registration_atoms(
      strand_a_atom=strand_2a_atom,
      strand_b_atom=strand_2b_atom,
      strand_a_position=strand_2a_position,
      strand_b_position=strand_2b_position,
      sense=sense,log=log)
    if not registration_pair_2.ok: self.ok=False

    if not self.ok:
      self.pairs_are_equivalent=None
    elif registration_pair_1.is_equivalent_to(registration_pair_2):
      self.pairs_are_equivalent=True
    else:
      self.pairs_are_equivalent=False

  def get_positions_and_atoms_for_one_strand(self,prev=True):
    if prev:
      atom_selection="( (%s) or (%s)) and name ca" %(
       self.strand_1a.as_atom_selections(), self.strand_2a.as_atom_selections())
    else:
      atom_selection="( (%s) or (%s)) and name ca" %(
       self.strand_1b.as_atom_selections(), self.strand_2b.as_atom_selections())
    asc=self.hierarchy.atom_selection_cache()
    sel=asc.selection(string = atom_selection)
    ph=self.hierarchy.deep_copy().select(sel)  # keep original hierarchy too
    # now ph is list of CA by residue.  Find the position of our N or O

    if not self.registration_1 or not self.registration_2:
      return False,None,None,None,None

    strand_1_position=self.get_registration_position(
      hierarchy=ph,
      registration=self.registration_1,
      prev=prev)

    strand_2_position=self.get_registration_position(
      hierarchy=ph,
      registration=self.registration_2,
      prev=prev)

    if prev:
      strand_1_atom=self.registration_1.prev_atom
      strand_2_atom=self.registration_2.prev_atom
    else:
      strand_1_atom=self.registration_1.cur_atom
      strand_2_atom=self.registration_2.cur_atom

    ok=None
    if strand_1_position is None and strand_2_position is not None:
      ok=False
    if strand_2_position is None and strand_1_position is not None:
      ok=False
    if strand_1_position is not None and strand_2_position is not None:
      ok=True

    return ok,strand_1_position,strand_2_position,strand_1_atom,strand_2_atom

  def get_registration_position(self,hierarchy=None,
      registration=None,prev=False):
    if not registration:
      return None

    if prev:
      resname=registration.prev_resname
      resseq=registration.prev_resseq
      chain_id=registration.prev_chain_id
      icode=registration.prev_icode
    else:
      resname=registration.cur_resname
      resseq=registration.cur_resseq
      chain_id=registration.cur_chain_id
      icode=registration.cur_icode
    # Assumes that there are no alternate conformations
    pos=-1
    for model in hierarchy.models()[:1]:
      for chain in model.chains()[:1]:
        for conformer in chain.conformers()[:1]:
          for residue in conformer.residues():
            pos+=1
            if residue.resname==resname and \
               residue.resseq.replace(" ","")==str(resseq) and \
               chain.id==chain_id and residue.icode==icode:
              return pos
    return None



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

  def count_residues(self,hierarchy=None):
    "Count residues in this secondary structure"

    if hierarchy is None:
      raise AssertionError,"Require hierarchy for count_residues"

    atom_selection=self.combine_atom_selections(self.as_atom_selections())
    if not atom_selection:
       return 0

    asc=hierarchy.atom_selection_cache()
    sel = asc.selection(string = atom_selection)
    ph=hierarchy.deep_copy().select(sel)
    return ph.overall_counts().n_residues

  def count_h_bonds(self,hierarchy=None,
       max_h_bond_length=None):
    "Count good and poor H-bonds in this secondary structure"

    if hierarchy is None:
      raise AssertionError,"Require hierarchy for count_h_bonds"

    atom_selection=self.combine_atom_selections(self.as_atom_selections())
    if not atom_selection:
       return 0,0

    asc=hierarchy.atom_selection_cache()
    sel = asc.selection(string = atom_selection)
    ph=hierarchy.deep_copy().select(sel)

    from mmtbx.secondary_structure.find_ss_from_ca import \
      find_secondary_structure
    from libtbx.utils import null_out
    fss=find_secondary_structure(hierarchy=ph,
      user_annotation_text=self.as_pdb_str(),
      force_secondary_structure_input=True,
      combine_annotations=False,
      max_h_bond_length=max_h_bond_length,out=null_out())
    return fss.number_of_good_h_bonds,fss.number_of_poor_h_bonds


  def is_similar_to(self,other=None,hierarchy=None,
      maximum_length_difference=None, minimum_overlap=None):
    "Return True if this annotation is similar to other."
    "Allow equivalent representations (different alignment atoms with)"
    "equivalent results, backwards order of strands in sheet, different "
    "lengths of strands optional."
    if self.is_same_as(other=other): return True  # always quicker

    if hierarchy is None:
      raise AssertionError,"Require hierarchy for is_similar_to"

    # Check for different objects
    if type(self) != type(other): return False

    self_atom_selections=self.as_atom_selections()
    other_atom_selections=other.as_atom_selections()
    if type(self_atom_selections)==type([1,2,3]):
       assert len(self_atom_selections)==1
       self_atom_selections=self_atom_selections[0]
    if type(other_atom_selections)==type([1,2,3]):
       assert len(other_atom_selections)==1
       other_atom_selections=other_atom_selections[0]
    return segments_are_similar(
        atom_selection_1=self_atom_selections,
        atom_selection_2=other_atom_selections,
        hierarchy=hierarchy,
        maximum_length_difference=maximum_length_difference,
        minimum_overlap=minimum_overlap)

  def combine_atom_selections(self,atom_selections_list,require_all=False):
    if require_all:
      joiner="and"
    else:
      joiner="or"

    all_list=[]
    for x in atom_selections_list:
      if type(x)==type([1,2,3]):
        all_list+=x
      else:
        all_list.append(x)
    if not all_list: return ""

    atom_selection="(%s)" %(all_list[0])
    for x in all_list[1:]:
      atom_selection+=" %s (%s) " %(joiner,x)
    return atom_selection

  def overlaps_with(self,other=None,hierarchy=None):
    assert hierarchy
    s1=self.as_atom_selections()
    s2=other.as_atom_selections()
  
    atom_selection=self.combine_atom_selections([s1,s2],require_all=True)
    if atom_selection:
      asc=hierarchy.atom_selection_cache()
      sel = asc.selection(string = atom_selection)
      ph=hierarchy.deep_copy().select(sel)  #ph = needed part of hierarchy
      if ph.overall_counts().n_residues>0:
        return True
    return False


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
      except ValueError, e:
        print >> log, "Bad SHEET records, were skipped:\n"
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

  def multiply_to_asu(self, ncs_copies_chain_names, n_copies):
    from iotbx.ncs.ncs_preprocess import format_num_as_str
    import copy
    # ncs_copies_chain_names = {'chain B_022':'CD' ...}
    def make_key(chain_id, n):
      return "chain %s_%s" % (chain_id, format_num_as_str(n))
    def get_new_chain_id(old_chain_id, n_copy):
      new_chain_id = ncs_copies_chain_names.get(
          make_key(old_chain_id,n_copy), None)
      if new_chain_id is None:
        raise Sorry("Something went wrong in mulitplying HELIX/SHEET records")
      return new_chain_id
    new_helices = []
    new_sheets = []
    new_h_serial = 1
    new_sheet_id = 1
    for n_copy in xrange(1, n_copies+1):
      for helix in self.helices:
        new_helix = copy.deepcopy(helix)
        new_helix.erase_hbond_list()
        new_helix.set_new_chain_ids(
            get_new_chain_id(new_helix.start_chain_id, n_copy))
        new_helix.set_new_serial(new_h_serial, adopt_as_id=True)
        new_h_serial += 1
        new_helices.append(new_helix)
      for sheet in self.sheets:
        new_sheet = copy.deepcopy(sheet)
        new_sheet.sheet_id = new_sheet_id
        new_sheet_id += 1
        new_sheet.erase_hbond_list()
        for strand in new_sheet.strands:
          strand.set_new_chain_ids(
              get_new_chain_id(strand.start_chain_id, n_copy))
        for reg in new_sheet.registrations:
          if reg is not None:
            reg.set_new_cur_chain_id(
                get_new_chain_id(reg.cur_chain_id, n_copy))
            reg.set_new_prev_chain_id(
                get_new_chain_id(reg.prev_chain_id, n_copy))
        new_sheets.append(new_sheet)
    self.helices = new_helices
    self.sheets = new_sheets

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

  def as_atom_selections(self, add_segid=None):
    selections = []
    for helix in self.helices :
      try :
        selections.extend(helix.as_atom_selections(add_segid=add_segid))
      except RuntimeError, e :
        pass
    for sheet in self.sheets :
      selections.extend(sheet.as_atom_selections(add_segid=add_segid))
    return selections

  def overall_selection(self,add_segid=None):
    selections = []
    for helix in self.helices:
      try :
        selections.extend(helix.as_atom_selections(add_segid=add_segid))
      except RuntimeError, e :
        pass
    for sheet in self.sheets:
      try:
        selections.extend(sheet.as_atom_selections(add_segid=add_segid))
      except RuntimeError, e :
        pass
    return "(" + ") or (".join(selections) + ")"

  def overall_helix_selection(self, add_segid=None):
    selections = []
    for helix in self.helices:
      try :
        selections.extend(helix.as_atom_selections(add_segid=add_segid))
      except RuntimeError, e :
        pass
    print "selections", selections
    return "(" + ") or (".join(selections) + ")"

  def overall_sheet_selection(self, add_segid=None):
    selections = []
    for sheet in self.sheets:
      try:
        selections.extend(sheet.as_atom_selections(add_segid=add_segid))
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

  def split_sheets(self):
    "Split all multi-strand sheets into 2-strand sheets"
    new_sheets=[]
    from copy import deepcopy
    for sheet in self.sheets:
      new_sheets+=sheet.split(starting_sheet_id_number=len(new_sheets)+1)
    return annotation(
      helices=deepcopy(self.helices),
      sheets=new_sheets)

  def merge_sheets(self):
    "Group 2-strand sheets into larger sheets if identical component strands"
    # Assumes that all the sheets are non-overlapping
    # First run sheet.split() on all sheets or split_sheets on the annotation
    #  as this requires 2-strand sheets (not 3 or more)
    from copy import deepcopy

    sheet_pointer_0={}
    sheet_pointer_1={}
    all_strands=[]
    for sheet in self.sheets:
      assert len(sheet.strands)==2 and len(sheet.registrations)==2
      strand_0_text=sheet.strands[0].as_atom_selections()
      strand_1_text=sheet.strands[1].as_atom_selections()
      sheet_pointer_0[strand_0_text]=sheet
      sheet_pointer_1[strand_1_text]=sheet
      if not sheet.strands[0] in all_strands:
        all_strands.append(sheet.strands[0])
      if not sheet.strands[1] in all_strands:
        all_strands.append(sheet.strands[1])

    used_strand_selections=[] # strands that have been used
    # Start with a strand that is in position 1 only (if any) and work through
    #   all strands in that sheet. Then on second iteration take all unused
    new_sheets=[]
    for iter in [0,1]:
      for strand in all_strands:
        key=strand.as_atom_selections()
        if key in used_strand_selections: continue
        if key in sheet_pointer_0.keys() and (
            iter==1 or not key in sheet_pointer_1.keys()):
          used_strand_selections.append(key)
          working_sheet=deepcopy(sheet_pointer_0.get(key))
          new_sheets.append(working_sheet)  # now we will extend this sheet

          second_strand=working_sheet.strands[1] # second strand
          next_key=second_strand.as_atom_selections()
          if next_key in used_strand_selections: continue

          while next_key: # points to next sheet to merge 2nd strand
            # find sheet where next_strand is the first strand
            next_sheet=sheet_pointer_0.get(next_key)
            if not next_sheet: break
            used_strand_selections.append(next_key)

            next_strand=deepcopy(next_sheet.strands[1])
            next_registration=deepcopy(next_sheet.registrations[1])

            working_sheet.add_strand(next_strand)
            working_sheet.add_registration(next_registration)

            next_key=next_strand.as_atom_selections()
            if next_key in used_strand_selections: break
    # Now renumber the new sheets
    sheet_number=0
    for sheet in new_sheets:
      sheet_number+=1
      sheet.n_strands=len(sheet.strands)
      sheet.sheet_id="%d" %(sheet_number)
      strand_number=0
      for strand in sheet.strands:
        strand_number+=1
        strand.sheet_id=sheet.sheet_id
        strand.strand_id=strand_number

    # Now new_sheets has strands arranged in sheets.
    new_annotation=deepcopy(self)
    new_annotation.sheets=new_sheets
    return new_annotation

  def combine_annotations(self,hierarchy=None,other=None,
    keep_self=None,
    max_h_bond_length=None,
    maximize_h_bonds=True,
    maximum_length_difference=None,
    minimum_overlap=None,
    out=sys.stdout):
    "Create new annotation that is combination of self and other"

    # Hierarchy is required in order to identify what parts of annotation are
    #   in common for self and other and to identify which are better

    # Default behavior (keep_self=None):
    # For each annotation, keep the version that has the most information;
    #   if equal, keep self.  Then add things from other that are not present in
    #   self. Information is number of H-bonds (default) or total residues in
    #   secondary_structure.

    # Optional behavior (keep_self=True):
    # Keep all of self and add things from other not present in self.

    self.keep_self=keep_self
    self.maximize_h_bonds=maximize_h_bonds
    self.maximum_length_difference=maximum_length_difference
    self.minimum_overlap=minimum_overlap
    self.max_h_bond_length=max_h_bond_length

    # Select just the part of hierarchy we will need
    assert hierarchy
    s1=self.as_atom_selections()
    s2=other.as_atom_selections()
    atom_selection=self.combine_atom_selections([s1,s2])
    if not atom_selection:
      return None # nothing to do
    asc=hierarchy.atom_selection_cache()
    sel = asc.selection(string = atom_selection)
    ph=hierarchy.deep_copy().select(sel)  #ph = needed part of hierarchy
    # Split all sheets into pairs
    a1=self.split_sheets()
    a2=other.split_sheets()

    # Find all pairs of overlapping annotations and all unique annotations in
    #  self and other

    #print >>out,"\nFinding matching and unique helices:"
    helices=self.get_unique_set(
       a1.helices,a2.helices,hierarchy=hierarchy,out=out)

    #print >>out,"\nFinding matching and unique sheets:"
    sheets=self.get_unique_set(a1.sheets,a2.sheets,hierarchy=hierarchy,out=out)

    new_annotation=annotation(sheets=sheets,helices=helices)
    new_annotation=new_annotation.merge_sheets()
    return new_annotation

  def score_pair(self,h1,h2,
      maximize_h_bonds=None,
      hierarchy=None,
      max_h_bond_length=None,
      poor_h_bond_weight=0.5,
      keep_self=None):

      assert h1 # h2 can be None
      if maximize_h_bonds:
        n_good_1,n_poor_1=h1.count_h_bonds(hierarchy=hierarchy,
          max_h_bond_length=self.max_h_bond_length)
        score_1=n_good_1+poor_h_bond_weight*n_poor_1
        if h2:
          n_good_2,n_poor_2=h2.count_h_bonds(hierarchy=hierarchy,
            max_h_bond_length=self.max_h_bond_length)
          score_2=n_good_2+poor_h_bond_weight*n_poor_2
        else:
          score_2=None
      elif keep_self:
        score_1=1
        score_2=0
      else: # take the one with more residues in secondary structure
        score_1=h1.count_residues(hierarchy=hierarchy)
        if h2:
          score_2=h2.count_residues(hierarchy=hierarchy)
        else:
          score_2=None
      return score_1,score_2

  def remove_overlapping_annotations(self,hierarchy=None,
     maximize_h_bonds=True,max_h_bond_length=None,
     maximum_length_difference=None,minimum_overlap=None):
    if not hierarchy:
      raise Sorry("Need hierarchy for remove_overlapping_annotations")
    self.maximum_length_difference=maximum_length_difference
    self.minimum_overlap=minimum_overlap
    self.maximize_h_bonds=maximize_h_bonds
    self.max_h_bond_length=max_h_bond_length

    self.sheets=self.split_sheets().sheets
    new_sheets=self.select_best_overlapping_annotations(hierarchy=hierarchy,
      sheet_or_helix_list=self.sheets)
    new_helices=self.select_best_overlapping_annotations(hierarchy=hierarchy,
      sheet_or_helix_list=self.helices)
    from copy import deepcopy
    new_annotation=annotation(
      helices=deepcopy(new_helices),
      sheets=deepcopy(new_sheets))
    return new_annotation.merge_sheets()

  def select_best_overlapping_annotations(self,hierarchy=None,
      sheet_or_helix_list=None):
    assert hierarchy
    score_list=[]
    for s1 in sheet_or_helix_list:
      score_1,score_2=self.score_pair(s1,None,
         maximize_h_bonds=self.maximize_h_bonds,
         hierarchy=hierarchy,
         max_h_bond_length=self.max_h_bond_length,
         keep_self=None)
      score_list.append([score_1,s1])
    score_list.sort()
    score_list.reverse()

    remove_list=[]
    keep_list=[]
    for score1,s1 in score_list: # remove things that overlap and lower score
      if not s1 in remove_list and not s1 in keep_list:
        keep_list.append(s1)
      for score2,s2 in score_list:
        if s2==s1 or s2 in keep_list or s2 in remove_list: continue
        if s1.overlaps_with(other=s2,hierarchy=hierarchy):
          remove_list.append(s2)
    return keep_list
 
  def get_unique_set(self,h1_list,h2_list,hierarchy=None,
     out=sys.stdout):
    # Sheets should be split before using get_unique_set
    pairs=[]
    overlapping_but_not_matching_pairs=[]
    unique_h1=[]
    unique_h2=[]
    used_h1=[]
    used_h2=[]
    for h1 in h1_list:
      for h2 in h2_list:

        if h2.is_similar_to(other=h1,hierarchy=hierarchy,
            maximum_length_difference=self.maximum_length_difference,
            minimum_overlap=self.minimum_overlap):
          if not h1 in used_h1: used_h1.append(h1)
          if not h2 in used_h2: used_h2.append(h2)
          pairs.append([h1,h2])
        elif h2.overlaps_with(other=h1,hierarchy=hierarchy):
          overlapping_but_not_matching_pairs.append([h1,h2])
          if not h1 in used_h1: used_h1.append(h1)
          if not h2 in used_h2: used_h2.append(h2)
    for h1 in h1_list:
      if not h1 in used_h1: unique_h1.append(h1)
    for h2 in h2_list:
      if not h2 in used_h2: unique_h2.append(h2)

    if pairs:
      #print >>out,"\nMatching pairs:"
      for [h1,h2] in pairs:

        score_1,score_2=self.score_pair(h1,h2,
          maximize_h_bonds=self.maximize_h_bonds,
          hierarchy=hierarchy,
          max_h_bond_length=self.max_h_bond_length,
          keep_self=self.keep_self)
        #print >>out,"SELF : %7.1f\n%s" %(score_1,h1.as_pdb_str())
        #print >>out,"OTHER: %7.1f\n%s" %(score_2,h2.as_pdb_str())

    if overlapping_but_not_matching_pairs:
      #print >>out,"\nOverlapping non-matching pairs:"
      for [h1,h2] in overlapping_but_not_matching_pairs:
        score_1,score_2=self.score_pair(h1,h2,
          maximize_h_bonds=self.maximize_h_bonds,
          hierarchy=hierarchy,
          max_h_bond_length=self.max_h_bond_length,
          keep_self=self.keep_self)
        #print >>out,"SELF : %7.1f\n%s" %(score_1,h1.as_pdb_str())
        #print >>out,"OTHER: %7.1f\n%s\n" %(score_2,h2.as_pdb_str())

    if unique_h1 or unique_h2:
      #print >>out,"\nNon-matching:"
      for h1 in unique_h1:
        score_1,score_2=self.score_pair(h1,None,
          maximize_h_bonds=self.maximize_h_bonds,
          hierarchy=hierarchy,
          max_h_bond_length=self.max_h_bond_length,
          keep_self=self.keep_self)
        #print >>out,"SELF : %7.1f\n%s" %(score_1,h1.as_pdb_str())
      for h2 in unique_h2:
        score_1,score_2=self.score_pair(h2,None,
          maximize_h_bonds=self.maximize_h_bonds,
          hierarchy=hierarchy,
          max_h_bond_length=self.max_h_bond_length,
          keep_self=self.keep_self)
        #print >>out,"OTHER: %7.1f\n%s" %(score_1,h2.as_pdb_str())

    unique=unique_h1+unique_h2

    # Now take the best (or self if desired) of pairs and add on unique
    final_list=[]
    for [h1,h2] in pairs+overlapping_but_not_matching_pairs:
      score_1,score_2=self.score_pair(h1,h2,
        maximize_h_bonds=self.maximize_h_bonds,
        hierarchy=hierarchy,
        max_h_bond_length=self.max_h_bond_length,
        keep_self=self.keep_self)

      #print "score 1 and 2: ",score_1,score_2
      if score_1>=score_2:
        final_list.append(h1)
      else:
        final_list.append(h2)

    for h in unique:
      final_list.append(h)

    return final_list

  def is_comparable_to(self,other=None):

    # Just see if self and other are remotely comparable

    if other is None: return False,None,None

    if type(self) != type(other):
      return False,None,None # Check for different objects

    # Check for completely different helices
    if len(self.helices) != len(other.helices):
      return False,None,None

    # Split all sheets so that everything is comparable
    a1=self.split_sheets()
    a2=other.split_sheets()
    # Check for completely different sheets
    if len(a1.sheets) != len(a2.sheets):
      return False,None,None

    return True,a1,a2

  def overlaps_with(self,other=None,hierarchy=None):
    "Returns True if any element of the annotation overlap"

    a1=self.split_sheets()
    a2=other.split_sheets()
    assert hierarchy
    for h1 in a1.helices:
      for h2 in a2.helices:
        if h1.overlaps_with(other=h2,hierarchy=hierarchy):
          return True

    for s1 in a1.sheets:
      for s2 in a2.sheets:
        if s1.overlaps_with(other=s2,hierarchy=hierarchy):
          return True

    return False

  def is_similar_to(self,other=None,hierarchy=None,
      maximum_length_difference=None, minimum_overlap=None):
    "Return True if this annotation is similar to other."
    "Allow equivalent representations (different alignment atoms with)"
    "equivalent results, backwards order of strands in sheet, different "
    "lengths of strands optional."

    if self.is_same_as(other=other): return True  # always quicker

    is_comparable,a1,a2=self.is_comparable_to(other=other)
    if not is_comparable: return False

    # Now a1 and a2 both have split sheets

    # Check if helices are comparable
    for h1 in a1.helices:
      found=False
      for h2 in a2.helices:
        if h1.is_similar_to(other=h2,hierarchy=hierarchy,
            minimum_overlap=minimum_overlap,
            maximum_length_difference=maximum_length_difference):
          found=True
          break
      if not found: return False
    # Matched all the helices

    # Check if sheets are comparable. Note all sheets have 2 strands
    for s1 in a1.sheets:
      found=False
      for s2 in a2.sheets:
        if s1.is_similar_to(other=s2,hierarchy=hierarchy,
            minimum_overlap=minimum_overlap,
            maximum_length_difference=maximum_length_difference):
          found=True
          break
      if not found: return False
    # Matched all the sheets
    return True

  def is_same_as(self,other=None):
    "Return True if this annotation is the same as other. Allows different"
    "representations of same sheet (as pairs, as sheet)"
    "Does not allow different order of sheet or backwards"
    "Does not allow equivalent but different specifications of the alignment"

    is_comparable,a1,a2=self.is_comparable_to(other=other)
    if not is_comparable: return False

    def sort_strings(h):
      sorted=[]
      for x in h:
        sorted.append(x.as_pdb_str(set_id_zero=True))
      sorted.sort()
      return sorted

    # Sort helices and compare
    a1_helices=sort_strings(a1.helices)
    a2_helices=sort_strings(a2.helices)
    if a1_helices != a2_helices: return False

    # Sort sheets and compare
    a1_sheets=sort_strings(a1.sheets)
    a2_sheets=sort_strings(a2.sheets)
    if a1_sheets != a2_sheets: return False

    # Everything is the same
    return True

  def count_h_bonds(self,hierarchy=None,max_h_bond_length=None):
    # go through ss and count good and poor H-bonds
    number_of_good_h_bonds=0
    number_of_poor_h_bonds=0
    split_sheets=self.split_sheets() # annotation with split sheets
    for h in self.helices:
      n_good,n_poor=h.count_h_bonds(hierarchy=hierarchy,
        max_h_bond_length=max_h_bond_length)
      number_of_good_h_bonds+=n_good
      number_of_poor_h_bonds+=n_poor
    for s in self.sheets:
      n_good,n_poor=s.count_h_bonds(hierarchy=hierarchy,
        max_h_bond_length=max_h_bond_length)
      number_of_good_h_bonds+=n_good
      number_of_poor_h_bonds+=n_poor
    return number_of_good_h_bonds,number_of_poor_h_bonds


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
      self.adopt_serial_as_id()
    else:
      assert isinstance(self.helix_id, str)
    if self.start_chain_id != self.end_chain_id:
      raise RuntimeError("Don't know how to deal with helices with multiple "+
        "chain IDs ('%s' vs. '%s')." % (self.start_chain_id, self.end_chain_id))

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
      start_resseq=(line[21:25]),
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
    if amide_isel is None or len(amide_isel) == 0:
      error_msg = "Error in helix definition.\n"
      error_msg += "String '%s' selected 0 atoms.\n" % sele_str
      error_msg += "Most likely the definition of helix does not match model.\n"
      raise Sorry(error_msg)
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
      start_resseq=start_atom.resseq,
      start_icode=start_atom.icode,
      end_resname=end_atom.resname,
      end_chain_id=end_atom.chain_id,
      end_resseq=end_atom.resseq,
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

  def set_new_serial(self, serial, adopt_as_id=False):
    self.serial = serial
    if adopt_as_id:
      self.adopt_serial_as_id()

  def adopt_serial_as_id(self):
    self.helix_id = "%s" % self.serial
    self.helix_id = self.helix_id[:3]

  def set_new_chain_ids(self, new_chain_id):
    self.start_chain_id = new_chain_id
    self.end_chain_id = new_chain_id

  def erase_hbond_list(self):
    self.hbond_list = []

  def as_restraint_group(self, log=sys.stdout, prefix_scope="",
      add_segid=None, show_hbonds=False):
    if self.start_chain_id != self.end_chain_id :
      print >> log, "Helix chain ID mismatch: starts in %s, ends in %s" % (
        self.start_chain_id, self.end_chain_id)
      return None
    sele = self.as_atom_selections(add_segid=add_segid)[0]
    if prefix_scope != "" and not prefix_scope.endswith(".") :
      prefix_scope += "."
    serial_and_id = ""
    if self.serial is not None and self.serial > 0:
      serial_and_id += "\n  serial_number = %d" % self.serial
    if self.helix_id is not None:
      serial_and_id += "\n  helix_identifier = %s" % self.helix_id
    hbond_restr = ""
    if show_hbonds:
      if self.get_n_defined_hbonds() > 0:
        for hb in self.hbond_list:
          hb_str = "\n  hbond {\n    donor = %s\n    acceptor = %s\n  }" % (
            hb[0], hb[1])
          hbond_restr += hb_str
    rg = """\
%sprotein.helix {%s
  selection = %s
  helix_type = %s%s
}""" % (prefix_scope, serial_and_id, sele,
        self.helix_class, hbond_restr)
    return rg

  def as_pdb_str (self, set_id_zero=False):
    def h_class_to_pdb_int(h_class):
      h_class_int = self.helix_class_to_int(h_class)
      if h_class_int == 0:
        return 1
      return h_class_int
    format = "HELIX  %3d %3s %3s%2s %4s%1s %3s%2s %4s%1s%2d%30s %5d"
    if set_id_zero:
      serial=0
      helix_id=0
    else:
      serial=self.serial
      if self.helix_id is None:
        helix_id=self.serial
      else:
        helix_id=self.helix_id[:3]
    out = format % (serial,
      helix_id,
      self.start_resname,
      self.start_chain_id, self.start_resseq, self.start_icode,
      self.end_resname, self.end_chain_id, self.end_resseq, self.end_icode,
      h_class_to_pdb_int(self.helix_class), self.comment, self.length)
    return out.strip()

  def as_atom_selections(self, add_segid=None):
    segid_extra = ""
    if add_segid is not None :
      segid_extra = "and segid '%s' " % add_segid
    resid_start = "%s%s" % (self.start_resseq, self.start_icode)
    resid_end = "%s%s" % (self.end_resseq, self.end_icode)
    sele = "chain '%s' %sand resid %s through %s" % (self.start_chain_id,
      segid_extra, resid_start, resid_end)
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
      assert 0, "Wrong helix_class creeped in object fields: %d" % (
          self.helix_class)

  def is_same_as(self,other=None):
    # return True if it self has the same values as other

    if not other: return False
    if type(self) != type(other): return False # Check for different objects

    for key in [ 'start_resname',
        'start_chain_id',
        'start_resseq',
        'start_icode',
        'end_resname',
        'end_chain_id',
        'end_resseq',
        'end_icode',
        'helix_class' ]:
      if getattr(self,key,None) != getattr(other,key,None):
        return False
    return True

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
    if self.start_chain_id != self.end_chain_id:
      raise RuntimeError("Don't know how to deal with helices with multiple "+
        "chain IDs ('%s' vs. '%s')." % (self.start_chain_id, self.end_chain_id))

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

  def set_new_chain_ids(self, new_chain_id):
    self.start_chain_id = new_chain_id
    self.end_chain_id = new_chain_id

  def as_atom_selections(self, add_segid=None):
    segid_extra = ""
    if add_segid is not None :
      segid_extra = "and segid '%s' " % add_segid
    resid_start = "%s%s" % (self.start_resseq, self.start_icode)
    resid_end = "%s%s" % (self.end_resseq, self.end_icode)
    sele = "chain '%s' %sand resid %s through %s" % (self.start_chain_id,
      segid_extra, resid_start, resid_end)
    return sele


  def is_same_as(self,other=None):
    # return True if it self has the same values as other
    # include sense (as it is relative to another strand)

    if not other: return False
    if type(self) != type(other): return False # Check for different objects

    key_list= [ 'start_resname',
        'start_chain_id',
        'start_resseq',
        'start_icode',
        'end_resname',
        'end_chain_id',
        'end_resseq',
        'end_icode',
        'sense']

    for key in key_list:
      if getattr(self,key,None) != getattr(other,key,None):
        return False
    return True

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
    if len(line.strip()) < 67:
      return None
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

  def set_new_cur_chain_id(self, new_chain_id):
    self.cur_chain_id = new_chain_id

  def set_new_prev_chain_id(self, new_chain_id):
    self.prev_chain_id = new_chain_id

  def as_atom_selections(self, add_segid=None):
    segid_extra = ""
    if add_segid is not None :
      segid_extra = "and segid '%s' " % add_segid
    sele_base = "chain '%s' %sand resid %s and name %s"
    resid_curr = "%d%s" % (self.cur_resseq,self.cur_icode)
    resid_prev = "%d%s" % (self.prev_resseq,
        self.prev_icode)
    sele_curr = sele_base % (self.cur_chain_id, segid_extra,
        resid_curr, self.cur_atom.strip())
    sele_prev = sele_base % (self.prev_chain_id,segid_extra,
        resid_prev, self.prev_atom.strip())
    return sele_curr, sele_prev

  def is_same_as(self,other=None):

    if not other: return False
    if type(self) != type(other): return False # Check for different objects

    for key in [ 'cur_atom',
        'cur_resname',
        'cur_chain_id',
        'cur_resseq',
        'cur_icode',
        'prev_atom',
        'prev_resname',
        'prev_chain_id',
        'prev_resseq',
        'prev_icode']:
      if getattr(self,key,None) != getattr(other,key,None):
        return False
    return True


  def reversed(self):
    # swap cur and prev
    from copy import deepcopy
    new_register=deepcopy(self)
    for key in [
        'atom',
        'resname',
        'chain_id',
        'resseq',
        'icode']:
      setattr(new_register,'prev_'+key,getattr(self,'cur_'+key))
      setattr(new_register,'cur_'+key,getattr(self,'prev_'+key))
    return new_register



class pdb_sheet(structure_base):
  def __init__ (self,
      sheet_id,
      n_strands,
      strands,
      registrations,
      hbond_list=[], # list of (donor, acceptor) selections
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
      reg = pdb_strand_register.from_pdb_record(r)
      # do we really need to stop on this?
      if sense == 0 and reg is not None:
        raise Sorry("Sense should be 1 or -1 for non-first strand:\n%s" % r)
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
    if amide_isel is None or len(amide_isel) == 0:
      error_msg = "Error in sheet definition.\n"
      error_msg += "String '%s' selected 0 atoms.\n" % sele_str
      error_msg += "Most likely the definition of sheet does not match model.\n"
      raise Sorry(error_msg)
    start_atom = atoms[amide_isel[0]]
    end_atom = atoms[amide_isel[-1]]
    first_strand = pdb_strand(
        sheet_id=sheet_id,
        strand_id=1,
        start_resname=start_atom.resname,
        start_chain_id=start_atom.chain_id,
        start_resseq=start_atom.resseq,
        start_icode=start_atom.icode,
        end_resname=end_atom.resname,
        end_chain_id=end_atom.chain_id,
        end_resseq=end_atom.resseq,
        end_icode=end_atom.icode,
        sense=0)
    strands = [first_strand]
    registrations = [None]
    for i, strand_param in enumerate(sheet_params.strand):
      sele_str = ("(%s) and (name N) and (altloc 'A' or altloc ' ')" %
                      strand_param.selection)
      amide_isel = isel(sele_str)
      if amide_isel is None or len(amide_isel) == 0:
        error_msg = "Error in helix definition.\n"
        error_msg += "String '%s' selected 0 atoms.\n" % sele_str
        error_msg += "Most likely the definition of helix does not match model.\n"
        raise Sorry(error_msg)
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
      reg_cur_atom = None
      if strand_param.bond_start_current is not None:
        reg_cur_sel = isel(strand_param.bond_start_current)
        if reg_cur_sel is None or len(reg_cur_sel) == 0:
          error_msg = "Error in sheet definition.\n"
          error_msg += "String '%s' selected 0 atoms.\n" % strand_param.bond_start_current
          error_msg += "Most likely the definition of sheet does not match model.\n"
          raise Sorry(error_msg)
        reg_cur_atom = atoms[reg_cur_sel[0]]
      if reg_cur_atom is None: # No current atom in registration
        pass
        # raise Sorry("This bond_start_current yields 0 atoms:\n %s" % strand_param.bond_start_current)
      reg_prev_atom = None
      if strand_param.bond_start_previous is not None:
        reg_prev_sel = isel(strand_param.bond_start_previous)
        if reg_prev_sel is None or len(reg_prev_sel) == 0:
          error_msg = "Error in sheet definition.\n"
          error_msg += "String '%s' selected 0 atoms.\n" % strand_param.bond_start_previous
          error_msg += "Most likely the definition of sheet does not match model.\n"
          raise Sorry(error_msg)
        reg_prev_atom = atoms[reg_prev_sel[0]]
      if reg_prev_atom is None: # No previous atom in registration
        pass
        # raise Sorry("This bond_start_previous yields 0 atoms:\n %s" % strand_param.bond_start_previous)
      reg = None
      if reg_cur_atom is not None and reg_prev_atom is not None:
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

  def erase_hbond_list(self):
    self.hbond_list = []

  def add_strand (self, strand) :
    self.strands.append(strand)

  def add_registration (self, registration) :
    self.registrations.append(registration)

  def as_atom_selections (self, add_segid=None):
    strand_selections = []
    for strand in self.strands:
      strand_selections.append(strand.as_atom_selections(add_segid=add_segid))
    return strand_selections

  def get_n_defined_hbonds(self):
    if self.hbond_list is not None:
      return len(self.hbond_list)
    return 0

  def as_pdb_str(self, strand_id=None, set_id_zero=False):
    assert len(self.strands) == len(self.registrations)
    lines = []
    for strand, reg in zip(self.strands, self.registrations) :
      format1 = "SHEET  %3d %3s%2d %3s%2s%4s%1s %3s%2s%4s%1s%2d"
      format2 = "%4s%3s%2s%4s%1s %4s%3s%2s%4s%1s"
      # print "STRAND, REG", strand, reg
      if set_id_zero:
        strand_id_for_print=0
        sheet_id_for_print=0
      else:
        strand_id_for_print=strand.strand_id
        sheet_id_for_print=strand.sheet_id

      line = format1 % (strand_id_for_print, sheet_id_for_print, self.n_strands,
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
      if strand_id is not None and strand_id == strand.strand_id:
        return line
      lines.append(line.strip())
    return "\n".join(lines)

  def as_restraint_group(self, log=sys.stdout, prefix_scope="",
      add_segid=None, show_hbonds=False):
    if len(self.strands) == 0 :
      return None
    selections = []
    senses = []
    reg_curr = []
    reg_prev = []
    for (strand,registration) in zip(self.strands, self.registrations) :
      sele = strand.as_atom_selections(add_segid=add_segid)
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
        curr, prev = registration.as_atom_selections(add_segid=add_segid)
        reg_curr.append(curr)
        reg_prev.append(prev)
      else :
        reg_curr.append(None)
        reg_prev.append(None)
    # print "selections", selections
    # print "reg_curr", reg_curr
    # print "reg_prev", reg_prev
    # STOP()
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
    hbond_restr = ""
    if show_hbonds:
      if self.get_n_defined_hbonds() > 0:
        for hb in self.hbond_list:
          hb_str = "\n  hbond {\n    donor = %s\n    acceptor = %s\n  }" % (
            hb[0], hb[1])
          hbond_restr += hb_str
    phil_str = """
%sprotein.sheet {
  sheet_id = "%s"
  first_strand = %s
%s%s
}""" % (prefix_scope, self.sheet_id, first_strand, "\n".join(strands), hbond_restr)
    # print "phil_str", phil_str
    return phil_str

  def split(self,starting_sheet_id_number=1):
    assert len(self.strands)==len(self.registrations)
    new_sheets=[]
    from copy import deepcopy
    for s1,s2,r2 in zip(
      self.strands[:-1],
      self.strands[1:],self.registrations[1:]):
      self_id="%d" %(len(new_sheets)+1)
      new_s1=deepcopy(s1)
      new_s1.sense=0
      new_s1.strand_id=1
      new_s1.sheet_id="%d" %(starting_sheet_id_number)
      new_s2=deepcopy(s2)
      new_s2.strand_id=2
      new_s2.sheet_id="%d" %(starting_sheet_id_number)
      new_r2=deepcopy(r2)
      new_sheet=pdb_sheet(
             sheet_id="%d" %(starting_sheet_id_number),
             n_strands=2,
             strands=[new_s1,new_s2],
             registrations=[None,new_r2])
      new_sheets.append(new_sheet)
      starting_sheet_id_number+=1
    return new_sheets

  def overlaps_with(self,other=None,hierarchy=None):
    "For sheet, overlaps if both members of any strand pair overlap"
    assert hierarchy

    if len(self.strands)<2: return False
    if len(other.strands)<2: return False
    self_sheets=[self]
    other_sheets=[other]
    if len(self.strands)>2 or len(other.strands)>2: # break down
      self_sheets=self.split()
      other_sheets=other.split()
    asc=hierarchy.atom_selection_cache()
    for self_sheet in self_sheets:
      for other_sheet in other_sheets: 
        assert len(self_sheet.strands)==2 
        assert len(other_sheet.strands)==2 
        
        s1a=self_sheet.strands[0].as_atom_selections()
        s1b=self_sheet.strands[1].as_atom_selections()
        s2a=other_sheet.strands[0].as_atom_selections()
        s2b=other_sheet.strands[1].as_atom_selections()
        for pair_1,pair_2 in [[[s1a,s2a],[s1b,s2b]],[[s1a,s2b],[s1b,s2a]] ]:
          atom_selection=self.combine_atom_selections(pair_1,require_all=True)
          if not atom_selection: continue
          sel = asc.selection(string = atom_selection)
          ph=hierarchy.deep_copy().select(sel)  #ph = needed part of hierarchy
          if ph.overall_counts().n_residues==0: continue

          atom_selection=self.combine_atom_selections(pair_2,require_all=True)
          if not atom_selection: continue
          sel = asc.selection(string = atom_selection)
          ph=hierarchy.deep_copy().select(sel)  #ph = needed part of hierarchy
          if ph.overall_counts().n_residues>0:
            return True  # both match
    return False

  def is_similar_to(self,other=None,hierarchy=None,
      maximum_length_difference=None, minimum_overlap=None,log=sys.stdout):
    # Two sheets are similar if their strands are similar and their
    # alignments are similar, or if reversing the order, they are similar
    if self.is_same_as(other=other): return True  # always quicker

    # Check if two strands are similar to the other two strands
    s1a=self.strands[0]
    s1b=self.strands[1]
    s2a=other.strands[0]
    s2b=other.strands[1]

    # both must be parallel or antiparallel
    if s1b.sense!=s2b.sense: return False

    match_forward=None
    if s1a.is_similar_to(other=s2a,hierarchy=hierarchy,
        maximum_length_difference=maximum_length_difference,
        minimum_overlap=minimum_overlap) and \
       s1b.is_similar_to(other=s2b,hierarchy=hierarchy,
        maximum_length_difference=maximum_length_difference,
        minimum_overlap=minimum_overlap):
      match_forward=True
    elif s1a.is_similar_to(other=s2b,hierarchy=hierarchy,
        maximum_length_difference=maximum_length_difference,
        minimum_overlap=minimum_overlap) and \
       s1b.is_similar_to(other=s2a,hierarchy=hierarchy,
        maximum_length_difference=maximum_length_difference,
        minimum_overlap=minimum_overlap):
      match_forward=False
    if match_forward is None:
      return False

    assert len(self.registrations)==len(other.registrations)
    reg1=self.registrations[1]
    reg2=other.registrations[1]

    if reg1 is None and reg2 is None:
      return True # seems ok
    elif reg1 is None or reg2 is None:
      return False # one has it and other does not

    # otherwise check the two for similarity 
    if not match_forward:  # swap other so they should be comparable
      xx=s2a
      s2a=s2b
      s2b=xx
      reg2=reg2.reversed()

    # Now check alignments


    # Find the registration atom positions
    ra=registration_atoms(
      hierarchy=hierarchy,
      strand_1a=s1a,strand_1b=s1b,strand_2a=s2a,strand_2b=s2b,
      registration_1=reg1,registration_2=reg2,
      sense=s1b.sense,log=log)

    return ra.pairs_are_equivalent

  def is_same_as(self,other=None):

    if not other: return False
    if type(self) != type(other): return False # Check for different objects

    if self.n_strands != other.n_strands: return False
    for s1,s2 in zip(self.strands,other.strands):
      if not s1.is_same_as(s2): return False
    for r1,r2 in zip(self.registrations,other.registrations):
      if r1 is None and r2 is None: continue # special case..can be None
      if r1 is None or r2 is None: return False
      if not r1.is_same_as(r2): return False
    return True

if __name__ == "__main__" :
  pass
