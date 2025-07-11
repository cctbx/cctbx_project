"""Tools for identification of secondary structure in a model"""
from __future__ import absolute_import, division, print_function
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

from libtbx import adopt_init_args
import sys
from iotbx.pdb.hybrid_36 import hy36encode, hy36decode
import iotbx.cif.model
import copy
from libtbx.utils import null_out
from six.moves import range
from six.moves import zip


def lists_have_comment_element(a,b):
  for x in a:
    if x in b:
      return True
  return False

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
      h1=hierarchy.select(sel, copy_atoms=True)  # detach from original hierarchy
      number_self=h1.overall_counts().n_residues
    except Exception as e:
      return False

    asc=hierarchy.atom_selection_cache()
    sel=asc.selection(string = atom_selection_2)
    try:
      h2=hierarchy.select(sel, copy_atoms=True)
      number_other=h2.overall_counts().n_residues
    except Exception as e:
      return False

    if maximum_length_difference is not None and \
        abs(number_self-number_other)>maximum_length_difference:
      return False

    asc=hierarchy.atom_selection_cache()
    atom_selection="(%s) and (%s)" %(atom_selection_1,atom_selection_2)
    sel=asc.selection(string = atom_selection)
    try:
      h12=hierarchy.select(sel, copy_atoms=True)
      number_both=h12.overall_counts().n_residues
    except Exception as e:
      return False

    if number_both<minimum_overlap:
      return False

    return True

def choose_correct_cif_record(cif_dict, key1, key2, mandatory=True):
  if key1 in cif_dict:
    val = cif_dict[key1]
    if val != '?' and val != '.':
      return val
  if mandatory:
    return cif_dict[key2]
  else:
    return cif_dict.get(key2, '.')

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
          print("Cannot interpret H-bonding of %s to %s " %(
            strand_a_atom,strand_b_atom), file=log)
          ok=False
          return

    if not sense in [-1,1]:
      print("Cannot interpret bonding with sense not -1 or 1:", file=log)
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
    ph=self.hierarchy.select(sel, copy_atoms=True)  # detach from original hierarchy
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
               residue.resseq.replace(" ","")==resseq.replace(" ", "") and \
               chain.id==chain_id and residue.icode==icode:
              return pos
    return None


def get_amide_isel(asc, ss_element_string_selection):
  # There are cases where SS element can be only in B conformation
  # and split in A/B conformations
  sele_str_main_conf = "(%s) and (name N) and (altloc 'A' or altloc ' ')" % ss_element_string_selection
  amide_isel_main_conf = asc.iselection(sele_str_main_conf)
  if len(amide_isel_main_conf) > 0:
    return amide_isel_main_conf
  sele_str_all = "(%s) and (name N)" % ss_element_string_selection
  amide_isel_all = asc.iselection(sele_str_all)
  if len(amide_isel_all) > 0:
    return amide_isel_all
  error_msg = "Error in helix definition. It will be skipped.\n"
  error_msg += "String '%s' selected 0 atoms.\n" % sele_str_main_conf
  error_msg += "String '%s' selected 0 atoms.\n" % sele_str_all
  error_msg += "Most likely the definition of SS element does not match model.\n"
  return error_msg


class structure_base(object):

  def as_pdb_str(self):
    return None

  def as_pdb_or_mmcif_str(self):
    return None

  def __str__(self):
    return self.as_pdb_or_mmcif_str()

  @staticmethod
  def convert_resseq(resseq):
    if isinstance(resseq, str):
      return "%4s" % resseq[:4]
    elif isinstance(resseq, int):
      return hy36encode(4, resseq)

  @staticmethod
  def convert_id(id):
    if isinstance(id, str):
      return "%3s" % id[:3]
    elif isinstance(id, int):
      return hy36encode(3, id)

  def get_start_resseq_as_int(self):
    if self.start_resseq is not None:
      if isinstance(self.start_resseq, str):
        return hy36decode(4, self.start_resseq)
      elif isinstance(self.start_resseq, int):
        return self.start_resseq
    return None

  def get_end_resseq_as_int(self):
    if self.end_resseq is not None:
      if isinstance(self.end_resseq, str):
        return hy36decode(4, self.end_resseq)
      elif isinstance(self.end_resseq, int):
        return self.end_resseq
    return None

  @staticmethod
  def id_as_int(id):
    if isinstance(id, str):
      return hy36decode(min(len(id), 3), id)
    elif isinstance(id, int):
      return id
    raise ValueError("String or int is needed")

  @staticmethod
  def icode_to_cif(icode):
    if icode == ' ':
      return '?'
    else:
      return icode

  @staticmethod
  def parse_chain_id(chars):
    # assert len(chars) == 2
    if chars == "":
      # most likely somebody stripped empty chain id
      return " "
    if len(chars) == 1:
      return chars
    if chars == "  ":
      return " "
    else:
      return chars.strip()

  @staticmethod
  def parse_cif_insertion_code(chars):
    if chars == '?' or chars == '.' or chars == '':
      return ' '
    return chars

  @staticmethod
  def filter_helix_records(lines):
    result = []
    if lines is None:
      return result
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
    if lines is None:
      return result
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
      raise AssertionError("Require hierarchy for count_residues")

    atom_selection=self.combine_atom_selections(self.as_atom_selections())
    if not atom_selection:
       return 0

    asc=hierarchy.atom_selection_cache()
    sel = asc.selection(string = atom_selection)
    ph=hierarchy.select(sel, copy_atoms=True)
    return ph.overall_counts().n_residues

  def present_in_hierarchy(self, hierarchy):
    start_present = False
    end_present = False
    start_resseq = self.start_resseq if isinstance(self.start_resseq, str) else self.convert_resseq(self.start_resseq)
    end_resseq = self.end_resseq if isinstance(self.end_resseq, str) else self.convert_resseq(self.end_resseq)
    if len(hierarchy.models()) == 0:
      return False
    for chain in hierarchy.models()[0].chains():
      if chain.id == self.start_chain_id:
        for rg in chain.residue_groups():
          resname = rg.atom_groups()[0].resname
          if rg.resseq == start_resseq and self.start_resname == resname:
            start_present = True
          if rg.resseq == end_resseq and self.end_resname == resname:
            end_present = True
          if start_present and end_present:
            return True
    return False

  def _change_residue_numbering_in_place_helper(self, renumbering_dictionary, se=["start", "end"]):
    """ Changes residue numbers and insertion codes in annotations.
    For cases where one needs to renumber hierarchy and keep annotations
    consisten with the new numbering.
    Not for cases where residues removed or added.
    Therefore residue name stays the same.

    Called from derived classes.

    Args:
        renumbering_dictionary (_type_): data structure to match old and new numbering:
        {'chain id':
          {(old resseq, icode):(new resseq, icode)}
        }
        se: labels for start and end. Also can be ["cur", "prev"] in sheet registrations
    """
    for start_end in se:
      chain_id_attr = "%s_chain_id" % start_end
      resseq_attr  = "%s_resseq" % start_end
      icode_attr = "%s_icode" % start_end
      chain_dic = renumbering_dictionary.get(getattr(self, chain_id_attr), None)
      if chain_dic is not None:
          resseq, icode = getattr(self, resseq_attr), getattr(self, icode_attr)
          new_resseq, new_icode = chain_dic.get((resseq, icode), (None, None))
          if new_resseq is not None and new_icode is not None:
              setattr(self, resseq_attr, new_resseq)
              setattr(self, icode_attr, new_icode)


  def count_h_bonds(self,hierarchy=None,
       max_h_bond_length=None,ss_by_chain=False):
    "Count good and poor H-bonds in this hierarchy"
    "This is generic version"
    "This uses the annotation in self"
    "Not appropriate for large structures unless you set ss_by_chain=True"
    "Use instead annotation.count_h_bonds in general."
    if hierarchy is None:
      raise AssertionError("Require hierarchy for count_h_bonds")

    atom_selection=self.combine_atom_selections(self.as_atom_selections())
    if not atom_selection:
       return 0,0

    asc=hierarchy.atom_selection_cache()
    sel = asc.selection(string = atom_selection)
    ph=hierarchy.select(sel, copy_atoms=True)

    from mmtbx.secondary_structure.find_ss_from_ca import \
      find_secondary_structure
    fss=find_secondary_structure(hierarchy=ph,
      user_annotation_text=self.as_pdb_or_mmcif_str(),
      force_secondary_structure_input=True,
      combine_annotations=False,
      ss_by_chain=ss_by_chain,
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
      raise AssertionError("Require hierarchy for is_similar_to")

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

  def has_chains_in_common(self,other=None):
    if lists_have_comment_element(
       self.chains_included(),
       other.chains_included()):
      return True
    else:
      return False

  def overlaps_with(self,other=None,hierarchy=None):
    assert hierarchy
    s1=self.as_atom_selections()
    s2=other.as_atom_selections()

    atom_selection=self.combine_atom_selections([s1,s2],require_all=True)
    if atom_selection:
      asc=hierarchy.atom_selection_cache()
      sel = asc.selection(string = atom_selection)
      ph=hierarchy.select(sel, copy_atoms=True)  #ph = needed part of hierarchy
      if ph.overall_counts().n_residues>0:
        return True
    return False



class annotation(structure_base):
  def __init__(self, helices=None, sheets=None):
    assert (not None in [helices, sheets])
    self.helices = helices
    self.sheets = sheets

  @classmethod
  def resseq_as_int(cls, resseq):
    if resseq is not None:
      if isinstance(resseq, str):
        return hy36decode(4, resseq)
      elif isinstance(resseq, int):
        return resseq
    return None

  @classmethod
  def from_records(cls, records=None, log=None):
    "Initialize annotation from pdb HELIX/SHEET records"
    helices = []
    sheets = []
    for line in cls.filter_helix_records(records):
      try:
        h = pdb_helix.from_pdb_record(line)
      except ValueError:
        print("Bad HELIX record, was skipped:\n%s" % line, file=log)
      else:
        helices.append(h)
    for sh_lines in cls.filter_and_split_sheet_records(records):
      try:
        sh = pdb_sheet.from_pdb_records(sh_lines)
      except ValueError as e:
        print("Bad SHEET records, were skipped:\n", file=log)
        for l in sh_lines:
          print("  %s" % l, file=log)
      else:
        sheets.append(sh)
    return cls(helices=helices, sheets=sheets)

  @classmethod
  def from_cif_block(cls, cif_block, log=None):
    "initialize annotation from cif HELIX/SHEET records"
    helices = []
    serial = 1
    helix_loop = cif_block.get_loop_or_row("_struct_conf")
    if helix_loop is not None:
      for helix_row in helix_loop.iterrows():
        try:
          h = pdb_helix.from_cif_row(helix_row, serial)
          if h is not None:
            helices.append(h)
            serial += 1
        except ValueError:
          print("Bad HELIX records!", file=log)
    sheets = []
    struct_sheet_loop = cif_block.get_loop_or_row("_struct_sheet")
    struct_sheet_order_loop = cif_block.get_loop_or_row("_struct_sheet_order")
    struct_sheet_range_loop = cif_block.get_loop_or_row("_struct_sheet_range")
    struct_sheet_hbond_loop = cif_block.get_loop_or_row("_pdbx_struct_sheet_hbond")
    if (struct_sheet_loop is not None
        and struct_sheet_order_loop is not None
        and struct_sheet_range_loop is not None
        and struct_sheet_hbond_loop is not None):
      # Check that keys are unique
      error_msg_list = []
      dups_sheet = struct_sheet_loop.check_key_is_unique(key_list=["_struct_sheet.id"])
      dups_order = struct_sheet_order_loop.check_key_is_unique(key_list=["_struct_sheet_order.sheet_id",
          "_struct_sheet_order.range_id_1", "_struct_sheet_order.range_id_2"])
      dups_range = struct_sheet_range_loop.check_key_is_unique(key_list=["_struct_sheet_range.sheet_id",
            "_struct_sheet_range.id"])
      dups_hbond = struct_sheet_hbond_loop.check_key_is_unique(key_list=["_pdbx_struct_sheet_hbond.sheet_id",
          "_pdbx_struct_sheet_hbond.range_id_1", "_pdbx_struct_sheet_hbond.range_id_2"])
      for e_str, dups in [("Duplication in _struct_sheet.id: %s", dups_sheet),
          ("Duplication in _struct_sheet_order: %s %s %s", dups_order),
          ("Duplication in _struct_sheet_range: %s %s", dups_range),
          ("Duplication in _pdbx_struct_sheet_hbond: %s %s %s", dups_hbond)]:
        for dup in dups:
          error_msg_list.append(e_str % dup)
      if len(error_msg_list) > 0:
        msg = "Error in sheet definitions:\n  " + "\n  ".join(error_msg_list)
        raise Sorry(msg)
      for sheet_row in struct_sheet_loop.iterrows():
        sheet_id = sheet_row['_struct_sheet.id']
        # we will count number_of_strands in from_cif_rows
        # number_of_strands = int(sheet_row['_struct_sheet.number_strands'])
        try:
          sh = pdb_sheet.from_cif_rows(sheet_id,
              struct_sheet_order_loop, struct_sheet_range_loop,
              struct_sheet_hbond_loop)
        except ValueError as e:
          print("Bad sheet records.\n", file=log)
        else:
          sheets.append(sh)
    return cls(helices=helices, sheets=sheets)

  @classmethod
  def from_phil(cls, phil_helices, phil_sheets, pdb_hierarchy,
      atom_selection_cache=None, log=None):
    helices = []
    sheets = []
    cache = atom_selection_cache
    if cache is None:
      cache = pdb_hierarchy.atom_selection_cache()
    h_atoms = pdb_hierarchy.atoms()
    for i, helix_param in enumerate(phil_helices):
      if helix_param.selection is not None:
        h = pdb_helix.from_phil_params(helix_param, pdb_hierarchy, h_atoms, cache, i, log)
        if h is not None:
          helices.append(h)
    for i, sheet_param in enumerate(phil_sheets):
      if sheet_param.first_strand is not None:
        sh = pdb_sheet.from_phil_params(sheet_param, pdb_hierarchy, h_atoms, cache, log)
        if sh is not None:
          sheets.append(sh)
    return cls(helices=helices, sheets=sheets)

  def deep_copy(self):
    return copy.deepcopy(self)

  def simple_elements(self):
    """ helices and strands"""
    all_elem_list = []
    all_elem_list += self.helices
    for sh in self.sheets:
      all_elem_list += sh.strands
    return all_elem_list

  def filter_annotation(
      self,
      hierarchy=None,
      asc=None,
      remove_short_annotations=True,
      remove_3_10_helices=True,
      remove_empty_annotations=True,
      concatenate_consecutive_helices=True,
      split_helices_with_prolines=True,
      filter_sheets_with_long_hbonds=True
      ):
    """ Returns filtered annotation"""
    result = self.deep_copy()
    if remove_short_annotations:
      result.remove_short_annotations()
    if asc is None and hierarchy is not None:
      asc = hierarchy.atom_selection_cache()
    if remove_3_10_helices:
      result.remove_3_10_helices()
    if hierarchy is not None:
      if asc is None:
        asc = hierarchy.atom_selection_cache()
      if remove_empty_annotations:
        result.remove_empty_annotations(hierarchy)
      if concatenate_consecutive_helices:
        result.concatenate_consecutive_helices(hierarchy, asc)
      if split_helices_with_prolines:
        result.split_helices_with_prolines(hierarchy, asc)
      if filter_sheets_with_long_hbonds:
        result.filter_sheets_with_long_hbonds(hierarchy, asc)
    return result

  def remove_empty_annotations(self, hierarchy):
    # returns annotation of deleted helices and sheets
    h_indeces_to_delete = []
    sh_indeces_to_delete = []
    for i, h in enumerate(self.helices):
      if not h.present_in_hierarchy(hierarchy):
        h_indeces_to_delete.append(i)
    for i, sh in enumerate(self.sheets):
      for st in sh.strands:
        if not st.present_in_hierarchy(hierarchy):
          if i not in sh_indeces_to_delete:
            sh_indeces_to_delete.append(i)
    deleted_helices = []
    deleted_sheets = []
    if len(h_indeces_to_delete) > 0:
      for i in reversed(h_indeces_to_delete):
        deleted_helices.append(self.helices[i])
        del self.helices[i]
    if len(sh_indeces_to_delete) > 0:
      for i in reversed(sh_indeces_to_delete):
        deleted_sheets.append(self.sheets[i])
        del self.sheets[i]
    return annotation(helices=deleted_helices, sheets=deleted_sheets)

  def remove_short_annotations(self,
      helix_min_len=5, sheet_min_len=3, keep_one_stranded_sheets=False):
    # returns nothing
    # Remove short annotations
    h_indeces_to_delete = []
    for i, h in enumerate(self.helices):
      if h.length < helix_min_len:
        h_indeces_to_delete.append(i)
    if len(h_indeces_to_delete) > 0:
      for i in reversed(h_indeces_to_delete):
        del self.helices[i]
    sh_indeces_to_delete = []
    for i, sh in enumerate(self.sheets):
      sh.remove_short_strands(size=sheet_min_len)
      if sh.n_strands < 2 and not keep_one_stranded_sheets:
        sh_indeces_to_delete.append(i)
    if len(sh_indeces_to_delete) > 0:
      for i in reversed(sh_indeces_to_delete):
        del self.sheets[i]

  def concatenate_consecutive_helices(self, hierarchy=None, asc=None):
    # also should divide them
    from mmtbx.command_line.angle import calculate_axes_and_angle_directional
    if asc is None and hierarchy is not None:
      asc = hierarchy.atom_selection_cache()
    concatenations_made = True
    while concatenations_made:
      concatenations_made = False
      new_helices = []
      if self.get_n_helices() < 2:
        return
      new_helices.append(self.helices[0])
      for i in range(1, len(self.helices)):
        # checking angle
        angle = 0
        if hierarchy is not None:
          h1_hierarchy = hierarchy.select(asc.selection(self.helices[i-1].as_atom_selections()[0]))
          h2_hierarchy = hierarchy.select(asc.selection(self.helices[i].as_atom_selections()[0]))
          xrs1 = h1_hierarchy.extract_xray_structure()
          xrs2 = h2_hierarchy.extract_xray_structure()
          if xrs1.sites_cart().size() < 2 or xrs2.sites_cart().size() < 2:
            new_helices.append(self.helices[i])
            continue
          a1, a2, angle = calculate_axes_and_angle_directional(
              xrs1,
              xrs2)
        # print "angle between '%s', '%s': %.2f" % (
        #     self.helices[i-1].as_pdb_str()[7:40],
        #     self.helices[i].as_pdb_str()[7:40],
        #     angle)
        if (new_helices[-1].end_chain_id == self.helices[i].start_chain_id and
            abs(new_helices[-1].get_end_resseq_as_int() - self.helices[i].get_start_resseq_as_int()) < 2):
            # helices are next to each other, so we will either
          if (angle < 15 and
              new_helices[-1].end_resname != "PRO" and
              self.helices[i].start_resname != "PRO"):
            # concatenate
            new_helices[-1].end_resname = self.helices[i].start_resname
            new_helices[-1].set_end_resseq(self.helices[i].get_end_resseq_as_int())
            new_helices[-1].end_icode = self.helices[i].start_icode
            new_helices[-1].length = (new_helices[-1].get_end_resseq_as_int() -
              new_helices[-1].get_start_resseq_as_int() + 1)
            concatenations_made = True
          else:
            # or separate them by moving borders by 1 residue
            h1_rg = [x for x in h1_hierarchy.residue_groups()][-2]
            h2_rgs = [x for x in h2_hierarchy.residue_groups()]
            if len(h2_rgs) > 3:
              # short one could creep here
              h2_rg = [x for x in h2_hierarchy.residue_groups()][1]

              new_helices[-1].end_resname = h1_rg.atom_groups()[0].resname
              new_helices[-1].set_end_resseq(h1_rg.resseq)
              new_helices[-1].end_icode = h1_rg.icode
              new_helices[-1].length -= 1

              self.helices[i].start_resname = h2_rg.atom_groups()[0].resname
              self.helices[i].set_start_resseq(h2_rg.resseq)
              self.helices[i].start_icode = h2_rg.icode
              self.helices[i].length -= 1

              new_helices.append(self.helices[i])

        else:
          new_helices.append(self.helices[i])
      self.helices = new_helices

  def split_helices_with_prolines(self, hierarchy, asc=None):
    # If hierarchy is present: break helices with PRO inside
    if asc is None:
      asc = hierarchy.atom_selection_cache()
    new_helices = []
    for i,h in enumerate(self.helices):
      selected_h = hierarchy.select(asc.selection(h.as_atom_selections()[0]))
      i_pro_res = []
      rgs = selected_h.only_model().only_chain().residue_groups()
      for j,rg in enumerate(rgs):
        if rg.atom_groups()[0].resname.strip() == "PRO":
          i_pro_res.append(j)
      for j in i_pro_res:
        if j > 0:
          h1 = h.deep_copy()
          h1.end_resname = rgs[j-1].atom_groups()[0].resname.strip()
          h1.end_resseq = rgs[j-1].resseq
          h1.end_icode = rgs[j-1].icode
          h1.length = j
          h1.erase_hbond_list()
          h.start_resname = rgs[j].atom_groups()[0].resname.strip()
          h.start_resseq = rgs[j].resseq
          h.start_icode = rgs[j].icode
          h.length = len(rgs) - j
          h.erase_hbond_list()
          if h1.length > 2:
            new_helices.append(h1)
      if h.length > 3:
        new_helices.append(h)
    for i, h in enumerate(new_helices):
      h.set_new_serial(serial=i+1, adopt_as_id=True)
    self.helices=new_helices

  def filter_sheets_with_long_hbonds(self, hierarchy, asc=None):
    """ Currently using from_ca method becuase ksdssp works wrong with
    partial model. In some cases it cannot detect chain break and extends
    strand e.g. from 191 to 226 (in 1ubf) when the strand ends at 194th residue."""
    assert hierarchy is not None, "Cannot measure bonds w/o hierarchy!"
    from mmtbx.secondary_structure.ss_validation import gather_ss_stats
    from mmtbx.secondary_structure import manager as ss_manager
    from mmtbx.secondary_structure import sec_str_master_phil
    if asc is None:
      asc = hierarchy.atom_selection_cache()
    ss_stats_obj = gather_ss_stats(pdb_h=hierarchy)
    new_sheets = []
    sh_indeces_to_delete = []
    asc = asc
    if asc is None:
      asc = hierarchy.atom_selection_cache()
    # prepare atom_selections for sheets
    sh_atom_selections = []
    for sh in self.sheets:
      atom_selections = sh.as_atom_selections()
      total_selection = " or ".join(["(%s)" % x for x in atom_selections])
      sel = asc.selection(total_selection)
      sh_atom_selections.append(sel)
    for i, sh in enumerate(self.sheets):
      if i not in sh_indeces_to_delete:
        sh_tuple = ([], [sh])
        (n_hbonds, n_bad_hbonds, n_mediocre_hbonds, hb_lens,
            n_outliers, n_wrong_region) = ss_stats_obj(sh_tuple)
        if n_bad_hbonds > 0:
          # SHEET is bad, work with it
          sh_indeces_to_delete.append(i)
          # check other sheets for intersection with this one
          for j in range(i+1, self.get_n_sheets()):
            intersect = sh_atom_selections[i].deep_copy()
            if j not in sh_indeces_to_delete: # don't check deleted sheets
              intersect &= sh_atom_selections[j]
              if not intersect.all_eq(False):
                sh_indeces_to_delete.append(j)
                sh_atom_selections[i] |= sh_atom_selections[j]
          # find new sheet structure for this sheet
          sh_hierarchy = hierarchy.select(sh_atom_selections[i])
          # sh_hierarchy.write_pdb_file(file_name="sheet_%d.pdb"%i)
          n_rgs = len(list(sh_hierarchy.residue_groups()))
          ss_def_pars = sec_str_master_phil.extract()
          ss_def_pars.secondary_structure.protein.search_method="from_ca"
          ss_m = ss_manager(
              pdb_hierarchy=sh_hierarchy,
              params=ss_def_pars.secondary_structure,
              log=null_out())
          fresh_sheets = ss_m.actual_sec_str.sheets
          # checking for bug occuring in 4a7h where one strand happens to be
          # too long
          # ksdssp_bug = False
          # for f_sh in fresh_sheets:
          #   for f_strand in f_sh.strands:
          #     if f_strand.get_end_resseq_as_int() - f_strand.get_start_resseq_as_int() > n_rgs:
          #       sh_hierarchy.write_pdb_file(file_name="ksdssp_failure.pdb")
          #       ksdssp_bug = True
          #       break
          #       # raise Sorry("It is 4a7h or ksdssp failed on another structure.")
          # if ksdssp_bug:
          #   ss_def_pars = sec_str_master_phil.extract()
          #   ss_def_pars.secondary_structure.protein.search_method="from_ca"
          #   ss_m = ss_manager(
          #       pdb_hierarchy=sh_hierarchy,
          #       params=ss_def_pars.secondary_structure,
          #       log=null_out())
          new_sheets += ss_m.actual_sec_str.sheets
    if len(sh_indeces_to_delete) > 0:
      for i in sorted(sh_indeces_to_delete, reverse=True):
        del self.sheets[i]
    self.sheets = self.sheets + new_sheets
    if len(new_sheets) > 0:
      self.reset_sheet_ids()
    self.remove_short_annotations()

  def reset_sheet_ids(self):
    import itertools, string
    ch = string.digits+string.ascii_uppercase
    ids = [''.join(x) for x in itertools.product(ch, repeat = 1)]
    if self.get_n_sheets() > len(ids):
      ids = [''.join(x) for x in itertools.product(ch, repeat = 2)]
    if self.get_n_sheets() > len(ids):
      ids = [''.join(x) for x in itertools.product(ch, repeat = 3)]
    if self.get_n_sheets() > len(ids):
      # Too many sheets for PDB format (max 3 chars for sheet id)
      raise Sorry("Too many sheets in annotations: %d" % self.get_n_sheets())
    for i, sh in enumerate(self.sheets):
      sh.sheet_id = ids[i]

  def renumber_helices(self):
    i = 1
    for helix in self.helices:
      helix.set_new_serial(i, adopt_as_id=True)
      i += 1

  def renumber_sheets(self):
    i = 1
    for sheet in self.sheets:
      sheet.sheet_id = self.convert_id(i)
      for strand in sheet.strands:
        strand.sheet_id = "%s" % self.convert_id(i)
      i += 1

  def renumber_helices_and_sheets(self):
    self.renumber_sheets()
    self.renumber_helices()

  def multiply_to_asu_2(self, chain_ids_dict):
    def get_new_chain_id(old_chain_id, n_copy):
      if n_copy == 0 and old_chain_id in chain_ids_dict:
        return old_chain_id
      return chain_ids_dict[old_chain_id][n_copy-1]
    n_copies = len(list(chain_ids_dict.values())[0]) + 1
    new_helices = []
    new_sheets = []
    new_h_serial = 0
    new_sheet_id = 0
    for n_copy in range(n_copies):
      for helix in self.helices:
        new_helix = copy.deepcopy(helix)
        new_helix.erase_hbond_list()
        try:
          new_h_serial += 1
          new_helix.set_new_chain_ids(
              get_new_chain_id(new_helix.start_chain_id, n_copy))
          new_helix.set_new_serial(new_h_serial, adopt_as_id=True)
          new_helices.append(new_helix)
        except KeyError:
          continue
      for sheet in self.sheets:
        new_sheet = copy.deepcopy(sheet)
        new_sheet_id += 1
        new_sheet.sheet_id = self.convert_id(new_sheet_id)
        new_sheet.erase_hbond_list()
        try:
          for strand in new_sheet.strands:
            strand.set_new_chain_ids(
                get_new_chain_id(strand.start_chain_id, n_copy))
            strand.sheet_id = "%s" % self.convert_id(new_sheet_id)
          for reg in new_sheet.registrations:
            if reg is not None:
              reg.set_new_cur_chain_id(
                  get_new_chain_id(reg.cur_chain_id, n_copy))
              reg.set_new_prev_chain_id(
                  get_new_chain_id(reg.prev_chain_id, n_copy))
        except KeyError:
          continue
        new_sheets.append(new_sheet)
    self.helices = new_helices
    self.sheets = new_sheets

  def remove_1hb_helices(self):
    filtered_helices = []
    for h in self.helices:
      if h.get_n_maximum_hbonds() <= 1:
        continue
      else:
        filtered_helices.append(h)
    self.helices = filtered_helices

  def remove_3_10_helices(self):
    filtered_helices = []
    for h in self.helices:
      if h.get_class_as_int() == 5:
        continue
      else:
        filtered_helices.append(h)
    self.helices = filtered_helices

  def change_residue_numbering_in_place(self, renumbering_dictionary):
    """ Changes residue numbers and insertion codes in annotations.
    For cases where one needs to renumber hierarchy and keep annotations
    consisten with the new numbering.
    Not for cases where residues removed or added.
    Therefore residue name stays the same.

    Args:
        renumbering_dictionary (_type_): data structure to match old and new numbering:
        {'chain id':
          {(old resseq, icode):(new resseq, icode)}
        }
    """
    for h in self.helices:
      h.change_residue_numbering_in_place(renumbering_dictionary)
    for sh in self.sheets:
      sh.change_residue_numbering_in_place(renumbering_dictionary)

  def as_cif_loops(self):
    """
    Returns list of loops needed to represent SS annotation. The first for
    helix, others for sheets. If there's no helix, there will be only sheet
    loops. Or empty list if there's nothing to output."""

    loops = []
    if self.get_n_helices() > 0:
      helix_info_prefix = '_struct_conf.'
      helix_info_cif_names = (
            'conf_type_id',
            'id',
            'pdbx_PDB_helix_id',
            'beg_label_comp_id',
            'beg_label_asym_id',
            'beg_label_seq_id',
            'pdbx_beg_PDB_ins_code',
            'end_label_comp_id',
            'end_label_asym_id',
            'end_label_seq_id',
            'pdbx_end_PDB_ins_code',
            'pdbx_PDB_helix_class',
            'details',
            'pdbx_PDB_helix_length')
      helix_loop = iotbx.cif.model.loop(header=(
          ["%s%s" % (helix_info_prefix, x) for x in helix_info_cif_names]))
      for h in self.helices:
        h_dict = h.as_cif_dict()
        row = []
        for cif_name in helix_info_cif_names:
          v = h_dict.get(cif_name, '?')
          row.append(v)
        helix_loop.add_row(row)
      loops.append(helix_loop)
      type_prefix = "_struct_conf_type."
      type_names = ["id", "criteria", "reference"]
      conf_type_ids = set(helix_loop["_struct_conf.conf_type_id"])
      type_loop = iotbx.cif.model.loop(header=(
          ["%s%s" % (type_prefix, x) for x in type_names]))
      for conf_type in conf_type_ids:
        type_loop.add_row([conf_type, "?", "?"])
      loops.append(type_loop)
    if self.get_n_sheets() > 0:
      struct_sheet_loop = iotbx.cif.model.loop(header=(
          '_struct_sheet.id',
          '_struct_sheet.type',
          '_struct_sheet.number_strands',
          '_struct_sheet.details'))
      struct_sheet_order_loop = iotbx.cif.model.loop(header=(
          '_struct_sheet_order.sheet_id',
          '_struct_sheet_order.range_id_1',
          '_struct_sheet_order.range_id_2',
          '_struct_sheet_order.offset',
          '_struct_sheet_order.sense'))
      struct_sheet_range_loop = iotbx.cif.model.loop(header=(
          '_struct_sheet_range.sheet_id',
          '_struct_sheet_range.id',
          '_struct_sheet_range.beg_label_comp_id',
          '_struct_sheet_range.beg_label_asym_id',
          '_struct_sheet_range.beg_label_seq_id',
          '_struct_sheet_range.pdbx_beg_PDB_ins_code',
          '_struct_sheet_range.end_label_comp_id',
          '_struct_sheet_range.end_label_asym_id',
          '_struct_sheet_range.end_label_seq_id',
          '_struct_sheet_range.pdbx_end_PDB_ins_code'))
      struct_sheet_hbond_loop = iotbx.cif.model.loop(header=(
          '_pdbx_struct_sheet_hbond.sheet_id',
          '_pdbx_struct_sheet_hbond.range_id_1',
          '_pdbx_struct_sheet_hbond.range_id_2',
          '_pdbx_struct_sheet_hbond.range_1_label_atom_id',
          '_pdbx_struct_sheet_hbond.range_1_label_comp_id',
          '_pdbx_struct_sheet_hbond.range_1_label_asym_id',
          '_pdbx_struct_sheet_hbond.range_1_label_seq_id',
          '_pdbx_struct_sheet_hbond.range_1_PDB_ins_code',
          '_pdbx_struct_sheet_hbond.range_2_label_atom_id',
          '_pdbx_struct_sheet_hbond.range_2_label_comp_id',
          '_pdbx_struct_sheet_hbond.range_2_label_asym_id',
          '_pdbx_struct_sheet_hbond.range_2_label_seq_id',
          '_pdbx_struct_sheet_hbond.range_2_PDB_ins_code'))
      for sh in self.sheets:
        sh_dict = sh.as_cif_dict()
        # parse it here and toss into loops
        struct_sheet_loop.add_row((
            sh_dict['_struct_sheet.id'],
            sh_dict['_struct_sheet.type'],
            sh_dict['_struct_sheet.number_strands'],
            sh_dict['_struct_sheet.details']))
        for struct_sheet_loop_row in zip(
            sh_dict['_struct_sheet_order.sheet_id'],
            sh_dict['_struct_sheet_order.range_id_1'],
            sh_dict['_struct_sheet_order.range_id_2'],
            sh_dict['_struct_sheet_order.offset'],
            sh_dict['_struct_sheet_order.sense']):
          struct_sheet_order_loop.add_row(struct_sheet_loop_row)
        for struct_sheet_range_row in zip(
            sh_dict['_struct_sheet_range.sheet_id'],
            sh_dict['_struct_sheet_range.id'],
            sh_dict['_struct_sheet_range.beg_label_comp_id'],
            sh_dict['_struct_sheet_range.beg_label_asym_id'],
            sh_dict['_struct_sheet_range.beg_label_seq_id'],
            sh_dict['_struct_sheet_range.pdbx_beg_PDB_ins_code'],
            sh_dict['_struct_sheet_range.end_label_comp_id'],
            sh_dict['_struct_sheet_range.end_label_asym_id'],
            sh_dict['_struct_sheet_range.end_label_seq_id'],
            sh_dict['_struct_sheet_range.pdbx_end_PDB_ins_code']):
          struct_sheet_range_loop.add_row(struct_sheet_range_row)
        for struct_sheet_hbond_row in zip(
            sh_dict['_pdbx_struct_sheet_hbond.sheet_id'],
            sh_dict['_pdbx_struct_sheet_hbond.range_id_1'],
            sh_dict['_pdbx_struct_sheet_hbond.range_id_2'],
            sh_dict['_pdbx_struct_sheet_hbond.range_1_label_atom_id'],
            sh_dict['_pdbx_struct_sheet_hbond.range_1_label_comp_id'],
            sh_dict['_pdbx_struct_sheet_hbond.range_1_label_asym_id'],
            sh_dict['_pdbx_struct_sheet_hbond.range_1_label_seq_id'],
            sh_dict['_pdbx_struct_sheet_hbond.range_1_PDB_ins_code'],
            sh_dict['_pdbx_struct_sheet_hbond.range_2_label_atom_id'],
            sh_dict['_pdbx_struct_sheet_hbond.range_2_label_comp_id'],
            sh_dict['_pdbx_struct_sheet_hbond.range_2_label_asym_id'],
            sh_dict['_pdbx_struct_sheet_hbond.range_2_label_seq_id'],
            sh_dict['_pdbx_struct_sheet_hbond.range_2_PDB_ins_code']):
          struct_sheet_hbond_loop.add_row(struct_sheet_hbond_row)
      loops.append(struct_sheet_loop)
      loops.append(struct_sheet_order_loop)
      loops.append(struct_sheet_range_loop)
      loops.append(struct_sheet_hbond_loop)
    return loops

  def as_mmcif_str(self, data_block_name=None):
    cif_object = iotbx.cif.model.cif()
    if data_block_name is None:
      data_block_name = "phenix"
    cif_object[data_block_name] = self.as_cif_block()
    from six.moves import cStringIO as StringIO
    f = StringIO()
    cif_object.show(out = f)
    return f.getvalue()

  def as_cif_block(self):
    cif_block = iotbx.cif.model.block()
    ss_cif_loops = self.as_cif_loops()
    for loop in ss_cif_loops:
      cif_block.add_loop(loop)
    return cif_block

  def fits_in_pdb_format(self):
    for helix in self.helices :
      if (not helix.fits_in_pdb_format()):  return False
    for sheet in self.sheets :
      if (not sheet.fits_in_pdb_format()):  return False
    return True

  def as_pdb_str(self):
    records = []
    for helix in self.helices :
      records.append(helix.as_pdb_str())
    for sheet in self.sheets :
      records.append(sheet.as_pdb_str())
    return "\n".join(records)

  def as_pdb_or_mmcif_str(self, target_format = 'pdb'):
    # Return str in target format if possible, otherwise in mmcif
    if target_format == 'pdb' and self.fits_in_pdb_format():
      return self.as_pdb_str()
    else:
      return self.as_mmcif_str()

  def as_restraint_groups(self, log=sys.stdout, prefix_scope="",
      add_segid=None):
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
      except RuntimeError as e :
        pass
    for sheet in self.sheets :
      selections.extend(sheet.as_atom_selections(add_segid=add_segid))
    return selections

  def overall_selection(self,add_segid=None,trim_ends_by=None):
    result = ""
    result = self.overall_helices_selection(add_segid=add_segid,trim_ends_by=trim_ends_by)
    s_s = self.overall_sheets_selection(add_segid=add_segid,trim_ends_by=trim_ends_by)
    if len(s_s) > 0 and len(result) > 0:
      result += " or %s" % s_s
      return result
    if len(result) == 0:
      return s_s
    else:
      return result

  def overall_helices_selection(self, add_segid=None,trim_ends_by=None):
    selections = []
    for helix in self.helices:
      try :
        selections.extend(helix.as_atom_selections(add_segid=add_segid,
         trim_ends_by=trim_ends_by))
      except RuntimeError as e :
        pass
    return "(" + ") or (".join(selections) + ")" if len(selections) > 0 else ""

  def overall_sheets_selection(self, add_segid=None,trim_ends_by=None):
    selections = []
    for sheet in self.sheets:
      try:
        selections.extend(sheet.as_atom_selections(add_segid=add_segid,
          trim_ends_by=trim_ends_by))
      except RuntimeError as e :
        pass
    return "(" + ") or (".join(selections) + ")" if len(selections) > 0 else ""


  def overall_helix_selection(self,
                              add_segid=None,
                              helix_types_selection=None,
                              n_hbonds_selection=None,
                              ):
    selections = []
    for helix in self.helices:
      if helix_types_selection is not None:
        if helix.helix_class not in helix_types_selection: continue
      if n_hbonds_selection is not None:
        if helix.get_n_maximum_hbonds()<n_hbonds_selection: continue
      try :
        selections.extend(helix.as_atom_selections(add_segid=add_segid))
      except RuntimeError as e :
        pass
    return "(" + ") or (".join(selections) + ")"

  def overall_sheet_selection(self, add_segid=None):
    selections = []
    for sheet in self.sheets:
      try:
        selections.extend(sheet.as_atom_selections(add_segid=add_segid))
      except RuntimeError as e :
        pass
    return "(" + ") or (".join(selections) + ")"

  def as_bond_selections(self):
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

  def is_empty(self):
    return self.get_n_helices() + self.get_n_sheets() == 0

  def get_n_helix_residues(self):
    result = 0
    for h in self.helices:
      result += h.length
    return result

  def get_n_sheet_residues(self):
    result = 0
    for sh in self.sheets:
      result += sh.get_approx_size()
    return result

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
    for sheet in self.sheets:
      new_sheets+=sheet.split(starting_sheet_id_number=len(new_sheets)+1)
    return annotation(
      helices=copy.deepcopy(self.helices),
      sheets=new_sheets)

  def merge_sheets(self):
    "Group 2-strand sheets into larger sheets if identical component strands"
    # Assumes that all the sheets are non-overlapping
    # First run sheet.split() on all sheets or split_sheets on the annotation
    #  as this requires 2-strand sheets (not 3 or more)

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
          working_sheet=copy.deepcopy(sheet_pointer_0.get(key))
          new_sheets.append(working_sheet)  # now we will extend this sheet

          second_strand=working_sheet.strands[1] # second strand
          next_key=second_strand.as_atom_selections()
          if next_key in used_strand_selections: continue

          while next_key: # points to next sheet to merge 2nd strand
            # find sheet where next_strand is the first strand
            next_sheet=sheet_pointer_0.get(next_key)
            if not next_sheet: break
            used_strand_selections.append(next_key)

            next_strand=copy.deepcopy(next_sheet.strands[1])
            next_registration=copy.deepcopy(next_sheet.registrations[1])

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
    new_annotation=copy.deepcopy(self)
    new_annotation.sheets=new_sheets
    return new_annotation

  def add_helices_and_sheets_simple(self, other_annot):
    self.helices += other_annot.helices
    self.sheets += other_annot.sheets
    self.renumber_helices_and_sheets()

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
    ph=hierarchy.select(sel, copy_atoms=True)  #ph = needed part of hierarchy
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
      rescore_if_zero_scores=True,
      poor_h_bond_weight=0.5,
      keep_self=None):

      assert h1 # h2 can be None
      score_1=None
      score_2=None
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

      if (score_1 is None and score_2 is None ) or \
         (rescore_if_zero_scores and score_1==0 and score_2==0): #nothing yet
        if keep_self:
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
    new_annotation=annotation(
      helices=copy.deepcopy(new_helices),
      sheets=copy.deepcopy(new_sheets))
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
    score_list = sorted(score_list, key = lambda s: s[0], reverse = True)

    remove_list=[]
    keep_list=[]
    for score1,s1 in score_list: # remove things that overlap and lower score
      if not s1 in remove_list and not s1 in keep_list:
        keep_list.append(s1)
      for score2,s2 in score_list:
        if s2==s1 or s2 in keep_list or s2 in remove_list: continue
        if not s1.has_chains_in_common(other=s2): continue
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
        if not h2.has_chains_in_common(other=h1):
          pass # nothting to do
        elif h2.is_similar_to(other=h1,hierarchy=hierarchy,
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
          rescore_if_zero_scores=True,
          keep_self=self.keep_self)
        #print >>out,"SELF : %7.1f\n%s" %(score_1,h1.as_pdb_str())
        #print >>out,"OTHER: %7.1f\n%s" %(score_2,h2.as_pdb_str())

    if overlapping_but_not_matching_pairs:
      #print >>out,"\nOverlapping non-matching pairs:"
      for [h1,h2] in overlapping_but_not_matching_pairs:
        score_1,score_2=self.score_pair(h1,h2,
          maximize_h_bonds=self.maximize_h_bonds,
          hierarchy=hierarchy,
          rescore_if_zero_scores=True,
          max_h_bond_length=self.max_h_bond_length,
          keep_self=self.keep_self)
        #print >>out,"SELF : %7.1f\n%s" %(score_1,h1.as_pdb_str())
        #print >>out,"OTHER: %7.1f\n%s\n" %(score_2,h2.as_pdb_str())

    if unique_h1 or unique_h2:
      #print >>out,"\nNon-matching:"
      for h1 in unique_h1:
        score_1,score_2=self.score_pair(h1,None,
          maximize_h_bonds=self.maximize_h_bonds,
          rescore_if_zero_scores=True,
          hierarchy=hierarchy,
          max_h_bond_length=self.max_h_bond_length,
          keep_self=self.keep_self)
        #print >>out,"SELF : %7.1f\n%s" %(score_1,h1.as_pdb_str())
      for h2 in unique_h2:
        score_1,score_2=self.score_pair(h2,None,
          maximize_h_bonds=self.maximize_h_bonds,
          hierarchy=hierarchy,
          rescore_if_zero_scores=True,
          max_h_bond_length=self.max_h_bond_length,
          keep_self=self.keep_self)
        #print >>out,"OTHER: %7.1f\n%s" %(score_1,h2.as_pdb_str())

    unique=unique_h1+unique_h2

    # Now take the best (or self if desired) of pairs and add on unique
    final_list=[]
    for [h1,h2] in pairs+overlapping_but_not_matching_pairs:
      score_1,score_2=self.score_pair(h1,h2,
        maximize_h_bonds=self.maximize_h_bonds,
        rescore_if_zero_scores=True,
        hierarchy=hierarchy,
        max_h_bond_length=self.max_h_bond_length,
        keep_self=self.keep_self)

      #print "score 1 and 2: ",score_1,score_2
      if score_1>=score_2:
        if not h1 in final_list:final_list.append(h1)
        if h2 in final_list:final_list.remove(h2)
      else:
        if not h2 in final_list:final_list.append(h2)
        if h1 in final_list:final_list.remove(h1)

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

  def chains_included(self):
    "Report list of all chains involved "
    chain_list=[]
    for x in self.helices+self.sheets:
      for c in x.chains_included():
        if not x in chain_list: chain_list.append(x)
    return chain_list

  def overlaps_with(self,other=None,hierarchy=None):
    "Returns True if any element of the annotation overlap"

    a1=self.split_sheets()
    a2=other.split_sheets()
    assert hierarchy
    for h1 in a1.helices:
      for h2 in a2.helices:
        if h1.has_chains_in_common(other=h2) and \
          h1.overlaps_with(other=h2,hierarchy=hierarchy):
          return True

    for s1 in a1.sheets:
      for s2 in a2.sheets:
        if s1.has_chains_in_common(other=s2) and \
           s1.overlaps_with(other=s2,hierarchy=hierarchy):
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
        sorted.append(x.as_pdb_str(set_id_zero=True, force_format = True))
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

class pdb_helix(structure_base):
  _helix_class_array = ['unknown','alpha', 'unknown', 'pi', 'unknown',
        '3_10', 'unknown', 'unknown', 'unknown', 'unknown', 'unknown']
  _cif_helix_classes = {'HELX_RH_AL_P': 1, 'HELX_RH_PI_P':3, 'HELX_RH_3T_P':5}

  def __init__(self,
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
    if (length <= 0):
      raise Sorry("Bad helix length. Check HELIX records.\n%s" % self.as_pdb_str())
    if isinstance(self.helix_class, int):
      self.helix_class = self._helix_class_array[helix_class]
    if self.helix_class not in self._helix_class_array:
      raise Sorry("Bad helix class: %s. Check HELIX records." % helix_class)

    self.start_chain_id = self.parse_chain_id(start_chain_id)
    self.end_chain_id = self.parse_chain_id(end_chain_id)
    if isinstance(self.helix_id, int):
      self.helix_id = "%s" % self.helix_id
      self.helix_id = self.helix_id[:3]
    elif self.helix_id is None or self.helix_id.strip()=="":
      self.adopt_serial_as_id()
    else:
      assert isinstance(self.helix_id, str)
    if self.start_chain_id != self.end_chain_id:
      raise Sorry("Don't know how to deal with helices with multiple "+
        "chain IDs ('%s' vs. '%s')." % (self.start_chain_id, self.end_chain_id))

  @classmethod
  def helix_class_to_int(cls, h_class):
    return cls._helix_class_array.index(h_class)

  @classmethod
  def helix_class_to_str(cls, h_class):
    return cls._helix_class_array[h_class]

  @classmethod
  def get_helix_class(cls, cif_row):
    conf_type_class = cif_row.get('_struct_conf.conf_type_id', None)
    if conf_type_class not in ['HELX_P', 'HELX_RH_AL_P', 'HELX_RH_PI_P', 'HELX_RH_3T_P']:
      return None
    if conf_type_class == 'HELX_P':
      pdbx_class = int(cif_row.get('_struct_conf.pdbx_PDB_helix_class', 0))
      return cls._helix_class_array[pdbx_class]
    else:
      return cls._helix_class_array[cls._cif_helix_classes[conf_type_class]]

  @classmethod
  def from_cif_row(cls, cif_row, serial):
    h_class = cls.get_helix_class(cif_row)
    if h_class is None:
      return None
    start_resname = choose_correct_cif_record(
        cif_row,
        '_struct_conf.beg_auth_comp_id',
        '_struct_conf.beg_label_comp_id')
    start_chain_id = cls.parse_chain_id(
        choose_correct_cif_record(
            cif_row,
            '_struct_conf.beg_auth_asym_id',
            '_struct_conf.beg_label_asym_id'))
    start_resseq = int(choose_correct_cif_record(
        cif_row,
        '_struct_conf.beg_auth_seq_id',
        '_struct_conf.beg_label_seq_id'))
    end_resname = choose_correct_cif_record(
        cif_row,
        '_struct_conf.end_auth_comp_id',
        '_struct_conf.end_label_comp_id')
    end_chain_id = cls.parse_chain_id(
        choose_correct_cif_record(
            cif_row,
            '_struct_conf.end_auth_asym_id',
            '_struct_conf.end_label_asym_id'))
    end_resseq = int(choose_correct_cif_record(
        cif_row,
        '_struct_conf.end_auth_seq_id',
        '_struct_conf.end_label_seq_id'))
    comment = ""
    if ('_struct_conf.details' in cif_row and
        cif_row['_struct_conf.details'] != '?' and
        cif_row['_struct_conf.details'] != '.'):
      comment = cif_row['_struct_conf.details']
    # this is not mandatory item in mmCIF, so we'll try to estimate it in case
    # it is absent. Since it is mmCIF, it should not contain hybrid36 notation,
    # so int() should be fine
    length = int(cif_row.get('_struct_conf.pdbx_PDB_helix_length', 0))
    if length == 0:
      length = end_resseq - start_resseq

    return cls(
      serial=serial,
      helix_id=cif_row.get('_struct_conf.pdbx_PDB_helix_id', None),
      start_resname=start_resname,
      start_chain_id=start_chain_id,
      start_resseq=cls.convert_resseq(start_resseq),
      start_icode=cls.parse_cif_insertion_code(cif_row.get('_struct_conf.pdbx_beg_PDB_ins_code', '.')),
      end_resname=end_resname,
      end_chain_id=end_chain_id,
      end_resseq=cls.convert_resseq(end_resseq),
      end_icode=cls.parse_cif_insertion_code(cif_row.get('_struct_conf.pdbx_end_PDB_ins_code', '.')),
      helix_class=h_class,
      helix_selection=None,
      comment=comment,
      length=length)

  @classmethod
  def from_pdb_record(cls, line):
    "Raises ValueError in case of corrupted HELIX record!!!"
    if len(line) < 76:
      line += " "*(80-len(line))
    return cls(
      serial=cls.convert_id(line[7:10]),
      helix_id=line[11:14].strip(),
      start_resname=line[15:18],
      start_chain_id=cls.parse_chain_id(line[18:20]),
      start_resseq=line[21:25],
      start_icode=line[25],
      end_resname=line[27:30],
      end_chain_id=cls.parse_chain_id(line[30:32]),
      end_resseq=line[33:37],
      end_icode=line[37],
      helix_class=cls.helix_class_to_str(int(line[38:40])),
      helix_selection=None,
      comment=line[40:70],
      length=int(line[71:76])) #string.atoi(line[71:76]))

  @classmethod
  def from_phil_params(cls, helix_params, pdb_hierarchy, h_atoms, cache, serial=0, log=None):
    if log is None:
      log = sys.stdout
    if helix_params.selection is None :
      print("Empty helix at serial %d." % (serial), file=log)
      # continue
    if helix_params.serial_number is not None:
      serial = helix_params.serial_number
    amide_isel = get_amide_isel(cache, helix_params.selection)
    if isinstance(amide_isel, str):
      print(amide_isel, file=log)
      return None
    start_atom = h_atoms[amide_isel[0]]
    end_atom = h_atoms[amide_isel[-1]]
    hbonds = []
    for hb in helix_params.hbond:
      if hb.donor is None:
        print("Donor selection in hbond cannot be None", file=log)
        continue
      if hb.acceptor is None:
        print("Acceptor selection in hbond cannot be None", file=log)
        continue
      hbonds.append((hb.donor, hb.acceptor))
    return cls(
      serial=serial,
      helix_id=helix_params.helix_identifier,
      start_resname=start_atom.parent().resname,
      start_chain_id=start_atom.parent().parent().parent().id,
      start_resseq=start_atom.parent().parent().resseq,
      start_icode=start_atom.parent().parent().icode,
      end_resname=end_atom.parent().resname,
      end_chain_id=end_atom.parent().parent().parent().id,
      end_resseq=end_atom.parent().parent().resseq,
      end_icode=end_atom.parent().parent().icode,
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

  def deep_copy(self):
    return copy.deepcopy(self)

  def get_class_as_int(self):
    return self.helix_class_to_int(self.helix_class)

  def get_class_as_str(self):
    return self.helix_class

  def set_start_resseq(self, resseq):
    self.start_resseq = self.convert_resseq(resseq)

  def set_end_resseq(self, resseq):
    self.end_resseq = self.convert_resseq(resseq)

  def set_new_serial(self, serial, adopt_as_id=False):
    self.serial = hy36encode(3, serial)
    if adopt_as_id:
      self.adopt_serial_as_id()

  def adopt_serial_as_id(self):
    self.helix_id = "%s" % self.serial

  def set_new_chain_ids(self, new_chain_id):
    self.start_chain_id = new_chain_id
    self.end_chain_id = new_chain_id

  def erase_hbond_list(self):
    self.hbond_list = []

  def as_restraint_group(self, log=sys.stdout, prefix_scope="",
      add_segid=None, show_hbonds=False):
    if self.start_chain_id != self.end_chain_id :
      print("Helix chain ID mismatch: starts in %s, ends in %s" % (
        self.start_chain_id, self.end_chain_id), file=log)
      return None
    sele = self.as_atom_selections(add_segid=add_segid)[0]
    if prefix_scope != "" and not prefix_scope.endswith("."):
      prefix_scope += "."
    serial_and_id = ""
    if self.serial is not None and self.id_as_int(self.serial) > 0:
      serial_and_id += "\n  serial_number = %s" % self.serial
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

  def comment_to_cif(self):
    stripped = self.comment.strip()
    if len(stripped) == 0:
      return "?"
    else:
      return stripped

  def as_cif_dict(self):
    """Returns dict. keys - cif field names, values - appropriate values."""
    result = {}
    # XXX This is most default type of helix.
    # Refer here for the others:
    # http://mmcif.wwpdb.org/dictionaries/mmcif_mdb.dic/Items/_struct_conf_type.id.html
    # Note, that PDB itself does not assign proper helix types in mmCIF although it has
    # no problems doing it for PDB format and pretty pictures on the web-site.
    result['conf_type_id'] = "HELX_P"
    result['id'] = self.serial if isinstance(self.serial, int) else self.serial.strip()
    result['pdbx_PDB_helix_id'] = self.helix_id
    result['beg_label_comp_id'] = self.start_resname
    result['beg_label_asym_id'] = self.start_chain_id
    result['beg_label_seq_id'] =  "%d" % self.get_start_resseq_as_int()
    result['pdbx_beg_PDB_ins_code'] = self.icode_to_cif(self.start_icode)
    result['end_label_comp_id'] = self.end_resname
    result['end_label_asym_id'] = self.end_chain_id
    result['end_label_seq_id'] =  "%d" % self.get_end_resseq_as_int()
    result['pdbx_end_PDB_ins_code'] = self.icode_to_cif(self.end_icode)
    result['pdbx_PDB_helix_class'] = self.helix_class_to_int(self.helix_class)
    result['details'] = self.comment_to_cif()
    result['pdbx_PDB_helix_length'] = self.length
    return result

  def as_pdb_or_mmcif_str(self, target_format = 'pdb'):
    # Return str in target format if possible, otherwise in mmcif
    if target_format == 'pdb' and self.fits_in_pdb_format():
      return self.as_pdb_str()
    else:
      return self.as_mmcif_str()

  def fits_in_pdb_format(self):
    if len(self.start_resname.strip()) > 3: return False
    if len(self.end_resname.strip()) > 3: return False
    if len(self.start_chain_id.strip()) > 2: return False
    if len(self.end_chain_id.strip()) > 2: return False
    return True

  def as_mmcif_str(self):
    ann = annotation(helices = [self], sheets = [])
    text = ann.as_mmcif_str()
    return text

  def as_pdb_str(self, set_id_zero=False, force_format = False):
    if (not force_format) and (not self.fits_in_pdb_format()):
      raise AssertionError(
       "Helix does not fit in PDB format. "+
       "Please fix code to use as_pdb_or_mmcif_str instead of as_pdb_str")
    def h_class_to_pdb_int(h_class):
      h_class_int = self.helix_class_to_int(h_class)
      if h_class_int == 0:
        return 1
      return h_class_int
    format = "HELIX  %3s %3s %3s%2s %4s%1s %3s%2s %4s%1s%2d%30s %5d"
    if set_id_zero:
      serial=0
      helix_id=0
    else:
      serial=self.serial
      if self.helix_id is None:
        helix_id=self.serial
      else:
        helix_id=self.helix_id[:3]
    out = format % (self.convert_id(serial),
      helix_id,
      self.start_resname,
      self.start_chain_id, self.start_resseq, self.start_icode,
      self.end_resname, self.end_chain_id, self.end_resseq, self.end_icode,
      h_class_to_pdb_int(self.helix_class), self.comment, self.length)
    return out.strip()

  def as_atom_selections(self, add_segid=None, trim_ends_by=None):
    segid_extra = ""
    if add_segid is not None :
      segid_extra = "and segid '%s' " % add_segid
    if trim_ends_by and \
     self.get_end_resseq_as_int()-self.get_start_resseq_as_int()>2*trim_ends_by:
      resid_start = "%s%s" % (self.convert_resseq(
         self.get_start_resseq_as_int()+trim_ends_by), self.start_icode)
      resid_end = "%s%s" % (self.convert_resseq(
         self.get_end_resseq_as_int()-trim_ends_by), self.end_icode)
    else:  # usual
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
    elif self.helix_class=='pi':
      return self.length-5
    elif self.helix_class=='unknown':
      return 0
    else:
      # Should never happen
      assert 0, "Wrong helix_class creeped in object fields: %d" % (
          self.helix_class)

  def chains_included(self):
    "Report list of all chains involved in this helix"
    chain_list=[self.start_chain_id]
    if self.start_chain_id != self.start_chain_id:
      chain_list.append(self.end_chain_id)
    return chain_list

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

  def change_residue_numbering_in_place(self, renumbering_dictionary):
    self._change_residue_numbering_in_place_helper(renumbering_dictionary)


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
  def __init__(self,
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
      sense):
    adopt_init_args(self, locals())
    # Python 3 prevents comparisons between str and int
    # assert (sheet_id > 0) and (strand_id > 0)
    assert (sheet_id is not None) and (strand_id is not None)
    if sense not in [-1, 0, 1]:
      raise Sorry("Bad sense in SHEET record: '%s'" % sense)
    self.start_chain_id = self.parse_chain_id(start_chain_id)
    self.end_chain_id = self.parse_chain_id(end_chain_id)
    self.set_start_resseq(self.start_resseq)
    self.set_end_resseq(self.end_resseq)
    if self.start_chain_id != self.end_chain_id:
      raise Sorry("Don't know how to deal with strands with different "+
        "chain IDs ('%s' vs. '%s')." % (self.start_chain_id, self.end_chain_id))
    self.approx_length = self.get_end_resseq_as_int() - self.get_start_resseq_as_int()

  @classmethod
  def from_pdb_record(cls, line):
    if len(line) < 76:
      line += " "*(80-len(line))
    return cls(sheet_id=line[11:14],
        strand_id=cls.convert_id(line[7:10]),
        start_resname=line[17:20],
        start_chain_id=cls.parse_chain_id(line[20:22]),
        start_resseq=line[22:26],
        start_icode=line[26],
        end_resname=line[28:31],
        end_chain_id=cls.parse_chain_id(line[31:33]),
        end_resseq=line[33:37],
        end_icode=line[37],
        sense=int(line[38:40]))

  @classmethod
  def from_cif_dict(cls, cif_dict, sense):
    start_resname = choose_correct_cif_record(
        cif_dict,
        '_struct_sheet_range.beg_auth_comp_id',
        '_struct_sheet_range.beg_label_comp_id')
    start_chain_id = choose_correct_cif_record(
        cif_dict,
        '_struct_sheet_range.beg_auth_asym_id',
        '_struct_sheet_range.beg_label_asym_id')
    start_resseq = int(choose_correct_cif_record(
        cif_dict,
        '_struct_sheet_range.beg_auth_seq_id',
        '_struct_sheet_range.beg_label_seq_id'))
    end_resname = choose_correct_cif_record(
        cif_dict,
        '_struct_sheet_range.end_auth_comp_id',
        '_struct_sheet_range.end_label_comp_id')
    end_chain_id = choose_correct_cif_record(
        cif_dict,
        '_struct_sheet_range.end_auth_asym_id',
        '_struct_sheet_range.end_label_asym_id')
    end_resseq = int(choose_correct_cif_record(
        cif_dict,
        '_struct_sheet_range.end_auth_seq_id',
        '_struct_sheet_range.end_label_seq_id'))
    return cls(
        sheet_id=cif_dict['_struct_sheet_range.sheet_id'],
        strand_id=int(cif_dict['_struct_sheet_range.id']),
        start_resname=start_resname,
        start_chain_id=cls.parse_chain_id(start_chain_id),
        start_resseq=cls.convert_resseq(start_resseq),
        start_icode=cls.parse_cif_insertion_code(
            cif_dict.get('_struct_sheet_range.pdbx_beg_PDB_ins_code','.')),
        end_resname=end_resname,
        end_chain_id=cls.parse_chain_id(end_chain_id),
        end_resseq=cls.convert_resseq(end_resseq),
        end_icode=cls.parse_cif_insertion_code(
            cif_dict.get('_struct_sheet_range.pdbx_end_PDB_ins_code', '.')),
        sense=sense,
      )
  def deep_copy(self):
    return copy.deepcopy(self)

  def set_new_chain_ids(self, new_chain_id):
    self.start_chain_id = new_chain_id
    self.end_chain_id = new_chain_id

  def change_residue_numbering_in_place(self, renumbering_dictionary):
    self._change_residue_numbering_in_place_helper(renumbering_dictionary)

  def sense_as_cif(self):
    if self.sense == 0:
      return '?'
    elif self.sense == 1:
      return "parallel"
    elif self.sense == -1:
      return "anti-parallel"
    else:
      raise Sorry("Invalid sense creeped in object: %s", self.sense)

  def set_start_resseq(self, resseq):
    self.start_resseq = self.convert_resseq(resseq)

  def set_end_resseq(self, resseq):
    self.end_resseq = self.convert_resseq(resseq)

  def as_atom_selections(self, add_segid=None,trim_ends_by=None):
    segid_extra = ""
    if add_segid is not None :
      segid_extra = "and segid '%s' " % add_segid
    if trim_ends_by and \
     self.get_end_resseq_as_int()-self.get_start_resseq_as_int()>2*trim_ends_by:
      resid_start = "%s%s" % (self.convert_resseq(
         self.get_start_resseq_as_int()+trim_ends_by), self.start_icode)
      resid_end = "%s%s" % (self.convert_resseq(
         self.get_end_resseq_as_int()-trim_ends_by), self.end_icode)
    else:  # usual
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
  def __init__(self,
      cur_atom,
      cur_resname,
      cur_chain_id,
      cur_resseq,
      cur_icode,
      prev_atom,
      prev_resname,
      prev_chain_id,
      prev_resseq,
      prev_icode):
    adopt_init_args(self, locals())

  @classmethod
  def from_pdb_record(cls, line):
    if len(line.strip()) < 67:
      return None
    return cls(cur_atom=line[41:45],
        cur_resname=line[45:48],
        cur_chain_id=cls.parse_chain_id(line[48:50]),
        cur_resseq=line[50:54],
        cur_icode=line[54],
        prev_atom=line[56:60],
        prev_resname=line[60:63],
        prev_chain_id=cls.parse_chain_id(line[63:65]),
        prev_resseq=line[65:69],
        prev_icode=line[69])

  @staticmethod
  def adjust_cif_atom_name(name):
    if len(name) == 1:
      return " %s  " % name
    elif len(name) == 2:
      return " %s " % name
    elif len(name) == 3:
      return "%s " % name
    return name

  @classmethod
  def from_cif_dict(cls, cif_dict):
    cur_atom = cls.adjust_cif_atom_name(
        choose_correct_cif_record(
            cif_dict,
            '_pdbx_struct_sheet_hbond.range_2_auth_atom_id',
            '_pdbx_struct_sheet_hbond.range_2_label_atom_id'))
    cur_resname = choose_correct_cif_record(
        cif_dict,
        '_pdbx_struct_sheet_hbond.range_2_auth_comp_id',
        '_pdbx_struct_sheet_hbond.range_2_label_comp_id',
        mandatory=False)
    cur_chain_id = choose_correct_cif_record(
        cif_dict,
        '_pdbx_struct_sheet_hbond.range_2_auth_asym_id',
        '_pdbx_struct_sheet_hbond.range_2_label_asym_id',
        mandatory=False)
    cur_resseq = int(choose_correct_cif_record(
        cif_dict,
        '_pdbx_struct_sheet_hbond.range_2_auth_seq_id',
        '_pdbx_struct_sheet_hbond.range_2_label_seq_id'))
    cur_icode = cls.parse_cif_insertion_code(
        cif_dict.get('_pdbx_struct_sheet_hbond.range_2_PDB_ins_code','.'))
    prev_atom = cls.adjust_cif_atom_name(
        choose_correct_cif_record(
          cif_dict,
          '_pdbx_struct_sheet_hbond.range_1_auth_atom_id',
          '_pdbx_struct_sheet_hbond.range_1_label_atom_id'))
    prev_resname = choose_correct_cif_record(
        cif_dict,
        '_pdbx_struct_sheet_hbond.range_1_auth_comp_id',
        '_pdbx_struct_sheet_hbond.range_1_label_comp_id',
        mandatory=False)
    prev_chain_id = choose_correct_cif_record(
        cif_dict,
        '_pdbx_struct_sheet_hbond.range_1_auth_asym_id',
        '_pdbx_struct_sheet_hbond.range_1_label_asym_id',
        mandatory=False)
    prev_resseq = int(choose_correct_cif_record(
        cif_dict,
        '_pdbx_struct_sheet_hbond.range_1_auth_seq_id',
        '_pdbx_struct_sheet_hbond.range_1_label_seq_id'))
    prev_icode = cls.parse_cif_insertion_code(
        cif_dict.get('_pdbx_struct_sheet_hbond.range_1_PDB_ins_code','.'))
    return cls(
        cur_atom=cur_atom,
        cur_resname=cur_resname,
        cur_chain_id=cls.parse_chain_id(cur_chain_id),
        cur_resseq=cls.convert_resseq(cur_resseq),
        cur_icode=cur_icode,
        prev_atom=prev_atom,
        prev_resname=prev_resname,
        prev_chain_id=cls.parse_chain_id(prev_chain_id),
        prev_resseq=cls.convert_resseq(prev_resseq),
        prev_icode=prev_icode)

  def get_cur_resseq_as_int(self):
    if self.cur_resseq is not None:
      return hy36decode(4, self.cur_resseq)
    return None

  def get_prev_resseq_as_int(self):
    if self.prev_resseq is not None:
      return hy36decode(4, self.prev_resseq)
    return None

  def set_cur_resseq(self, resseq):
    self.cur_resseq = self.convert_resseq(resseq)

  def set_prev_resseq(self, resseq):
    self.prev_resseq = self.convert_resseq(resseq)

  def set_new_cur_chain_id(self, new_chain_id):
    self.cur_chain_id = new_chain_id

  def set_new_prev_chain_id(self, new_chain_id):
    self.prev_chain_id = new_chain_id

  def as_atom_selections(self, add_segid=None):
    segid_extra = ""
    if add_segid is not None :
      segid_extra = "and segid '%s' " % add_segid
    sele_base = "chain '%s' %sand resid %s and name %s"
    resid_curr = "%s%s" % (self.cur_resseq,self.cur_icode)
    resid_prev = "%s%s" % (self.prev_resseq,
        self.prev_icode)
    sele_curr = sele_base % (self.cur_chain_id, segid_extra,
        resid_curr, self.cur_atom.strip())
    sele_prev = sele_base % (self.prev_chain_id,segid_extra,
        resid_prev, self.prev_atom.strip())
    return sele_curr, sele_prev

  def change_residue_numbering_in_place(self, renumbering_dictionary):
    self._change_residue_numbering_in_place_helper(renumbering_dictionary, se=["cur", "prev"])

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
    new_register=copy.deepcopy(self)
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
  def __init__(self,
      sheet_id,
      n_strands,
      strands,
      registrations,
      hbond_list=[], # list of (donor, acceptor) selections
      ):
    adopt_init_args(self, locals())
    if isinstance(self.sheet_id, int):
      self.sheet_id = hy36encode(3, self.sheet_id)
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
  def from_cif_rows(cls, sheet_id,
            struct_sheet_order_loop, struct_sheet_range_loop,
            struct_sheet_hbond_loop):
    strands = []
    registrations = []
    number_of_strands = 0
    # counting number of strands
    number_of_strands = struct_sheet_range_loop['_struct_sheet_range.sheet_id'].\
        count(sheet_id)
    for i in range(1, number_of_strands+1):
      f_rows = struct_sheet_range_loop.find_row(kv_dict={
          '_struct_sheet_range.sheet_id' : sheet_id,
          '_struct_sheet_range.id' : str(i),
        })
      if i == 1:
        # the first strand, no sense, no hbond
        strands.append(pdb_strand.from_cif_dict(f_rows[0], 0))
        registrations.append(None)
      else:
        # all the rest
        res_range = None
        registration = None
        sense = None
        if len(f_rows) == 1:
          res_range = f_rows[0]
        sense_rows = struct_sheet_order_loop.find_row(kv_dict = {
            '_struct_sheet_order.sheet_id':sheet_id,
            '_struct_sheet_order.range_id_1':str(i-1),
            '_struct_sheet_order.range_id_2':str(i)
          })
        if len(sense_rows) >= 1:
          sense = 0
          str_sense = sense_rows[0].get('_struct_sheet_order.sense', '.')
          if str_sense == 'parallel':
            sense = 1
          elif str_sense == 'anti-parallel':
            sense = -1
        registration_rows = struct_sheet_hbond_loop.find_row(kv_dict = {
              '_pdbx_struct_sheet_hbond.sheet_id':sheet_id,
              '_pdbx_struct_sheet_hbond.range_id_1':str(i-1),
              '_pdbx_struct_sheet_hbond.range_id_2':str(i)
          })
        if len(registration_rows) >= 1:
          registration = registration_rows[0]
        if (res_range is not None and
            registration is not None and sense is not None):
          strands.append(pdb_strand.from_cif_dict(res_range, sense))
          registrations.append(pdb_strand_register.from_cif_dict(registration))
    return cls(sheet_id=sheet_id,
               n_strands=number_of_strands,
               strands=strands,
               registrations=registrations)

  @classmethod
  def from_phil_params(cls, sheet_params, pdb_hierarchy, h_atoms, cache, log=None):
    if log is None:
      log = sys.stdout
    if sheet_params.first_strand is None:
      raise Sorry("Empty first strand selection")
    sheet_id="1"
    if sheet_params.sheet_id is not None:
      sheet_id =  "%3s" % sheet_params.sheet_id[:3]
    n_strands = len(sheet_params.strand) + 1
    amide_isel = get_amide_isel(cache, sheet_params.first_strand)
    if isinstance(amide_isel, str):
      print(amide_isel, file=log)
      return None
    start_atom = h_atoms[amide_isel[0]]
    end_atom = h_atoms[amide_isel[-1]]
    first_strand = pdb_strand(
        sheet_id=sheet_id,
        strand_id=1,
        start_resname=start_atom.parent().resname,
        start_chain_id=start_atom.parent().parent().parent().id,
        start_resseq=start_atom.parent().parent().resseq,
        start_icode=start_atom.parent().parent().icode,
        end_resname=end_atom.parent().resname,
        end_chain_id=end_atom.parent().parent().parent().id,
        end_resseq=end_atom.parent().parent().resseq,
        end_icode=end_atom.parent().parent().icode,
        sense=0)
    strands = [first_strand]
    registrations = [None]
    for i, strand_param in enumerate(sheet_params.strand):
      amide_isel = get_amide_isel(cache, strand_param.selection)
      if isinstance(amide_isel, str):
        print(amide_isel, file=log)
        return None
      start_atom = h_atoms[amide_isel[0]]
      end_atom = h_atoms[amide_isel[-1]]
      sense = cls.sense_to_int(strand_param.sense)
      strand = pdb_strand(
          sheet_id=sheet_id,
          strand_id=i+2,
          start_resname=start_atom.parent().resname,
          start_chain_id=start_atom.parent().parent().parent().id,
          start_resseq=start_atom.parent().parent().resseq,
          start_icode=start_atom.parent().parent().icode,
          end_resname=end_atom.parent().resname,
          end_chain_id=end_atom.parent().parent().parent().id,
          end_resseq=end_atom.parent().parent().resseq,
          end_icode=end_atom.parent().parent().icode,
          sense=sense)
      reg_cur_atom = None
      if strand_param.bond_start_current is not None:
        reg_cur_sel = cache.iselection(strand_param.bond_start_current)
        if reg_cur_sel is None or len(reg_cur_sel) == 0:
          error_msg = "Error in sheet definition. Whole sheet will be skipped.\n"
          error_msg += "String '%s' selected 0 atoms.\n" % strand_param.bond_start_current
          error_msg += "Most likely the definition of sheet does not match model.\n"
          print(error_msg, file=log)
          return None
        reg_cur_atom = h_atoms[reg_cur_sel[0]]
      if reg_cur_atom is None: # No current atom in registration
        pass
        # raise Sorry("This bond_start_current yields 0 atoms:\n %s" % strand_param.bond_start_current)
      reg_prev_atom = None
      if strand_param.bond_start_previous is not None:
        reg_prev_sel = cache.iselection(strand_param.bond_start_previous)
        if reg_prev_sel is None or len(reg_prev_sel) == 0:
          error_msg = "Error in sheet definition. Whole sheet will be skipped.\n"
          error_msg += "String '%s' selected 0 atoms.\n" % strand_param.bond_start_previous
          error_msg += "Most likely the definition of sheet does not match model.\n"
          print(error_msg, file=log)
          return None
        reg_prev_atom = h_atoms[reg_prev_sel[0]]
      if reg_prev_atom is None: # No previous atom in registration
        pass
        # raise Sorry("This bond_start_previous yields 0 atoms:\n %s" % strand_param.bond_start_previous)
      reg = None
      if reg_cur_atom is not None and reg_prev_atom is not None:
        reg = pdb_strand_register(
            cur_atom=reg_cur_atom.name,
            cur_resname=reg_cur_atom.parent().resname,
            cur_chain_id=reg_cur_atom.parent().parent().parent().id,
            cur_resseq=reg_cur_atom.parent().parent().resseq,
            cur_icode=reg_cur_atom.parent().parent().icode,
            prev_atom=reg_prev_atom.name,
            prev_resname=reg_prev_atom.parent().resname,
            prev_chain_id=reg_prev_atom.parent().parent().parent().id,
            prev_resseq=reg_prev_atom.parent().parent().resseq,
            prev_icode=reg_prev_atom.parent().parent().icode)
      strands.append(strand)
      registrations.append(reg)
    hbonds = []
    for hb in sheet_params.hbond:
      if hb.donor is None:
        print("Donor selection in hbond cannot be None", file=log)
        continue
      if hb.acceptor is None:
        print("Acceptor selection in hbond cannot be None", file=log)
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

  def get_approx_size(self):
    s = 0
    for strand in self.strands:
      s += strand.approx_length
    return s

  def remove_short_strands(self, size=3):
    strand_indices_to_delete = []
    # need to make sure we are not deleting strands from the middle and
    # brake registrations...
    cont = True
    i = 0
    while cont and i < self.n_strands:
      if self.strands[i].approx_length < size:
        if i not in strand_indices_to_delete:
          strand_indices_to_delete.append(i)
      else:
        cont = False
      i += 1
    cont = True
    i = len(self.strands) - 1
    while cont and i >= 0:
      if self.strands[i].approx_length < size:
        if i not in strand_indices_to_delete:
          strand_indices_to_delete.append(i)
      else:
        cont = False
      i -= 1
    if len(strand_indices_to_delete) == self.n_strands:
      self.sheet_id=0
      self.n_strands=0
      self.strands=[]
      self.registrations=[]
    else:
      for i in reversed(sorted(strand_indices_to_delete)):
        del self.strands[i]
        del self.registrations[i]
        self.n_strands -= 1
      self.registrations[0] = None
      self.strands[0].sense = 0
      if len(strand_indices_to_delete) > 0:
        self.renumber_strands()
        self.erase_hbond_list()

  def change_residue_numbering_in_place(self, renumbering_dictionary):
    self.erase_hbond_list()
    for st in self.strands:
      st.change_residue_numbering_in_place(renumbering_dictionary)
    for reg in self.registrations:
      if reg is not None:
        reg.change_residue_numbering_in_place(renumbering_dictionary)

  def deep_copy(self):
    return copy.deepcopy(self)

  def renumber_strands(self):
    i = 1
    for s in self.strands:
      s.strand_id = i
      i += 1

  def erase_hbond_list(self):
    self.hbond_list = []

  def add_strand(self, strand):
    self.strands.append(strand)

  def add_registration(self, registration):
    self.registrations.append(registration)

  def as_atom_selections(self, trim_ends_by=None, add_segid=None):
    strand_selections = []
    for strand in self.strands:
      strand_selections.append(strand.as_atom_selections(add_segid=add_segid,
       trim_ends_by=trim_ends_by))
    return strand_selections

  def get_n_defined_hbonds(self):
    if self.hbond_list is not None:
      return len(self.hbond_list)
    return 0

  def as_pdb_or_mmcif_str(self, target_format = 'pdb'):
    # Return str in target format if possible, otherwise in mmcif
    if target_format == 'pdb' and self.fits_in_pdb_format():
      return self.as_pdb_str()
    else:
      return self.as_mmcif_str()

  def fits_in_pdb_format(self):
    for strand in self.strands:
      if len(strand.start_resname.strip()) > 3: return False
      if len(strand.end_resname.strip()) > 3: return False
      if len(strand.start_chain_id.strip()) > 2: return False
      if len(strand.end_chain_id.strip()) > 2: return False
    return True

  def as_mmcif_str(self):
    ann = annotation(helices = [], sheets = [self])
    text = ann.as_mmcif_str()
    return text

  def as_pdb_str(self, strand_id=None, set_id_zero=False, force_format = False):
    if (not force_format) and (not self.fits_in_pdb_format()):
      raise AssertionError("Sheet does not fit in PDB format"+
       "Please fix code to use as_pdb_or_mmcif_str instead of as_pdb_str")
    assert len(self.strands) == len(self.registrations)
    lines = []
    for strand, reg in zip(self.strands, self.registrations):
      format1 = "SHEET  %3s %3s%2d %3s%2s%4s%1s %3s%2s%4s%1s%2d"
      format2 = "%4s%3s%2s%4s%1s %4s%3s%2s%4s%1s"
      # print "STRAND, REG", strand, reg
      if set_id_zero:
        strand_id_for_print=0
        sheet_id_for_print=self.convert_id(0)
      else:
        strand_id_for_print=strand.strand_id
        sheet_id_for_print=self.convert_id(strand.sheet_id)

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

  def as_cif_dict(self):
    """Returns dict. keys - cif field names, values - appropriate values,
    lists where needed."""
    sheet_id_to_output = self.sheet_id if isinstance(self.sheet_id, int) else self.sheet_id.strip()
    result = {}
    result['_struct_sheet.id'] = sheet_id_to_output
    result['_struct_sheet.type'] = '?'
    result['_struct_sheet.number_strands'] = self.n_strands
    result['_struct_sheet.details'] = '?'

    result['_struct_sheet_order.sheet_id'] = []
    result['_struct_sheet_order.range_id_1'] = []
    result['_struct_sheet_order.range_id_2'] = []
    result['_struct_sheet_order.offset'] = []
    result['_struct_sheet_order.sense'] = []

    result['_struct_sheet_range.sheet_id'] = []
    result['_struct_sheet_range.id'] = []
    result['_struct_sheet_range.beg_label_comp_id'] = []
    result['_struct_sheet_range.beg_label_asym_id'] = []
    result['_struct_sheet_range.beg_label_seq_id'] = []
    result['_struct_sheet_range.pdbx_beg_PDB_ins_code'] = []
    result['_struct_sheet_range.end_label_comp_id'] = []
    result['_struct_sheet_range.end_label_asym_id'] = []
    result['_struct_sheet_range.end_label_seq_id'] = []
    result['_struct_sheet_range.pdbx_end_PDB_ins_code'] = []

    result['_pdbx_struct_sheet_hbond.sheet_id'] = []
    result['_pdbx_struct_sheet_hbond.range_id_1'] = []
    result['_pdbx_struct_sheet_hbond.range_id_2'] = []
    result['_pdbx_struct_sheet_hbond.range_1_label_atom_id'] = []
    result['_pdbx_struct_sheet_hbond.range_1_label_comp_id'] = []
    result['_pdbx_struct_sheet_hbond.range_1_label_asym_id'] = []
    result['_pdbx_struct_sheet_hbond.range_1_label_seq_id'] = []
    result['_pdbx_struct_sheet_hbond.range_1_PDB_ins_code'] = []
    result['_pdbx_struct_sheet_hbond.range_2_label_atom_id'] = []
    result['_pdbx_struct_sheet_hbond.range_2_label_comp_id'] = []
    result['_pdbx_struct_sheet_hbond.range_2_label_asym_id'] = []
    result['_pdbx_struct_sheet_hbond.range_2_label_seq_id'] = []
    result['_pdbx_struct_sheet_hbond.range_2_PDB_ins_code'] = []

    for i, strand, registration in zip(range(self.n_strands), self.strands, self.registrations):
      # _struct_sheet_order
      if i != 0:
        result['_struct_sheet_order.sheet_id'].append(sheet_id_to_output)
        result['_struct_sheet_order.range_id_1'].append(i)
        result['_struct_sheet_order.range_id_2'].append(i+1)
        result['_struct_sheet_order.offset'].append('?')
        result['_struct_sheet_order.sense'].append(strand.sense_as_cif())

      if strand is not None: # should always be True
        result['_struct_sheet_range.sheet_id'].append(sheet_id_to_output)
        result['_struct_sheet_range.id'].append(i+1)
        result['_struct_sheet_range.beg_label_comp_id'].append(strand.start_resname)
        result['_struct_sheet_range.beg_label_asym_id'].append(strand.start_chain_id)
        result['_struct_sheet_range.beg_label_seq_id'].append(annotation.resseq_as_int(strand.start_resseq))
        result['_struct_sheet_range.pdbx_beg_PDB_ins_code'].append(self.icode_to_cif(strand.start_icode))
        result['_struct_sheet_range.end_label_comp_id'].append(strand.end_resname)
        result['_struct_sheet_range.end_label_asym_id'].append(strand.end_chain_id)
        result['_struct_sheet_range.end_label_seq_id'].append(annotation.resseq_as_int(strand.end_resseq))
        result['_struct_sheet_range.pdbx_end_PDB_ins_code'].append(self.icode_to_cif(strand.end_icode))

      if registration is not None:
        result['_pdbx_struct_sheet_hbond.sheet_id'].append(sheet_id_to_output)
        result['_pdbx_struct_sheet_hbond.range_id_1'].append(i)
        result['_pdbx_struct_sheet_hbond.range_id_2'].append(i+1)
        result['_pdbx_struct_sheet_hbond.range_1_label_atom_id'].append(registration.prev_atom.strip())
        result['_pdbx_struct_sheet_hbond.range_1_label_comp_id'].append(registration.prev_resname)
        result['_pdbx_struct_sheet_hbond.range_1_label_asym_id'].append(registration.prev_chain_id)
        result['_pdbx_struct_sheet_hbond.range_1_label_seq_id'].append(annotation.resseq_as_int(registration.prev_resseq))
        result['_pdbx_struct_sheet_hbond.range_1_PDB_ins_code'].append(self.icode_to_cif(registration.prev_icode))
        result['_pdbx_struct_sheet_hbond.range_2_label_atom_id'].append(registration.cur_atom.strip())
        result['_pdbx_struct_sheet_hbond.range_2_label_comp_id'].append(registration.cur_resname)
        result['_pdbx_struct_sheet_hbond.range_2_label_asym_id'].append(registration.cur_chain_id)
        result['_pdbx_struct_sheet_hbond.range_2_label_seq_id'].append(annotation.resseq_as_int(registration.cur_resseq))
        result['_pdbx_struct_sheet_hbond.range_2_PDB_ins_code'].append(self.icode_to_cif(registration.cur_icode))
    return result


  def as_restraint_group(self, log=sys.stdout, prefix_scope="",
      add_segid=None, show_hbonds=False):
    if len(self.strands) == 0 :
      return None
    selections = []
    senses = []
    reg_curr = []
    reg_prev = []
    for (strand,registration) in zip(self.strands, self.registrations):
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
    for (sele, sense, curr, prev) in zip(selections,senses,reg_curr,reg_prev):
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
    if prefix_scope != "" and not prefix_scope.endswith("."):
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
    for s1,s2,r2 in zip(
      self.strands[:-1],
      self.strands[1:],self.registrations[1:]):
      self_id="%d" %(len(new_sheets)+1)
      new_s1=copy.deepcopy(s1)
      new_s1.sense=0
      new_s1.strand_id=1
      new_s1.sheet_id="%d" %(starting_sheet_id_number)
      new_s2=copy.deepcopy(s2)
      new_s2.strand_id=2
      new_s2.sheet_id="%d" %(starting_sheet_id_number)
      new_r2=copy.deepcopy(r2)
      new_sheet=pdb_sheet(
             sheet_id="%d" %(starting_sheet_id_number),
             n_strands=2,
             strands=[new_s1,new_s2],
             registrations=[None,new_r2])
      new_sheets.append(new_sheet)
      starting_sheet_id_number+=1
    return new_sheets

  def chains_included(self):
    "Report list of all chains involved in this sheet"
    chain_list=[]
    for strand in self.strands:
      if not strand.start_chain_id in chain_list:
        chain_list.append(strand.start_chain_id)
      if strand.start_chain_id != strand.start_chain_id and \
         not strand.end_chain_id in chain_list:
          chain_list.append(strand.end_chain_id)
    return chain_list

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
          ph=hierarchy.select(sel, copy_atoms=True)  #ph = needed part of hierarchy
          if ph.overall_counts().n_residues==0: continue

          atom_selection=self.combine_atom_selections(pair_2,require_all=True)
          if not atom_selection: continue
          sel = asc.selection(string = atom_selection)
          ph=hierarchy.select(sel, copy_atoms=True)  #ph = needed part of hierarchy
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
