"""
Tools for reading, writing, and manipulating PDB-formatted files and
for managing their data as a PDB hierarchy.
"""
from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex

import boost_adaptbx.boost.python as bp
from six.moves import zip
ext = bp.import_ext("iotbx_pdb_ext")
from iotbx_pdb_ext import *

import iotbx.pdb.records
import iotbx.pdb.hierarchy
from iotbx.pdb.experiment_type import experiment_type
from scitbx import matrix

from iotbx.pdb.atom_name_interpretation import \
  interpreters as protein_atom_name_interpreters
import scitbx.array_family.shared # import dependency
import scitbx.stl.set
from libtbx import smart_open
from libtbx.str_utils import show_string
from libtbx.utils import plural_s, hashlib_md5, date_and_time, to_bytes, Sorry
from libtbx import Auto
from six.moves import cStringIO as StringIO
import sys
import calendar
import six
import os
op = os.path

def construct_special_position_settings(
      crystal_symmetry,
      special_position_settings=None,
      weak_symmetry=False,
      min_distance_sym_equiv=0.5,
      u_star_tolerance=0):
  """Construct the special-position settings for a crystal symmetry"""
  #crystal_symmetry = crystal_symmetry(
  #  crystal_symmetry=special_position_settings,
  #  weak_symmetry=weak_symmetry)
  if (crystal_symmetry is None): return None
  if (special_position_settings is not None):
    min_distance_sym_equiv=special_position_settings.min_distance_sym_equiv()
    u_star_tolerance = special_position_settings.u_star_tolerance()
  return crystal_symmetry.special_position_settings(
    min_distance_sym_equiv=min_distance_sym_equiv,
    u_star_tolerance=u_star_tolerance)

def is_pdb_file(file_name):
  """Return True if this is a PDB file"""
  for known_binary_extension in ['mtz', 'ccp4', 'mrc', 'pickle', 'pkl']:
    if file_name.endswith(known_binary_extension):
      return False
  with smart_open.for_reading(file_name=file_name) as f:
    pdb_raw_records = f.read().splitlines()
  for pdb_str in pdb_raw_records:
    if (pdb_str.startswith("CRYST1")):
      try: cryst1 = iotbx.pdb.records.cryst1(pdb_str=pdb_str)
      except iotbx.pdb.records.FormatError: continue
      if (cryst1.ucparams is not None and cryst1.sgroup is not None):
        return True
    elif (   pdb_str.startswith("ATOM  ")
          or pdb_str.startswith("HETATM")): # PDB OK
      try: pdb_inp = ext.input(
        source_info=None, lines=flex.std_string([pdb_str]))
      except KeyboardInterrupt: raise
      except Exception: continue
      if (pdb_inp.atoms().size() == 1):
        atom = pdb_inp.atoms()[0]
        if (atom.name != "    "):
          return True
  return False

def is_pdb_mmcif_file(file_name):
  """Return True if this is an mmCIF file"""
  try:
    cif_model = iotbx.cif.reader(file_path=file_name).model()
    cif_block = cif_model.values()[0]
    if "_atom_site" in cif_block:
      return True
  except Exception as e:
    return False

def systematic_chain_ids():
  """Return a list of possible 2-character chain IDS"""
  import string
  u, l, d = string.ascii_uppercase, string.ascii_lowercase, string.digits
  _ = result = list(u)
  _.extend(l)
  _.extend(d)
  def xy(first, second):
    for f in first:
      for s in second:
        _.append(f+s)
  a = u+l+d
  xy(u, a)
  xy(l, a)
  xy(d, a)
  return result

cns_dna_rna_residue_names = {
  "ADE": "A",
  "CYT": "C",
  "GUA": "G",
  "THY": "T",
  "URI": "U"
}

mon_lib_dna_rna_cif = set(["AD", "A", "CD", "C", "GD", "G", "TD", "U"])

rna_dna_reference_residue_names = {
  "A": "?A",
  "C": "?C",
  "G": "?G",
  "U": "U",
  "T": "DT",
  "+A": "?A",
  "+C": "?C",
  "+G": "?G",
  "+U": "U",
  "+T": "DT",
  "DA": "DA",
  "DC": "DC",
  "DG": "DG",
  "DT": "DT",
  "ADE": "?A",
  "CYT": "?C",
  "GUA": "?G",
  "URI": "U",
  "THY": "DT",
  # "AR": "A",
  # "CR": "C",
  # "GR": "G",
  # "UR": "U",
  "AD": "DA",
  "CD": "DC",
  "GD": "DG",
  "TD": "DT"
}

def get_one_letter_rna_dna_name(resname):
  """ Truncate the residue name to one letter. Return the letter or None."""
  if resname is None:
    return None
  assert isinstance(resname, str)
  result = resname.strip()
  result = result.strip("r").strip("R")
  result = result.strip("d").strip("D")
  result = result.strip("+")
  if len(result) == 1:
    return result.upper()
  else:
    result = cns_dna_rna_residue_names.get(result, None)
    return result

def rna_dna_reference_residue_name(common_name):
  """Return the standard RNA or DNA reference name for this common name"""
  return rna_dna_reference_residue_names.get(common_name.strip().upper())

rna_dna_atom_names_reference_to_mon_lib_translation_dict = {
  " C1'": "C1*",
  " C2 ": "C2",
  " C2'": "C2*",
  " C3'": "C3*",
  " C4 ": "C4",
  " C4'": "C4*",
  " C5 ": "C5",
  " C5'": "C5*",
  " C6 ": "C6",
  " C7 ": "C5M",
  " C8 ": "C8",
  " H1 ": "H1",
  " H1'": "H1*",
  " H2 ": "H2",
# " H2'": special case: rna: "H2*", dna: "H2*1"
  " H21": "H21",
  " H22": "H22",
  " H3 ": "H3",
  " H3'": "H3*",
  " H4'": "H4*",
  " H41": "H41",
  " H42": "H42",
  " H5 ": "H5",
  " H5'": "H5*1",
  " H6 ": "H6",
  " H61": "H61",
  " H62": "H62",
  " H71": "H5M1",
  " H72": "H5M2",
  " H73": "H5M3",
  " H8 ": "H8",
  " N1 ": "N1",
  " N2 ": "N2",
  " N3 ": "N3",
  " N4 ": "N4",
  " N6 ": "N6",
  " N7 ": "N7",
  " N9 ": "N9",
  " O2 ": "O2",
  " O2'": "O2*",
  " O3'": "O3*",
  " O4 ": "O4",
  " O4'": "O4*",
  " O5'": "O5*",
  " O6 ": "O6",
  " OP1": "O1P",
  " OP2": "O2P",
  " OP3": "O3T",
  " P  ": "P",
  "H2''": "H2*2",
  "H5''": "H5*2",
  "HO2'": "HO2*",
  "HO3'": "HO3*",
  "HO5'": "HO5*",
  "HOP3": "HOP3" # added to monomer library
}

protein_atom_names_backbone = [
  " N  ",
  " CA ",
  " C  ",
  " O  ",
  " OXT",
  " H  ",
  " D  ",
  " HXT",
  " DXT",
  " HA ",
  " HA2",
  " HA3",
  " DA ",
  " DA2",
  " DA3",
  " H1 ",
  " H2 ",
  " H3 ",
  " D1 ",
  " D2 ",
  " D3 ",
]

rna_dna_atom_names_backbone_reference_set = set([
  " C1'",
  " C2'",
  " C3'",
  " C4'",
  " C5'",
  " H1'",
  " H2'",
  " H3'",
  " H4'",
  " H5'",
  " O2'",
  " O3'",
  " O4'",
  " O5'",
  " OP1",
  " OP2",
  " OP3",
  " P  ",
  "H2''",
  "H5''",
  "HO2'",
  "HO3'",
  "HO5'",
  "HOP3"
])

rna_dna_atom_name_aliases = [
  ("1D2", " H21", "G DG"),
  ("1D2'", " H2'", "ANY"),
  ("1D2*", " H2'", "ANY"),
  ("1D4", " H41", "C DC"),
  ("1D5'", " H5'", "ANY"),
  ("1D5*", " H5'", "ANY"),
  ("1D5M", " H71", "DT"),
  ("1D6", " H61", "A DA"),
  ("1H2", " H21", "G DG"),
  ("1H2'", " H2'", "ANY"),
  ("1H2*", " H2'", "ANY"),
  ("1H4", " H41", "C DC"),
  ("1H5'", " H5'", "ANY"),
  ("1H5*", " H5'", "ANY"),
  ("1H5M", " H71", "DT"),
  ("1H6", " H61", "A DA"),
  ("2D2", " H22", "G DG"),
  ("2D2'", "H2''", "DA DC DG DT"),
  ("2D2*", "H2''", "DA DC DG DT"),
  ("2D4", " H42", "C DC"),
  ("2D5'", "H5''", "ANY"),
  ("2D5*", "H5''", "ANY"),
  ("2D5M", " H72", "DT"),
  ("2D6", " H62", "A DA"),
  ("2DO'", "HO2'", "A C G U"),
  ("2DO*", "HO2'", "A C G U"),
  ("2DOP", "HOP2", "ANY"),
  ("2H2", " H22", "G DG"),
  ("2H2'", "H2''", "DA DC DG DT"),
  ("2H2*", "H2''", "DA DC DG DT"),
  ("2H4", " H42", "C DC"),
  ("2H5'", "H5''", "ANY"),
  ("2H5*", "H5''", "ANY"),
  ("2H5M", " H72", "DT"),
  ("2H6", " H62", "A DA"),
  ("2HO'", "HO2'", "A C G U"),
  ("2HO*", "HO2'", "A C G U"),
  ("2HOP", "HOP2", "ANY"),
  ("3D5M", " H73", "DT"),
  ("3DOP", "HOP3", "ANY"),
  ("3H5M", " H73", "DT"),
  ("3HOP", "HOP3", "ANY"),
  ("C1'", " C1'", "ANY"),
  ("C1*", " C1'", "ANY"),
  ("C2", " C2 ", "ANY"),
  ("C2'", " C2'", "ANY"),
  ("C2*", " C2'", "ANY"),
  ("C3'", " C3'", "ANY"),
  ("C3*", " C3'", "ANY"),
  ("C4", " C4 ", "ANY"),
  ("C4'", " C4'", "ANY"),
  ("C4*", " C4'", "ANY"),
  ("C5", " C5 ", "ANY"),
  ("C5'", " C5'", "ANY"),
  ("C5*", " C5'", "ANY"),
  ("C5M", " C7 ", "DT"),
  ("C6", " C6 ", "ANY"),
  ("C7", " C7 ", "DT"),
  ("C8", " C8 ", "A G DA DG"),
  ("D1", " H1 ", "G DG"),
  ("D1'", " H1'", "ANY"),
  ("D1*", " H1'", "ANY"),
  ("D2", " H2 ", "A DA"),
  ("D2'", " H2'", "ANY"),
  ("D2*", " H2'", "ANY"),
  ("D2''", "H2''", "DA DC DG DT"),
  ("D2'1", " H2'", "ANY"),
  ("D2*1", " H2'", "ANY"),
  ("D2'2", "H2''", "DA DC DG DT"),
  ("D2*2", "H2''", "DA DC DG DT"),
  ("D21", " H21", "G DG"),
  ("D22", " H22", "G DG"),
  ("D3", " H3 ", "U DT"),
  ("D3'", " H3'", "ANY"),
  ("D3*", " H3'", "ANY"),
  ("D3T", "HO3'", "ANY"),
  ("D4'", " H4'", "ANY"),
  ("D4*", " H4'", "ANY"),
  ("D41", " H41", "C DC"),
  ("D42", " H42", "C DC"),
  ("D5", " H5 ", "C U DC"),
  ("D5'", " H5'", "ANY"),
  ("D5*", "HO5'", "ANY"),
  ("D5''", "H5''", "ANY"),
  ("D5'1", " H5'", "ANY"),
  ("D5*1", " H5'", "ANY"),
  ("D5'2", "H5''", "ANY"),
  ("D5*2", "H5''", "ANY"),
  ("D5M1", " H71", "DT"),
  ("D5M2", " H72", "DT"),
  ("D5M3", " H73", "DT"),
  ("D5T", "HO5'", "ANY"),
  ("D6", " H6 ", "C U DC DT"),
  ("D61", " H61", "A DA"),
  ("D62", " H62", "A DA"),
  ("D71", " H71", "DT"),
  ("D72", " H72", "DT"),
  ("D73", " H73", "DT"),
  ("D8", " H8 ", "A G DA DG"),
  ("DO2'", "HO2'", "A C G U"),
  ("DO2*", "HO2'", "A C G U"),
  ("H1", " H1 ", "G DG"),
  ("H1'", " H1'", "ANY"),
  ("H1*", " H1'", "ANY"),
  ("H2", " H2 ", "A DA"),
  ("H2'", " H2'", "ANY"),
  ("H2*", " H2'", "ANY"),
  ("H2''", "H2''", "DA DC DG DT"),
  ("H2'1", " H2'", "ANY"),
  ("H2*1", " H2'", "ANY"),
  ("H2'2", "H2''", "DA DC DG DT"),
  ("H2*2", "H2''", "DA DC DG DT"),
  ("H21", " H21", "G DG"),
  ("H22", " H22", "G DG"),
  ("H3", " H3 ", "U DT"),
  ("H3'", " H3'", "ANY"),
  ("H3*", " H3'", "ANY"),
  ("H3T", "HO3'", "ANY"),
  ("H4'", " H4'", "ANY"),
  ("H4*", " H4'", "ANY"),
  ("H41", " H41", "C DC"),
  ("H42", " H42", "C DC"),
  ("H5", " H5 ", "C U DC"),
  ("H5'", " H5'", "ANY"),
  ("H5*", "HO5'", "ANY"),
  ("H5''", "H5''", "ANY"),
  ("H5'1", " H5'", "ANY"),
  ("H5*1", " H5'", "ANY"),
  ("H5'2", "H5''", "ANY"),
  ("H5*2", "H5''", "ANY"),
  ("H5M1", " H71", "DT"),
  ("H5M2", " H72", "DT"),
  ("H5M3", " H73", "DT"),
  ("H5T", "HO5'", "ANY"),
  ("H6", " H6 ", "C U DC DT"),
  ("H61", " H61", "A DA"),
  ("H62", " H62", "A DA"),
  ("H71", " H71", "DT"),
  ("H72", " H72", "DT"),
  ("H73", " H73", "DT"),
  ("H8", " H8 ", "A G DA DG"),
  ("HO2'", "HO2'", "A C G U"),
  ("HO2*", "HO2'", "A C G U"),
  ("HO3'", "HO3'", "ANY"),
  ("HO3*", "HO3'", "ANY"),
  ("HO5'", "HO5'", "ANY"),
  ("HO5*", "HO5'", "ANY"),
  ("HOP2", "HOP2", "ANY"),
  ("HOP3", "HOP3", "ANY"),
  ("N1", " N1 ", "ANY"),
  ("N2", " N2 ", "G DG"),
  ("N3", " N3 ", "ANY"),
  ("N4", " N4 ", "C DC"),
  ("N6", " N6 ", "A DA"),
  ("N7", " N7 ", "A G DA DG"),
  ("N9", " N9 ", "A G DA DG"),
  ("O1P", " OP1", "ANY"),
  ("O2", " O2 ", "C U DC DT"),
  ("O2'", " O2'", "A C G U"),
  ("O2*", " O2'", "A C G U"),
  ("O2P", " OP2", "ANY"),
  ("O3'", " O3'", "ANY"),
  ("O3*", " O3'", "ANY"),
  ("O3P", " OP3", "ANY"),
  ("O3T", " OP3", "ANY"),
  ("O4", " O4 ", "U DT"),
  ("O4'", " O4'", "ANY"),
  ("O4*", " O4'", "ANY"),
  ("O5'", " O5'", "ANY"),
  ("O5*", " O5'", "ANY"),
  ("O5T", " OP3", "ANY"),
  ("O6", " O6 ", "G DG"),
  ("OP1", " OP1", "ANY"),
  ("OP2", " OP2", "ANY"),
  ("OP3", " OP3", "ANY"),
  ("P", " P  ", "ANY")]

def __rna_dna_atom_names_backbone_aliases():
  result = {}
  for a,r,f in rna_dna_atom_name_aliases:
    if (r in rna_dna_atom_names_backbone_reference_set):
      result[a] = r
  return result
rna_dna_atom_names_backbone_aliases = __rna_dna_atom_names_backbone_aliases()

class rna_dna_atom_names_interpretation(object):

  def __init__(self, residue_name, atom_names):
    if (residue_name == "T"):
      residue_name = "DT"
    else:
      assert residue_name in ["?A", "?C", "?G",
                              "A", "C", "G", "U",
                              "DA", "DC", "DG", "DT", "T"]
    self.residue_name = residue_name
    self.atom_names = atom_names
    rna_dna_atom_names_interpretation_core(self)

  def unexpected_atom_names(self):
    result = []
    for atom_name,info in zip(self.atom_names, self.infos):
      if (info.reference_name is None):
        result.append(atom_name)
    return result

  def mon_lib_names(self):
    result = []
    for info in self.infos:
      rn = info.reference_name
      if (rn is None):
        result.append(None)
      else:
        mn = rna_dna_atom_names_reference_to_mon_lib_translation_dict.get(rn)
        if (mn is not None):
          result.append(mn)
        elif (rn == " H2'"):
          if (self.residue_name.startswith("D")):
            result.append("H2*1")
          else:
            result.append("H2*")
        else:
          assert rn == "HOP3", '%s!="HOP3"' % rn # only atom not covered by monomer library
          result.append(None)
    return result

class residue_name_plus_atom_names_interpreter(object):
  """Create an object that has standard values of work_residue_name and
     atom_name_interpretation"""
  def __init__(self,
        residue_name,
        atom_names,
        translate_cns_dna_rna_residue_names=None,
        return_mon_lib_dna_name=False):
    work_residue_name = residue_name.strip().upper()
    if (len(work_residue_name) == 0):
      self.work_residue_name = None
      self.atom_name_interpretation = None
      return
    from iotbx.pdb.amino_acid_codes import three_letter_l_given_three_letter_d
    l_aa_rn = three_letter_l_given_three_letter_d.get(work_residue_name)
    if (l_aa_rn is None):
      d_aa_rn = None
    else:
      d_aa_rn = work_residue_name
      work_residue_name = l_aa_rn
    protein_interpreter = protein_atom_name_interpreters.get(
      work_residue_name)
    atom_name_interpretation = None
    if (protein_interpreter is not None):
      atom_name_interpretation = protein_interpreter.match_atom_names(
        atom_names=atom_names)
      if (atom_name_interpretation is not None):
        atom_name_interpretation.d_aa_residue_name = d_aa_rn
    else:
      assert d_aa_rn is None
      if (    translate_cns_dna_rna_residue_names is not None
          and not translate_cns_dna_rna_residue_names
          and work_residue_name in cns_dna_rna_residue_names):
        rna_dna_ref_residue_name = None
      else:
        rna_dna_ref_residue_name = rna_dna_reference_residue_name(
          common_name=work_residue_name)
      if (rna_dna_ref_residue_name is not None):
        atom_name_interpretation = rna_dna_atom_names_interpretation(
          residue_name=rna_dna_ref_residue_name,
          atom_names=atom_names)
        if (atom_name_interpretation.n_unexpected != 0):
          if (    len(atom_names) == 1
              and work_residue_name in mon_lib_dna_rna_cif):
            self.work_residue_name = None
            self.atom_name_interpretation = None
            return
          if (    translate_cns_dna_rna_residue_names is None
              and work_residue_name in cns_dna_rna_residue_names):
            atom_name_interpretation = None
        if (atom_name_interpretation is not None):
          work_residue_name = atom_name_interpretation.residue_name
          if (return_mon_lib_dna_name):
            work_residue_name = {
              "A": "A",
              "C": "C",
              "G": "G",
              "U": "U",
              "DA": "AD",
              "DC": "CD",
              "DG": "GD",
              "DT": "TD"}[work_residue_name]
    self.work_residue_name = work_residue_name
    self.atom_name_interpretation = atom_name_interpretation

#
# check out merge_files_and_check_for_overlap
class combine_unique_pdb_files(object):

  def __init__(self, file_names):
    self.file_name_registry = {}
    self.md5_registry = {}
    self.unique_file_names = []
    self.raw_records = []
    self.raw_text_block_list = []
    for file_name in file_names:
      if (file_name in self.file_name_registry):
        self.file_name_registry[file_name] += 1
      else:
        self.file_name_registry[file_name] = 1
        with smart_open.for_reading(file_name=file_name) as f:
          text = f.read()
          r = [s.expandtabs().rstrip() for s in text.splitlines()]
        m = hashlib_md5()
        m.update(to_bytes("\n".join(r), codec='utf8'))
        m = m.hexdigest()
        l = self.md5_registry.get(m)
        if (l is not None):
          l.append(file_name)
        else:
          self.md5_registry[m] = [file_name]
          self.unique_file_names.append(file_name)
          self.raw_records.extend(r)
          self.raw_text_block_list.append(text)

  def report_non_unique(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    n_ignored = 0
    for file_name in sorted(self.file_name_registry.keys()):
      n = self.file_name_registry[file_name]
      if (n != 1):
        print(prefix+"INFO: PDB file name appears %d times: %s" % (
          n, show_string(file_name)), file=out)
        n_ignored += (n-1)
    if (n_ignored != 0):
      print(prefix+"  %d repeated file name%s ignored." % \
        plural_s(n=n_ignored), file=out)
    n_identical = 0
    for file_names in sorted(self.md5_registry.values()):
      if (len(file_names) != 1):
        print(prefix+"INFO: PDB files with identical content:", file=out)
        for file_name in file_names:
          print(prefix+"  %s" % show_string(file_name), file=out)
        n_identical += len(file_names)-1
    if (n_identical != 0):
      print(prefix+"%d file%s with repeated content ignored." % \
        plural_s(n=n_identical), file=out)
    if (n_ignored != 0 or n_identical != 0):
      print(prefix.rstrip(), file=out)

class header_date(object):

  def __init__(self, field):
    "Expected format: DD-MMM-YY"
    self.dd = None
    self.mmm = None
    self.yy = None
    self.yyyy = None
    if (len(field) != 9): return
    if (field.count("-") != 2): return
    if (field[2] != "-" or field[6] != "-"): return
    dd, mmm, yy = field.split("-")
    try: self.dd = int(dd)
    except ValueError: pass
    else:
      if (self.dd < 1 or self.dd > 31): self.dd = None
    if (mmm.upper() in [
          "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
          "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]):
      self.mmm = mmm.upper()
    try: self.yy = int(yy)
    except ValueError: pass
    else:
      if (self.yy < 0 or self.yy > 99): self.yy = None
    if self.yy is not None:
      if self.yy < 60: # I hope by 2060 no one uses the PDB format seriously!
        self.yyyy = 2000 + self.yy
      else:
        self.yyyy = 1900 + self.yy

  def is_fully_defined(self):
    return self.dd is not None \
       and self.mmm is not None \
       and self.yy is not None \
       and self.yyyy is not None

def header_year(record):
  """Return year from header record"""
  if (record.startswith("HEADER")):
    date = header_date(field=record[50:59])
    if (date.is_fully_defined()): return date.yyyy
    fields = record.split()
    fields.reverse()
    for field in fields:
      date = header_date(field=field)
      if (date.is_fully_defined()): return date.yyyy
  return None

class Please_pass_string_or_None(object): pass

class pdb_input_from_any(object):

  def __init__(self,
               file_name=None,
               source_info=Please_pass_string_or_None,
               lines=None,
               raise_sorry_if_format_error=False):
    self.file_format = None
    content = None
    from iotbx.pdb.mmcif import cif_input
    mmcif_exts = ('.cif', '.mmcif')
    if file_name is not None and file_name.strip(".gz").endswith(mmcif_exts):
      file_inputs = (cif_input, pdb_input)
    else:
      file_inputs = (pdb_input, cif_input)
    exc_info = None
    for file_input in file_inputs:
      try:
        content = file_input(
          file_name=file_name,
          source_info=source_info,
          lines=lines,
          raise_sorry_if_format_error=raise_sorry_if_format_error)
      except Exception as e:
        # store the first error encountered and re-raise later if can't
        # interpret as any file type
        if exc_info is None: exc_info = sys.exc_info()
        continue
      else: exc_info = None
      if file_input is pdb_input:
        # XXX nasty hack:
        #   pdb_input only raises an error if there are lines starting with
        #   "ATOM  " or "HETATM" and it subsequently fails to interpret these
        #   lines as ATOM/HETATM records
        #
        # XXX This hack fails to recognize /net/cci/pdb_mirror/mmcif/of/2of6.cif.gz
        # Reason: atom coordinates luckily can be parsed (no letters there)
        #
        n_unknown_records = content.unknown_section().size()
        n_records = sum(content.record_type_counts().values())
        n_blank_records = content.record_type_counts().get('      ', 0)
        # n_records > 0 for empty input, to get empty pdb_input object instead
        # of Sorry from cif_input
        if (((n_unknown_records == n_records) or
              (n_unknown_records == (n_records - n_blank_records)
                  and n_unknown_records > 0))
            and n_records>0):
          continue
        # Additional check that solves 2of6:
        # if the first non-comment non-empty line contains data_ this is mmCIF
        if lines is not None and len(lines)>0:
          len_lines = len(lines)
          i = 0
          while i < len_lines and (
            lines[i].strip().startswith('#') or len(lines[i].strip()) == 0):
            i += 1
          if i < len_lines and lines[i][:5].strip() == 'data_':
            continue
        self.file_format = "pdb"
      else :
        self.file_format = "cif"
      break
    if exc_info is not None:
      six.reraise(exc_info[0], exc_info[1], exc_info[2])
    if content is None:
      raise Sorry("Could not interpret input as any file type.")
    self._file_content = content

  def file_content(self):
    return self._file_content

def pdb_input(
    file_name=None,
    source_info=Please_pass_string_or_None,
    lines=None,
    raise_sorry_if_format_error=False):
  """Read in a model file and return a pdb_input object. Normally use input() instead"""
  if (file_name is not None):
    try :
      with smart_open.for_reading(file_name, gzip_mode='rt') as f:
        lines = f.read()
      return ext.input(
        source_info="file " + str(file_name), # XXX unicode hack - dangerous
        lines=flex.split_lines(lines))
    except ValueError as e :
      if (raise_sorry_if_format_error):
        raise Sorry("Format error in %s:\n%s" % (str(file_name), str(e)))
      else :
        raise
  assert source_info is not Please_pass_string_or_None
  if (isinstance(lines, str)):
    lines = flex.split_lines(lines)
  elif (isinstance(lines, (list, tuple))):
    lines = flex.std_string(lines)
  try :
    return ext.input(source_info=source_info, lines=lines)
  except ValueError as e :
    if (raise_sorry_if_format_error):
      raise Sorry("Format error:\n%s" % str(e))
    else :
      raise

def input(
    file_name=None,
    source_info=Please_pass_string_or_None,
    lines=None,
    raise_sorry_if_format_error=False):
  """
  Main input method for both PDB and mmCIF files; will automatically determine
  the actual format and return the appropriate data type.

  Parameters
  ----------
  file_name: path to PDB or mmCIF file
  source_info: string describing source of input (e.g. file name)
  lines: flex.std_string array of input lines
  raise_sorry_if_format_error: re-raise any low-level parser errors as a
    libtbx.utils.Sorry exception instance for clean user feedback

  Returns
  -------
  An object representing the result of parsing, including an array of atom
  objects; the actual class will differ depending on the input format.  Much of
  the API will be the same in either case.
  """
  return pdb_input_from_any(
    file_name=file_name,
    source_info=source_info,
    lines=lines,
    raise_sorry_if_format_error=raise_sorry_if_format_error).file_content()

default_atom_names_scattering_type_const = ["PEAK", "SITE"]

input_sections = (
  "unknown_section",
  "title_section",
  "remark_section",
  "primary_structure_section",
  "heterogen_section",
  "secondary_structure_section",
  "connectivity_annotation_section",
  "miscellaneous_features_section",
  "crystallographic_section",
  "connectivity_section",
  "bookkeeping_section")


class pdb_input_mixin(object):

  def label_to_auth_asym_id_dictionary(self):
    """ Only avaliable for cif_input
    """
    return None

  def deposition_date(self, us_style=True):
    """
    Placeholder to match mmCIF functionality. Probably could parse
    REVDAT.
    """
    result = None
    for line in self.title_section():
      if(line.startswith("HEADER ")):
        date = header_date(field=line[50:59])
        if(date.is_fully_defined()):
          dd = str(date.dd).strip()
          if(len(dd)==1): dd = "0"+dd
          result = "%s-%s-%s"%(dd, str(date.mmm), str(date.yyyy))
          if(us_style):
            months = dict((v.upper(),k) for k,v in enumerate(calendar.month_abbr))
            m=str(months[str(date.mmm).upper()])
            if(len(m)==1): m = "0"+m
            result = "%s-%s-%s"%(str(date.yyyy), m, dd)
    return result

  def special_position_settings(self,
        special_position_settings=None,
        weak_symmetry=False,
        min_distance_sym_equiv=0.5,
        u_star_tolerance=0):
    crystal_symmetry = self.crystal_symmetry(
      crystal_symmetry=special_position_settings,
      weak_symmetry=weak_symmetry)
    return construct_special_position_settings(
      crystal_symmetry = crystal_symmetry,
      special_position_settings=special_position_settings,
      weak_symmetry=weak_symmetry,
      min_distance_sym_equiv=min_distance_sym_equiv,
      u_star_tolerance=u_star_tolerance)

  def as_pdb_string(self,
        crystal_symmetry=Auto,
        cryst1_z=Auto,
        write_scale_records=True,
        append_end=False,
        atom_hetatm=True,
        sigatm=True,
        anisou=True,
        siguij=True,
        cstringio=None,
        link_records=Auto,
        return_cstringio=Auto):
    """
    Generate standard PDB format.  Will use built-in crystal symmetry if
    available.
    """
    if (cstringio is None):
      cstringio = StringIO()
      if (return_cstringio is Auto):
        return_cstringio = False
    elif (return_cstringio is Auto):
      return_cstringio = True
    if 0:
      if (link_records is Auto):
        print(format_link_records(self.get_link_records()), file=cstringio)
      elif (link_records is not None):
        print(format_link_records(link_records), file=cstringio)
    if (crystal_symmetry is Auto):
      crystal_symmetry = self.crystal_symmetry()
    if (cryst1_z is Auto):
      cryst1_z = self.extract_cryst1_z_columns()
    if (crystal_symmetry is not None or cryst1_z is not None):
      print(format_cryst1_and_scale_records(
        crystal_symmetry=crystal_symmetry,
        cryst1_z=cryst1_z,
        write_scale_records=write_scale_records), file=cstringio)

    py3out = self._as_pdb_string_cstringio(
      cstringio=cstringio,
      append_end=append_end,
      atom_hetatm=atom_hetatm,
      sigatm=sigatm,
      anisou=anisou,
      siguij=siguij)
    if six.PY3:
      cstringio.write( py3out)
    if (return_cstringio):
      return cstringio
    return cstringio.getvalue()

  def write_pdb_file(self,
        file_name,
        open_append=False,
        crystal_symmetry=Auto,
        cryst1_z=Auto,
        write_scale_records=True,
        append_end=False,
        atom_hetatm=True,
        sigatm=True,
        anisou=True,
        siguij=True):
    if (crystal_symmetry is Auto):
      crystal_symmetry = self.crystal_symmetry()
    if (cryst1_z is Auto):
      cryst1_z = self.extract_cryst1_z_columns()
    if (crystal_symmetry is not None or cryst1_z is not None):
      if (open_append): mode = "a"
      else:             mode = "w"
      with open(file_name, mode) as f:
        print(format_cryst1_and_scale_records(
          crystal_symmetry=crystal_symmetry,
          cryst1_z=cryst1_z,
          write_scale_records=write_scale_records), file=f)
      open_append = True
    self._write_pdb_file(
      file_name=file_name,
      open_append=open_append,
      append_end=append_end,
      atom_hetatm=atom_hetatm,
      sigatm=sigatm,
      anisou=anisou,
      siguij=siguij)

  def xray_structure_simple(self,
        crystal_symmetry=None,
        weak_symmetry=False,
        cryst1_substitution_buffer_layer=None,
        unit_cube_pseudo_crystal=False,
        fractional_coordinates=False,
        use_scale_matrix_if_available=True,
        min_distance_sym_equiv=0.5,
        non_unit_occupancy_implies_min_distance_sym_equiv_zero=True,
        scattering_type_exact=False,
        enable_scattering_type_unknown=False,
        atom_names_scattering_type_const
          =default_atom_names_scattering_type_const):
    """
    Create a single cctbx.xray.structure object from the atom records, using
    only the first model found.
    """
    if(crystal_symmetry is not None): self._scale_matrix = None
    return self.xray_structures_simple(
      one_structure_for_each_model=False,
      crystal_symmetry=crystal_symmetry,
      weak_symmetry=weak_symmetry,
      cryst1_substitution_buffer_layer=cryst1_substitution_buffer_layer,
      unit_cube_pseudo_crystal=unit_cube_pseudo_crystal,
      fractional_coordinates=fractional_coordinates,
      use_scale_matrix_if_available=use_scale_matrix_if_available,
      min_distance_sym_equiv=min_distance_sym_equiv,
      non_unit_occupancy_implies_min_distance_sym_equiv_zero=
        non_unit_occupancy_implies_min_distance_sym_equiv_zero,
      scattering_type_exact=scattering_type_exact,
      enable_scattering_type_unknown=enable_scattering_type_unknown,
      atom_names_scattering_type_const=atom_names_scattering_type_const)[0]

  def xray_structures_simple(self,
        one_structure_for_each_model=True,
        crystal_symmetry=None,
        weak_symmetry=False,
        cryst1_substitution_buffer_layer=None,
        unit_cube_pseudo_crystal=False,
        fractional_coordinates=False,
        min_distance_sym_equiv=0.5,
        non_unit_occupancy_implies_min_distance_sym_equiv_zero=True,
        use_scale_matrix_if_available=True,
        scattering_type_exact=False,
        enable_scattering_type_unknown=False,
        atom_names_scattering_type_const
          =default_atom_names_scattering_type_const):
    """
    Create a list of cctbx.xray.structure objects, one per model in the
    input file.  Note that for most single-model structures (i.e. nearly all
    crystal structures), this will be a single-item list.
    """
    from cctbx import xray
    from cctbx import crystal
    from cctbx import uctbx
    if (unit_cube_pseudo_crystal):
      assert crystal_symmetry is None and cryst1_substitution_buffer_layer is None
      crystal_symmetry = crystal.symmetry(
        unit_cell=(1,1,1,90,90,90),
        space_group_symbol="P1")
    else:
      crystal_symmetry = self.crystal_symmetry(
        crystal_symmetry=crystal_symmetry,
        weak_symmetry=weak_symmetry)
      if (crystal_symmetry is None):
        crystal_symmetry = crystal.symmetry()
      if (crystal_symmetry.unit_cell() is None):
        crystal_symmetry = crystal_symmetry.customized_copy(
          unit_cell=uctbx.non_crystallographic_unit_cell(
            sites_cart=self.atoms().extract_xyz(),
            buffer_layer=cryst1_substitution_buffer_layer))
      if (crystal_symmetry.space_group_info() is None):
        crystal_symmetry = crystal_symmetry.cell_equivalent_p1()
    unit_cell = crystal_symmetry.unit_cell()
    scale_r = (0,0,0,0,0,0,0,0,0)
    scale_t = (0,0,0)
    if (not unit_cube_pseudo_crystal):
      if (use_scale_matrix_if_available):
        scale_matrix = self.scale_matrix()
        if (scale_matrix is not None):
          # Avoid subtle inconsistencies due to rounding errors.
          # 1.e-6 is the precision of the values on the SCALE records.
          if (max([abs(s-f) for s,f in zip(
                     scale_matrix[0],
                     unit_cell.fractionalization_matrix())]) < 1.e-6):
            if (scale_matrix[1] != [0,0,0]):
              scale_matrix[0] = unit_cell.fractionalization_matrix()
            else:
              scale_matrix = None
      else:
        scale_matrix = None
      if (scale_matrix is not None):
        scale_r = scale_matrix[0]
        scale_t = scale_matrix[1]
    result = []
    if (atom_names_scattering_type_const is None):
      atom_names_scattering_type_const = []
    loop = xray_structures_simple_extension(
      one_structure_for_each_model,
      unit_cube_pseudo_crystal,
      fractional_coordinates,
      scattering_type_exact,
      enable_scattering_type_unknown,
      self.atoms_with_labels(),
      self.model_indices(),
      scitbx.stl.set.stl_string(atom_names_scattering_type_const),
      unit_cell,
      scale_r,
      scale_t)
    special_position_settings = crystal_symmetry.special_position_settings(
      min_distance_sym_equiv=min_distance_sym_equiv)
    try :
      while (next(loop)):
        result.append(xray.structure(
          special_position_settings=special_position_settings,
          scatterers=loop.scatterers,
          non_unit_occupancy_implies_min_distance_sym_equiv_zero=
            non_unit_occupancy_implies_min_distance_sym_equiv_zero))
    except ValueError as e :
      raise Sorry(str(e))
    return result

bp.inject(ext.input, pdb_input_mixin)
@bp.inject_into(ext.input)
class _():

  """
  This class parses PDB format, including non-ATOM records.  Atom objects will
  be created as part of the parsing, but the full PDB hierarchy object requires
  calling the construct_hierarchy() method.
  """

  def __getinitargs__(self):
    lines = flex.std_string()
    for section in input_sections[:-2]:
      lines.extend(getattr(self, section)())
    pdb_string = StringIO()

    py3out = self._as_pdb_string_cstringio(  # NOTE py3out is None in python 2
      cstringio=pdb_string,
      append_end=False,
      atom_hetatm=True,
      sigatm=True,
      anisou=True,
      siguij=True)
    if six.PY3:
      pdb_string.write(py3out)
    lines.extend(flex.split_lines(pdb_string.getvalue()))
    for section in input_sections[-2:]:
      lines.extend(getattr(self, section)())
    return ("pickle", lines)

  def file_type(self):
    return "pdb"

  def sequence_from_SEQRES(self):
    from iotbx.pdb import amino_acid_codes
    d = {}
    ps = self.primary_structure_section()
    for l in ps:
      l = l.strip()
      ls = l.split()
      if(l.startswith("SEQRES")):
        kw, i_seq, chid, rns = ls[0], ls[1], ls[2], ls[4:]
        d.setdefault(chid, []).extend(rns)
    result = []
    ott = amino_acid_codes.one_letter_given_three_letter
    for k, vs in zip(d.keys(), d.values()): # FIXME use iteritems?
      result.append(">chain %s"%k)
      result.append("".join([ott.get(v,"?") for v in vs]))
    return "\n".join(result)

  def extract_header_year(self):
    for line in self.title_section():
      if (line.startswith("HEADER ")):
        return header_year(line)
    return None

  def extract_authors(self):
    trigger = "AUTHOR"
    result = []
    def is_number(s):
      try:
        float(s)
        return True
      except ValueError: return False
    # extract and put into one string...
    lt = ""
    cntr=0
    for l in self.title_section():
      if(l.startswith(trigger)):
        l_=l.strip().replace(trigger,"").strip()
        if(is_number(l_[0])): l_ = l_.replace(l_[0],"").strip()
        #
        l1 = lt.split(",")
        l1 = l1[len(l1)-1]
        l2 = l_.split(",")[0]
        if((l1.count(".")>0 and l2.count(".")==0) or l1.endswith("-")): j=""
        else: j = ","
        #
        lt = lt + j + l_
        cntr+=1
    # ...then analyze
    l_=lt.strip().replace(trigger,"").strip()
    if(is_number(l_[0])): l_ = l_.replace(l_[0],"").strip()
    if(not (l_.startswith(trigger) or is_number(l_))):
      l_ = l_.split(",")
      for l__ in l_:
        l__ = "".join([x.strip() for x in l__])
        if(len(l__)>0):
          l__ = l__.split(".")
          l__.sort()
          l__ = [x.upper() for x in l__]
          l__ = ".".join(l__)
          if(l__[0].isalpha() and l__[len(l__)-1].isalpha()):
            result.append(l__)
    return result

  def extract_remark_iii_records(self, iii):
    result = []
    pattern = "REMARK %3d " % iii
    for line in self.remark_section():
      if (line.startswith(pattern)):
        result.append(line)
    return result

  def extract_secondary_structure(self, log=None):
    from iotbx.pdb import secondary_structure
    records = self.secondary_structure_section()
    return secondary_structure.annotation.from_records(records, log)

  def extract_LINK_records(self):
    '''
    Collect link records from PDB file
    '''
    result = []
    for line in self.connectivity_annotation_section():
      if (line.startswith('LINK') or line.startswith('link')):
        result.append(line)
    return result

  def crystal_symmetry_from_cryst1(self):
    from iotbx.pdb import cryst1_interpretation
    for line in self.crystallographic_section():
      if (line.startswith("CRYST1")):
        return cryst1_interpretation.crystal_symmetry(cryst1_record=line)
    return None

  def extract_cryst1_z_columns(self):
    for line in self.crystallographic_section():
      if (line.startswith("CRYST1")):
        result = line[66:]
        if (len(result) < 4): result += " " * (4-len(result))
        return result
    return None

  def _crystal_symmetry_from_cns_remark_sg(self):
    from iotbx.cns import pdb_remarks
    for line in self.remark_section():
      if (line.startswith("REMARK sg=")):
        crystal_symmetry = pdb_remarks.extract_symmetry(pdb_record=line)
        if (crystal_symmetry is not None):
          return crystal_symmetry
    return None

  def crystal_symmetry(self,
        crystal_symmetry=None,
        weak_symmetry=False):
    self_symmetry = self.crystal_symmetry_from_cryst1()
    if (self_symmetry is None):
      self_symmetry = self._crystal_symmetry_from_cns_remark_sg()
    if (crystal_symmetry is None):
      return self_symmetry
    if (self_symmetry is None):
      return crystal_symmetry
    return self_symmetry.join_symmetry(
      other_symmetry=crystal_symmetry,
      force=not weak_symmetry)

  def scale_matrix(self):
    if (not hasattr(self, "_scale_matrix")):
      source_info = self.source_info()
      if (len(source_info) > 0): source_info = " (%s)" % source_info
      self._scale_matrix = [[None]*9,[None]*3]
      done_set = set()
      done_list = []
      for line in self.crystallographic_section():
        if (line.startswith("SCALE") and line[5:6] in ["1", "2", "3"]):
          r = read_scale_record(line=line, source_info=source_info)
          if (r.n not in done_set):
            for i_col,v in enumerate(r.r):
              self._scale_matrix[0][(r.n-1)*3+i_col] = v
            self._scale_matrix[1][r.n-1] = r.t
            done_set.add(r.n)
          done_list.append(r.n)
      if (len(done_list) == 0):
        self._scale_matrix = None
      elif (sorted(done_list[:3]) != [1,2,3]):
        raise ValueError(
          "Improper set of PDB SCALE records%s" % source_info)
    return self._scale_matrix

  def process_BIOMT_records(self):
    import iotbx.mtrix_biomt
    return iotbx.mtrix_biomt.process_BIOMT_records_pdb(
      lines = self.extract_remark_iii_records(350))

  def process_MTRIX_records(self):
    import iotbx.mtrix_biomt
    return iotbx.mtrix_biomt.process_MTRIX_records_pdb(
      lines=self.crystallographic_section())

  def get_r_rfree_sigma(self, file_name=None):
    from iotbx.pdb import extract_rfactors_resolutions_sigma
    remark_2_and_3_records = self.extract_remark_iii_records(2)
    remark_2_and_3_records.extend(self.extract_remark_iii_records(3))
    return extract_rfactors_resolutions_sigma.get_r_rfree_sigma(
      remark_2_and_3_records, file_name)

  def resolution(self):
    return self.get_r_rfree_sigma().resolution

  def get_program_name(self):
    remark_3_lines = self.extract_remark_iii_records(3)
    result = None
    for line in remark_3_lines:
      line = line.strip()
      result = iotbx.pdb.remark_3_interpretation.get_program(st = line)
      if(result is not None): return result
    if(result is not None):
      result = "_".join(result.split())
    return result

  def get_solvent_content(self):
    remark_280_lines = self.extract_remark_iii_records(280)
    mc = []
    for remark in remark_280_lines:
      remark = remark.upper()
      if(remark.count("SOLVENT")==1 and
         remark.count("CONTENT")==1):
        try:
          mc.append(remark.split()[6])
        except Exception:
          try:
            mc.append(remark[remark.index(":")+1:])
          except Exception:
            mc.append(remark)
    result = None
    if(len(mc) == 1):
      try: result = float(mc[0])
      except IndexError: pass
      except ValueError: pass
    return result

  def get_matthews_coeff(self):
    remark_280_lines = self.extract_remark_iii_records(280)
    mc = []
    for remark in remark_280_lines:
      remark = remark.upper()
      if(remark.count("MATTHEWS")==1 and
         remark.count("COEFFICIENT")==1):
        try:
          mc.append(remark.split()[6])
        except Exception:
          try:
            mc.append(remark[remark.index(":")+1:])
          except Exception:
            mc.append(remark)
    result = None
    if(len(mc) == 1):
      try: result = float(mc[0])
      except IndexError: pass
      except ValueError: pass
    return result

  def extract_tls_params(self, hierarchy):
    import iotbx.pdb.remark_3_interpretation
    remark_3_records = self.extract_remark_iii_records(3)
    chain_ids = []
    for model in hierarchy.models():
      for chain in model.chains():
        chain_ids.append(chain.id)
    return iotbx.pdb.remark_3_interpretation.extract_tls_parameters(
      remark_3_records = remark_3_records,
      pdb_hierarchy    = hierarchy,
      chain_ids        = chain_ids)

  def extract_f_model_core_constants(self):

    remark_3_records = self.extract_remark_iii_records(3)
    return remark_3_interpretation.extract_f_model_core_constants(remark_3_records)

  def extract_wavelength(self, first_only=True):
    for line in self.remark_section():
      # XXX this will miss multi-line records!
      if (line.startswith("REMARK 200  WAVELENGTH OR RANGE")):
        fields = line.split(":")
        assert (len(fields) == 2)
        subfields = fields[1].replace(";", ",").strip().split(",")
        wavelengths = []
        for field in subfields :
          if (field.strip() == ""):
            continue
          elif (field.strip() == "NULL"):
            wavelengths.append(None)
          else :
            try :
              wavelengths.append(float(field.strip()))
            except ValueError :
              wavelengths.append(None)
        if (first_only):
          if (len(wavelengths) > 0):
            return wavelengths[0]
          return None
        return wavelengths
    return None

  def get_experiment_type(self):
    for line in self.title_section():
      if (line.startswith("EXPDTA")):
        return experiment_type(iotbx.pdb.records.expdta(line).keywords)
    return experiment_type([])

  def extract_connectivity(self):
    """
    Parse CONECT records and extract the indices of bonded atoms.  Returns
    a scitbx.array_family.shared.stl_set_unsigned object corresponding to the
    atoms array, with each element being the list of indices of bonded atoms
    (if any).  If no CONECT records are found, returns None.

    Note that the ordering of atoms may be altered by construct_hierarchy(), so
    this method should probably be called after the hierarchy is made.
    """
    lines = self.connectivity_section()
    if (len(lines) == 0):
      return None
    from scitbx.array_family import shared
    bonds = shared.stl_set_unsigned(len(self.atoms()), [])
    serial_ref_hash = {}
    for i_seq, atom in enumerate(self.atoms()):
      serial = atom.serial.strip()
      serial_ref_hash[serial] = i_seq
    for line in lines :
      assert (line.startswith("CONECT"))
      record = iotbx.pdb.records.conect(line)
      i_seq = serial_ref_hash[record.serial.strip()]
      assert (record.serial_numbers_bonded_atoms.count('') != 4)
      for j_seq_str in record.serial_numbers_bonded_atoms :
        if (j_seq_str != ''):
          bonds[i_seq].append(serial_ref_hash[j_seq_str.strip()])
    return bonds

  def get_link_records(self):
    for atom1 in self.atoms_with_labels():
      if not atom1.hetero: continue
      for atom2 in self.atoms_with_labels():
        if not atom2.hetero: continue
        if atom1.resname==atom2.resname: continue
        yield atom1, atom2, "1555", "1555"
        break
      break

  def get_restraints_used(self):
    return {'CDL' : self.used_cdl_restraints(),
            'omega' : self.used_omega_restraints(),
            'Amber' : self.used_amber_restraints(),
           }

  def _used_what_restraints(self, what):
    rc = False
    for line in self.remark_section():
      if line.startswith("REMARK   3") and (what in line):
        rc = True
        break
    return rc

  def used_cdl_restraints(self):
    return self._used_what_restraints('CDL')

  def used_omega_cdl_restraints(self):
    return self._used_what_restraints('omega-cdl')

  def used_amber_restraints(self):
    return self._used_what_restraints('Amber')

# Table of structures split into multiple PDB files.
# Assembled manually.
# Based on 46377 PDB files as of Tuesday Oct 02, 2007
#   noticed in passing: misleading REMARK 400 in 1VSA (1vs9 and 2i1c
#   don't exist)
# Updated 2009-04-07, based on 56751 PDB files, using SPLIT records.
pdb_codes_fragment_files = """\
1crp 1crr
1f49 1gho
1gix 1giy
1j4z 1kpo
1jgo 1jgp 1jgq
1jyy 1jyz
1jz0 1jz1
1otz 1p0t
1pns 1pnu
1pnx 1pny
1s1h 1s1i
1ti2 1vld
1ti4 1vle
1ti6 1vlf
1utf 1utv
1voq 1vor 1vos 1vou 1vov 1vow 1vox 1voy 1voz 1vp0
1vs5 1vs6 1vs7 1vs8
1vsa 2ow8
1vsp 2qnh
1we3 1wf4
1yl3 1yl4
2avy 2aw4 2aw7 2awb
2b64 2b66
2b9m 2b9n
2b9o 2b9p
2bld 2bvi
2gy9 2gya
2gyb 2gyc
2hgi 2hgj
2hgp 2hgq
2hgr 2hgu
2i2p 2i2t 2i2u 2i2v
2j00 2j01 2j02 2j03
2jl5 2jl6 2jl7 2jl8
2qal 2qam 2qan 2qao
2qb9 2qba 2qbb 2qbc
2qbd 2qbe 2qbf 2qbg
2qbh 2qbi 2qbj 2qbk
2qou 2qov 2qow 2qox
2qoy 2qoz 2qp0 2qp1
2uv9 2uva
2uvb 2uvc
2v46 2v47 2v48 2v49
2vhm 2vhn 2vho 2vhp
2z4k 2z4l 2z4m 2z4n
2zkq 2zkr
2zuo 2zv4 2zv5
3bz1 3bz2
3d5a 3d5b 3d5c 3d5d
3df1 3df2 3df3 3df4
3f1e 3f1f 3f1g 3f1h
"""

class join_fragment_files(object):
  def __init__(self, file_names):
    info = flex.std_string()
    info.append("REMARK JOINED FRAGMENT FILES (iotbx.pdb)")
    info.append("REMARK " + date_and_time())
    roots = []
    z = None
    from cctbx import crystal
    self.crystal_symmetry = crystal.symmetry()
    z_warning = 'REMARK ' \
      'Warning: CRYST1 Z field (columns 67-70) is not an integer: "%-4.4s"'
    for file_name in file_names:
      pdb_inp = iotbx.pdb.input(file_name=file_name)
      z_ = pdb_inp.extract_cryst1_z_columns()
      if (z_ is not None):
        z_ = z_.strip()
        if (z_ != "") : z_ = int(z_)
        else : z_ = None
      else: z_ = None
      if z is None:
        z = z_
      else:
        z = max(z, z_)
      cs = pdb_inp.crystal_symmetry()
      info.append("REMARK %s" % show_string(file_name))
      if cs is not None and cs.unit_cell() is not None:
        info.append("REMARK %s" % iotbx.pdb.format_cryst1_record(cs, z=z_))
        self.crystal_symmetry = self.crystal_symmetry.join_symmetry(
          cs, force=True)
      roots.append(pdb_inp.construct_hierarchy(sort_atoms=False))
    if self.crystal_symmetry.unit_cell() is not None:
      info.append(iotbx.pdb.format_cryst1_record(
        crystal_symmetry=self.crystal_symmetry, z=z))
    result = iotbx.pdb.hierarchy.join_roots(roots=roots)
    result.info.extend(info)
    result.reset_i_seq_if_necessary()
    self.joined = result

def merge_files_and_check_for_overlap(file_names, output_file,
    site_clash_cutoff=0.5, log=sys.stdout):
  """ Merge models and write composite model"""
  assert len(file_names) > 0
  merged_records = combine_unique_pdb_files(file_names)
  warnings = StringIO()
  merged_records.report_non_unique(out=warnings)
  merged_hierarchy = join_fragment_files(file_names).joined
  f = open(output_file, "w")
  f.write(merged_hierarchy.as_pdb_string())
  f.close()
  n_clashes = quick_clash_check(output_file,
    site_clash_cutoff=site_clash_cutoff,
    out=log)
  return n_clashes

def quick_clash_check(file_name, site_clash_cutoff=0.5, out=sys.stdout,
    show_outliers=5):
  """Carry out clash check on a model file"""
  pdb_inp = input(file_name=file_name)
  pdb_atoms = pdb_inp.atoms_with_labels()
  xray_structure = pdb_inp.xray_structure_simple(
    cryst1_substitution_buffer_layer=10,
    enable_scattering_type_unknown=True)
  sites_frac = xray_structure.sites_frac()
  unit_cell = xray_structure.unit_cell()
  pair_asu_table = xray_structure.pair_asu_table(
    distance_cutoff=site_clash_cutoff)
  pair_sym_table = pair_asu_table.extract_pair_sym_table()
  atom_pairs = pair_sym_table.simple_edge_list()
  return len(atom_pairs)

standard_rhombohedral_space_group_symbols = [
"R 3 :H",
"R 3 :R",
"R -3 :H",
"R -3 :R",
"R 3 2 :H",
"R 3 2 :R",
"R 3 m :H",
"R 3 m :R",
"R 3 c :H",
"R 3 c :R",
"R -3 m :H",
"R -3 m :R",
"R -3 c :H",
"R -3 c :R"]
if ("set" in __builtins__):
  standard_rhombohedral_space_group_symbols = set(
    standard_rhombohedral_space_group_symbols)

def format_cryst1_sgroup(space_group_info):
  """Format space_group for CRYST1 record"""
  result = space_group_info.type().lookup_symbol()
  if (result in standard_rhombohedral_space_group_symbols):
    result = result[-1] + result[1:-3]
  def compress(s):
    if (len(s) > 11): return s.replace(" ", "")
    return s
  result = compress(result)
  if (len(result) > 11 and not space_group_info.group().is_centric()):
    from iotbx.mtz.extract_from_symmetry_lib import ccp4_symbol
    alt = ccp4_symbol(
      space_group_info=space_group_info,
      lib_name="syminfo.lib",
      require_at_least_one_lib=False)
    if (alt is not None and alt != result.replace(" ", "")):
      result = compress(alt)
  return result

def format_cryst1_record(crystal_symmetry, z=None):
  """Format CRYST1 record from crystal_symmetry"""
  # CRYST1
  #  7 - 15       Real(9.3)      a             a (Angstroms).
  # 16 - 24       Real(9.3)      b             b (Angstroms).
  # 25 - 33       Real(9.3)      c             c (Angstroms).
  # 34 - 40       Real(7.2)      alpha         alpha (degrees).
  # 41 - 47       Real(7.2)      beta          beta (degrees).
  # 48 - 54       Real(7.2)      gamma         gamma (degrees).
  # 56 - 66       LString        sGroup        Space group.
  # 67 - 70       Integer        z             Z value.
  if (z is None): z = ""
  else: z = str(z)
  return ("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11.11s%4.4s" % (
    crystal_symmetry.unit_cell().parameters()
    + (format_cryst1_sgroup(
         space_group_info=crystal_symmetry.space_group_info()),
       z))).rstrip()

def format_scale_records(unit_cell=None,
                         fractionalization_matrix=None,
                         u=[0,0,0]):
  """Format SCALE records from unit_cell"""
  #  1 -  6       Record name    "SCALEn"       n=1, 2, or 3
  # 11 - 20       Real(10.6)     s[n][1]        Sn1
  # 21 - 30       Real(10.6)     s[n][2]        Sn2
  # 31 - 40       Real(10.6)     s[n][3]        Sn3
  # 46 - 55       Real(10.5)     u[n]           Un
  assert [unit_cell, fractionalization_matrix].count(None) == 1
  if (unit_cell is not None):
    f = unit_cell.fractionalization_matrix()
  else:
    assert len(fractionalization_matrix) == 9
    f = fractionalization_matrix
  assert len(u) == 3
  return (("SCALE1    %10.6f%10.6f%10.6f     %10.5f\n"
           "SCALE2    %10.6f%10.6f%10.6f     %10.5f\n"
           "SCALE3    %10.6f%10.6f%10.6f     %10.5f") % (
    f[0], f[1], f[2], u[0],
    f[3], f[4], f[5], u[1],
    f[6], f[7], f[8], u[2])).replace(" -0.000000", "  0.000000")

def format_cryst1_and_scale_records(
      crystal_symmetry=None,
      cryst1_z=None,
      write_scale_records=True,
      scale_fractionalization_matrix=None,
      scale_u=[0,0,0]):
  """Format CRYST1 and SCALE records from crystal_symmetry"""
  from cctbx import crystal
  from cctbx import sgtbx
  from cctbx import uctbx
  if (crystal_symmetry is None):
    unit_cell = None
    space_group_info = None
  elif (isinstance(crystal_symmetry, crystal.symmetry)):
    unit_cell = crystal_symmetry.unit_cell()
    space_group_info = crystal_symmetry.space_group_info()
  elif (isinstance(crystal_symmetry, uctbx.ext.unit_cell)):
    unit_cell = crystal_symmetry
    space_group_info = None
  elif (isinstance(crystal_symmetry, (list, tuple))):
    assert len(crystal_symmetry) == 6 # unit cell parameters
    unit_cell = uctbx.unit_cell(crystal_symmetry)
    space_group_info = None
  else:
    raise ValueError("invalid crystal_symmetry object")
  if (unit_cell is None):
    if (scale_fractionalization_matrix is None):
      unit_cell = uctbx.unit_cell((1,1,1,90,90,90))
    else:
      unit_cell = uctbx.unit_cell(
        orthogonalization_matrix=matrix.sqr(
          scale_fractionalization_matrix).inverse())
  if (space_group_info is None):
    space_group_info = sgtbx.space_group_info(symbol="P 1")
  result = format_cryst1_record(
    crystal_symmetry=crystal.symmetry(
      unit_cell=unit_cell, space_group_info=space_group_info),
    z=cryst1_z)
  if (write_scale_records):
    if (scale_fractionalization_matrix is None):
      scale_fractionalization_matrix = unit_cell.fractionalization_matrix()
    result += "\n" + format_scale_records(
      fractionalization_matrix=scale_fractionalization_matrix,
      u=scale_u)
  return result

def format_link_records(link_list):
  """Format LINK records"""

  """
COLUMNS         DATA TYPE      FIELD           DEFINITION
-----------------------------------------------------------------------------
 1 -  6         Record name    "LINK  "
13 - 16         Atom           name1           Atom name.
17              Character      altLoc1         Alternate location indicator.
18 - 20         Residue name   resName1        Residue  name.
22              Character      chainID1        Chain identifier.
23 - 26         Integer        resSeq1         Residue sequence number.
27              AChar          iCode1          Insertion code.
43 - 46         Atom           name2           Atom name.
47              Character      altLoc2         Alternate location indicator.
48 - 50         Residue name   resName2        Residue name.
52              Character      chainID2        Chain identifier.
53 - 56         Integer        resSeq2         Residue sequence number.
57              AChar          iCode2          Insertion code.
60 - 65         SymOP          sym1            Symmetry operator atom 1.
67 - 72         SymOP          sym2            Symmetry operator atom 2.
74 - 78         Real(5.2)      Length          Link distance
"""
  test = """
         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
LINK         O   GLY A  49                NA    NA A6001     1555   1555  2.98
"""
  def _format_link_atom(atom):
    result = "%4s%s%-3s %s%4s%s" % (atom.name,
                                   atom.altloc,
                                   atom.resname,
                                   atom.chain_id,
                                   atom.resseq,
                                   atom.icode,
                                   )
    return result
  result = ""
  for atom1, atom2, sym_op1, sym_op2 in link_list:
    result += "LINK        "
    result += _format_link_atom(atom1)
    result += " "*15
    result += _format_link_atom(atom2)
    result += " "
    result += "   1555"
    result += "   1555"
    result += " %5.2f" % 2.9
  return result

class read_scale_record(object):

  __slots__ = ["n", "r", "t"]

  def __init__(O, line, source_info=""):
    try: O.n = int(line[5:6])
    except ValueError: O.n = None
    if (O.n not in [1,2,3]):
      raise RuntimeError(
        "Unknown PDB record %s%s" % (show_string(line[:6]), source_info))
    values = []
    for i in [10,20,30,45]:
      fld = line[i:i+10]
      if (len(fld.strip()) == 0):
        value = 0
      else:
        try: value = float(fld)
        except ValueError:
          raise RuntimeError(
            "Not a floating-point value, PDB record %s%s:\n" % (
              show_string(line[:6]), source_info)
            + "  " + line + "\n"
            + "  %s%s" % (" "*i, "^"*10))
      values.append(value)
    O.r, O.t = values[:3], values[3]

def resseq_decode(s):
  """Convert from hybrid-36 to integer number """
  try: return hy36decode(width=4, s="%4s" % s)
  except ValueError:
    raise ValueError('invalid residue sequence number: "%4s"' % s)

def resseq_encode(value):
  """Convert from integer number to hybrid-36"""
  return hy36encode(width=4, value=value)

def encode_serial_number(width, value):
  """Convert from serial number to hybrid-36"""

  if (isinstance(value, str)):
    assert len(value) <= width
    return value
  if (isinstance(value, int)):
    return hy36encode(width=width, value=value)
  raise RuntimeError("serial number value must be str or int.")

def make_atom_with_labels(
      result=None,
      xyz=None,
      sigxyz=None,
      occ=None,
      sigocc=None,
      b=None,
      sigb=None,
      uij=None,
      siguij=None,
      hetero=None,
      serial=None,
      name=None,
      segid=None,
      element=None,
      charge=None,
      model_id=None,
      chain_id=None,
      resseq=None,
      icode=None,
      altloc=None,
      resname=None):
  """Make an atom_with_labels-like object """
  if (result is None):
    result = hierarchy.atom_with_labels()
  else :
    assert type(result).__name__ == 'atom_with_labels'
  if (xyz is not None): result.xyz = xyz
  if (sigxyz is not None): result.sigxyz = sigxyz
  if (occ is not None): result.occ = occ
  if (sigocc is not None): result.sigocc = sigocc
  if (b is not None): result.b = b
  if (sigb is not None): result.sigb = sigb
  if (uij is not None): result.uij = uij
  if (siguij is not None): result.siguij = siguij
  if (hetero is not None): result.hetero = hetero
  if (serial is not None): result.serial = serial
  if (name is not None): result.name = name
  if (segid is not None): result.segid = segid
  if (element is not None): result.element = element
  if (charge is not None): result.charge = charge
  if (model_id is not None): result.model_id = model_id
  if (chain_id is not None): result.chain_id = chain_id
  if (resseq is not None): result.resseq = resseq
  if (icode is not None): result.icode = icode
  if (altloc is not None): result.altloc = altloc
  if (resname is not None): result.resname = resname
  return result

def get_file_summary(pdb_in, hierarchy=None):
  """Summarize model file"""
  if (hierarchy is None):
    hierarchy = pdb_in.construct_hierarchy()
  counts = hierarchy.overall_counts()
  chain_ids = []
  for id in sorted(counts.chain_ids.keys()):
    if (id == " "):
      chain_ids.append("' '")
    else :
      chain_ids.append(id)
  info_list = [
    ("Number of atoms", counts.n_atoms),
    ("Number of chains", counts.n_chains),
    ("Chain IDs", ", ".join(chain_ids)),
    ("Alternate conformations", counts.n_alt_conf),
  ]
  if (counts.n_models > 1):
    info_list.insert(0, ("Number of models", counts.n_models))
  cl = counts.resname_classes
  if ("common_amino_acid" in cl):
    info_list.append(("Amino acid residues", cl['common_amino_acid']))
  if ("common_nucleic_acid" in cl) or ("ccp4_mon_lib_rna_dna" in cl):
    n_atoms = cl.get("common_nucleic_acid", 0) + \
              cl.get("ccp4_mon_lib_rna_dna", 0)
    info_list.append(("Nucleic acid residues", n_atoms))
  if ("common_water" in cl):
    info_list.append(("Water molecules", cl['common_water']))
  if ("common_element" in cl):
    names = []
    for name in counts.resnames :
      if (iotbx.pdb.common_residue_names_get_class(name)=="common_element"):
        names.append(name)
    value = "%d (%s)" % (cl['common_element'], ", ".join(names))
    info_list.append(("Elemental ions", value))
  if ("common_small_molecule" in cl) or ("other" in cl):
    names = []
    for name in counts.resnames :
      res_class = iotbx.pdb.common_residue_names_get_class(name)
      if (res_class in ["common_small_molecule", "other"]):
        names.append(name)
    n_atoms = cl.get("common_small_molecule", 0) + cl.get("other", 0)
    value = "%d (%s)" % (n_atoms, ", ".join(names))
    info_list.append(("Other molecules", value))
  atoms = hierarchy.atoms()
  b_factors = atoms.extract_b()
  mean_b = flex.mean(b_factors)
  min_b = flex.min(b_factors)
  max_b = flex.max(b_factors)
  info_list.append(("Mean isotropic B-factor", "%.2f (range: %.2f - %.2f)" %
    (mean_b, min_b, max_b)))
  if (min_b <= 0):
    n_bad_adp = (b_factors <= 0).count(True)
    info_list.append(("Atoms with iso. B <= 0", "%d ***" % n_bad_adp))
  occ = atoms.extract_occ()
  if (flex.min(occ) <= 0):
    n_zero_occ = (occ <= 0).count(True)
    info_list.append(("Atoms with zero ocupancy", "%d ***" % n_zero_occ))
  symm = pdb_in.crystal_symmetry()
  if (symm is not None):
    space_group = symm.space_group_info()
    info_list.append(("Space group", str(space_group)))
    unit_cell = symm.unit_cell()
    if (unit_cell is not None):
      uc_str = " ".join([ "%g" % x for x in unit_cell.parameters() ])
      info_list.append(("Unit cell", uc_str))
  return info_list

def show_file_summary(pdb_in, hierarchy=None, out=None):
  """Display summary of model file"""
  if (out is None):
    out = sys.stdout
  info = get_file_summary(pdb_in, hierarchy)
  label_width = max([ len(l) for l,v in info ]) + 2
  format = "%%-%ds %%s" % label_width
  for label, value in info :
    print(format % (label + ":", str(value)), file=out)
  return info
