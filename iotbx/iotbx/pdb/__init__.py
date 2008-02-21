from cctbx.array_family import flex

import boost.python
ext = boost.python.import_ext("iotbx_pdb_ext")
from iotbx_pdb_ext import *

from iotbx.pdb.atom_name_interpretation import \
  interpreters as protein_atom_name_interpreters
import scitbx.array_family.shared
import scitbx.stl.set
from libtbx import smart_open
from libtbx.math_utils import iround
from libtbx.str_utils import show_string, show_sorted_by_counts
from libtbx.utils import plural_s, Sorry
from libtbx.itertbx import count
from libtbx import dict_with_default_0
from cStringIO import StringIO
import md5
import sys

cns_dna_rna_residue_names = {
  "ADE": "A",
  "CYT": "C",
  "GUA": "G",
  "THY": "T",
  "URI": "U"
}

mon_lib_dna_rna_cif = ["AD", "AR", "CD", "CR", "GD", "GR", "TD", "UR"]
if ("set" in __builtins__):
  mon_lib_dna_rna_cif = set(mon_lib_dna_rna_cif)

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
  "AR": "A",
  "CR": "C",
  "GR": "G",
  "UR": "U",
  "AD": "DA",
  "CD": "DC",
  "GD": "DG",
  "TD": "DT"
}

def rna_dna_reference_residue_name(common_name):
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
  "HO5'": "HO5*"
# "HOP3": not covered by monomer library
}

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
          assert rn == "HOP3" # only atom not covered by monomer library
          result.append(None)
    return result

class residue_name_plus_atom_names_interpreter(object):

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
    protein_interpreter = protein_atom_name_interpreters.get(work_residue_name)
    atom_name_interpretation = None
    if (protein_interpreter is not None):
      atom_name_interpretation = protein_interpreter.match_atom_names(
        atom_names=atom_names)
    else:
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
              "A": "AR",
              "C": "CR",
              "G": "GR",
              "U": "UR",
              "DA": "AD",
              "DC": "CD",
              "DG": "GD",
              "DT": "TD"}[work_residue_name]
    self.work_residue_name = work_residue_name
    self.atom_name_interpretation = atom_name_interpretation

class combine_unique_pdb_files(object):

  def __init__(self, file_names):
    self.file_name_registry = {}
    self.md5_registry = {}
    self.unique_file_names = []
    self.raw_records = []
    for file_name in file_names:
      if (file_name in self.file_name_registry):
        self.file_name_registry[file_name] += 1
      else:
        self.file_name_registry[file_name] = 1
        r = [s.expandtabs().rstrip()
          for s in smart_open.for_reading(
            file_name=file_name).read().splitlines()]
        m = md5.new()
        m.update("\n".join(r))
        m = m.hexdigest()
        l = self.md5_registry.get(m)
        if (l is not None):
          l.append(file_name)
        else:
          self.md5_registry[m] = [file_name]
          self.unique_file_names.append(file_name)
          self.raw_records.extend(r)

  def report_non_unique(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    n_ignored = 0
    for file_name in sorted(self.file_name_registry.keys()):
      n = self.file_name_registry[file_name]
      if (n != 1):
        print >> out, prefix+"INFO: PDB file name appears %d times: %s" % (
          n, show_string(file_name))
        n_ignored += (n-1)
    if (n_ignored != 0):
      print >> out, prefix+"  %d repeated file name%s ignored." % \
        plural_s(n=n_ignored)
    n_identical = 0
    for file_names in self.md5_registry.values():
      if (len(file_names) != 1):
        print >> out, prefix+"INFO: PDB files with identical content:"
        for file_name in file_names:
          print >> out, prefix+"  %s" % show_string(file_name)
        n_identical += len(file_names)-1
    if (n_identical != 0):
      print >> out, prefix+"%d file%s with repeated content ignored." % \
        plural_s(n=n_identical)
    if (n_ignored != 0 or n_identical != 0):
      print >> out, prefix.rstrip()

class Please_pass_string_or_None(object): pass

def input(
    file_name=None,
    source_info=Please_pass_string_or_None,
    lines=None):
  if (file_name is not None):
    return ext.input(
      source_info="file " + file_name,
      lines=flex.split_lines(smart_open.for_reading(file_name).read()))
  assert source_info is not Please_pass_string_or_None
  return ext.input(source_info=source_info, lines=lines)

def show_summary(
      file_name=None,
      source_info=Please_pass_string_or_None,
      lines=None,
      level_id=None,
      level_id_exception=ValueError,
      duplicate_max_show=10,
      out=None,
      prefix=""):
  pdb_inp = input(file_name=file_name, source_info=source_info, lines=lines)
  print >> out, prefix+"source info:", pdb_inp.source_info()
  hierarchy = pdb_inp.construct_hierarchy()
  overall_counts = hierarchy.overall_counts()
  assert overall_counts.n_atoms == pdb_inp.input_atom_labels_list().size()
  fmt = "%%%dd" % len(str(overall_counts.n_atoms))
  print >> out, prefix+"  total number of:"
  print >> out, prefix+"    models:    ", fmt % overall_counts.n_models
  print >> out, prefix+"    chains:    ", fmt % overall_counts.n_chains
  print >> out, prefix+"    alt. conf.:", fmt % (
    overall_counts.n_conformers - overall_counts.n_chains)
  print >> out, prefix+"    residues:  ", fmt % overall_counts.n_residues
  print >> out, prefix+"    atoms:     ", fmt % overall_counts.n_atoms
  #
  c = pdb_inp.atom_element_counts()
  print >> out, prefix+"  number of atom element types: %d"%len(c)
  print >> out, prefix+"  histogram of atom element frequency:"
  show_sorted_by_counts(c.items(), out=out, prefix=prefix+"    ")
  #
  c = overall_counts.residue_name_classes
  print >> out, prefix+"  residue name classes:"
  show_sorted_by_counts(c.items(), out=out, prefix=prefix+"    ")
  #
  c = overall_counts.chain_ids
  print >> out, prefix+"  number of chain ids: %d"% len(c)
  print >> out, prefix+"  histogram of chain id frequency:"
  show_sorted_by_counts(c.items(), out=out, prefix=prefix+"    ")
  #
  c = overall_counts.conformer_ids
  print >> out, prefix+"  number of conformer ids: %d"%len(c)
  print >> out, prefix+"  histogram of conformer id frequency:"
  show_sorted_by_counts(c.items(), out=out, prefix=prefix+"    ")
  #
  c = overall_counts.residue_names
  print >> out, prefix+"  number of residue names: %d"%len(c)
  print >> out, prefix+"  histogram of residue name frequency:"
  annotation_appearance = {
    "common_amino_acid": None,
    "common_rna_dna": None,
    "common_water": "   common water",
    "common_small_molecule": "   common small molecule",
    "common_element": "   common element",
    "other": "   other"
  }
  show_sorted_by_counts(c.items(), out=out, prefix=prefix+"    ",
    annotations=[
      annotation_appearance[common_residue_names_get_class(name=name)]
        for name in c.keys()])
  #
  msg = pdb_inp.have_duplicate_atom_labels_message(
    max_show=duplicate_max_show,
    prefix=prefix+"  ")
  if (msg is not None): print >> out, msg
  #
  msg = pdb_inp.have_altloc_mix_message(prefix=prefix+"  ")
  if (msg is None):
    msg = pdb_inp.have_blank_altloc_message(prefix=prefix+"  ")
  if (msg is not None): print >> out, msg
  #
  if (level_id is not None):
    print >> out, prefix+"  hierarchy level of detail = %s" % level_id
    hierarchy.show(
      out=out,
      prefix=prefix+"    ",
      level_id=level_id,
      level_id_exception=level_id_exception)
  #
  return pdb_inp, hierarchy

class _atom(boost.python.injector, ext.atom):

  def format_atom_record(self,
        serial,
        input_atom_labels=None,
        record_name=None,
        name=None,
        altloc=None,
        resname=None,
        chain=None,
        resseq=None,
        icode=None,
        xyz=None,
        occ=None,
        b=None,
        segid=None,
        element=None,
        charge=None):
    if (altloc is None):
      if (input_atom_labels is None):
        altloc = " "
      else:
        altloc = input_atom_labels.altloc()
    if (resname is None):
      if (input_atom_labels is None):
        resname = "DUM"
      else:
        resname = input_atom_labels.resname()
    if (chain is None):
      if (input_atom_labels is None):
        chain = " "
      else:
        chain = input_atom_labels.chain()
    if (resseq is None):
      if (input_atom_labels is None):
        resseq = 1
      else:
        resseq = input_atom_labels.resseq()
    if (icode is None):
      if (input_atom_labels is None):
        icode = " "
      else:
        icode = input_atom_labels.icode()
    if (record_name is None):
      if (self.hetero): record_name = "HETATM"
      else:             record_name = "ATOM"
    if (name is None): name = self.name
    if (xyz is None): xyz = self.xyz
    if (occ is None): occ = self.occ
    if (b is None): b = self.b
    if (segid is None): segid = self.segid
    if (element is None): element = self.element
    if (charge is None): charge = self.charge
    return format_atom_record(
      record_name=record_name,
      serial=serial,
      name=name,
      altLoc=altloc,
      resName=resname,
      chainID=chain,
      resSeq=resseq,
      iCode=icode,
      site=xyz,
      occupancy=occ,
      tempFactor=b,
      segID=segid,
      element=element,
      charge=charge)

  def format_anisou_record(self,
        serial,
        input_atom_labels=None,
        name=None,
        altloc=None,
        resname=None,
        chain=None,
        resseq=None,
        icode=None,
        uij=None,
        segid=None,
        element=None,
        charge=None):
    if (altloc is None):
      if (input_atom_labels is None):
        altloc = " "
      else:
        altloc = input_atom_labels.altloc()
    if (resname is None):
      if (input_atom_labels is None):
        resname = "DUM"
      else:
        resname = input_atom_labels.resname()
    if (chain is None):
      if (input_atom_labels is None):
        chain = " "
      else:
        chain = input_atom_labels.chain()
    if (resseq is None):
      if (input_atom_labels is None):
        resseq = 1
      else:
        resseq = input_atom_labels.resseq()
    if (icode is None):
      if (input_atom_labels is None):
        icode = " "
      else:
        icode = input_atom_labels.resseq()
    if (name is None): name = self.name
    if (uij is None): uij = self.uij
    if (segid is None): segid = self.segid
    if (element is None): element = self.element
    if (charge is None): charge = self.charge
    return format_anisou_record(
      serial=serial,
      name=name,
      altLoc=altloc,
      resName=resname,
      chainID=chain,
      resSeq=resseq,
      iCode=icode,
      u_cart=uij,
      segID=segid,
      element=element,
      charge=charge)

  def format_record_group(self,
        serial,
        input_atom_labels=None,
        name=None,
        altloc=None,
        resname=None,
        chain=None,
        resseq=None,
        icode=None,
        xyz=None,
        occ=None,
        b=None,
        uij=None,
        segid=None,
        element=None,
        charge=None):
     result = [self.format_atom_record(
       serial=serial,
       input_atom_labels=input_atom_labels,
       name=name,
       altloc=altloc,
       resname=resname,
       chain=chain,
       resseq=resseq,
       icode=icode,
       xyz=xyz,
       occ=occ,
       b=b,
       segid=segid,
       element=element,
       charge=charge)]
     if (uij is None): uij = self.uij
     if (uij != (-1,-1,-1,-1,-1,-1)):
       result.append(self.format_anisou_record(
         serial=serial,
         input_atom_labels=input_atom_labels,
         name=name,
         altloc=altloc,
         resname=resname,
         chain=chain,
         resseq=resseq,
         icode=icode,
         uij=uij,
         segid=segid,
         element=element,
         charge=charge))
     return result

class _residue(boost.python.injector, ext.residue):

  def residue_name_plus_atom_names_interpreter(self,
        translate_cns_dna_rna_residue_names=None,
        return_mon_lib_dna_name=False):
    return residue_name_plus_atom_names_interpreter(
      residue_name=self.name,
      atom_names=self.atom_names(),
      translate_cns_dna_rna_residue_names=translate_cns_dna_rna_residue_names,
      return_mon_lib_dna_name=return_mon_lib_dna_name)

class _conformer(boost.python.injector, ext.conformer):

  def format_fasta(self, max_line_length=79):
    rn_seq = []
    residue_classes = dict_with_default_0()
    for residue in self.residues():
      rnpani = residue.residue_name_plus_atom_names_interpreter()
      rn = rnpani.work_residue_name
      rn_seq.append(rn)
      if (rn is None):
        c = None
      else:
        c = common_residue_names_get_class(name=rn)
      residue_classes[c] += 1
    seq = []
    n_aa = residue_classes["common_amino_acid"]
    n_na = residue_classes["common_rna_dna"]
    if (n_aa > n_na):
      from iotbx.pdb.amino_acid_codes import one_letter_given_three_letter
      aa_3_as_1 = one_letter_given_three_letter.get
      for rn in rn_seq:
        seq.append(aa_3_as_1(rn, "X"))
    elif (n_na != 0):
      for rn in rn_seq:
        seq.append({
          "A": "A",
          "C": "C",
          "G": "G",
          "U": "U",
          "DA": "A",
          "DC": "C",
          "DG": "G",
          "DT": "T"}.get(rn, "N"))
    n = len(seq)
    if (n == 0): return None
    comment = [">"]
    p = self.parent()
    if (p is not None):
      comment.append('chain "%2s"' % p.id)
    comment.append('conformer "%s"' % self.id)
    result = [" ".join(comment)]
    i = 0
    while True:
      j = min(n, i+max_line_length)
      if (j == i): break
      result.append("".join(seq[i:j]))
      i = j
    return result

class _chain(boost.python.injector, ext.chain):

  def combine_conformers(self):
    result = []
    conformers = [conformer.residues() for conformer in self.conformers()]
    pivots = [0] * len(conformers)
    while True:
      pivot_memory_ids = {}
      for i_conf,i_piv in enumerate(pivots):
        residues = conformers[i_conf]
        if (i_piv == len(residues)): continue
        pivot_memory_ids[residues[i_piv].memory_id()] = i_conf
      if (len(pivot_memory_ids) == 0):
        break
      def find_i_conf_next():
        for i_conf,i_piv in enumerate(pivots):
          if (i_piv == len(conformers[i_conf])): continue
          def can_be_next():
            result = set([i_conf])
            for atom in conformers[i_conf][i_piv].atoms():
              if (atom.is_alternative()): continue
              for parent in atom.parents():
                j_conf = pivot_memory_ids.get(parent.memory_id())
                if (j_conf is None):
                  return None
                result.add(j_conf)
            return result
          i_confs = can_be_next()
          if (i_confs is not None):
            return i_confs
        raise RuntimeError
      i_confs = sorted(find_i_conf_next())
      flag = True
      for i_conf in i_confs:
        i_residue = pivots[i_conf]
        residue = conformers[i_conf][i_residue]
        suppress_chain_break = (i_residue == 0)
        for atom in residue.atoms():
          if (flag or atom.is_alternative()):
            if (not suppress_chain_break and not residue.link_to_previous):
              result.append(None)
            result.append(atom)
            suppress_chain_break = True
        flag = False
        pivots[i_conf] += 1
    return result

hierarchy_level_ids = ["model", "chain", "conformer", "residue", "atom"]

class _hierarchy(boost.python.injector, ext.hierarchy):

  def show(self,
        out=None,
        prefix="",
        level_id=None,
        level_id_exception=ValueError):
    if (level_id == None): level_id = "atom"
    try: level_no = hierarchy_level_ids.index(level_id)
    except ValueError:
      raise level_id_exception('Unknown level_id="%s"' % level_id)
    if (out is None): out = sys.stdout
    for model in self.models():
      chains = model.chains()
      print >> out, prefix+"model id=%d" % model.id, \
        "#chains=%d" % len(chains)
      if (level_no == 0): continue
      assert model.parent().memory_id() == self.memory_id()
      for chain in chains:
        conformers = chain.conformers()
        print >> out, prefix+'  chain id="%s"' % chain.id, \
          "#conformers=%d" % len(conformers)
        if (level_no == 1): continue
        assert chain.parent().memory_id() == model.memory_id()
        for conformer in conformers:
          residues = conformer.residues()
          print >> out, prefix+'    conformer id="%s"' % conformer.id, \
            "#residues=%d" % len(residues), \
            "#atoms=%d" % conformer.number_of_atoms(),
          n_alt = conformer.number_of_alternative_atoms()
          if (n_alt != 0):
            print >> out, "#alt. atoms=%d" % n_alt,
          print >> out
          if (level_no == 2): continue
          assert conformer.parent().memory_id() == chain.memory_id()
          suppress_chain_break = True
          for residue in residues:
            if (not residue.link_to_previous and not suppress_chain_break):
              print >> out, prefix+"      ### chain break ###"
            suppress_chain_break = False
            atoms = residue.atoms()
            print >> out, prefix+'      residue name="%s"' % residue.name, \
              'seq="%s"' % residue.seq, 'icode="%s"' % residue.icode, \
              "#atoms=%d" % len(atoms)
            if (level_no == 3): continue
            assert residue.parent().memory_id() == conformer.memory_id()
            for atom in atoms:
              if (atom.is_alternative()):
                mark = "&"
              else:
                mark = " "
              print >> out, prefix+'       %s "%s"' % (mark, atom.name)
              for parent in atom.parents():
                if (parent.memory_id() == residue.memory_id()):
                  break
              else:
                raise RuntimeError(
                  "parent residue not in list of atom parents")

  def as_str(self,
        prefix="",
        level_id=None,
        level_id_exception=ValueError):
    out = StringIO()
    self.show(
      out=out,
      prefix=prefix,
      level_id=level_id,
      level_id_exception=level_id_exception)
    return out.getvalue()

  def as_pdb_records(self, pdb_inp, append_end=False):
    result = []
    iall = pdb_inp.input_atom_labels_list()
    pdb_inp.reset_atom_tmp()
    models = self.models()
    for model in models:
      if (len(models) != 1):
        result.append("MODEL %7d" % model.id)
      atom_serial = 0
      for chain in model.chains():
        for atom in chain.combine_conformers():
          if (atom is None):
            result.append("BREAK")
            continue
          atom_serial += 1
          result.append(atom.format_atom_record(
            serial=atom_serial,
            input_atom_labels=iall[atom.tmp]))
        result.append("TER")
      if (len(models) != 1):
        result.append("ENDMDL")
    if (append_end):
      result.append("END")
    return result

  def as_pdb_string(self, pdb_inp, append_end=False):
    return "\n".join(self.as_pdb_records(
      pdb_inp=pdb_inp, append_end=append_end))+"\n"

default_atom_names_scattering_type_const = ["PEAK", "SITE"]

def _one_of(n):
  if (n != 1): return " (one of %d)" % n
  return ""

class _input(boost.python.injector, ext.input):

  def have_blank_altloc_message(self, prefix=""):
    n = self.number_of_alternative_groups_with_blank_altloc()
    if (n == 0):
      return None
    result = [
      prefix + "alternative group with blank altloc%s:" % _one_of(n=n)]
    iall = self.input_atom_labels_list()
    for i_seq in self.i_seqs_alternative_group_with_blank_altloc():
      result.append(prefix + "  " + iall[i_seq].pdb_format())
    return "\n".join(result)

  def raise_blank_altloc_if_necessary(self):
    msg = self.have_blank_altloc_message()
    if (msg is not None): raise Sorry(msg)

  def have_altloc_mix_message(self, prefix=""):
    n = self.number_of_chains_with_altloc_mix()
    if (n == 0):
      return None
    iall = self.input_atom_labels_list()
    result = [
      prefix + "mix of alternative groups with and without blank altlocs"
        " (in %d chain%s):" % plural_s(n)]
    result.append(self.have_blank_altloc_message(prefix=prefix+"  "))
    result.append(prefix + "  alternative group without blank altloc%s:" %
      _one_of(n=self.number_of_alternative_groups_without_blank_altloc()))
    for i_seq in self.i_seqs_alternative_group_without_blank_altloc():
      result.append(prefix + "    " + iall[i_seq].pdb_format())
    return "\n".join(result)

  def raise_altloc_mix_if_necessary(self):
    msg = self.have_altloc_mix_message()
    if (msg is not None): raise Sorry(msg)

  def have_duplicate_atom_labels_message(self, max_show=10, prefix=""):
    dup = self.find_duplicate_atom_labels()
    if (dup.size() == 0):
      return None
    result = [
      prefix + "number of groups of duplicate atom labels:  %3d" % dup.size(),
      prefix + "  total number of affected atoms:           %3d" %
        sum([i_seqs.size() for i_seqs in dup])]
    if (max_show > 0):
      iall = self.input_atom_labels_list()
      for i_seqs in dup[:max_show]:
        prfx = "  group "
        for i_seq in i_seqs:
          result.append(prefix + prfx + iall[i_seq].pdb_format())
          prfx = "        "
      if (dup.size() > max_show):
        result.append(prefix + "  ... %d remaining group%s not shown" %
          plural_s(dup.size()-max_show))
    return "\n".join(result)

  def raise_duplicate_atom_labels_if_necessary(self, max_show=10):
    msg = self.have_duplicate_atom_labels_message(max_show=max_show)
    if (msg is not None): raise Sorry(msg)

  def crystal_symmetry_from_cryst1(self):
    from iotbx.pdb import cryst1_interpretation
    for line in self.crystallographic_section():
      if (line.startswith("CRYST1")):
        return cryst1_interpretation.crystal_symmetry(cryst1_record=line)
    return None

  def scale_matrix(self):
    if (not hasattr(self, "_scale_matrix")):
      source_info = self.source_info()
      if (len(source_info) > 0): source_info = " (%s)" % source_info
      self._scale_matrix = [[None]*9,[None]*3]
      done = {}
      for line in self.crystallographic_section():
        if (line.startswith("SCALE") and line[5:6] in ["1", "2", "3"]):
          r = read_scale_record(line=line, source_info=source_info)
          for i_col,v in enumerate([r.sn1, r.sn2, r.sn3]):
            self._scale_matrix[0][(r.n-1)*3+i_col] = v
          self._scale_matrix[1][r.n-1] = r.un
          done[r.n] = None
      done = done.keys()
      done.sort()
      if (len(done) == 0):
        self._scale_matrix = None
      elif (done != [1,2,3]):
        raise RuntimeError(
          "Incomplete set of PDB SCALE records%s" % source_info)
    return self._scale_matrix

  def xray_structure_simple(self,
        crystal_symmetry=None,
        weak_symmetry=False,
        unit_cube_pseudo_crystal=False,
        fractional_coordinates=False,
        use_scale_matrix_if_available=True,
        min_distance_sym_equiv=0.5,
        scattering_type_exact=False,
        enable_scattering_type_unknown=False,
        atom_names_scattering_type_const
          =default_atom_names_scattering_type_const):
    return self.xray_structures_simple(
      one_structure_for_each_model=False,
      crystal_symmetry=crystal_symmetry,
      weak_symmetry=weak_symmetry,
      unit_cube_pseudo_crystal=unit_cube_pseudo_crystal,
      fractional_coordinates=fractional_coordinates,
      use_scale_matrix_if_available=use_scale_matrix_if_available,
      min_distance_sym_equiv=min_distance_sym_equiv,
      scattering_type_exact=scattering_type_exact,
      enable_scattering_type_unknown=enable_scattering_type_unknown,
      atom_names_scattering_type_const=atom_names_scattering_type_const)[0]

  def xray_structures_simple(self,
        one_structure_for_each_model=True,
        crystal_symmetry=None,
        weak_symmetry=False,
        unit_cube_pseudo_crystal=False,
        fractional_coordinates=False,
        min_distance_sym_equiv=0.5,
        use_scale_matrix_if_available=True,
        scattering_type_exact=False,
        enable_scattering_type_unknown=False,
        atom_names_scattering_type_const
          =default_atom_names_scattering_type_const):
    from cctbx import xray
    from cctbx import crystal
    assert crystal_symmetry is None or not unit_cube_pseudo_crystal
    if (not unit_cube_pseudo_crystal):
      cryst1_symmetry = self.crystal_symmetry_from_cryst1()
      if (crystal_symmetry is None):
        crystal_symmetry = cryst1_symmetry
      elif (cryst1_symmetry is not None):
        crystal_symmetry = cryst1_symmetry.join_symmetry(
          other_symmetry=crystal_symmetry,
          force=not weak_symmetry)
    if (crystal_symmetry is None
        or (    crystal_symmetry.unit_cell() is None
            and crystal_symmetry.space_group_info() is None)):
      unit_cube_pseudo_crystal = True
      crystal_symmetry = crystal.symmetry(
        unit_cell=(1,1,1,90,90,90),
        space_group_symbol="P1")
    assert crystal_symmetry.unit_cell() is not None
    assert crystal_symmetry.space_group_info() is not None
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
      self,
      one_structure_for_each_model,
      unit_cube_pseudo_crystal,
      fractional_coordinates,
      scattering_type_exact,
      enable_scattering_type_unknown,
      scitbx.stl.set.stl_string(atom_names_scattering_type_const),
      unit_cell,
      scale_r,
      scale_t)
    special_position_settings = crystal_symmetry.special_position_settings(
      min_distance_sym_equiv=min_distance_sym_equiv)
    while (loop.next()):
      result.append(xray.structure(
        special_position_settings=special_position_settings,
        scatterers=loop.scatterers))
    return result

def atom_labels_pdb_format(model, chain, conformer, residue, atom):
  result = []
  if (model.id != 0):
    result.append('model="%4d"' % model.id)
  result.append('"%4s%1s%3s%2s%4s%1s"' % (
    atom.name,
    conformer.id,
    residue.name,
    chain.id,
    residue.seq,
    residue.icode))
  if (atom.segid != "    "):
    result.append('segid="%s"' % atom.segid)
  return " ".join(result)

def format_cryst1_record(crystal_symmetry, z=None):
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
    + (str(crystal_symmetry.space_group_info()), z))).rstrip()

def format_scale_records(unit_cell=None,
                         fractionalization_matrix=None,
                         u=[0,0,0]):
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
  return ("SCALE1    %10.6f%10.6f%10.6f     %10.5f\n"
          "SCALE2    %10.6f%10.6f%10.6f     %10.5f\n"
          "SCALE3    %10.6f%10.6f%10.6f     %10.5f") % (
    f[0], f[1], f[2], u[0],
    f[3], f[4], f[5], u[1],
    f[6], f[7], f[8], u[2])

class read_scale_record:

  def __init__(self, line, source_info=""):
    try: self.n = int(line[5:6])
    except ValueError: self.n = None
    if (self.n not in [1,2,3]):
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
    self.sn1, self.sn2, self.sn3, self.un = values

def resseq_decode(s):
  return hy36decode(width=4, s="%4s" % s)

def resseq_encode(value):
  return hy36encode(width=4, value=value)

def encode_serial_number(width, value):
  if (isinstance(value, str)):
    assert len(value) <= width
    return value
  if (isinstance(value, int)):
    return hy36encode(width=width, value=value)
  raise RuntimeError("serial number value must be str or int.")

def format_atom_record(record_name="ATOM",
                       serial=0,
                       name=" C  ",
                       altLoc=" ",
                       resName="DUM",
                       chainID=" ",
                       resSeq=1,
                       iCode=" ",
                       site=(0,0,0),
                       occupancy=1,
                       tempFactor=0,
                       segID="    ",
                       element="  ",
                       charge="  "):
  # ATOM
  #  7 - 11  Integer       serial        Atom serial number.
  # 13 - 16  Atom          name          Atom name.
  # 17       Character     altLoc        Alternate location indicator.
  # 18 - 20  Residue name  resName       Residue name.
  # 21 - 22                chainID       Chain identifier.
  # 23 - 26  Integer       resSeq        Residue sequence number.
  # 27       AChar         iCode         Code for insertion of residues.
  # 31 - 38  Real(8.3)     x             Orthogonal coordinates for X in
  #                                      Angstroms.
  # 39 - 46  Real(8.3)     y             Orthogonal coordinates for Y in
  #                                      Angstroms.
  # 47 - 54  Real(8.3)     z             Orthogonal coordinates for Z in
  #                                      Angstroms.
  # 55 - 60  Real(6.2)     occupancy     Occupancy.
  # 61 - 66  Real(6.2)     tempFactor    Temperature factor.
  # 73 - 76  LString(4)    segID         Segment identifier, left-justified.
  # 77 - 78  LString(2)    element       Element symbol, right-justified.
  # 79 - 80  LString(2)    charge        Charge on the atom.
  return ((
    "%-6.6s%5.5s %-4.4s%1.1s%-3.3s%2.2s%4.4s%1.1s"
    "   %8.3f%8.3f%8.3f%6.2f%6.2f    "
    "  %-4.4s%2.2s%2.2s") % (
      record_name,
      encode_serial_number(width=5, value=serial),
      name, altLoc, resName, chainID,
      encode_serial_number(width=4, value=resSeq), iCode,
      site[0], site[1], site[2], occupancy, tempFactor,
      segID, element, charge)).rstrip()

def format_anisou_record(
      serial=0,
      name=" C  ",
      altLoc=" ",
      resName="DUM",
      chainID=" ",
      resSeq=1,
      iCode=" ",
      u_cart=(0,0,0,0,0,0),
      segID="    ",
      element="  ",
      charge="  "):
  # ANISOU
  #  7 - 11  Integer       serial        Atom serial number.
  # 13 - 16  Atom          name          Atom name.
  # 17       Character     altLoc        Alternate location indicator.
  # 18 - 20  Residue name  resName       Residue name.
  # 21 - 22                chainID       Chain identifier.
  # 23 - 26  Integer       resSeq        Residue sequence number.
  # 27       AChar         iCode         Code for insertion of residues.
  # 29 - 35  Integer       u[0][0]       U(1,1)
  # 36 - 42  Integer       u[1][1]       U(2,2)
  # 43 - 49  Integer       u[2][2]       U(3,3)
  # 50 - 56  Integer       u[0][1]       U(1,2)
  # 57 - 63  Integer       u[0][2]       U(1,3)
  # 64 - 70  Integer       u[1][2]       U(2,3)
  # 73 - 76  LString(4)    segID         Segment identifier, left-justified.
  # 77 - 78  LString(2)    element       Element symbol, right-justified.
  # 79 - 80  LString(2)    charge        Charge on the atom.
  return ((
    "%-6.6s%5.5s %-4.4s%1.1s%-3.3s%2.2s%4.4s%1.1s"
    " %7d%7d%7d%7d%7d%7d"
    "  %-4.4s%2.2s%2.2s") % ((
      "ANISOU",
      encode_serial_number(width=5, value=serial),
      name, altLoc, resName, chainID,
      encode_serial_number(width=4, value=resSeq), iCode)
    + tuple([iround(u*10000) for u in u_cart])
    + (segID, element, charge))).rstrip()

def format_ter_record(serial=0,
                      resName="DUM",
                      chainID=" ",
                      resSeq=1,
                      iCode=" "):
  #  7 - 11  Integer         serial     Serial number.
  # 18 - 20  Residue name    resName    Residue name.
  # 21 - 22  Character       chainID    Chain identifier.
  # 23 - 26  Integer         resSeq     Residue sequence number.
  # 27       AChar           iCode      Insertion code.
  return ("%-6.6s%5.5s      %-3.3s%2.2s%4.4s%1.1s" % (
    "TER",
    encode_serial_number(width=5, value=serial), resName, chainID,
    encode_serial_number(width=4, value=resSeq), iCode)).rstrip()
