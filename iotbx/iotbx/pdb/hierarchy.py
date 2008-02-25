import boost.python
ext = boost.python.import_ext("iotbx_pdb_hierarchy_ext")
from iotbx_pdb_hierarchy_ext import *

from libtbx import dict_with_default_0
from cStringIO import StringIO
import sys

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
    from iotbx.pdb import format_atom_record
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
    from iotbx.pdb import format_anisou_record
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
    from iotbx.pdb import residue_name_plus_atom_names_interpreter
    return residue_name_plus_atom_names_interpreter(
      residue_name=self.name,
      atom_names=self.atom_names(),
      translate_cns_dna_rna_residue_names=translate_cns_dna_rna_residue_names,
      return_mon_lib_dna_name=return_mon_lib_dna_name)

class _conformer(boost.python.injector, ext.conformer):

  def format_fasta(self, max_line_length=79):
    from iotbx.pdb import common_residue_names_get_class
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

level_ids = ["model", "chain", "conformer", "residue", "atom"]

class _root(boost.python.injector, ext.root):

  def show(self,
        out=None,
        prefix="",
        level_id=None,
        level_id_exception=ValueError):
    if (level_id == None): level_id = "atom"
    try: level_no = level_ids.index(level_id)
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

  def as_pdb_records(self, append_end=False):
    result = []
    models = self.models()
    for model in models:
      if (len(models) != 1):
        result.append("MODEL %7d" % model.id)
      atom_serial = 0
      for chain in model.chains():
        atom_serial = chain.append_atom_records(
          pdb_records=result,
          atom_serial=atom_serial)
        result.append("TER")
      if (len(models) != 1):
        result.append("ENDMDL")
    if (append_end):
      result.append("END")
    return result

  def as_pdb_string(self, append_end=False):
    return "\n".join(self.as_pdb_records(append_end=append_end))+"\n"
