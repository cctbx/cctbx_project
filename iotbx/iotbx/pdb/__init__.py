import cctbx.array_family.flex

import boost.python
ext = boost.python.import_ext("iotbx_pdb_ext")
from iotbx_pdb_ext import *

from iotbx.pdb.xray_structure import from_pdb as as_xray_structure
from scitbx.python_utils.math_utils import iround
from libtbx.str_utils import show_string, show_sorted_by_counts
import sys

def show_summary(
      file_name=None,
      source_info=None,
      lines=None,
      level_id=None,
      level_id_exception=ValueError,
      out=None,
      prefix=""):
  if (file_name is not None):
    assert source_info is None and lines is None
    pdb_inp = input(file_name=file_name)
  else:
    assert file_name is None
    pdb_inp = input(source_info=source_info, lines=lines)
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
  dup = pdb_inp.find_duplicate_atom_labels()
  if (dup.size() > 0):
    print >> out, prefix+"  number of groups of duplicate atom lables:  %3d" \
      % dup.size()
    print >> out, prefix+"    total number of affected atoms:           %3d" \
      % sum([i_seqs.size() for i_seqs in dup])
    iall = pdb_inp.input_atom_labels_list()
    for i_seqs in dup[:10]:
      prfx = "    group"
      for i_seq in i_seqs:
        print >> out, prefix+prfx, iall[i_seq].pdb_format()
        prfx = "         "
    if (dup.size() > 10):
      print >> out, prefix+"    ... remaining %d groups not shown" % (
        dup.size()-10)
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
              "seq=%4d" % residue.seq, 'icode="%s"' % residue.icode, \
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
  # 22       Character     chainID       Chain identifier.
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
    "%-6.6s%5d %-4.4s%1.1s%-4.4s%1.1s%4d%1.1s"
    "   %8.3f%8.3f%8.3f%6.2f%6.2f    "
    "  %-4.4s%2.2s%2.2s") % (
      record_name,
      serial%100000, name, altLoc, resName, chainID, resSeq%10000, iCode,
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
  # 22       Character     chainID       Chain identifier.
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
    "%-6.6s%5d %-4.4s%1.1s%-3.3s %1.1s%4d%1.1s"
    " %7d%7d%7d%7d%7d%7d"
    "  %-4.4s%2.2s%2.2s") % ((
      "ANISOU",
      serial%100000, name, altLoc, resName, chainID, resSeq%10000, iCode)
    + tuple([iround(u*10000) for u in u_cart])
    + (segID, element, charge))).rstrip()

def format_ter_record(serial=0,
                      resName="DUM",
                      chainID=" ",
                      resSeq=1,
                      iCode=" "):
  #  7 - 11  Integer         serial     Serial number.
  # 18 - 20  Residue name    resName    Residue name.
  # 22       Character       chainID    Chain identifier.
  # 23 - 26  Integer         resSeq     Residue sequence number.
  # 27       AChar           iCode      Insertion code.
  return ("%-6.6s%5d      %-3.3s %1.1s%4d%1.1s" % (
    "TER", serial, resName, chainID, resSeq, iCode)).rstrip()
