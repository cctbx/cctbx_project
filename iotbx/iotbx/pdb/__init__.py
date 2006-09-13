from cctbx.array_family import flex

import boost.python
ext = boost.python.import_ext("iotbx_pdb_ext")
from iotbx_pdb_ext import *

from iotbx.pdb.xray_structure import from_pdb as as_xray_structure
from scitbx import matrix
from scitbx.python_utils.math_utils import iround
from libtbx.str_utils import show_string, show_sorted_by_counts
from libtbx.itertbx import count
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
  msg = pdb_inp.have_altloc_mix_message(prefix=prefix+"  ")
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

default_atom_names_scattering_type_const = ["PEAK", "SITE"]

class _input(boost.python.injector, ext.input):

  def have_altloc_mix_message(self, prefix=""):
    if (self.number_of_chains_with_altloc_mix() == 0):
      return None
    iall = self.input_atom_labels_list()
    result = [
      prefix + "mix of alternative groups with and without blank altlocs:"]
    result.append(prefix + "  alternative group with blank altloc:")
    for i_seq in self.i_seqs_alternative_group_with_blank_altloc():
      result.append(prefix + "    " + iall[i_seq].pdb_format())
    result.append(prefix + "  alternative group without blank altloc:")
    for i_seq in self.i_seqs_alternative_group_without_blank_altloc():
      result.append(prefix + "    " + iall[i_seq].pdb_format())
    return "\n".join(result)

  def raise_altloc_mix_if_necessary(self):
    msg = self.have_altloc_mix_message()
    if (msg is not None): raise RuntimeError(msg)

  def crystal_symmetry_from_cryst1(self):
    from iotbx.pdb import cryst1_interpretation
    for line in self.crystallographic_section():
      if (line.startswith("CRYST1")):
        return cryst1_interpretation.crystal_symmetry(cryst1_record=line)

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
        unit_cube_pseudo_crystal=False,
        use_scale_matrix_if_available=True,
        scattering_type_exact=False,
        enable_scattering_type_unknown=False,
        atom_element_q_substitute=None,
        atom_names_scattering_type_const
          =default_atom_names_scattering_type_const):
    return self.xray_structures_simple(
      one_structure_for_each_model=False,
      crystal_symmetry=crystal_symmetry,
      unit_cube_pseudo_crystal=unit_cube_pseudo_crystal,
      use_scale_matrix_if_available=use_scale_matrix_if_available,
      scattering_type_exact=scattering_type_exact,
      enable_scattering_type_unknown=enable_scattering_type_unknown,
      atom_element_q_substitute=atom_element_q_substitute,
      atom_names_scattering_type_const=atom_names_scattering_type_const)[0]

  def xray_structures_simple(self,
        one_structure_for_each_model=True,
        crystal_symmetry=None,
        unit_cube_pseudo_crystal=False,
        use_scale_matrix_if_available=True,
        scattering_type_exact=False,
        enable_scattering_type_unknown=False,
        atom_element_q_substitute=None,
        atom_names_scattering_type_const
          =default_atom_names_scattering_type_const):
    from cctbx import xray
    from cctbx import crystal
    from cctbx import eltbx
    import cctbx.eltbx.xray_scattering
    from cctbx import adptbx
    assert crystal_symmetry is None or not unit_cube_pseudo_crystal
    if (crystal_symmetry is None and not unit_cube_pseudo_crystal):
      crystal_symmetry = self.crystal_symmetry_from_cryst1()
    if (crystal_symmetry is None
        or (    crystal_symmetry.unit_cell() is None
            and crystal_symmetry.space_group_info() is None)):
      unit_cube_pseudo_crystal = True
      crystal_symmetry = crystal.symmetry(
        unit_cell=(1,1,1,90,90,90),
        space_group_symbol="P1")
    if (not unit_cube_pseudo_crystal):
      uc = crystal_symmetry.unit_cell()
      if (use_scale_matrix_if_available):
        scale_matrix = self.scale_matrix()
      else:
        scale_matrix = None
      if (scale_matrix is None):
        frac = uc.fractionalize
      else:
        scale_r = matrix.sqr(scale_matrix[0])
        scale_t = matrix.col(scale_matrix[1])
    result = []
    scatterers = flex.xray_scatterer()
    model_indices = list(self.model_indices())
    model_indices.reverse()
    model_index = model_indices.pop()
    for i_atom,labels,atom in zip(count(),
                                  self.input_atom_labels_list(),
                                  self.atoms()):
      if (i_atom == model_index):
        if (one_structure_for_each_model):
          result.append(xray.structure(
            crystal_symmetry=crystal_symmetry,
            scatterers=scatterers))
          scatterers = flex.xray_scatterer()
        model_index = model_indices.pop()
      label = labels.pdb_format()
      if (unit_cube_pseudo_crystal):
        site = atom.xyz
      elif (scale_matrix is None):
        site = frac(atom.xyz)
      else:
        site = (scale_r * matrix.col(atom.xyz) + scale_t).elems
      if (atom.uij == (-1,-1,-1,-1,-1,-1)):
        u = adptbx.b_as_u(atom.b)
      elif (unit_cube_pseudo_crystal):
        u = atom.uij
      else:
        u = adptbx.u_cart_as_u_star(uc, atom.uij)
      if (atom_element_q_substitute is not None
          and atom.element == " Q"):
        scattering_type = atom_element_q_substitute
      elif (atom_names_scattering_type_const is not None
            and atom.name in atom_names_scattering_type_const):
        scattering_type = "const"
      else:
        chemical_element = atom.determine_chemical_element_simple()
        if (chemical_element is None
            and not enable_scattering_type_unknown):
          raise RuntimeError(
            'Unknown chemical element type: PDB ATOM %s element="%s"\n'
              % (label, atom.element)
            + "  To resolve this problem, specify a chemical element type in\n"
            + '  columns 77-78 of the PDB file, right justified (e.g. " C").')
        if (atom.charge != "  " and atom.charge[0] == " "):
          if (not enable_scattering_type_unknown):
            raise RuntimeError(
              'Unknown charge: PDB ATOM %s element="%s" charge="%s"'
                % (label, atom.element, atom.charge))
          chemical_element = None
        if (chemical_element is not None):
          chemical_element += atom.charge
          scattering_type = eltbx.xray_scattering.get_standard_label(
            label=chemical_element, exact=scattering_type_exact, optional=True)
        else:
          scattering_type = None
        if (scattering_type is None):
          if (not enable_scattering_type_unknown):
            raise RuntimeError(
              'Unknown scattering type: PDB ATOM %s element="%s" charge="%s"'
                % (label, atom.element, atom.charge))
          else:
            scattering_type = "unknown"
      scatterers.append(xray.scatterer(
        label=label,
        site=site,
        u=u,
        occupancy=atom.occ,
        scattering_type=scattering_type))
    result.append(xray.structure(
      crystal_symmetry=crystal_symmetry,
      scatterers=scatterers))
    return result

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
