from __future__ import division
import boost.python
ext = boost.python.import_ext("iotbx_pdb_hierarchy_ext")
from iotbx_pdb_hierarchy_ext import *

from cctbx.array_family import flex
from libtbx.str_utils import show_sorted_by_counts
from libtbx.utils import Sorry, plural_s, null_out
from libtbx import Auto, dict_with_default_0
from cStringIO import StringIO
import sys

class pickle_import_trigger(object): pass

level_ids = ["model", "chain", "residue_group", "atom_group", "atom"]

def _show_residue_group(rg, out, prefix):
  atoms = rg.atoms()
  if (atoms.size() == 0):
    ch = rg.parent()
    if (ch is None): ch = "  "
    else:            ch = "%s" % ch.id
    print >> out, prefix+'empty: "%s%s"' % (ch, rg.resid())
  else:
    def show_atom(atom):
      print >> out, prefix+'"%s"' % atom.format_atom_record(
        replace_floats_with=".*.")
    if (atoms.size() <= 3):
      for atom in atoms: show_atom(atom)
    else:
      show_atom(atoms[0])
      print >> out, prefix+'... %d atom%s not shown' % plural_s(
        atoms.size()-2)
      show_atom(atoms[-1])

class overall_counts(object):

  def __init__(self):
    self._errors = None
    self._warnings = None

  def show(self,
        out=None,
        prefix="",
        flag_errors=True,
        flag_warnings=True,
        residue_groups_max_show=10,
        duplicate_atom_labels_max_show=10):
    from iotbx.pdb import common_residue_names_get_class
    if (out is None): out = sys.stdout
    self._errors = []
    self._warnings = []
    def add_err(msg):
      if (flag_errors): print >> out, prefix+msg
      self._errors.append(msg.strip())
    def add_warn(msg):
      if (flag_warnings): print >> out, prefix+msg
      self._warnings.append(msg.strip())
    fmt = "%%%dd" % len(str(self.n_atoms))
    print >> out, prefix+"total number of:"
    if (self.n_duplicate_model_ids != 0):
      add_err("  ### ERROR: duplicate model ids ###")
    if (self.n_empty_models != 0):
      add_warn("  ### WARNING: empty model ###")
    print >> out, prefix+"  models:    ", fmt % self.n_models,
    infos = []
    if (self.n_duplicate_model_ids != 0):
      infos.append("%d with duplicate model id%s" % plural_s(
        self.n_duplicate_model_ids))
    if (self.n_empty_models != 0):
      infos.append("%d empty" % self.n_empty_models)
    if (len(infos) != 0): print >> out, "(%s)" % "; ".join(infos),
    print >> out
    if (self.n_duplicate_chain_ids != 0):
      add_warn("  ### WARNING: duplicate chain ids ###")
    if (self.n_empty_chains != 0):
      add_warn("  ### WARNING: empty chain ###")
    print >> out, prefix+"  chains:    ", fmt % self.n_chains,
    infos = []
    if (self.n_duplicate_chain_ids != 0):
      infos.append("%d with duplicate chain id%s" % plural_s(
        self.n_duplicate_chain_ids))
    if (self.n_empty_chains != 0):
      infos.append("%d empty" % self.n_empty_chains)
    if (self.n_explicit_chain_breaks != 0):
      infos.append("%d explicit chain break%s" % plural_s(
        self.n_explicit_chain_breaks))
    if (len(infos) != 0): print >> out, "(%s)" % "; ".join(infos),
    print >> out
    print >> out, prefix+"  alt. conf.:", fmt % self.n_alt_conf
    print >> out, prefix+"  residues:  ", fmt % (
      self.n_residues + self.n_residue_groups + self.n_empty_residue_groups),
    if (self.n_residue_groups != 0):
      print >> out, "(%d with mixed residue names)" % self.n_residue_groups,
    print >> out
    if (self.n_duplicate_atom_labels != 0):
      add_err("  ### ERROR: duplicate atom labels ###")
    print >> out, prefix+"  atoms:     ", fmt % self.n_atoms,
    if (self.n_duplicate_atom_labels != 0):
      print >> out, "(%d with duplicate labels)" %self.n_duplicate_atom_labels,
    print >> out
    print >> out, prefix+"  anisou:    ", fmt % self.n_anisou
    if (self.n_empty_residue_groups != 0):
      add_warn("  ### WARNING: empty residue_group ###")
      print >> out, prefix+"  empty residue_groups:", \
        fmt % self.n_empty_residue_groups
    if (self.n_empty_atom_groups != 0):
      add_warn("  ### WARNING: empty atom_group ###")
      print >> out, prefix+"  empty atom_groups:", \
        fmt % self.n_empty_atom_groups
    #
    c = self.element_charge_types
    print >> out, prefix+"number of atom element+charge types:", len(c)
    if (len(c) != 0):
      print >> out, prefix+"histogram of atom element+charge frequency:"
      show_sorted_by_counts(c.items(), out=out, prefix=prefix+"  ")
    #
    c = self.resname_classes
    print >> out, prefix+"residue name classes:",
    if (len(c) == 0): print >> out, None,
    print >> out
    show_sorted_by_counts(c.items(), out=out, prefix=prefix+"  ")
    #
    c = self.chain_ids
    print >> out, prefix+"number of chain ids: %d" % len(c)
    if (len(c) != 0):
      print >> out, prefix+"histogram of chain id frequency:"
      show_sorted_by_counts(c.items(), out=out, prefix=prefix+"  ")
    #
    c = self.alt_conf_ids
    print >> out, prefix+"number of alt. conf. ids: %d" % len(c)
    if (len(c) != 0):
      print >> out, prefix+"histogram of alt. conf. id frequency:"
      show_sorted_by_counts(c.items(), out=out, prefix=prefix+"  ")
      #
      fmt = "%%%dd" % len(str(max(
        self.n_alt_conf_none,
        self.n_alt_conf_pure,
        self.n_alt_conf_proper,
        self.n_alt_conf_improper)))
      print >> out, prefix+"residue alt. conf. situations:"
      print >> out, prefix+"  pure main conf.:    ", fmt%self.n_alt_conf_none
      print >> out, prefix+"  pure alt. conf.:    ", fmt%self.n_alt_conf_pure
      print >> out, prefix+"  proper alt. conf.:  ", fmt%self.n_alt_conf_proper
      if (self.n_alt_conf_improper != 0):
        add_err("  ### ERROR: improper alt. conf. ###")
      print >> out, prefix+"  improper alt. conf.:", \
        fmt % self.n_alt_conf_improper
      self.show_chains_with_mix_of_proper_and_improper_alt_conf(
        out=out, prefix=prefix)
    #
    c = self.resnames
    print >> out, prefix+"number of residue names: %d" % len(c)
    if (len(c) != 0):
      print >> out, prefix+"histogram of residue name frequency:"
      annotation_appearance = {
        "common_amino_acid": None,
        "common_rna_dna": None,
        "common_water": "   common water",
        "common_small_molecule": "   common small molecule",
        "common_element": "   common element",
        "other": "   other"
      }
      show_sorted_by_counts(c.items(), out=out, prefix=prefix+"  ",
        annotations=[
          annotation_appearance[common_residue_names_get_class(name=name)]
            for name in c.keys()])
    #
    if (len(self.consecutive_residue_groups_with_same_resid) != 0):
      add_warn("### WARNING: consecutive residue_groups with same resid ###")
    self.show_consecutive_residue_groups_with_same_resid(
      out=out, prefix=prefix, max_show=residue_groups_max_show)
    #
    if (len(self.residue_groups_with_multiple_resnames_using_same_altloc)!= 0):
      add_err("### ERROR: residue group with multiple resnames using"
        " same altloc ###")
      self.show_residue_groups_with_multiple_resnames_using_same_altloc(
        out=out, prefix=prefix, max_show=residue_groups_max_show)
    #
    self.show_duplicate_atom_labels(
      out=out, prefix=prefix, max_show=duplicate_atom_labels_max_show)

  def as_str(self,
        prefix="",
        residue_groups_max_show=10,
        duplicate_atom_labels_max_show=10):
    out = StringIO()
    self.show(
      out=out,
      prefix=prefix,
      residue_groups_max_show=residue_groups_max_show,
      duplicate_atom_labels_max_show=duplicate_atom_labels_max_show)
    return out.getvalue()

  def errors(self):
    if (self._errors is None): self.show(out=null_out())
    return self._errors

  def warnings(self):
    if (self._warnings is None): self.show(out=null_out())
    return self._warnings

  def errors_and_warnings(self):
    return self.errors() + self.warnings()

  def show_improper_alt_conf(self, out=None, prefix=""):
    if (self.n_alt_conf_improper == 0): return
    if (out is None): out = sys.stdout
    for residue_group,label in [(self.alt_conf_proper, "proper"),
                                (self.alt_conf_improper, "improper")]:
      if (residue_group is None): continue
      print >> out, prefix+"residue with %s altloc" % label
      for ag in residue_group.atom_groups():
        for atom in ag.atoms():
          print >> out, prefix+'  "%s"' % atom.format_atom_record(
            replace_floats_with=".*.")

  def raise_improper_alt_conf_if_necessary(self):
    sio = StringIO()
    self.show_improper_alt_conf(out=sio)
    msg = sio.getvalue()
    if (len(msg) != 0): raise Sorry(msg.rstrip())

  def show_chains_with_mix_of_proper_and_improper_alt_conf(self,
        out=None,
        prefix=""):
    if (out is None): out = sys.stdout
    n = self.n_chains_with_mix_of_proper_and_improper_alt_conf
    print >> out, \
      prefix+"chains with mix of proper and improper alt. conf.:", n
    if (n != 0): prefix = prefix + "  "
    self.show_improper_alt_conf(out=out, prefix=prefix)

  def raise_chains_with_mix_of_proper_and_improper_alt_conf_if_necessary(self):
    if (self.n_chains_with_mix_of_proper_and_improper_alt_conf == 0):
      return
    sio = StringIO()
    self.show_chains_with_mix_of_proper_and_improper_alt_conf(out=sio)
    raise Sorry(sio.getvalue().rstrip())

  def show_consecutive_residue_groups_with_same_resid(self,
        out=None,
        prefix="",
        max_show=10):
    cons = self.consecutive_residue_groups_with_same_resid
    if (len(cons) == 0): return
    if (out is None): out = sys.stdout
    print >> out, \
      prefix+"number of consecutive residue groups with same resid: %d" % \
        len(cons)
    if (max_show is None): max_show = len(cons)
    elif (max_show <= 0): return
    delim = prefix+"  "+"-"*42
    prev_rg = None
    for rgs in cons[:max_show]:
      for next,rg in zip(["", "next "], rgs):
        if (    prev_rg is not None
            and prev_rg.memory_id() == rg.memory_id()): continue
        elif (next == "" and prev_rg is not None):
          print >> out, delim
        prev_rg = rg
        print >> out, prefix+"  %sresidue group:" % next
        _show_residue_group(rg=rg, out=out, prefix=prefix+"    ")
    if (len(cons) > max_show):
      print >> out, delim
      print >> out, prefix + "  ... %d remaining instance%s not shown" % \
        plural_s(len(cons)-max_show)

  def show_residue_groups_with_multiple_resnames_using_same_altloc(self,
        out=None,
        prefix="",
        max_show=10):
    rgs = self.residue_groups_with_multiple_resnames_using_same_altloc
    if (len(rgs) == 0): return
    print >> out, prefix+"residue groups with multiple resnames using" \
      " same altloc:", len(rgs)
    if (max_show is None): max_show = len(cons)
    elif (max_show <= 0): return
    for rg in rgs[:max_show]:
      print >> out, prefix+"  residue group:"
      _show_residue_group(rg=rg, out=out, prefix=prefix+"    ")
    if (len(rgs) > max_show):
      print >> out, prefix + "  ... %d remaining instance%s not shown" % \
        plural_s(len(rgs)-max_show)

  def \
    raise_residue_groups_with_multiple_resnames_using_same_altloc_if_necessary(
        self, max_show=10):
    sio = StringIO()
    self.show_residue_groups_with_multiple_resnames_using_same_altloc(
      out=sio, max_show=max_show)
    msg = sio.getvalue()
    if (len(msg) != 0): raise Sorry(msg.rstrip())

  def show_duplicate_atom_labels(self, out=None, prefix="", max_show=10):
    dup = self.duplicate_atom_labels
    if (len(dup) == 0): return
    if (out is None): out = sys.stdout
    fmt = "%%%dd" % len(str(self.n_duplicate_atom_labels))
    print >> out, prefix+"number of groups of duplicate atom labels:", \
      fmt % len(dup)
    print >> out, prefix+"  total number of affected atoms:         ", \
      fmt % self.n_duplicate_atom_labels
    if (max_show is None): max_show = len(dup)
    elif (max_show <= 0): return
    for atoms in dup[:max_show]:
      prfx = "  group "
      for atom in atoms:
        print >> out, prefix+prfx+'"%s"' % atom.format_atom_record(
          replace_floats_with=".*.")
        prfx = "        "
    if (len(dup) > max_show):
      print >> out, prefix+"  ... %d remaining group%s not shown" % \
        plural_s(len(dup)-max_show)

  def raise_duplicate_atom_labels_if_necessary(self, max_show=10):
    sio = StringIO()
    self.show_duplicate_atom_labels(out=sio, max_show=max_show)
    msg = sio.getvalue()
    if (len(msg) != 0): raise Sorry(msg.rstrip())

class __hash_eq_mixin(object):

  def __hash__(self):
    return hash(self.memory_id())

  def __eq__(self, other):
    if (isinstance(other, self.__class__)):
      return (self.memory_id() == other.memory_id())
    return False

  def __ne__(self, other):
    return not ( self == other )

class _(boost.python.injector, ext.root, __hash_eq_mixin):

  def __getstate__(self):
    version = 2
    pdb_string = StringIO()
    self._as_pdb_string_cstringio(
      cstringio=pdb_string,
      append_end=True,
      interleaved_conf=0,
      atoms_reset_serial_first_value=None,
      atom_hetatm=True,
      sigatm=True,
      anisou=True,
      siguij=True)
    return (version, pickle_import_trigger(), self.info, pdb_string.getvalue())

  def __setstate__(self, state):
    assert len(state) >= 3
    version = state[0]
    if   (version == 1): assert len(state) == 3
    elif (version == 2): assert len(state) == 4
    else: raise RuntimeError("Unknown version of pickled state.")
    self.info = state[-2]
    import iotbx.pdb
    models = iotbx.pdb.input(
      source_info="pickle",
      lines=flex.split_lines(state[-1])).construct_hierarchy().models()
    self.pre_allocate_models(number_of_additional_models=len(models))
    for model in models:
      self.append_model(model=model)

  def chains(self):
    for model in self.models():
      for chain in model.chains():
        yield chain

  def residue_groups(self):
    for model in self.models():
      for chain in model.chains():
        for rg in chain.residue_groups():
          yield rg

  def atom_groups(self):
    for model in self.models():
      for chain in model.chains():
        for rg in chain.residue_groups():
          for ag in rg.atom_groups():
            yield ag

  def only_model(self):
    assert self.models_size() == 1
    return self.models()[0]

  def only_chain(self):
    return self.only_model().only_chain()

  def only_residue_group(self):
    return self.only_chain().only_residue_group()

  def only_conformer(self):
    return self.only_chain().only_conformer()

  def only_atom_group(self):
    return self.only_residue_group().only_atom_group()

  def only_residue(self):
    return self.only_conformer().only_residue()

  def only_atom(self):
    return self.only_atom_group().only_atom()

  def overall_counts(self):
    result = overall_counts()
    self.get_overall_counts(result)
    return result

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
    if (self.models_size() == 0):
      print >> out, prefix+'### WARNING: empty hierarchy ###'
    model_ids = dict_with_default_0()
    for model in self.models():
      model_ids[model.id] += 1
    for model in self.models():
      chains = model.chains()
      if (model_ids[model.id] != 1):
        s = "  ### ERROR: duplicate model id ###"
      else: s = ""
      print >> out, prefix+'model id="%s"' % model.id, \
        "#chains=%d%s" % (len(chains), s)
      if (level_no == 0): continue
      if (model.chains_size() == 0):
        print >> out, prefix+'  ### WARNING: empty model ###'
      model_chain_ids = dict_with_default_0()
      for chain in chains:
        model_chain_ids[chain.id] += 1
      for chain in chains:
        rgs = chain.residue_groups()
        if (model_chain_ids[chain.id] != 1):
          s = "  ### WARNING: duplicate chain id ###"
        else: s = ""
        print >> out, prefix+'  chain id="%s"' % chain.id, \
          "#residue_groups=%d%s" % (len(rgs), s)
        if (level_no == 1): continue
        if (chain.residue_groups_size() == 0):
          print >> out, prefix+'    ### WARNING: empty chain ###'
        suppress_chain_break = True
        prev_resid = ""
        for rg in rgs:
          if (not rg.link_to_previous and not suppress_chain_break):
            print >> out, prefix+"    ### chain break ###"
          suppress_chain_break = False
          ags = rg.atom_groups()
          resnames = set()
          for ag in rg.atom_groups():
            resnames.add(ag.resname)
          infos = []
          if (len(resnames) > 1): infos.append("with mixed residue names")
          resid = rg.resid()
          if (prev_resid == resid): infos.append("same as previous resid")
          prev_resid = resid
          if (len(infos) != 0): s = "  ### Info: %s ###" % "; ".join(infos)
          else: s = ""
          print >> out, prefix+'    resid="%s"' % resid, \
            "#atom_groups=%d%s" % (len(ags), s)
          if (level_no == 2): continue
          if (rg.atom_groups_size() == 0):
            print >> out, prefix+'      ### WARNING: empty residue_group ###'
          for ag in ags:
            atoms = ag.atoms()
            print >> out, prefix+'      altloc="%s"' % ag.altloc, \
              'resname="%s"' % ag.resname, \
              "#atoms=%d" % len(atoms)
            if (level_no == 3): continue
            if (ag.atoms_size() == 0):
              print >> out, prefix+'        ### WARNING: empty atom_group ###'
            for atom in atoms:
              print >> out, prefix+'        "%s"' % atom.name

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

  def as_pdb_string(self,
        crystal_symmetry=None,
        cryst1_z=None,
        write_scale_records=True,
        append_end=False,
        interleaved_conf=0,
        atoms_reset_serial_first_value=None,
        atom_hetatm=True,
        sigatm=True,
        anisou=True,
        siguij=True,
        cstringio=None,
        return_cstringio=Auto):
    if (cstringio is None):
      cstringio = StringIO()
      if (return_cstringio is Auto):
        return_cstringio = False
    elif (return_cstringio is Auto):
      return_cstringio = True
    if (crystal_symmetry is not None or cryst1_z is not None):
      from iotbx.pdb import format_cryst1_and_scale_records
      print >> cstringio, format_cryst1_and_scale_records(
        crystal_symmetry=crystal_symmetry,
        cryst1_z=cryst1_z,
        write_scale_records=write_scale_records)
    self._as_pdb_string_cstringio(
      cstringio=cstringio,
      append_end=append_end,
      interleaved_conf=interleaved_conf,
      atoms_reset_serial_first_value=atoms_reset_serial_first_value,
      atom_hetatm=atom_hetatm,
      sigatm=sigatm,
      anisou=anisou,
      siguij=siguij)
    if (return_cstringio):
      return cstringio
    return cstringio.getvalue()

  def extract_xray_structure(self, crystal_symmetry=None) :
    import iotbx.pdb
    pdb_str = self.as_pdb_string(crystal_symmetry=crystal_symmetry)
    pdb_inp = iotbx.pdb.input(
      source_info="pdb_hierarchy",
      lines=flex.split_lines(pdb_str))
    return pdb_inp.xray_structure_simple()

  def adopt_xray_structure(self, xray_structure):
    from iotbx.pdb import common_residue_names_get_class as gc
    from cctbx import adptbx
    if(self.atoms().size() != xray_structure.scatterers().size()):
      raise RuntimeError("Incompatible size of hierarchy and scatterers array.")
    awl = self.atoms_with_labels()
    scatterers = xray_structure.scatterers()
    uc = xray_structure.unit_cell()
    orth = uc.orthogonalize
    def set_attr(sc, a):
      a.set_xyz(new_xyz=orth(sc.site))
      a.set_occ(new_occ=sc.occupancy)
      a.set_element(sc.scattering_type[:2])
      a.set_charge(sc.scattering_type[1:3])
      if(sc.u_iso != -1):
        a.set_b(new_b=adptbx.u_as_b(sc.u_iso))
      if(sc.u_star != (-1.0, -1.0, -1.0, -1.0, -1.0, -1.0)):
        a.set_uij(new_uij = adptbx.u_star_as_u_cart(uc,sc.u_star))
    for sc, a in zip(scatterers, awl):
      id_str = a.id_str()
      resname_from_sc = id_str[10:13]
      cl1 = gc(resname_from_sc)
      cl2 = gc(a.resname)
      if([cl1,cl2].count("common_water")==2 or
         (id_str.strip()=="" and cl2=="common_water")):
        assert sc.scattering_type.strip().lower() in ["o","h","d"]
        set_attr(sc=sc, a=a)
      else:
        # XXX may be fix it when creating IS ? or make another special case?
        if(sc.scattering_type[:2] != "IS"):
          assert sc.label == a.id_str()
        set_attr(sc=sc, a=a)

  def write_pdb_file(self,
        file_name,
        open_append=False,
        crystal_symmetry=None,
        cryst1_z=None,
        write_scale_records=True,
        append_end=False,
        interleaved_conf=0,
        atoms_reset_serial_first_value=None,
        atom_hetatm=True,
        sigatm=True,
        anisou=True,
        siguij=True):
    if (crystal_symmetry is not None or cryst1_z is not None):
      from iotbx.pdb import format_cryst1_and_scale_records
      if (open_append): mode = "ab"
      else:             mode = "wb"
      print >> open(file_name, mode), format_cryst1_and_scale_records(
        crystal_symmetry=crystal_symmetry,
        cryst1_z=cryst1_z,
        write_scale_records=write_scale_records)
      open_append = True
    self._write_pdb_file(
      file_name=file_name,
      open_append=open_append,
      append_end=append_end,
      interleaved_conf=interleaved_conf,
      atoms_reset_serial_first_value=atoms_reset_serial_first_value,
      atom_hetatm=atom_hetatm,
      sigatm=sigatm,
      anisou=anisou,
      siguij=siguij)

  def atoms_with_labels(self):
    for model in self.models():
      for chain in model.chains():
        is_first_in_chain = True
        for rg in chain.residue_groups():
          is_first_after_break = not (is_first_in_chain or rg.link_to_previous)
          for ag in rg.atom_groups():
            for atom in ag.atoms():
              yield atom_with_labels(
                atom=atom,
                model_id=model.id,
                chain_id=chain.id,
                resseq=rg.resseq,
                icode=rg.icode,
                altloc=ag.altloc,
                resname=ag.resname,
                is_first_in_chain=is_first_in_chain,
                is_first_after_break=is_first_after_break)
              is_first_in_chain = False
              is_first_after_break = False

  def get_conformer_indices (self) :
    n_seq = self.atoms_size()
    conformer_indices = flex.size_t(n_seq, 0)
    altloc_indices = self.altloc_indices()
    if ("" in altloc_indices): p = 0
    else:                      p = 1
    altlocs = sorted(altloc_indices.keys())
    for i,altloc in enumerate(altlocs):
      if (altloc == ""): continue
      conformer_indices.set_selected(altloc_indices[altloc], i+p)
    return conformer_indices

  def transfer_chains_from_other(self, other):
    from iotbx.pdb import hy36encode
    i_model = 0
    other_models = other.models()
    for md,other_md in zip(self.models(), other_models):
      i_model += 1
      md.id = hy36encode(width=4, value=i_model)
      md.transfer_chains_from_other(other=other_md)
    msz, omsz = self.models_size(), other.models_size()
    if (omsz > msz):
      for other_md in other_models[msz:]:
        i_model += 1
        md = model(id = hy36encode(width=4, value=i_model))
        md.transfer_chains_from_other(other=other_md)
        self.append_model(model=md)

  def atom_selection_cache(self):
    from iotbx.pdb.atom_selection import cache
    return cache(root=self)

  def occupancy_groups_simple(self, common_residue_name_class_only=None,
                              ignore_hydrogens=True):
    if(ignore_hydrogens):
      sentinel = self.atoms().reset_tmp_for_occupancy_groups_simple()
    else:
      sentinel = self.atoms().reset_tmp(first_value=0, increment=1)
    result = []
    for chain in self.chains():
      result.extend(chain.occupancy_groups_simple(
        common_residue_name_class_only=common_residue_name_class_only))
    del sentinel
    return result

  def distance_based_simple_two_way_bond_sets(self,
        fallback_expected_bond_length=1.4,
        fallback_search_max_distance=2.5) :
    from cctbx.crystal import distance_based_connectivity
    atoms = self.atoms().deep_copy() # XXX potential bottleneck
    atoms.set_chemical_element_simple_if_necessary()
    sites_cart = atoms.extract_xyz()
    elements = atoms.extract_element()
    conformer_indices = self.get_conformer_indices()
    return distance_based_connectivity.build_simple_two_way_bond_sets(
      sites_cart=sites_cart,
      elements=elements,
      conformer_indices=conformer_indices,
      fallback_expected_bond_length=fallback_expected_bond_length,
      fallback_search_max_distance=fallback_search_max_distance)

class _(boost.python.injector, ext.model, __hash_eq_mixin):

  def residue_groups(self):
    for chain in self.chains():
      for rg in chain.residue_groups():
        yield rg

  def atom_groups(self):
    for chain in self.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          yield ag

  def only_chain(self):
    assert self.chains_size() == 1
    return self.chains()[0]

  def only_residue_group(self):
    return self.only_chain().only_residue_group()

  def only_conformer(self):
    return self.only_chain().only_conformer()

  def only_atom_group(self):
    return self.only_residue_group().only_atom_group()

  def only_residue(self):
    return self.only_conformer().only_residue()

  def only_atom(self):
    return self.only_atom_group().only_atom()

class _(boost.python.injector, ext.chain, __hash_eq_mixin):

  def atom_groups(self):
    for rg in self.residue_groups():
      for ag in rg.atom_groups():
        yield ag

  def only_residue_group(self):
    assert self.residue_groups_size() == 1
    return self.residue_groups()[0]

  def only_conformer(self):
    conformers = self.conformers()
    assert len(conformers) == 1
    return conformers[0]

  def only_atom_group(self):
    return self.only_residue_group().only_atom_group()

  def only_residue(self):
    return self.only_conformer().only_residue()

  def only_atom(self):
    return self.only_atom_group().only_atom()

  def residues(self):
    return self.only_conformer().residues()

  def occupancy_groups_simple(self, common_residue_name_class_only=None):
    result = []
    residue_groups = self.residue_groups()
    n_rg = len(residue_groups)
    done = [False] * n_rg
    def process_range(i_begin, i_end):
      isolated_var_occ = []
      groups = {}
      for i_rg in xrange(i_begin, i_end):
        done[i_rg] = True
        for ag in residue_groups[i_rg].atom_groups():
          altloc = ag.altloc
          if (altloc == ""):
            for atom in ag.atoms():
              if (atom.tmp < 0): continue
              if (atom.occ > 0 and atom.occ < 1):
                isolated_var_occ.append(atom.tmp)
          else:
            group = []
            for atom in ag.atoms():
              if (atom.tmp < 0): continue
              group.append(atom.tmp)
            if (len(group) != 0):
              groups.setdefault(altloc, []).extend(group)
      groups = groups.values()
      if (len(groups) != 0):
        for group in groups: group.sort()
        def group_cmp(a, b): return cmp(a[0], b[0])
        groups.sort(group_cmp)
        result.append(groups)
      for i in isolated_var_occ:
        result.append([[i]])
    for i_begin,i_end in self.find_pure_altloc_ranges(
          common_residue_name_class_only=common_residue_name_class_only):
      process_range(i_begin, i_end)
    for i_rg in xrange(n_rg):
      if (done[i_rg]): continue
      process_range(i_rg, i_rg+1)
    def groups_cmp(a, b):
      return cmp(a[0][0], b[0][0])
    result.sort(groups_cmp)
    return result

class _(boost.python.injector, ext.residue_group, __hash_eq_mixin):

  def only_atom_group(self):
    assert self.atom_groups_size() == 1
    return self.atom_groups()[0]

  def only_atom(self):
    return self.only_atom_group().only_atom()

class _(boost.python.injector, ext.atom_group, __hash_eq_mixin):

  def only_atom(self):
    assert self.atoms_size() == 1
    return self.atoms()[0]

class _(boost.python.injector, ext.atom, __hash_eq_mixin):

  def is_in_same_conformer_as(self, other):
    ag_i = self.parent(optional=False)
    ag_j = other.parent(optional=False)
    altloc_i = ag_i.altloc
    altloc_j = ag_j.altloc
    if (    len(altloc_i) != 0
        and len(altloc_j) != 0
        and altloc_i != altloc_j):
      return False
    def p3(ag):
      return ag.parent(optional=False) \
               .parent(optional=False) \
               .parent(optional=False)
    model_i = p3(ag_i)
    model_j = p3(ag_j)
    return model_i.memory_id() == model_j.memory_id()

  def set_element_and_charge_from_scattering_type_if_necessary(self,
        scattering_type):
    from cctbx.eltbx.xray_scattering \
      import get_element_and_charge_symbols \
        as gec
    sct_e, sct_c = gec(scattering_type=scattering_type, exact=False)
    pdb_ec = self.element.strip() + self.charge.strip()
    if (len(pdb_ec) != 0):
      if (sct_e == "" and sct_c == ""):
        return False
      pdb_e, pdb_c = gec(scattering_type=pdb_ec, exact=False)
      if (    pdb_e == sct_e
          and pdb_c == sct_c):
        return False
    self.element = "%2s" % sct_e.upper()
    self.charge = "%-2s" % sct_c
    return True

class _(boost.python.injector, ext.conformer):

  def only_residue(self):
    residues = self.residues()
    assert len(residues) == 1
    return residues[0]

  def only_atom(self):
    return self.only_residue().only_atom()

  def get_residue_names_and_classes (self) :
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
    return (rn_seq, residue_classes)

  def is_protein (self, min_content=0.8) :
    rn_seq, residue_classes = self.get_residue_names_and_classes()
    n_aa = residue_classes["common_amino_acid"]
    n_na = residue_classes["common_rna_dna"]
    if ((n_aa > n_na) and ((n_aa / len(rn_seq)) >= min_content)) :
      return True
    return False

  def is_na (self, min_content=0.8) :
    rn_seq, residue_classes = self.get_residue_names_and_classes()
    n_aa = residue_classes["common_amino_acid"]
    n_na = residue_classes["common_rna_dna"]
    if ((n_na > n_aa) and ((n_na / len(rn_seq)) >= min_content)) :
      return True
    return False

  def as_sequence (self) :
    rn_seq, residue_classes = self.get_residue_names_and_classes()
    n_aa = residue_classes["common_amino_acid"]
    n_na = residue_classes["common_rna_dna"]
    seq = []
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
    return seq

  def format_fasta(self, max_line_length=79):
    seq = self.as_sequence()
    n = len(seq)
    if (n == 0): return None
    comment = [">"]
    p = self.parent()
    if (p is not None):
      comment.append('chain "%2s"' % p.id)
    comment.append('conformer "%s"' % self.altloc)
    result = [" ".join(comment)]
    i = 0
    while True:
      j = min(n, i+max_line_length)
      if (j == i): break
      result.append("".join(seq[i:j]))
      i = j
    return result

  def as_padded_sequence (self, missing_char='X', skip_insertions=False,
      pad=True) :
    seq = self.as_sequence()
    padded_seq = []
    last_resseq = 0
    last_icode = " "
    i = 0
    for i, residue in enumerate(self.residues()) :
      if (skip_insertions) and (residue.icode != " ") :
        continue
      resseq = residue.resseq_as_int()
      if (pad) and (resseq > (last_resseq + 1)) :
        for x in range(resseq - last_resseq - 1) :
          padded_seq.append(missing_char)
      last_resseq = resseq
      padded_seq.append(seq[i])
    return "".join(padded_seq)

  def as_sec_str_sequence (self, helix_sele, sheet_sele, missing_char='X') :
    ss_seq = []
    last_resseq = 0
    for i, residue in enumerate(self.residues()) :
      resseq = residue.resseq_as_int()
      if resseq < 1 :
        continue
      if resseq > (last_resseq + 1) :
        for x in range(resseq - last_resseq - 1) :
          ss_seq.append(missing_char)
      found = False
      for atom in residue.atoms() :
        if helix_sele[atom.i_seq] :
          ss_seq.append('H')
          found = True
          break
        elif sheet_sele[atom.i_seq] :
          ss_seq.append('S')
          found = True
          break
      if not found :
        ss_seq.append('L')
      last_resseq = resseq
    return "".join(ss_seq)

class _(boost.python.injector, ext.residue):

  def __getinitargs__(self):
    result_root = self.root()
    if (result_root is None):
      orig_conformer = self.parent()
      assert orig_conformer is not None
      orig_chain = orig_conformer.parent()
      assert orig_chain is not None
      orig_model = orig_chain.parent()
      assert orig_model is not None
      result_atom_group = atom_group(
        altloc=orig_conformer.altloc, resname=self.resname)
      result_residue_group = residue_group(
        resseq=self.resseq, icode=self.icode)
      result_chain = chain(id=orig_chain.id)
      result_model = model(id=orig_model.id)
      result_root = root()
      result_root.append_model(result_model)
      result_model.append_chain(result_chain)
      result_chain.append_residue_group(result_residue_group)
      result_residue_group.append_atom_group(result_atom_group)
      for atom in self.atoms():
        result_atom_group.append_atom(atom.detached_copy())
    return (result_root,)

  def standalone_copy(self):
    return residue(root=self.__getinitargs__()[0])

  def only_atom(self):
    assert self.atoms_size() == 1
    return self.atoms()[0]

  def residue_name_plus_atom_names_interpreter(self,
        translate_cns_dna_rna_residue_names=None,
        return_mon_lib_dna_name=False):
    from iotbx.pdb import residue_name_plus_atom_names_interpreter
    return residue_name_plus_atom_names_interpreter(
      residue_name=self.resname,
      atom_names=[atom.name for atom in self.atoms()],
      translate_cns_dna_rna_residue_names=translate_cns_dna_rna_residue_names,
      return_mon_lib_dna_name=return_mon_lib_dna_name)

class input(object):

  def __init__(self, file_name=None, pdb_string=None, source_info=Auto):
    assert [file_name, pdb_string].count(None) == 1
    import iotbx.pdb
    if (file_name is not None):
      assert source_info is Auto
      self.input = iotbx.pdb.input(file_name=file_name)
    else:
      if (source_info is Auto): source_info = "string"
      self.input = iotbx.pdb.input(
        source_info=source_info, lines=flex.split_lines(pdb_string))
    self.hierarchy = self.input.construct_hierarchy()

  def __getinitargs__(self):
    from pickle import PicklingError
    raise PicklingError

  def hierarchy_to_input_atom_permutation(self):
    h_atoms = self.hierarchy.atoms()
    sentinel = h_atoms.reset_tmp(first_value=0, increment=1)
    return self.input.atoms().extract_tmp_as_size_t()

  def input_to_hierarchy_atom_permutation(self):
    i_atoms = self.input.atoms()
    sentinel = i_atoms.reset_tmp(first_value=0, increment=1)
    return self.hierarchy.atoms().extract_tmp_as_size_t()

  def xray_structure_simple (self, *args, **kwds) :
    perm = self.input_to_hierarchy_atom_permutation()
    xrs = self.input.xray_structure_simple(*args, **kwds)
    return xrs.select(perm)

class show_summary(input):

  def __init__(self,
        file_name=None,
        pdb_string=None,
        out=None,
        prefix="",
        flag_errors=True,
        flag_warnings=True,
        residue_groups_max_show=10,
        duplicate_atom_labels_max_show=10,
        level_id=None,
        level_id_exception=ValueError):
    input.__init__(self, file_name=file_name, pdb_string=pdb_string)
    print >> out, prefix+self.input.source_info()
    self.overall_counts = self.hierarchy.overall_counts()
    self.overall_counts.show(
      out=out,
      prefix=prefix+"  ",
      residue_groups_max_show=residue_groups_max_show,
      duplicate_atom_labels_max_show=duplicate_atom_labels_max_show)
    if (level_id is not None):
      self.hierarchy.show(
        out=out,
        prefix=prefix+"  ",
        level_id=level_id,
        level_id_exception=level_id_exception)

def append_chain_id_suffixes(roots, suffixes=Auto):
  if (suffixes is Auto):
    suffixes="123456789" \
             "ABCDEFGHIJKLMNOPQRSTUVWXYZ" \
             "abcdefghijklmnopqrstuvwxyz"
  assert len(roots) <= len(suffixes)
  for root,suffix in zip(roots, suffixes):
    for model in root.models():
      for chain in model.chains():
        assert len(chain.id) == 1
        chain.id += suffix

def join_roots(roots, chain_id_suffixes=Auto):
  if (chain_id_suffixes is not None):
    append_chain_id_suffixes(roots=roots, suffixes=chain_id_suffixes)
  result = root()
  for rt in roots:
    result.transfer_chains_from_other(other=rt)
  return result

# XXX: Nat's utility functions
def new_hierarchy_from_chain (chain) :
  import iotbx.pdb.hierarchy
  hierarchy = iotbx.pdb.hierarchy.root()
  model = iotbx.pdb.hierarchy.model()
  model.append_chain(chain.detached_copy())
  hierarchy.append_model(model)
  return hierarchy

# XXX: this will only replace the *first* chain it finds with an identical
# model/ID combination, for situations where water molecules are given the
# same ID as the nearest macromolecule.
def find_and_replace_chains (original_hierarchy, partial_hierarchy,
    log=sys.stdout) :
  for original_model in original_hierarchy.models() :
    for partial_model in partial_hierarchy.models() :
      if original_model.id == partial_model.id :
        #print >> log, "    found model '%s'" % partial_model.id
        i = 0
        while i < len(original_model.chains()) :
          original_chain = original_model.chains()[i]
          j = 0
          while j < len(partial_model.chains()) :
            partial_chain = partial_model.chains()[j]
            if original_chain.id == partial_chain.id :
              #print >> log, "      found chain '%s' at index %d" % (
              #  partial_chain.id, i)
              original_model.remove_chain(i)
              original_model.insert_chain(i, partial_chain.detached_copy())
              partial_model.remove_chain(j)
              break
            j += 1
          i += 1

# used for reporting build results in phenix
def get_residue_and_fragment_count (pdb_file=None, pdb_hierarchy=None) :
  from iotbx.pdb import amino_acid_codes
  import iotbx.pdb
  from libtbx import smart_open
  covalent_residues = amino_acid_codes.one_letter_given_three_letter.keys()
  covalent_residues.extend(iotbx.pdb.cns_dna_rna_residue_names.keys())
  if (pdb_file is not None) :
    raw_records = flex.std_string()
    f = smart_open.for_reading(file_name=pdb_file)
    raw_records.extend(flex.split_lines(f.read()))
    pdb_in = iotbx.pdb.input(source_info=pdb_file, lines=raw_records)
    pdb_hierarchy = pdb_in.construct_hierarchy()
  assert (pdb_hierarchy is not None)
  models = pdb_hierarchy.models()
  if len(models) == 0 :
    return (0, 0, 0)
  chains = models[0].chains()
  if len(chains) == 0 :
    return (0, 0, 0)
  n_res = 0
  n_frag = 0
  n_h2o = 0
  for chain in chains :
    i = -999
    for res in chain.conformers()[0].residues() :
      if res.resname.strip() in covalent_residues :
        n_res += 1
        resseq = res.resseq_as_int()
        if resseq > (i + 1) :
          n_frag += 1
        i = resseq
      elif res.resname.strip() in ['HOH','WAT','H2O'] :
        n_h2o += 1
  return (n_res, n_frag, n_h2o)
