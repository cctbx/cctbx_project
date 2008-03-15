from __future__ import generators

import boost.python
ext = boost.python.import_ext("iotbx_pdb_hierarchy_v2_ext")
from iotbx_pdb_hierarchy_v2_ext import *

from libtbx.str_utils import show_sorted_by_counts
from libtbx.utils import Sorry, plural_s
from cStringIO import StringIO
import sys

level_ids = ["model", "chain", "residue_group", "atom_group", "atom"]

class overall_counts(object):

  def show(self,
        out=None,
        prefix="",
        consecutive_residue_groups_max_show=10,
        duplicate_atom_labels_max_show=10):
    from iotbx.pdb import common_residue_names_get_class
    if (out is None): out = sys.stdout
    fmt = "%%%dd" % len(str(self.n_atoms))
    print >> out, prefix+"total number of:"
    print >> out, prefix+"  models:    ", fmt % self.n_models,
    if (self.n_duplicate_model_ids != 0):
      print >> out, "(%d with duplicate model id%s)" % plural_s(
        self.n_duplicate_model_ids),
    print >> out
    print >> out, prefix+"  chains:    ", fmt % self.n_chains,
    problems = []
    if (self.n_duplicate_chain_ids != 0):
      problems.append("%d with duplicate chain id%s" % plural_s(
        self.n_duplicate_chain_ids))
    if (self.n_explicit_chain_breaks != 0):
      problems.append("%d explicit chain break%s" % plural_s(
        self.n_explicit_chain_breaks))
    if (len(problems) != 0):
      print >> out, "(%s)" % "; ".join(problems),
    print >> out
    print >> out, prefix+"  alt. conf.:", fmt % self.n_alt_conf
    print >> out, prefix+"  residues:  ", fmt % (
      self.n_residues + self.n_residue_groups),
    if (self.n_residue_groups != 0):
      print >> out, "(%d with mixed residue names)" % self.n_residue_groups,
    print >> out
    print >> out, prefix+"  atoms:     ", fmt % self.n_atoms,
    if (self.n_duplicate_atom_labels != 0):
      print >> out, "(%d with duplicate labels)" %self.n_duplicate_atom_labels,
    print >> out
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
      print >> out, prefix+"  improper alt. conf.:", \
        fmt % self.n_alt_conf_improper
      if (self.n_alt_conf_improper != 0):
        for residue_group,label in [(self.alt_conf_proper, "proper"),
                                    (self.alt_conf_improper, "improper")]:
          if (residue_group is None): continue
          print >> out, prefix+"residue with %s altloc" % label
          for ag in residue_group.atom_groups():
            for atom in ag.atoms():
              print >> out, prefix+"  "+atom.format_atom_record()[:27].rstrip()
      print >> out, \
        prefix+"chains with mix of proper and improper alt. conf.:", \
        self.n_chains_with_mix_of_proper_and_improper_alt_conf
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
    self.show_consecutive_residue_groups_with_same_resid(
      out=out, prefix=prefix, max_show=consecutive_residue_groups_max_show)
    #
    msg = self.have_duplicate_atom_labels_message(
      max_show=duplicate_atom_labels_max_show,
      prefix=prefix)
    if (msg is not None): print >> out, msg

  def as_str(self,
        prefix="",
        consecutive_residue_groups_max_show=10,
        duplicate_atom_labels_max_show=10):
    out = StringIO()
    self.show(
      out=out,
      prefix=prefix,
      consecutive_residue_groups_max_show=consecutive_residue_groups_max_show,
      duplicate_atom_labels_max_show=duplicate_atom_labels_max_show)
    return out.getvalue()

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
    if (max_show <= 0): return
    delim = prefix+"  "+"-"*31
    prev_rg = None
    for rgs in cons[:max_show]:
      for next,rg in zip(["", "next "], rgs):
        if (rg is prev_rg): continue
        elif (next == "" and prev_rg is not None):
          print >> out, delim
        prev_rg = rg
        print >> out, prefix+"  %sresidue group:" % next
        atoms = rg.atoms()
        if (atoms.size() == 0):
          ch = rg.parent()
          if (ch is None): ch = "  "
          else:            ch = "%s" % ch.id
          print >> out, prefix+'    empty: "%s%s"' % (ch, rg.resid())
        else:
          def show_atom(atom):
            print >> out, prefix+'    "%s"' % atom.format_atom_record()[:27]
          if (atoms.size() <= 3):
            for atom in atoms: show_atom(atom)
          else:
            show_atom(atoms[0])
            print >> out, prefix+'    ... %d atom%s not shown' % plural_s(
              atoms.size()-2)
            show_atom(atoms[-1])
    if (len(cons) > max_show):
      print >> out, delim
      print >> out, prefix + "  ... %d remaining instance%s not shown" % \
        plural_s(len(cons)-max_show)

  def have_duplicate_atom_labels_message(self, max_show=10, prefix=""):
    dup = self.duplicate_atom_labels
    if (len(dup) == 0):
      return None
    fmt = "%%%dd" % len(str(self.n_duplicate_atom_labels))
    result = [
      prefix + "number of groups of duplicate atom labels: " + fmt % len(dup),
      prefix + "  total number of affected atoms:          " + fmt %
        self.n_duplicate_atom_labels]
    if (max_show > 0):
      for atoms in dup[:max_show]:
        prfx = "  group "
        for atom in atoms:
          result.append(
            prefix + prfx + '"' + atom.format_atom_record()[:27] + '"')
          prfx = "        "
      if (len(dup) > max_show):
        result.append(prefix + "  ... %d remaining group%s not shown" %
          plural_s(len(dup)-max_show))
    return "\n".join(result)

  def raise_duplicate_atom_labels_if_necessary(self, max_show=10):
    msg = self.have_duplicate_atom_labels_message(max_show=max_show)
    if (msg is not None): raise Sorry(msg)

class _root(boost.python.injector, ext.root):

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

  def only_atom_group(self):
    return self.only_residue_group().only_atom_group()

  def only_atom(self):
    return self.only_atom_group().only_atom()

  def overall_counts(self):
    from iotbx.pdb import common_residue_names_get_class
    from libtbx import dict_with_default_0
    blank_altloc_char = " "
    n_models = self.models_size()
    n_residues = 0
    n_residue_groups = 0
    n_explicit_chain_breaks = 0
    chain_ids = dict_with_default_0()
    alt_conf_ids = dict_with_default_0()
    resnames = dict_with_default_0()
    element_charge_types = dict_with_default_0()
    n_alt_conf_none = 0
    n_alt_conf_pure = 0
    n_alt_conf_proper = 0
    n_alt_conf_improper = 0
    alt_conf_proper = None
    alt_conf_improper = None
    consecutive_residue_groups_with_same_resid = []
    n_chains_with_mix_of_proper_and_improper_alt_conf = 0
    n_duplicate_model_ids = 0
    n_duplicate_chain_ids = 0
    n_duplicate_atom_labels = 0
    duplicate_atom_labels = []
    model_ids = dict_with_default_0()
    atoms = self.atoms()
    atoms.reset_tmp()
    for model in self.models():
      model_ids[model.id] += 1
      model_chain_ids = dict_with_default_0()
      model_atom_labels_i_seqs = {}
      prev_rg = None
      for chain in model.chains():
        model_chain_ids[chain.id] += 1
        chain_ids[chain.id] += 1
        chain_altlocs = {} # FUTURE: set
        chain_alt_conf_proper = None
        chain_alt_conf_improper = None
        suppress_chain_break = True
        for rg in chain.residue_groups():
          if (not rg.link_to_previous and not suppress_chain_break):
            n_explicit_chain_breaks += 1
          suppress_chain_break = False
          have_main_conf = False
          have_blank_altloc = False
          rg_altlocs = {} # FUTURE: set
          rg_resnames = {} # FUTURE: set
          for ag in rg.atom_groups():
            if (ag.altloc == ""):
              have_main_conf = True
            else:
              if (ag.altloc == blank_altloc_char):
                have_blank_altloc = True
              rg_altlocs[ag.altloc] = None
            rg_resnames[ag.resname] = None
            for atom in ag.atoms():
              model_atom_labels_i_seqs.setdefault(
                atom.pdb_label_columns(), []).append(atom.tmp)
              element_charge_types[
                "%2s%-2s" % (atom.element, atom.charge)] += 1
              rg_last_atom = atom
          if (have_blank_altloc):
            n_alt_conf_improper += 1
            if (chain_alt_conf_improper is None):
              chain_alt_conf_improper = rg
          elif (have_main_conf):
            if (len(rg_altlocs) == 0):
              n_alt_conf_none += 1
            else:
              n_alt_conf_proper += 1
              if (chain_alt_conf_proper is None):
                chain_alt_conf_proper = rg
          elif (len(rg_altlocs) != 0):
            n_alt_conf_pure += 1
          chain_altlocs.update(rg_altlocs)
          if (len(rg_resnames) == 1):
            n_residues += 1
          else:
            n_residue_groups += 1
          for resname in rg_resnames:
            resnames[resname] += 1
          if (    prev_rg is not None
              and prev_rg.resid() == rg.resid()):
            consecutive_residue_groups_with_same_resid.append((prev_rg, rg))
          prev_rg = rg
        for altloc in chain_altlocs:
          alt_conf_ids[altloc] += 1
        if (    chain_alt_conf_proper is not None
            and chain_alt_conf_improper is not None):
          if (   alt_conf_proper is None
              or alt_conf_improper is None):
            alt_conf_proper = chain_alt_conf_proper
            alt_conf_improper = chain_alt_conf_improper
        else:
          if (alt_conf_proper is None):
            alt_conf_proper = chain_alt_conf_proper
          if (alt_conf_improper is None):
            alt_conf_improper = chain_alt_conf_improper
      for chain_id,count in model_chain_ids.items():
        if (count != 1): n_duplicate_chain_ids += count
      model_duplicate_atom_labels = []
      for i_seqs in model_atom_labels_i_seqs.values():
        if (len(i_seqs) != 1):
          model_duplicate_atom_labels.append(i_seqs)
      model_duplicate_atom_labels.sort()
      for i_seqs in model_duplicate_atom_labels:
        n_duplicate_atom_labels += len(i_seqs)
        duplicate_atom_labels.append([atoms[i_seq] for i_seq in i_seqs])
    for model_id,count in model_ids.items():
      if (count != 1): n_duplicate_model_ids += count
    resname_classes = dict_with_default_0()
    for resname,count in resnames.items():
      resname_classes[common_residue_names_get_class(name=resname)] += count
    result = overall_counts()
    result.root = self
    result.n_duplicate_model_ids = n_duplicate_model_ids
    result.n_duplicate_chain_ids = n_duplicate_chain_ids
    result.n_duplicate_atom_labels = n_duplicate_atom_labels
    result.duplicate_atom_labels = duplicate_atom_labels
    result.n_models = n_models
    result.n_chains = sum(chain_ids.values())
    result.n_alt_conf = sum(alt_conf_ids.values())
    result.n_residues = n_residues
    result.n_residue_groups = n_residue_groups
    result.n_explicit_chain_breaks = n_explicit_chain_breaks
    result.n_atoms = sum(element_charge_types.values())
    result.chain_ids = chain_ids
    result.alt_conf_ids = alt_conf_ids
    result.resnames = resnames
    result.resname_classes = resname_classes
    result.element_charge_types = element_charge_types
    result.n_alt_conf_none = n_alt_conf_none
    result.n_alt_conf_pure = n_alt_conf_pure
    result.n_alt_conf_proper = n_alt_conf_proper
    result.n_alt_conf_improper = n_alt_conf_improper
    result.alt_conf_proper = alt_conf_proper
    result.alt_conf_improper = alt_conf_improper
    result.consecutive_residue_groups_with_same_resid \
         = consecutive_residue_groups_with_same_resid
    result.n_chains_with_mix_of_proper_and_improper_alt_conf \
         = n_chains_with_mix_of_proper_and_improper_alt_conf
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
    for model in self.models():
      chains = model.chains()
      print >> out, prefix+'model id="%s"' % model.id, \
        "#chains=%d" % len(chains)
      if (level_no == 0): continue
      for chain in chains:
        rgs = chain.residue_groups()
        print >> out, prefix+'  chain id="%s"' % chain.id, \
          "#residue_groups=%d" % len(rgs)
        if (level_no == 1): continue
        suppress_chain_break = True
        for rg in rgs:
          if (not rg.link_to_previous and not suppress_chain_break):
            print >> out, prefix+"    ### chain break ###"
          suppress_chain_break = False
          ags = rg.atom_groups()
          resnames = {} # FUTURE: set
          for ag in rg.atom_groups():
            resnames[ag.resname] = None
          if (len(resnames) > 1): s = " *** with mixed residue names ***"
          else: s = ""
          print >> out, prefix+'    resid="%s"' % rg.resid(), \
            "#atom_groups=%d%s" % (len(ags), s)
          if (level_no == 2): continue
          for ag in ags:
            atoms = ag.atoms()
            print >> out, prefix+'      altloc="%s"' % ag.altloc, \
              'resname="%s"' % ag.resname, \
              "#atoms=%d" % len(atoms)
            if (level_no == 3): continue
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

  def as_pdb_records(self, append_end=False):
    result = []
    models = self.models()
    for model in models:
      if (len(models) != 1):
        result.append("MODEL %7s" % model.id)
      for chain in model.chains():
        atom_serial = chain.append_atom_records(
          pdb_records=result)
        result.append("TER")
      if (len(models) != 1):
        result.append("ENDMDL")
    if (append_end):
      result.append("END")
    return result

  def as_pdb_string(self, append_end=False):
    return "\n".join(self.as_pdb_records(append_end=append_end))+"\n"

  def occupancy_groups_simple(self, common_residue_name_class_only=None):
    self.atoms().reset_tmp_for_occupancy_groups_simple()
    result = []
    for chain in self.chains():
      result.extend(chain.occupancy_groups_simple(
        common_residue_name_class_only=common_residue_name_class_only))
    return result

class _model(boost.python.injector, ext.model):

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

  def only_atom_group(self):
    return self.only_residue_group().only_atom_group()

  def only_atom(self):
    return self.only_atom_group().only_atom()

class _chain(boost.python.injector, ext.chain):

  def atom_groups(self):
    for rg in self.residue_groups():
      for ag in rg.atom_groups():
        yield ag

  def only_residue_group(self):
    assert self.residue_groups_size() == 1
    return self.residue_groups()[0]

  def only_atom_group(self):
    return self.only_residue_group().only_atom_group()

  def only_atom(self):
    return self.only_atom_group().only_atom()

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

class _residue_group(boost.python.injector, ext.residue_group):

  def only_atom_group(self):
    assert self.atom_groups_size() == 1
    return self.atom_groups()[0]

  def only_atom(self):
    return self.only_atom_group().only_atom()

class _atom_group(boost.python.injector, ext.atom_group):

  def only_atom(self):
    assert self.atoms_size() == 1
    return self.atoms()[0]
