from __future__ import generators

import boost.python
ext = boost.python.import_ext("iotbx_pdb_hierarchy_v2_ext")
from iotbx_pdb_hierarchy_v2_ext import *

from cStringIO import StringIO

level_ids = ["model", "chain", "residue_group", "atom_group", "atom"]

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
          print >> out, prefix+'    resid="%s"' % rg.resid(), \
            "#atom_groups=%d" % len(ags)
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
