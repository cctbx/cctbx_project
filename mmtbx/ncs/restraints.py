from __future__ import division
from scitbx.array_family import flex

import boost.python
ext = boost.python.import_ext("mmtbx_ncs_restraints_ext")
from mmtbx_ncs_restraints_ext import *

from cctbx import adptbx
import scitbx.restraints
from scitbx.math import superpose
from scitbx import matrix
import scitbx.stl
from libtbx.str_utils import show_string
from libtbx.utils import Sorry
from libtbx import adopt_init_args
from itertools import count
import sys

class selection_properties(object):

  def __init__(self, string, iselection):
    self.string = string
    self.iselection = iselection

# XXX rename
def match_ordered_chain_members(a, b):
  result = []
  na = len(a)
  if (na == 0): return result
  nb = len(b)
  if (nb == 0): return result
  ia = 0
  ib = 0
  while True:
    ja = ia
    jb = ib
    while True:
      if (a[ia] == b[jb]):
        result.append((ia, jb))
        ib = jb
        break
      elif (jb != ib and b[ib] == a[ja]):
        result.append((ja, ib))
        ia = ja
        break
      else:
        ja += 1
        if (ja == na):
          while True:
            jb += 1
            if (jb == nb): break
            if (a[ia] == b[jb]):
              result.append((ia, jb))
              ib = jb
              break
          break
        jb += 1
        if (jb == nb):
          while True:
            if (b[ib] == a[ja]):
              result.append((ja, ib))
              ia = ja
              break
            ja += 1
            if (ja == na): break
          break
    ia += 1
    if (ia == na): break
    ib += 1
    if (ib == nb): break
  return result

class resconfid(object):

  def __init__(self, residue_group):
    self.residue_group = residue_group
    self.resid = residue_group.resid()
    confs = {}
    for ag in residue_group.atom_groups():
      altloc = ag.altloc
      assert not altloc in confs
      confs[altloc] = ag.resname
    self.confs = confs

  def __eq__(self, other):
    if (self.resid != other.resid): return False
    return self.confs_eq(other=other)

  def confs_eq(self, other):
    for altloc,resname in self.confs.items():
      other_resname = other.confs.get(altloc)
      if (other_resname is None): continue
      if (other_resname != resname): return False
    return True

  def __repr__(self): # XXX remove
    return str(self.resid) + ":" + str(self.confs)

def all_confs_eq(resconfids_i, resconfids_j):
  if (len(resconfids_i) != len(resconfids_j)): return False
  for i,j in zip(resconfids_i, resconfids_j):
    if (not i.confs_eq(other=j)): return False
  return True

def get_resconfids(model):
  result = []
  for rg in model.residue_groups():
    result.append(resconfid(residue_group=rg))
  return result

def match_atoms(
      special_position_indices,
      special_position_warnings_only,
      pdb_atoms,
      selection_strings,
      registry,
      reference_model,
      j_ncs,
      other_model,
      log):
  r_resconfids = get_resconfids(model=reference_model)
  o_resconfids = get_resconfids(model=other_model)
  if (all_confs_eq(resconfids_i=r_resconfids, resconfids_j=o_resconfids)):
    pairs_list = [(i,i) for i in xrange(len(r_resconfids))]
  else:
    pairs_list = match_ordered_chain_members(a=r_resconfids, b=o_resconfids)
  for i,j in pairs_list:
    r_rg = r_resconfids[i].residue_group
    o_rg = o_resconfids[j].residue_group
    r_ag_dict = {}
    for r_ag in r_rg.atom_groups():
      confid = r_ag.confid()
      assert not confid in r_ag_dict
      r_ag_dict[confid] = r_ag
    for o_ag in o_rg.atom_groups():
      r_ag = r_ag_dict.get(o_ag.confid())
      if (r_ag is None): continue
      r_atom_dict = {}
      for r_atom in r_ag.atoms():
        atom_name = r_atom.name
        assert atom_name not in r_atom_dict # duplicate atom names
        r_atom_dict[atom_name] = r_atom
      for o_atom in o_ag.atoms():
        r_atom = r_atom_dict.get(o_atom.name, o_atom)
        if (r_atom is o_atom): continue
        assert r_atom is not None # duplicate atom names
        r_atom_dict[r_atom.name] = None
        msg = None
        for i_ncs,atom in [(0,r_atom), (j_ncs,o_atom)]:
          if (atom.i_seq in special_position_indices):
            msg = (
                "NCS selection includes an atom on a special position:\n"
              + "  Selection: %s\n" % show_string(selection_strings[i_ncs])
              + '    %s' % atom.quote())
            if (not special_position_warnings_only):
              raise Sorry(msg)
            elif (log is not None):
              print >> log, "WARNING:", msg
        if (msg is not None): continue
        i_seq_i = r_atom.i_seq
        i_seq_j = o_atom.i_seq
        if (i_seq_i == i_seq_j):
          raise Sorry("NCS selections restrain atom to itself:\n"
            + "  Reference selection: %s\n" % show_string(selection_strings[0])
            + "      Other selection: %s\n" % show_string(selection_strings[j_ncs])
            + '    %s' % r_atom.quote())
        stat, i_diag = registry.enter(
          i_seq=i_seq_i, j_seq=i_seq_j, j_ncs=j_ncs)
        if (stat == 1):
          assert i_diag != j_ncs
          raise Sorry(
            "Two different NCS operators applied to same pair of atoms:\n"
            + "       Reference selection: %s\n" % show_string(
                selection_strings[0])
            + "  Previous other selection: %s\n" % show_string(
                selection_strings[i_diag])
            + "   Current other selection: %s\n" % show_string(
                selection_strings[j_ncs])
            + '    %s\n' % r_atom.quote()
            + '    %s' % o_atom.quote())
        if (stat == 2):
          raise Sorry(
            "Current reference atom previously other atom:\n"
            + "  Current reference selection: %s\n" % show_string(
                selection_strings[0])
            + "     Previous other selection: %s\n" % show_string(
                selection_strings[i_diag])
            + '    %s\n' % r_atom.quote())
        if (stat == 3):
          raise Sorry(
            "Two other atoms mapped to same reference atom:\n"
            + "  Reference selection: %s\n" % show_string(
                selection_strings[0])
            + "      Other selection: %s\n" % show_string(
                selection_strings[j_ncs])
            + '    Reference: %s\n' % r_atom.quote()
            + '      Other 1: %s\n' % pdb_atoms[i_diag].quote()
            + '      Other 2: %s' % o_atom.quote())
        assert stat == 0

class pair_lists_generator(object):

  def __init__(self,
        processed_pdb,
        reference_selection_string,
        selection_strings,
        special_position_warnings_only,
        log):
    self.processed_pdb = processed_pdb
    del processed_pdb
    self.reference_selection_string = reference_selection_string
    del reference_selection_string
    self.selection_strings = selection_strings
    del selection_strings
    if (self.reference_selection_string is not None):
      self.selection_strings = [self.reference_selection_string] \
                             + list(self.selection_strings)
    n_ncs = len(self.selection_strings)
    if (n_ncs < 2):
      raise Sorry("Only one NCS restraints selection: %s\n"
        "  At least two selections are required." %
          show_string(self.selection_strings[0]))
    # shortcuts
    all_chain_proxies = self.processed_pdb.all_chain_proxies
    self.n_seq = all_chain_proxies.pdb_atoms.size()
    #
    self.selection_properties = []
    for selection_string in self.selection_strings:
      selection = all_chain_proxies.selection(string=selection_string)
      iselection = selection.iselection()
      if (iselection.size() == 0):
        raise Sorry("Empty NCS restraints selection: %s" % (
          show_string(selection_string)))
      self.selection_properties.append(selection_properties(
        string=selection_string,
        iselection=iselection))
    #
    self.registry = pair_registry(n_seq=self.n_seq, n_ncs=n_ncs)
    special_position_indices = scitbx.stl.set.unsigned(iter(
      all_chain_proxies.site_symmetry_table().special_position_indices()))
    pdb_hierarchy = all_chain_proxies.pdb_hierarchy
    reference_hierarchy = pdb_hierarchy.select(
      atom_selection=self.selection_properties[0].iselection)
    for j_ncs in xrange(1,n_ncs):
      other_hierarchy = pdb_hierarchy.select(
        atom_selection=self.selection_properties[j_ncs].iselection)
      r_models = reference_hierarchy.models()
      o_models = other_hierarchy.models()
      if (len(r_models) != len(o_models)):
        raise Sorry(
          "NCS restraints selections yield different number of PDB MODELs:\n"
          + "  reference selection: %s\n" %
            show_string(self.selection_strings[0])
          + "      other selection: %s\n" %
            show_string(self.selection_strings[j_ncs])
          + "  number of reference MODELs: %d\n" % len(r_models)
          + "  number of     other MODELs: %d" % len(o_models))
      for r_model,o_model in zip(r_models, o_models):
        if (r_model.id != o_model.id):
          raise Sorry(
            "NCS restraints selections lead to PDB MODEL mismatches:\n"
            + "  reference selection: %s\n" %
              show_string(self.selection_strings[0])
            + "      other selection: %s\n" %
              show_string(self.selection_strings[j_ncs])
            + "  reference MODEL id: %s\n" % show_string(r_model.id)
            + "      other MODEL id: %s" % show_string(o_model.id))
        match_atoms(
          special_position_indices=special_position_indices,
          special_position_warnings_only=special_position_warnings_only,
          pdb_atoms=all_chain_proxies.pdb_atoms,
          selection_strings=self.selection_strings,
          registry=self.registry,
          reference_model=r_model,
          j_ncs=j_ncs,
          other_model=o_model,
          log=log)
    #
    for i_pair,pair in enumerate(self.registry.selection_pairs()):
      if (pair[0].size() < 2):
        detail = ["do not produce any pairs",
                  "produce only one pair"][pair[0].size()]
        raise Sorry("\n".join([
          "NCS restraints selections %s of matching atoms:" % detail,
          "  Reference selection: %s" % show_string(self.selection_strings[0]),
          "      Other selection: %s" % show_string(
            self.selection_strings[i_pair+1])]))

class group(object):

  def from_atom_selections(
        processed_pdb,
        reference_selection_string,
        selection_strings,
        coordinate_sigma,
        b_factor_weight,
        special_position_warnings_only,
        u_average_min=1.e-6,
        log=None):
    g = pair_lists_generator(
      processed_pdb=processed_pdb,
      reference_selection_string=reference_selection_string,
      selection_strings=selection_strings,
      special_position_warnings_only=special_position_warnings_only,
      log=log)
    return group(
      selection_strings=g.selection_strings,
      registry=g.registry,
      coordinate_sigma=coordinate_sigma,
      b_factor_weight=b_factor_weight,
      u_average_min=u_average_min)
  from_atom_selections = staticmethod(from_atom_selections)

  def __init__(self,
        selection_strings,
        registry,
        coordinate_sigma,
        b_factor_weight,
        u_average_min):
      adopt_init_args(self, locals())
      self.selection_pairs = registry.selection_pairs()

  def register_additional_isolated_sites(self, number):
    self.registry.register_additional_isolated_sites(number=number)

  def select(self, iselection):
    return group(
      selection_strings=self.selection_strings,
      registry=self.registry.proxy_select(iselection=iselection),
      coordinate_sigma=self.coordinate_sigma,
      b_factor_weight=self.b_factor_weight,
      u_average_min=self.u_average_min)

  def operators(self, sites_cart):
    return _operators(group=self, sites_cart=sites_cart)

  def energies_adp_iso(self,
        u_isos,
        average_power,
        compute_gradients=True,
        gradients=None):
    return _energies_adp_iso(
      group=self,
      u_isos=u_isos,
      average_power=average_power,
      compute_gradients=compute_gradients,
      gradients=gradients)

class _energies_adp_iso(scitbx.restraints.energies):

  def __init__(self,
        group,
        u_isos,
        average_power,
        compute_gradients,
        gradients):
    scitbx.restraints.energies.__init__(self,
      compute_gradients=compute_gradients,
      gradients=gradients,
      gradients_size=u_isos.size(),
      gradients_factory=flex.double,
      normalization=False)
    self.group = group
    self.u_isos = u_isos
    self.average_power = average_power
    # XXX registry should be condensed to map of std::vectors
    #     then all the code below and adp_iso_residual_sum could
    #     become a method of the reduced registry
    max_index = 0
    for pair in self.group.selection_pairs:
      max_index = max(max_index, pair[0][-1])
    u_isos_sum = flex.double(max_index+1, 0)
    u_isos_count = flex.size_t(max_index+1, 0)
    for pair in self.group.selection_pairs:
      u_isos_sum.set_selected(pair[0], u_isos.select(pair[0]))
      u_isos_count.set_selected(pair[0], 1)
    for pair in self.group.selection_pairs:
      u_isos_sum.add_selected(pair[0], u_isos.select(pair[1]))
      u_isos_count.set_selected(pair[0], u_isos_count.select(pair[0])+1)
    sel = u_isos_count == 0
    u_isos_count.set_selected(sel, 1)
    u_isos_average = u_isos_sum / u_isos_count.as_double()
    u_isos_count.set_selected(sel, 0)
    sel = (~sel).iselection()
    self.rms_with_respect_to_average = flex.double()
    def residual_contribution(u_isos_current, u_isos_average):
      diff = u_isos_current - u_isos_average
      self.rms_with_respect_to_average.append(adptbx.u_as_b(
        flex.mean_sq(diff)**0.5))
      self.number_of_restraints += diff.size()
    residual_contribution(
      u_isos_current=u_isos.select(sel),
      u_isos_average=u_isos_average.select(sel))
    for pair in self.group.selection_pairs:
      residual_contribution(
        u_isos_current=u_isos.select(pair[1]),
        u_isos_average=u_isos_average.select(pair[0]))
    self.weight = self.group.b_factor_weight
    self.residual_sum = self.group.registry.adp_iso_residual_sum(
      weight=self.weight,
      average_power=self.average_power,
      u_isos=u_isos,
      u_average_min=self.group.u_average_min,
      gradients=self.gradients)
    self.finalize_target_and_gradients()
    self.u_isos_count = u_isos_count
    self.u_isos_average = u_isos_average

  def show_differences_to_average(self, site_labels, out=None, prefix=""):
    if (out is None): out = sys.stdout
    assert len(site_labels) == self.u_isos.size()
    max_label_size = 1
    for label in site_labels:
      max_label_size = max(max_label_size, len(label))
    fmt = "  %%%ds: %%7.2f - %%7.2f = %%8.4f" % max_label_size
    def show_selection(i_ncs, pair):
      print >> out, prefix + "NCS selection:", \
        show_string(self.group.selection_strings[i_ncs])
      print >> out, prefix + " "*(max_label_size+2) \
        + "    B-iso   NCS ave  Difference"
      u_isos_current = self.u_isos.select(pair[1])
      u_isos_average = self.u_isos_average.select(pair[0])
      for i,c,a in zip(pair[1], u_isos_current, u_isos_average):
        c = adptbx.u_as_b(c)
        a = adptbx.u_as_b(a)
        print >> out, prefix + fmt % (site_labels[i], c, a, c-a)
    sel = (self.u_isos_count != 0).iselection()
    show_selection(i_ncs=0, pair=[sel, sel])
    for i_ncs,pair in zip(count(1), self.group.selection_pairs):
      show_selection(i_ncs=i_ncs, pair=pair)

class _operators(object):

  def __init__(self, group, sites_cart):
    self.group = group
    self.matrices = []
    self.rms = []
    for pair in self.group.selection_pairs:
      superposition = superpose.least_squares_fit(
        reference_sites=sites_cart.select(pair[0]),
        other_sites=sites_cart.select(pair[1]))
      rtmx = matrix.rt((superposition.r, superposition.t))
      self.matrices.append(rtmx)
      x = sites_cart.select(pair[0])
      y = rtmx * sites_cart.select(pair[1])
      d_sq = (x-y).dot()
      self.rms.append(flex.mean(d_sq)**0.5)

  def show(self,
        sites_cart,
        n_slots_difference_histogram=6,
        out=None,
        prefix=""):
    if (out is None): out = sys.stdout
    selection_strings = self.group.selection_strings
    for i_op,pair,mx,rms in zip(
          count(1),
          self.group.selection_pairs,
          self.matrices,
          self.rms):
      print >> out, prefix + "NCS operator %d:" % i_op
      print >> out, prefix + "  Reference selection:", \
        show_string(selection_strings[0])
      print >> out, prefix + "      Other selection:", \
        show_string(selection_strings[i_op])
      print >> out, prefix + "  Number of atom pairs:", len(pair[0])
      print >> out, mx.r.mathematica_form(
        label="Rotation", format="%.6g", one_row_per_line=True,
        prefix=prefix+"  ")
      print >> out, mx.t.mathematica_form(
        label="Translation", format="%.6g", prefix=prefix+"  ")
      x = sites_cart.select(pair[0])
      y = mx * sites_cart.select(pair[1])
      d_sq = (x-y).dot()
      if (n_slots_difference_histogram is not None):
        print >> out, prefix + "  Histogram of differences:"
        diff_histogram = flex.histogram(
          data=flex.sqrt(d_sq), n_slots=n_slots_difference_histogram)
        diff_histogram.show(
          f=out, prefix=prefix+"    ", format_cutoffs="%8.6f")
      print >> out, \
        prefix + "  RMS difference with respect to the reference: %8.6f" %(rms)

  def energies_sites(self,
        sites_cart,
        compute_gradients=True,
        gradients=None,
        sites_average=None):
    return _energies_sites(
      operators=self,
      sites_cart=sites_cart,
      compute_gradients=compute_gradients,
      gradients=gradients,
      sites_average=sites_average)

class _energies_sites(scitbx.restraints.energies):

  def __init__(self,
        operators,
        sites_cart,
        compute_gradients,
        gradients,
        sites_average):
    scitbx.restraints.energies.__init__(self,
      compute_gradients=compute_gradients,
      gradients=gradients,
      gradients_size=sites_cart.size(),
      gradients_factory=flex.vec3_double,
      normalization=False)
    self.operators = operators
    self.sites_cart = sites_cart
    selection_pairs = self.operators.group.selection_pairs
    max_index = 0
    for pair in selection_pairs:
      max_index = max(max_index, pair[0][-1])
    sites_sum = flex.vec3_double(max_index+1, (0,0,0))
    sites_count = flex.size_t(max_index+1, 0)
    for pair in selection_pairs:
      sites_sum.set_selected(pair[0], sites_cart.select(pair[0]))
      sites_count.set_selected(pair[0], 1)
    for pair,op in zip(selection_pairs, self.operators.matrices):
      sites_sum.add_selected(pair[0], op*sites_cart.select(pair[1]))
      sites_count.set_selected(pair[0], sites_count.select(pair[0])+1)
    sel = sites_count == 0
    sites_count.set_selected(sel, 1)
    if (sites_average is not None):
      assert sites_average.size() == sites_count.size()
    else:
      sites_average = sites_sum / sites_count.as_double()
    sites_count.set_selected(sel, 0)
    sel = (~sel).iselection()
    assert self.operators.group.coordinate_sigma > 0
    self.weight = 1/self.operators.group.coordinate_sigma**2
    self.rms_with_respect_to_average = flex.double()
    if (self.gradients is not None):
      self.gradients.add_selected(
        sel,
        self.residual_contribution(
          sites_current=sites_cart.select(sel),
          sites_average=sites_average.select(sel)))
      for pair,op in zip(selection_pairs, self.operators.matrices):
        self.gradients.add_selected(
          pair[1],
          self.residual_contribution(
            sites_current=sites_cart.select(pair[1]),
            sites_average=op.inverse() * sites_average.select(pair[0])))
    else:
      self.residual_contribution(
        sites_current=sites_cart.select(sel),
        sites_average=sites_average.select(sel))
      for pair,op in zip(selection_pairs, self.operators.matrices):
        self.residual_contribution(
          sites_current=sites_cart.select(pair[1]),
          sites_average=op.inverse() * sites_average.select(pair[0]))
    self.finalize_target_and_gradients()
    self.sites_count = sites_count
    self.sites_average = sites_average

  def residual_contribution(self, sites_current, sites_average):
    diff = sites_current - sites_average
    self.rms_with_respect_to_average.append(flex.mean(diff.dot())**0.5)
    self.number_of_restraints += diff.size()
    self.residual_sum += self.weight * diff.sum_sq()
    if (self.gradients is not None):
      return (2 * self.weight) * diff

  def show_distances_to_average(self,
        site_labels,
        excessive_distance_limit=None,
        out=None,
        prefix=""):
    if (out is None): out = sys.stdout
    assert len(site_labels) == self.sites_cart.size()
    max_label_size = 1
    for label in site_labels:
      max_label_size = max(max_label_size, len(label))
    fmt = "  %%%ds: %%8.4f" % max_label_size
    def show_selection(i_ncs, pair, op):
      print >> out, prefix + "NCS selection:", \
        show_string(self.operators.group.selection_strings[i_ncs])
      print >> out, prefix + " "*(max_label_size+2) \
        + "  Distance to NCS average"
      sites_current = self.sites_cart.select(pair[1])
      if (op is None):
        sites_average = self.sites_average.select(pair[0])
      else:
        sites_average = op.inverse() * self.sites_average.select(pair[0])
      n_excessive = 0
      for i,c,a in zip(pair[1], sites_current, sites_average):
        abs_diff = abs(matrix.col(c) - matrix.col(a))
        print >> out, prefix + fmt % (site_labels[i], abs_diff),
        if (excessive_distance_limit is not None
            and abs_diff >= excessive_distance_limit):
          print >> out, "EXCESSIVE",
          n_excessive += 1
        print >> out
      return n_excessive
    sel = (self.sites_count != 0).iselection()
    n_excessive = show_selection(i_ncs=0, pair=[sel, sel], op=None)
    for i_ncs,pair,op in zip(count(1),
                             self.operators.group.selection_pairs,
                             self.operators.matrices):
      n_excessive += show_selection(i_ncs=i_ncs, pair=pair, op=op)
    return n_excessive

class groups(object):

  def __init__(self, members=None):
    if (members is None):
      self.members = []
    else:
      self.members = members
    self.operators = None

  def register_additional_isolated_sites(self, number):
    for group in self.members:
      group.register_additional_isolated_sites(number=number)

  def select(self, iselection):
    members = []
    for group in self.members:
      members.append(group.select(iselection=iselection))
    return groups(members=members)

  def energies_adp_iso(self,
        u_isos,
        average_power,
        compute_gradients=True,
        gradients=None,
        normalization=False):
    result = scitbx.restraints.energies(
      compute_gradients=compute_gradients,
      gradients=gradients,
      gradients_size=u_isos.size(),
      gradients_factory=flex.double,
      normalization=normalization)
    result.rms_with_respect_to_averages = []
    for group in self.members:
      if (    group.b_factor_weight is not None
          and group.b_factor_weight > 0):
        contribution = group.energies_adp_iso(
          u_isos=u_isos,
          average_power=average_power,
          compute_gradients=compute_gradients,
          gradients=result.gradients)
        result += contribution
        result.rms_with_respect_to_averages.append(
          contribution.rms_with_respect_to_average)
      else:
        result.rms_with_respect_to_averages.append(None)
    result.finalize_target_and_gradients()
    return result

  def show_adp_iso_differences_to_average(self,
         u_isos,
         site_labels,
         out=None,
         prefix=""):
    for i_group,group in enumerate(self.members):
      print >> out, prefix + "NCS restraint group %d:" % (i_group+1)
      if (    group.b_factor_weight is not None
          and group.b_factor_weight > 0):
        energies_adp_iso = group.energies_adp_iso(
          u_isos=u_isos,
          average_power=1,
          compute_gradients=False)
        print >> out, prefix + "  weight: %.6g" % energies_adp_iso.weight
        energies_adp_iso.show_differences_to_average(
          site_labels=site_labels, out=out, prefix=prefix+"  ")
      else:
        print >> out, \
          prefix+"  b_factor_weight: %s  =>  restraints disabled" % (
            str(group.b_factor_weight))

  def compute_operators(self, sites_cart):
    self.operators = []
    for group in self.members:
      self.operators.append(group.operators(sites_cart=sites_cart))

  def energies_sites(self,
        sites_cart,
        compute_gradients=True,
        gradients=None,
        normalization=False):
    self.compute_operators(sites_cart=sites_cart)
    result = scitbx.restraints.energies(
      compute_gradients=compute_gradients,
      gradients=gradients,
      gradients_size=sites_cart.size(),
      gradients_factory=flex.vec3_double,
      normalization=normalization)
    result.rms_with_respect_to_averages = []
    for operators in self.operators:
      if (    operators.group.coordinate_sigma is not None
          and operators.group.coordinate_sigma > 0):
        contribution = operators.energies_sites(
          sites_cart=sites_cart,
          compute_gradients=compute_gradients,
          gradients=result.gradients)
        result += contribution
        result.rms_with_respect_to_averages.append(
          contribution.rms_with_respect_to_average)
      else:
        result.rms_with_respect_to_averages.append(None)
    result.finalize_target_and_gradients()
    return result

  def show_operators(self, sites_cart, out=None, prefix=""):
    for i_group,group in enumerate(self.members):
      print >> out, prefix + "NCS restraint group %d:" % (i_group+1)
      ncs_operators = group.operators(sites_cart=sites_cart)
      ncs_operators.show(sites_cart=sites_cart, out=out, prefix=prefix+"  ")

  def extract_ncs_groups(self, sites_cart):
    result = []
    for group in self.members:
      ncs_operators = group.operators(sites_cart=sites_cart)
      result.append(ncs_operators)
    return result

  def show_sites_distances_to_average(self,
         sites_cart,
         site_labels,
         excessive_distance_limit=None,
         out=None,
         prefix=""):
    n_excessive = 0
    for i_group,group in enumerate(self.members):
      print >> out, prefix + "NCS restraint group %d:" % (i_group+1)
      if (    group.coordinate_sigma is not None
          and group.coordinate_sigma > 0):
        print >> out, prefix + "  coordinate_sigma: %.6g" % (
          group.coordinate_sigma)
        operators = group.operators(sites_cart=sites_cart)
        energies_sites = operators.energies_sites(
          sites_cart=sites_cart,
          compute_gradients=False)
        print >> out, prefix + "  weight:  %.6g" % energies_sites.weight
        n_excessive += energies_sites.show_distances_to_average(
          site_labels=site_labels,
          excessive_distance_limit=excessive_distance_limit,
          out=out,
          prefix=prefix+"  ")
      else:
        print >> out, \
          prefix+"  coordinate_sigma: %s  =>  restraints disabled" % (
            str(group.coordinate_sigma))
    return n_excessive

  def selection_restrained(self, n_seq=None):
    if (n_seq is None):
      n_seq = -1
      for group in self.members:
        for pair in group.selection_pairs:
          for sel in pair:
            n_seq = max(n_seq, flex.max(sel))
      n_seq += 1
    result = flex.bool(n_seq, False)
    for group in self.members:
      for pair in group.selection_pairs:
        for sel in pair:
          result.set_selected(sel, True)
    return result
