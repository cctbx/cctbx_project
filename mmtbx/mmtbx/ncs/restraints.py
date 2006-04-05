from __future__ import generators
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
from libtbx.itertbx import count
from libtbx.utils import Sorry, buffered_indentor
import sys

class selection_properties(object):

  def __init__(self, string, iselection):
    self.string = string
    self.iselection = iselection

class mm_active_i_seqs(object):

  def __init__(self, mm, active_i_seqs):
    self.mm = mm
    self.active_i_seqs = active_i_seqs

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

class conformer_selection_properties(object):

  def __init__(self, selection_properties, mma_with_selected_atoms_by_chain):
    self.string = selection_properties.string
    self.iselection = selection_properties.iselection
    self.mma_with_selected_atoms_by_chain = mma_with_selected_atoms_by_chain
    self._residue_names_by_chain = None
    self._residue_ids_by_chain = None

  def residue_names_by_chain(self):
    if (self._residue_names_by_chain is None):
      self._residue_names_by_chain = []
      for mma_with_selected_atoms in self.mma_with_selected_atoms_by_chain:
        residue_names = []
        for mma in mma_with_selected_atoms:
          residue_names.append(mma.mm.residue_name.split("%")[0])
        self._residue_names_by_chain.append(residue_names)
    return self._residue_names_by_chain

  def residue_ids_by_chain(self, atom_attributes_list=None):
    if (self._residue_ids_by_chain is None):
      assert atom_attributes_list is not None
      self._residue_ids_by_chain = []
      for mma_with_selected_atoms in self.mma_with_selected_atoms_by_chain:
        residue_ids = []
        for mma in mma_with_selected_atoms:
          i_seq = mma.active_i_seqs[0]
          residue_ids.append(atom_attributes_list[i_seq].residue_id())
        self._residue_ids_by_chain.append(residue_ids)
    return self._residue_ids_by_chain

  def match_residue_ids(self, other, atom_attributes_list):
    result = []
    self_rnbc = self.residue_names_by_chain()
    other_rnbc = other.residue_names_by_chain()
    for i,an,bn in zip(count(), self_rnbc, other_rnbc):
      if (an == bn):
        result.append([(i,i) for i in xrange(len(an))])
      else:
        ai = self.residue_ids_by_chain(atom_attributes_list)[i]
        bi = other.residue_ids_by_chain(atom_attributes_list)[i]
        result.append(match_ordered_chain_members(a=ai, b=bi))
    return result

  def match_atoms(self,
        other,
        pairs_list,
        j_ncs,
        atom_attributes_list,
        special_position_indices,
        registry,
        selection_strings):
    for mma_i,mma_j,pairs in zip(
          self.mma_with_selected_atoms_by_chain,
          other.mma_with_selected_atoms_by_chain,
          pairs_list):
      for i_mm_i,i_mm_j in pairs:
        other_map = {}
        for i_seq_j in mma_j[i_mm_j].active_i_seqs:
          key = atom_attributes_list[i_seq_j].name
          assert not other_map.has_key(key)
          other_map[key] = i_seq_j
        for i_seq_i in mma_i[i_mm_i].active_i_seqs:
          key = atom_attributes_list[i_seq_i].name
          i_seq_j = other_map.get(key, None)
          if (i_seq_j is not None):
            assert i_seq_j >= 0
            for which,i_seq in [(self, i_seq_i), (other, i_seq_j)]:
              if (i_seq in special_position_indices):
                raise Sorry(
                  "NCS selection includes an atom in a special position:\n"
                + "  Selection: %s\n" % show_string(which.string)
                + "    Atom: %s" % atom_attributes_list[i_seq].pdb_format())
            if (i_seq_i == i_seq_j):
              raise Sorry("NCS selections restrain atom to itself:\n"
                + "  Reference selection: %s\n" % show_string(self.string)
                + "      Other selection: %s\n" % show_string(other.string)
                + "    Atom: %s" % atom_attributes_list[i_seq_i].pdb_format())
            registry_enter = registry.enter(
              i_seq=i_seq_i, j_seq=i_seq_j, j_ncs=j_ncs)
            if (registry_enter < 0):
              raise Sorry(
                "Two different NCS operators applied to same pair of atoms:\n"
                + "       Reference selection: %s\n" % show_string(
                    self.string)
                + "  Previous other selection: %s\n" % show_string(
                    selection_strings[-registry_enter])
                + "   Current other selection: %s\n" % show_string(
                    other.string)
                + "    Atom 1: %s\n" %
                    atom_attributes_list[i_seq_i].pdb_format()
                + "    Atom 2: %s" %
                    atom_attributes_list[i_seq_j].pdb_format())
            other_map[key] = -1

class pair_lists_generator(object):

  def __init__(self,
        processed_pdb,
        reference_selection_string,
        selection_strings):
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
    assert n_ncs > 0
    if (n_ncs < 2):
      raise Sorry("Only one NCS restraints selection: %s\n"
        "  At least two selections are required." %
          show_string(self.selection_strings[0]))
    # shortcuts
    all_chain_proxies = self.processed_pdb.all_chain_proxies
    atom_attributes_list = all_chain_proxies.stage_1.atom_attributes_list
    self.n_seq = len(all_chain_proxies.stage_1.atom_attributes_list)
    #
    self.selection_properties = []
    for selection_string in self.selection_strings:
      iselection = all_chain_proxies.iselection(string=selection_string)
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
    for model in all_chain_proxies.processed_models:
      for conformer in model.conformers:
        conformer_selection_properties = \
          self.get_conformer_selection_properties(conformer)
        reference = conformer_selection_properties[0]
        for j_ncs in xrange(1,n_ncs):
          other = conformer_selection_properties[j_ncs]
          reference.match_atoms(
            other=other,
            pairs_list=reference.match_residue_ids(
              other=other, atom_attributes_list=atom_attributes_list),
            j_ncs=j_ncs,
            atom_attributes_list=atom_attributes_list,
            special_position_indices=special_position_indices,
            registry=self.registry,
            selection_strings=self.selection_strings)
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

  def get_conformer_selection_properties(self, conformer):
    result = []
    for selection_properties in self.selection_properties:
      conformer_selection = flex.intersection(
        size=self.n_seq,
        iselections=[selection_properties.iselection, conformer.iselection])
      mma_with_selected_atoms_by_chain = []
      for chain in conformer.chains:
        mma_with_selected_atoms = []
        for mm in chain.monomer_mapping_summaries():
          active_i_seqs = conformer_selection.filter_indices(
            mm.all_associated_i_seqs())
          if (active_i_seqs.size() > 0):
            mma_with_selected_atoms.append(mm_active_i_seqs(mm, active_i_seqs))
        if (len(mma_with_selected_atoms) > 0):
          mma_with_selected_atoms_by_chain.append(mma_with_selected_atoms)
      result.append(conformer_selection_properties(
        selection_properties=selection_properties,
        mma_with_selected_atoms_by_chain=mma_with_selected_atoms_by_chain))
    return result

class group(object):

  def __init__(self,
        processed_pdb,
        reference_selection_string,
        selection_strings,
        coordinate_sigma,
        b_factor_weight,
        u_average_min=1.e-6):
    self.processed_pdb = processed_pdb
    g = pair_lists_generator(
      processed_pdb=processed_pdb,
      reference_selection_string=reference_selection_string,
      selection_strings=selection_strings)
    self.registry = g.registry
    self.selection_strings = g.selection_strings
    self.selection_pairs = g.registry.selection_pairs()
    self.coordinate_sigma = coordinate_sigma
    self.b_factor_weight = b_factor_weight
    self.u_average_min = u_average_min

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

  def show_differences_to_average(self, site_labels=None, out=None, prefix=""):
    if (out is None): out = sys.stdout
    if (site_labels is None):
      site_labels = [
        atom.pdb_format() for atom in
          self.group.processed_pdb.all_chain_proxies.stage_1
            .atom_attributes_list]
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
    sites_cart = self._sites_cart(sites_cart)
    self.matrices = []
    for pair in self.group.selection_pairs:
      superposition = superpose.least_squares_fit(
        reference_sites=sites_cart.select(pair[0]),
        other_sites=sites_cart.select(pair[1]))
      self.matrices.append(matrix.rt((superposition.r, superposition.t)))

  def show(self,
        sites_cart,
        n_slots_difference_histogram=6,
        out=None,
        prefix=""):
    if (out is None): out = sys.stdout
    selection_strings = self.group.selection_strings
    sites_cart = self._sites_cart(sites_cart)
    for i_op,pair,mx in zip(
          count(1),
          self.group.selection_pairs,
          self.matrices):
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
        prefix + "  RMS difference with respect to the reference: %8.6f" % (
          flex.mean(d_sq)**0.5)

  def _sites_cart(self, sites_cart):
    stage_1 = self.group.processed_pdb.all_chain_proxies.stage_1
    if (sites_cart is None):
      sites_cart = stage_1.get_sites_cart()
    assert sites_cart.size() == len(stage_1.atom_attributes_list)
    return sites_cart

  def energies_sites(self,
        sites_cart,
        compute_gradients=True,
        gradients=None,
        sites_average=None):
    return _energies_sites(
      operators=self,
      sites_cart=self._sites_cart(sites_cart),
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

  def show_distances_to_average(self, site_labels=None, out=None, prefix=""):
    if (out is None): out = sys.stdout
    if (site_labels is None):
      site_labels = [
        atom.pdb_format() for atom in
          self.operators.group.processed_pdb.all_chain_proxies.stage_1
            .atom_attributes_list]
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
      for i,c,a in zip(pair[1], sites_current, sites_average):
        diff = matrix.col(c) - matrix.col(a)
        print >> out, prefix + fmt % (site_labels[i], abs(diff))
    sel = (self.sites_count != 0).iselection()
    show_selection(i_ncs=0, pair=[sel, sel], op=None)
    for i_ncs,pair,op in zip(count(1),
                             self.operators.group.selection_pairs,
                             self.operators.matrices):
      show_selection(i_ncs=i_ncs, pair=pair, op=op)

class groups(object):

  def __init__(self, members=None):
    if (members is None):
      self.members = []
    else:
      self.members = members
    self.operators = None

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
         site_labels=None,
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
        lock_operators=False,
        normalization=False):
    if (not lock_operators):
      self.compute_operators(sites_cart=sites_cart)
    else:
      assert self.operators is not None
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

  def show_sites_distances_to_average(self,
         sites_cart,
         site_labels=None,
         out=None,
         prefix=""):
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
        energies_sites.show_distances_to_average(
          site_labels=site_labels, out=out, prefix=prefix+"  ")
      else:
        print >> out, \
          prefix+"  coordinate_sigma: %s  =>  restraints disabled" % (
            str(group.coordinate_sigma))

  def show_unrestrained_atoms(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    all_chain_proxies = None
    for group in self.members:
      if (all_chain_proxies is None):
        all_chain_proxies = group.processed_pdb.all_chain_proxies
        atom_attributes_list = all_chain_proxies.stage_1.atom_attributes_list
        to_display = flex.bool(len(atom_attributes_list), True)
      else:
        assert all_chain_proxies is group.processed_pdb.all_chain_proxies
      for pair in group.selection_pairs:
        for sel in pair:
          to_display.set_selected(sel, False)
    buffer_main = buffered_indentor(file_object=out, indent=prefix)
    print >> buffer_main, "Atoms without NCS restraints:"
    for i_model,model in enumerate(all_chain_proxies.processed_models):
      buffer_model = buffer_main.shift_right()
      print >> buffer_model, "Model %d," % (i_model+1), \
        "PDB serial number:", model.serial
      for i_conformer,conformer in enumerate(model.conformers):
        buffer_conformer = buffer_model.shift_right()
        print >> buffer_conformer, "Conformer %d," % (i_conformer+1), \
          "PDB altLoc:", show_string(conformer.altLoc)
        for i_chain,chain in enumerate(conformer.chains):
          buffer_chain = buffer_conformer.shift_right()
          print >> buffer_chain, "Chain %d," % (i_chain+1), \
            "PDB chainID: %s," % show_string(chain.chainID), \
            "segID: %s" % show_string(chain.segID)
          for mm in chain.monomer_mapping_summaries():
            i_seqs = mm.all_associated_i_seqs()
            if (to_display.select(i_seqs).count(True) > 0):
              buffer_residue = buffer_chain.shift_right()
              print >> buffer_residue, \
                '"' + atom_attributes_list[i_seqs[0]].pdb_format()[6:]
              for i_seq in i_seqs:
                buffer_atom = buffer_residue.shift_right()
                if (to_display[i_seq]):
                  to_display[i_seq] = False
                  print >> buffer_atom, \
                    atom_attributes_list[i_seq].pdb_format()[:6]+'"'
                  buffer_atom.write_buffer()
