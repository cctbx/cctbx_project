from __future__ import division
from scitbx.array_family import flex

import boost.python
ext = boost.python.import_ext("mmtbx_ncs_restraints_ext")
from mmtbx_ncs_restraints_ext import *

from cctbx import adptbx
import scitbx.restraints
from scitbx.math import superpose
from scitbx import matrix
from libtbx.str_utils import show_string
from libtbx.utils import Sorry
from libtbx import adopt_init_args
from itertools import count
import sys
from libtbx.str_utils import line_breaker

class cartesian_ncs_manager(object):
  def __init__(self, model, ncs_params, ext_groups=None):
    # create bunch of group objects
    if ext_groups is not None:
      self.groups_obj = ext_groups
    else:
      many_groups = []
      ncs_obj = model.get_ncs_obj()
      ncs_groups_selection_string_list = ncs_obj.get_array_of_selections()
      ncs_restraints_group_list = ncs_obj.get_ncs_restraints_group_list()
      for i_gr, gr in enumerate(ncs_restraints_group_list):
        n_copies = gr.get_number_of_copies()
        registry = pair_registry(n_seq=model.get_number_of_atoms(), n_ncs=n_copies+1)
        for i_copy, c in enumerate(gr.copies):
          for i_seq, j_seq in zip(gr.master_iselection, c.iselection):
            stat, i_diag = registry.enter(
                i_seq=i_seq, j_seq=j_seq, j_ncs=i_copy+1)
        for i_pair,pair in enumerate(registry.selection_pairs()):
          if (pair[0].size() < 2):
            detail = ["do not produce any pairs",
                      "produce only one pair"][pair[0].size()]
            raise Sorry("\n".join([
              "NCS restraints selections %s of matching atoms:" % detail,
              "  Reference selection: %s" % show_string(self.selection_strings[0]),
              "      Other selection: %s" % show_string(
                self.selection_strings[i_pair+1])]))
        g = group(selection_strings=ncs_groups_selection_string_list[i_gr],
            registry=registry,
            coordinate_sigma=ncs_params.coordinate_sigma, # XXX GLOBAL
            b_factor_weight=ncs_params.b_factor_weight, # XXX GLOBAL
            u_average_min=1.e-6,)
        many_groups.append(g)
      self.groups_obj = groups(members=many_groups)

  def select(self, iselection):
    return cartesian_ncs_manager(
        model=None,
        ncs_params=None,
        ext_groups=self.groups_obj.select(iselection))

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
    for group in self.groups_obj.members:
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
    for i_group,group in enumerate(self.groups_obj.members):
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

  def get_n_groups(self):
    return len(self.groups_obj.members)

  def register_additional_isolated_sites(self, number):
    for group in self.groups_obj.members:
      group.register_additional_isolated_sites(number=number)

  def compute_operators(self, sites_cart):
    self.groups_obj.operators = []
    for group in self.groups_obj.members:
      self.groups_obj.operators.append(group.operators(sites_cart=sites_cart))

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
    for operators in self.groups_obj.operators:
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
    for i_group,group in enumerate(self.groups_obj.members):
      print >> out, prefix + "NCS restraint group %d:" % (i_group+1)
      ncs_operators = group.operators(sites_cart=sites_cart)
      ncs_operators.show(sites_cart=sites_cart, out=out, prefix=prefix+"  ")

  def extract_ncs_groups(self, sites_cart):
    result = []
    for group in self.groups_obj.members:
      ncs_operators = group.operators(sites_cart=sites_cart)
      result.append(ncs_operators)
    return result

  def as_pdb(self, sites_cart, out):
    result = out
    ncs_groups = self.extract_ncs_groups(sites_cart=sites_cart)
    pr = "REMARK   3  "
    print >> result, pr+"NCS DETAILS."
    print >> result, pr+" NUMBER OF NCS GROUPS : %-6d"%len(ncs_groups)
    for i_group, ncs_group in enumerate(ncs_groups):
      print >>result,pr+" NCS GROUP : %-6d"%(i_group+1)
      selection_strings = ncs_group.group.selection_strings
      for i_op,pair,mx,rms in zip(
          count(1),
          ncs_group.group.selection_pairs,
          ncs_group.matrices,
          ncs_group.rms):
        print >> result,pr+"  NCS OPERATOR : %-d" % i_op
        lines = line_breaker(selection_strings[0], width=34)
        for i_line, line in enumerate(lines):
          if(i_line == 0):
            print >> result, pr+"   REFERENCE SELECTION: %s"%line
          else:
            print >> result, pr+"                      : %s"%line
        lines = line_breaker(selection_strings[i_op], width=34)
        for i_line, line in enumerate(lines):
          if(i_line == 0):
            print >> result, pr+"   SELECTION          : %s"%line
          else:
            print >> result, pr+"                      : %s"%line
        print >> result,pr+"   ATOM PAIRS NUMBER  : %-d" % len(pair[0])
        print >> result,pr+"   RMSD               : %-10.3f" % rms
    return result.getvalue()

  def as_cif_block(self, loops, cif_block, sites_cart):
    if cif_block is None:
      cif_block = iotbx.cif.model.block()
    (ncs_ens_loop, ncs_dom_loop, ncs_dom_lim_loop, ncs_oper_loop,
        ncs_ens_gen_loop) = loops

    oper_id = 0
    ncs_groups = self.extract_ncs_groups(sites_cart=sites_cart)
    if ncs_groups is not None:
      for i_group, ncs_group in enumerate(ncs_groups):
        ncs_ens_loop.add_row((i_group+1, "?"))
        selection_strings = ncs_group.group.selection_strings
        matrices = ncs_group.matrices
        rms = ncs_group.rms
        pair_count = len(ncs_group.group.selection_pairs[0])
        for i_domain, domain_selection in enumerate(selection_strings):
          ncs_dom_loop.add_row((i_domain+1, i_group+1, "?"))
          # XXX TODO: export individual sequence ranges from selection
          ncs_dom_lim_loop.add_row(
            (i_group+1, i_domain+1, "?", "?", "?", "?", domain_selection))
          if i_domain > 0:
            rt_mx = ncs_group.matrices[i_domain-1]
            oper_id += 1
            row = [oper_id, "given"]
            row.extend(rt_mx.r)
            row.extend(rt_mx.t)
            row.append("?")
            ncs_oper_loop.add_row(row)
            ncs_ens_gen_loop.add_row((1, i_domain+1, i_group+1, oper_id))
    cif_block.add_loop(ncs_ens_loop)
    cif_block.add_loop(ncs_dom_loop)
    cif_block.add_loop(ncs_dom_lim_loop)
    if ncs_groups is not None:
      cif_block.add_loop(ncs_oper_loop)
      cif_block.add_loop(ncs_ens_gen_loop)
    return cif_block

  def show_sites_distances_to_average(self,
         sites_cart,
         site_labels,
         excessive_distance_limit=None,
         out=None,
         prefix=""):
    n_excessive = 0
    for i_group,group in enumerate(self.groups_obj.members):
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
      for group in self.groups_obj.members:
        for pair in group.selection_pairs:
          for sel in pair:
            n_seq = max(n_seq, flex.max(sel))
      n_seq += 1
    result = flex.bool(n_seq, False)
    for group in self.groups_obj.members:
      for pair in group.selection_pairs:
        for sel in pair:
          result.set_selected(sel, True)
    return result

class group(object):
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
    if(not isinstance(iselection, flex.size_t)):
      iselection = iselection.iselection()
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

  def select(self, iselection):
    members = []
    for group in self.members:
      members.append(group.select(iselection=iselection))
    return groups(members=members)
