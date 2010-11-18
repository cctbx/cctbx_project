import cctbx.crystal.direct_space_asu # import dependency
from cctbx.array_family import flex
import scitbx.array_family.shared # import dependency
import cctbx.geometry # import dependency

from libtbx.test_utils import approx_equal
from libtbx.str_utils import show_string

import boost.python
ext = boost.python.import_ext("cctbx_geometry_restraints_ext")
from cctbx_geometry_restraints_ext import *

import scitbx.stl.map
from stdlib import math
import sys

nonbonded_radius_table = scitbx.stl.map.stl_string_double

nonbonded_distance_table = scitbx.stl.map.stl_string_stl_map_stl_string_double
nonbonded_distance_dict = scitbx.stl.map.stl_string_double

def weight_as_sigma(weight):
  if (weight <= 0): return 0
  return 1/weight**0.5

def sigma_as_weight(sigma):
  if (sigma <= 0): return 0
  return 1/sigma**2

def angle_delta_deg(angle_1, angle_2, periodicity=1):
  half_period = 180./max(1,periodicity)
  d = math.fmod(angle_2-angle_1, 2*half_period)
  if   (d < -half_period): d += 2*half_period
  elif (d >  half_period): d -= 2*half_period
  return d

class proxy_registry_process_result(object):

  def __init__(self,
        tabulated_proxy=None,
        is_new=False,
        is_conflicting=False,
        conflict_source_labels=None):
    self.tabulated_proxy = tabulated_proxy
    self.is_new = is_new
    self.is_conflicting = is_conflicting
    self.conflict_source_labels = conflict_source_labels

class proxy_registry_base(object):

  def __init__(self, proxies, strict_conflict_handling):
    self.proxies = proxies
    self.strict_conflict_handling = strict_conflict_handling
    self.n_resolved_conflicts = 0
    self.counts = flex.size_t()
    self.discard_table()

  def initialize_table(self):
    self.table = {}
    self.source_labels = flex.std_string()
    self.source_n_expected_atoms = flex.int()

  def discard_table(self):
    self.table = None
    self.source_labels = None
    self.source_n_expected_atoms = None

  def append_custom_proxy(self, proxy):
    assert self.table is None
    self.proxies.append(proxy)
    self.counts.append(1)

  def _append_proxy(self, source_info, proxy, process_result):
    self.proxies.append(proxy)
    self.counts.append(1)
    self.source_labels.append(source_info.labels())
    self.source_n_expected_atoms.append(source_info.n_expected_atoms())
    process_result.tabulated_proxy = proxy
    process_result.is_new = True

  def _handle_conflict(self,
        source_info,
        proxy,
        i_list,
        process_result):
    source_n_expected_atoms = source_info.n_expected_atoms()
    if (self.strict_conflict_handling
        or self.source_n_expected_atoms[i_list]
           == source_n_expected_atoms):
      process_result.is_conflicting = True
      process_result.conflict_source_labels = (self.source_labels[i_list],
                                               source_info.labels())
    else:
      self.n_resolved_conflicts += 1
      if (self.source_n_expected_atoms[i_list] < source_n_expected_atoms):
        self.proxies[i_list] = proxy
        self.source_labels[i_list] += ", " + source_info.labels()
        self.source_n_expected_atoms[i_list] = source_n_expected_atoms
        process_result.tabulated_proxy = proxy

class bond_simple_proxy_registry(proxy_registry_base):

  def __init__(self, n_seq, strict_conflict_handling):
    proxy_registry_base.__init__(self,
      proxies=shared_bond_simple_proxy(),
      strict_conflict_handling=strict_conflict_handling)
    self.n_seq = n_seq

  def initialize_table(self):
    proxy_registry_base.initialize_table(self)
    self.table = [{} for i in xrange(self.n_seq)]

  def process(self, source_info, proxy, tolerance=1.e-6):
    result = proxy_registry_process_result()
    proxy = proxy.sort_i_seqs()
    if (not self.table[proxy.i_seqs[0]].has_key(proxy.i_seqs[1])):
      self.table[proxy.i_seqs[0]][proxy.i_seqs[1]] = self.proxies.size()
      self._append_proxy(
        source_info=source_info,
        proxy=proxy,
        process_result=result)
    else:
      i_list = self.table[proxy.i_seqs[0]][proxy.i_seqs[1]]
      result.tabulated_proxy = self.proxies[i_list]
      if (   abs(result.tabulated_proxy.distance_ideal
                                - proxy.distance_ideal) > tolerance
          or abs(result.tabulated_proxy.weight - proxy.weight)
               > tolerance):
        self._handle_conflict(
          source_info=source_info,
          proxy=proxy,
          i_list=i_list,
          process_result=result)
      if (not result.is_conflicting):
        self.counts[i_list] += 1
    return result

class angle_proxy_registry(proxy_registry_base):

  def __init__(self, strict_conflict_handling):
    proxy_registry_base.__init__(self,
      proxies=shared_angle_proxy(),
      strict_conflict_handling=strict_conflict_handling)

  def process(self, source_info, proxy, tolerance=1.e-6):
    result = proxy_registry_process_result()
    proxy = proxy.sort_i_seqs()
    tab_i_seq_1 = self.table.setdefault(proxy.i_seqs[1], {})
    i_seqs_0_2 = (proxy.i_seqs[0], proxy.i_seqs[2])
    if (not tab_i_seq_1.has_key(i_seqs_0_2)):
      tab_i_seq_1[i_seqs_0_2] = self.proxies.size()
      self._append_proxy(
        source_info=source_info,
        proxy=proxy,
        process_result=result)
    else:
      i_list = tab_i_seq_1[i_seqs_0_2]
      result.tabulated_proxy = self.proxies[i_list]
      if (   abs(angle_delta_deg(result.tabulated_proxy.angle_ideal,
                                 proxy.angle_ideal)) > tolerance
          or abs(result.tabulated_proxy.weight - proxy.weight)
               > tolerance):
        self._handle_conflict(
          source_info=source_info,
          proxy=proxy,
          i_list=i_list,
          process_result=result)
      if (not result.is_conflicting):
        self.counts[i_list] += 1
    return result

  def lookup_i_proxy(self, i_seqs):
    "See also: cctbx::geometry_restraints::angle_proxy::sort_i_seqs()"
    (i0, i1, i2) = i_seqs
    tab_i_seq_1 = self.table.get(i1)
    if (tab_i_seq_1 is None):
      return None
    if (i0 > i2):
      return tab_i_seq_1.get((i2, i0))
    return tab_i_seq_1.get((i0, i2))

class dihedral_proxy_registry(proxy_registry_base):

  def __init__(self, strict_conflict_handling):
    proxy_registry_base.__init__(self,
      proxies=shared_dihedral_proxy(),
      strict_conflict_handling=strict_conflict_handling)

  def process(self, source_info, proxy, tolerance=1.e-6):
    result = proxy_registry_process_result()
    proxy = proxy.sort_i_seqs()
    tab_i_seq_0 = self.table.setdefault(proxy.i_seqs[0], {})
    i_seqs_1_2_3 = (proxy.i_seqs[1], proxy.i_seqs[2], proxy.i_seqs[3])
    if (not tab_i_seq_0.has_key(i_seqs_1_2_3)):
      tab_i_seq_0[i_seqs_1_2_3] = self.proxies.size()
      self._append_proxy(
        source_info=source_info,
        proxy=proxy,
        process_result=result)
    else:
      i_list = tab_i_seq_0[i_seqs_1_2_3]
      result.tabulated_proxy = self.proxies[i_list]
      if (   abs(angle_delta_deg(result.tabulated_proxy.angle_ideal,
                                 proxy.angle_ideal,
                                 proxy.periodicity)) > tolerance
          or abs(result.tabulated_proxy.weight - proxy.weight)
               > tolerance
          or result.tabulated_proxy.periodicity != proxy.periodicity):
        self._handle_conflict(
          source_info=source_info,
          proxy=proxy,
          i_list=i_list,
          process_result=result)
      if (not result.is_conflicting):
        self.counts[i_list] += 1
    return result

  def lookup_i_proxy(self, i_seqs):
    "See also: cctbx::geometry_restraints::dihedral_proxy::sort_i_seqs()"
    (i0, i1, i2, i3) = i_seqs
    angle_sign = 1
    if (i0 > i3):
      i0, i3 = i3, i0
      angle_sign *= -1
    if (i1 > i2):
      i1, i2 = i2, i1
      angle_sign *= -1
    tab_i_seq_0 = self.table.get(i0)
    if (tab_i_seq_0 is None):
      return (None, None)
    return (tab_i_seq_0.get((i1, i2, i3)), angle_sign)

class chirality_proxy_registry(proxy_registry_base):

  def __init__(self, strict_conflict_handling):
    proxy_registry_base.__init__(self,
      proxies=shared_chirality_proxy(),
      strict_conflict_handling=strict_conflict_handling)

  def process(self, source_info, proxy, tolerance=1.e-6):
    result = proxy_registry_process_result()
    proxy = proxy.sort_i_seqs()
    tab_i_seq_0 = self.table.setdefault(proxy.i_seqs[0], {})
    i_seqs_1_2_3 = (proxy.i_seqs[1], proxy.i_seqs[2], proxy.i_seqs[3])
    if (not tab_i_seq_0.has_key(i_seqs_1_2_3)):
      tab_i_seq_0[i_seqs_1_2_3] = self.proxies.size()
      self._append_proxy(
        source_info=source_info,
        proxy=proxy,
        process_result=result)
    else:
      i_list = tab_i_seq_0[i_seqs_1_2_3]
      result.tabulated_proxy = self.proxies[i_list]
      if (   abs(result.tabulated_proxy.volume_ideal - proxy.volume_ideal)
               > tolerance
          or abs(result.tabulated_proxy.weight - proxy.weight)
               > tolerance):
        self._handle_conflict(
          source_info=source_info,
          proxy=proxy,
          i_list=i_list,
          process_result=result)
      if (not result.is_conflicting):
        self.counts[i_list] += 1
    return result

class planarity_proxy_registry(proxy_registry_base):

  def __init__(self, strict_conflict_handling):
    proxy_registry_base.__init__(self,
      proxies=shared_planarity_proxy(),
      strict_conflict_handling=strict_conflict_handling)

  def process(self, source_info, proxy, tolerance=1.e-6):
    assert proxy.i_seqs.size() > 0
    result = proxy_registry_process_result()
    proxy = proxy.sort_i_seqs()
    tab_i_seq_0 = self.table.setdefault(proxy.i_seqs[0], {})
    i_seqs_1_up = tuple(proxy.i_seqs[1:])
    if (not tab_i_seq_0.has_key(i_seqs_1_up)):
      tab_i_seq_0[i_seqs_1_up] = self.proxies.size()
      self._append_proxy(
        source_info=source_info,
        proxy=proxy,
        process_result=result)
    else:
      i_list = tab_i_seq_0[i_seqs_1_up]
      result.tabulated_proxy = self.proxies[i_list]
      if (not approx_equal(result.tabulated_proxy.weights, proxy.weights,
                           eps=tolerance)):
        self._handle_conflict(
          source_info=source_info,
          proxy=proxy,
          i_list=i_list,
          process_result=result)
      if (not result.is_conflicting):
        self.counts[i_list] += 1
    return result

class _prolsq_repulsion_function(
        boost.python.injector, prolsq_repulsion_function):

  def customized_copy(O, c_rep=None, k_rep=None, irexp=None, rexp=None):
    if (c_rep is None): c_rep = O.c_rep
    if (k_rep is None): k_rep = O.k_rep
    if (irexp is None): irexp = O.irexp
    if (rexp is None): rexp = O.rexp
    return prolsq_repulsion_function(
      c_rep=c_rep, k_rep=k_rep, irexp=irexp, rexp=rexp)

def _bond_show_sorted_impl(self,
                           by_value,
                           sites_cart,
                           site_labels=None,
                           unit_cell=None,
                           f=None,
                           prefix="",
                           max_items=None):
  if unit_cell is None:
    sorted_table, n_not_shown = self.get_sorted(
      by_value=by_value,
      sites_cart=sites_cart,
      site_labels=site_labels,
      max_items=max_items)
  else:
    sorted_table, n_not_shown = self.get_sorted(
      by_value=by_value,
      sites_cart=sites_cart,
      unit_cell=unit_cell,
      site_labels=site_labels,
      max_items=max_items)
  if sorted_table is None :
    return
  if (f is None): f = sys.stdout
  print >> f, "%sSorted by %s:" % (prefix, by_value)
  for restraint_info in sorted_table :
    (labels, distance_ideal, distance_model, slack, delta, sigma, weight,
     residual, sym_op_j, rt_mx) = restraint_info
    s = "bond"
    for label in labels :
      print >> f, "%s%4s %s" % (prefix, s, label)
      s = ""
    if (slack == 0):
      l = ""
      v = ""
    else:
      l = "  slack"
      v = " %6.3f" % slack
    print >> f, "%s  ideal  model%s  delta    sigma   weight residual%s" % (
      prefix, l, sym_op_j)
    print >> f, "%s  %5.3f %6.3f%s %6.3f %6.2e %6.2e %6.2e" % (
      prefix, distance_ideal, distance_model, v, delta,
      sigma, weight, residual),
    if (rt_mx is not None):
      print >> f, rt_mx,
    print >> f
  if (n_not_shown != 0):
    print >> f, prefix + "... (remaining %d not shown)" % n_not_shown

class _bond_simple_proxy(boost.python.injector, shared_bond_simple_proxy):

  def get_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        unit_cell=None,
        max_items=None):
    assert by_value in ["residual", "delta"]
    assert site_labels is None or len(site_labels) == sites_cart.size()
    if (self.size() == 0): return None, None
    if (max_items is not None and max_items <= 0): return None, None
    if (by_value == "residual"):
      data_to_sort = self.residuals(sites_cart=sites_cart, unit_cell=unit_cell)
    elif (by_value == "delta"):
      data_to_sort = flex.abs(
        self.deltas(sites_cart=sites_cart, unit_cell=unit_cell))
    else:
      raise AssertionError
    i_proxies_sorted = flex.sort_permutation(data=data_to_sort, reverse=True)
    if (max_items is not None):
      i_proxies_sorted = i_proxies_sorted[:max_items]
    smallest_distance_model = None
    sorted_table = []
    for i_proxy in i_proxies_sorted:
      proxy = self[i_proxy]
      i_seq,j_seq = proxy.i_seqs
      rt_mx = proxy.rt_mx_ji
      if rt_mx is None: sym_op_j = ""
      else: sym_op_j = " sym.op."
      if unit_cell is None:
        restraint = bond(
          sites_cart=sites_cart,
          proxy=proxy)
      else:
        restraint = bond(
          unit_cell,
          sites_cart=sites_cart,
          proxy=proxy)
      labels = []
      for i in [i_seq, j_seq]:
        if (site_labels is None): l = str(i)
        else:                     l = site_labels[i]
        labels.append(l)
      sorted_table.append(
        (labels, restraint.distance_ideal, restraint.distance_model,
         restraint.slack, restraint.delta,
         weight_as_sigma(weight=restraint.weight), restraint.weight,
         restraint.residual(), sym_op_j, rt_mx))
      if (smallest_distance_model is None
          or smallest_distance_model > restraint.distance_model):
        smallest_distance_model = restraint.distance_model
    n_not_shown = data_to_sort.size() - i_proxies_sorted.size()
    return sorted_table, n_not_shown

  def show_sorted(self,
                  by_value,
                  sites_cart,
                  site_labels=None,
                  unit_cell=None,
                  f=None,
                  prefix="",
                  max_items=None):
    if f is None: f = sys.stdout
    print >> f, "%sBond restraints: %d" % (prefix, self.size())
    _bond_show_sorted_impl(self, by_value,
                           sites_cart=sites_cart,
                           site_labels=site_labels,
                           unit_cell=unit_cell,
                           f=f,
                           prefix=prefix,
                           max_items=max_items)

  def deltas(self, sites_cart, unit_cell=None):
    if unit_cell is None:
      return bond_deltas(sites_cart=sites_cart, proxies=self)
    else:
      return bond_deltas(
        unit_cell=unit_cell, sites_cart=sites_cart, proxies=self)

  def residuals(self, sites_cart, unit_cell=None):
    if unit_cell is None:
      return bond_residuals(sites_cart=sites_cart, proxies=self)
    else:
      return bond_residuals(
        unit_cell=unit_cell, sites_cart=sites_cart, proxies=self)


class _bond_sorted_asu_proxies(boost.python.injector, bond_sorted_asu_proxies):

  def show_histogram_of_model_distances(self,
        sites_cart,
        n_slots=5,
        cutoff_warn_small=0.5,
        cutoff_warn_large=5,
        cutoff_warn_extreme=20,
        f=None,
        prefix=""):
    if (self.n_total() == 0): return None
    if (f is None): f = sys.stdout
    print >> f, "%sHistogram of bond lengths:" % prefix
    histogram = flex.histogram(
      data=bond_distances_model(
        sites_cart=sites_cart,
        sorted_asu_proxies=self),
      n_slots=n_slots)
    low_cutoff = histogram.data_min()
    for i,n in enumerate(histogram.slots()):
      high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
      print >> f, "%s  %8.2f - %8.2f: %d" % (
        prefix, low_cutoff, high_cutoff, n)
      low_cutoff = high_cutoff
    if (cutoff_warn_small is not None
        and histogram.data_min() < cutoff_warn_small):
      print >> f, "%sWarning: very small bond lengths." % prefix
    if (cutoff_warn_extreme is not None
        and histogram.data_max() > cutoff_warn_extreme):
      print >> f, "%sWarning: extremely large bond lengths." % prefix
    elif (cutoff_warn_large is not None
          and histogram.data_max() > cutoff_warn_large):
      print >> f, "%sWarning: very large bond lengths." % prefix
    return histogram

  def deltas(self, sites_cart):
    return bond_deltas(sites_cart=sites_cart, sorted_asu_proxies=self)

  def residuals(self, sites_cart):
    return bond_residuals(sites_cart=sites_cart, sorted_asu_proxies=self)

  def show_histogram_of_deltas(self,
        sites_cart,
        n_slots=5,
        f=None,
        prefix=""):
    if (self.n_total() == 0): return
    if (f is None): f = sys.stdout
    print >> f, "%sHistogram of bond deltas:" % prefix
    histogram = flex.histogram(
      data=flex.abs(self.deltas(sites_cart=sites_cart)),
      n_slots=n_slots)
    low_cutoff = histogram.data_min()
    for i,n in enumerate(histogram.slots()):
      high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
      print >> f, "%s  %8.3f - %8.3f: %d" % (
        prefix, low_cutoff, high_cutoff, n)
      low_cutoff = high_cutoff
    return histogram

  def get_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        max_items=None):
    assert by_value in ["residual", "delta"]
    assert site_labels is None or len(site_labels) == sites_cart.size()
    if (self.n_total() == 0): return None, None
    if (max_items is not None and max_items <= 0): return None, None
    if (by_value == "residual"):
      data_to_sort = self.residuals(sites_cart=sites_cart)
    elif (by_value == "delta"):
      data_to_sort = flex.abs(self.deltas(sites_cart=sites_cart))
    else:
      raise AssertionError
    i_proxies_sorted = flex.sort_permutation(data=data_to_sort, reverse=True)
    if (max_items is not None):
      i_proxies_sorted = i_proxies_sorted[:max_items]
    if (self.asu.size() == 0):
      asu_mappings = None
    else:
      asu_mappings = self.asu_mappings()
    smallest_distance_model = None
    n_simple = self.simple.size()
    sorted_table = []
    for i_proxy in i_proxies_sorted:
      if (i_proxy < n_simple):
        proxy = self.simple[i_proxy]
        i_seq,j_seq = proxy.i_seqs
        rt_mx = None
        sym_op_j = ""
        restraint = bond(
          sites_cart=sites_cart,
          proxy=proxy)
      else:
        proxy = self.asu[i_proxy-n_simple]
        i_seq,j_seq = proxy.i_seq,proxy.j_seq
        rt_mx = asu_mappings.get_rt_mx_ji(pair=proxy)
        sym_op_j = " sym.op."
        restraint = bond(
          sites_cart=sites_cart,
          asu_mappings=asu_mappings,
          proxy=proxy)
      labels = []
      for i in [i_seq, j_seq]:
        if (site_labels is None): l = str(i)
        else:                     l = site_labels[i]
        labels.append(l)
      sorted_table.append(
        (labels, restraint.distance_ideal, restraint.distance_model,
         restraint.slack, restraint.delta,
         weight_as_sigma(weight=restraint.weight), restraint.weight,
         restraint.residual(), sym_op_j, rt_mx))
      if (smallest_distance_model is None
          or smallest_distance_model > restraint.distance_model):
        smallest_distance_model = restraint.distance_model
    n_not_shown = data_to_sort.size() - i_proxies_sorted.size()
    return sorted_table, n_not_shown

  def show_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        f=None,
        prefix="",
        max_items=None):
    if f is None: f = sys.stdout
    print >> f, "%sBond restraints: %d" % (prefix, self.n_total())
    _bond_show_sorted_impl(self, by_value,
                           sites_cart=sites_cart,
                           site_labels=site_labels,
                           f=f,
                           prefix=prefix,
                           max_items=max_items)

class _nonbonded_sorted_asu_proxies(boost.python.injector,
        nonbonded_sorted_asu_proxies):

  def deltas(self, sites_cart):
    return nonbonded_deltas(sites_cart=sites_cart, sorted_asu_proxies=self)

  def show_histogram_of_model_distances(self,
        sites_cart,
        n_slots=5,
        cutoff_warn_small=1.0,
        f=None,
        prefix=""):
    if (self.n_total() == 0): return None
    if (f is None): f = sys.stdout
    print >> f, "%sHistogram of nonbonded interaction distances:" % prefix
    histogram = flex.histogram(
      data=self.deltas(sites_cart=sites_cart),
      n_slots=n_slots)
    low_cutoff = histogram.data_min()
    for i,n in enumerate(histogram.slots()):
      high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
      print >> f, "%s  %8.2f - %8.2f: %d" % (
        prefix, low_cutoff, high_cutoff, n)
      low_cutoff = high_cutoff
    if (cutoff_warn_small is not None
        and histogram.data_min() < cutoff_warn_small):
      print >> f, "%sWarning: very small nonbonded interaction distances." % (
        prefix)
    return histogram

  def show_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        f=None,
        prefix="",
        max_items=None,
        suppress_model_minus_vdw_greater_than=0.2,
        but_show_all_model_up_to=3.5):
    assert by_value in ["delta"]
    if (f is None): f = sys.stdout
    deltas = nonbonded_deltas(sites_cart=sites_cart, sorted_asu_proxies=self)
    print >> f, "%sNonbonded interactions: %d" % (prefix, deltas.size())
    if (deltas.size() == 0): return
    if (max_items is not None and max_items <= 0): return
    i_proxies_sorted = flex.sort_permutation(data=deltas)
    if (max_items is not None):
      i_proxies_sorted = i_proxies_sorted[:max_items]
    if (self.asu.size() == 0):
      asu_mappings = None
    else:
      asu_mappings = self.asu_mappings()
    print >> f, "%sSorted by model distance:" % prefix
    n_simple = self.simple.size()
    for i_proxy in i_proxies_sorted:
      if (i_proxy < n_simple):
        proxy = self.simple[i_proxy]
        i_seq,j_seq = proxy.i_seqs
        rt_mx = None
        sym_op_j = ""
      else:
        proxy = self.asu[i_proxy-n_simple]
        i_seq,j_seq = proxy.i_seq,proxy.j_seq
        rt_mx = asu_mappings.get_rt_mx_ji(pair=proxy)
        sym_op_j = " sym.op."
      def suppress():
        m, v = deltas[i_proxy], proxy.vdw_distance
        if (suppress_model_minus_vdw_greater_than is None): return False
        if (m-v <= suppress_model_minus_vdw_greater_than): return False
        if (but_show_all_model_up_to is None): return True
        if (m <= but_show_all_model_up_to): return False
        return True
      if (suppress()): continue
      s = "nonbonded"
      for i in [i_seq, j_seq]:
        if (site_labels is None): l = str(i)
        else:                     l = site_labels[i]
        print >> f, "%s%9s %s" % (prefix, s, l)
        s = ""
      print >> f, "%s   model   vdw%s" % (prefix, sym_op_j)
      print >> f, "%s  %6.3f %5.3f" % (
        prefix, deltas[i_proxy], proxy.vdw_distance),
      if (rt_mx is not None):
        print >> f, rt_mx,
      print >> f
    n_not_shown = deltas.size() - i_proxies_sorted.size()
    if (n_not_shown != 0):
      print >> f, prefix + "... (remaining %d not shown)" % n_not_shown

class _angle(boost.python.injector, angle):

  def _show_sorted_item(O, f, prefix):
    print >> f, "%s    ideal   model   delta" \
      "    sigma   weight residual" % prefix
    print >> f, "%s  %7.2f %7.2f %7.2f %6.2e %6.2e %6.2e" % (
      prefix,
      O.angle_ideal, O.angle_model, O.delta,
      weight_as_sigma(weight=O.weight), O.weight, O.residual())

  def _get_sorted_item (O) :
    return [O.angle_ideal, O.angle_model, O.delta,
            weight_as_sigma(weight=O.weight), O.weight, O.residual()]

class _shared_angle_proxy(boost.python.injector, shared_angle_proxy):

  def deltas(self, sites_cart, unit_cell=None):
    if unit_cell is None:
      return angle_deltas(sites_cart=sites_cart, proxies=self)
    else:
      return angle_deltas(
        unit_cell=unit_cell, sites_cart=sites_cart, proxies=self)

  def residuals(self, sites_cart, unit_cell=None):
    if unit_cell is None:
      return angle_residuals(sites_cart=sites_cart, proxies=self)
    else:
      return angle_residuals(
        unit_cell=unit_cell, sites_cart=sites_cart, proxies=self)

  def show_histogram_of_deltas(self,
        sites_cart,
        unit_cell=None,
        n_slots=5,
        f=None,
        prefix=""):
    return _show_histogram_of_deltas_impl(O=self,
      proxy_label="bond angle",
      format_cutoffs="%8.2f",
      unit_cell=unit_cell,
      sites_cart=sites_cart, n_slots=n_slots, f=f, prefix=prefix)

  def show_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        unit_cell=None,
        f=None,
        prefix="",
        max_items=None):
    _show_sorted_impl(O=self,
        proxy_type=angle,
        proxy_label="Bond angle",
        item_label="angle",
        by_value=by_value, unit_cell=unit_cell, sites_cart=sites_cart,
        site_labels=site_labels, f=f, prefix=prefix, max_items=max_items)

  def get_sorted (self,
        by_value,
        sites_cart,
        site_labels=None,
        unit_cell=None,
        max_items=None):
    return _get_sorted_impl(O=self,
        proxy_type=angle,
        by_value=by_value, unit_cell=unit_cell, sites_cart=sites_cart,
        site_labels=site_labels, max_items=max_items,
        get_restraints_only=False)

class _dihedral(boost.python.injector, dihedral):

  def _show_sorted_item(O, f, prefix):
    print >> f, "%s    ideal   model   delta" \
      " %        s    sigma   weight residual" % (
        prefix, {False: "sinusoidal", True: " harmonic "}[O.periodicity <= 0])
    print >> f, "%s  %7.2f %7.2f %7.2f %5d      %6.2e %6.2e %6.2e" % (
      prefix,
      O.angle_ideal, O.angle_model, O.delta, O.periodicity,
      weight_as_sigma(weight=O.weight), O.weight, O.residual())

  def _get_sorted_item (O) :
    return [O.angle_ideal, O.angle_model, O.delta, O.periodicity,
            weight_as_sigma(weight=O.weight), O.weight, O.residual()]

class _shared_dihedral_proxy(boost.python.injector, shared_dihedral_proxy):

  def deltas(self, sites_cart, unit_cell=None):
    if unit_cell is None:
      return dihedral_deltas(sites_cart=sites_cart, proxies=self)
    else:
      return dihedral_deltas(
        unit_cell=unit_cell, sites_cart=sites_cart, proxies=self)

  def residuals(self, sites_cart, unit_cell=None):
    if unit_cell is None:
      return dihedral_residuals(sites_cart=sites_cart, proxies=self)
    else:
      return dihedral_residuals(
        unit_cell=unit_cell, sites_cart=sites_cart, proxies=self)

  def show_histogram_of_deltas(self,
        sites_cart,
        unit_cell=None,
        n_slots=5,
        f=None,
        prefix=""):
    return _show_histogram_of_deltas_impl(O=self,
      proxy_label="dihedral angle",
      format_cutoffs="%8.2f",
      unit_cell=unit_cell,
      sites_cart=sites_cart, n_slots=n_slots, f=f, prefix=prefix)

  def show_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        unit_cell=None,
        f=None,
        prefix="",
        max_items=None,
        is_reference=False):
    if is_reference:
      proxy_label = "Reference dihedral angle"
    else:
      proxy_label = "Dihedral angle"
    _show_sorted_impl(O=self,
        proxy_type=dihedral,
        proxy_label=proxy_label,
        item_label="dihedral",
        by_value=by_value, unit_cell=unit_cell, sites_cart=sites_cart,
        site_labels=site_labels, f=f, prefix=prefix, max_items=max_items)

  def get_sorted (self,
        by_value,
        sites_cart,
        site_labels=None,
        unit_cell=None,
        max_items=None):
    return _get_sorted_impl(O=self,
        proxy_type=dihedral,
        by_value=by_value, unit_cell=unit_cell, sites_cart=sites_cart,
        site_labels=site_labels, max_items=max_items,
        get_restraints_only=False)

class _chirality(boost.python.injector, chirality):

  def _show_sorted_item(O, f, prefix):
    print >> f, "%s  both_signs  ideal   model" \
      "   delta    sigma   weight residual" % prefix
    print >> f, "%s    %-5s   %7.2f %7.2f %7.2f %6.2e %6.2e %6.2e" % (
      prefix,
      str(O.both_signs), O.volume_ideal, O.volume_model, O.delta,
      weight_as_sigma(weight=O.weight), O.weight, O.residual())

  def _get_sorted_item (O) :
    return [str(O.both_signs), O.volume_ideal, O.volume_model, O.delta,
      weight_as_sigma(weight=O.weight), O.weight, O.residual()]

class _shared_chirality_proxy(boost.python.injector, shared_chirality_proxy):

  def deltas(self, sites_cart):
    return chirality_deltas(sites_cart=sites_cart, proxies=self)

  def residuals(self, sites_cart):
    return chirality_residuals(sites_cart=sites_cart, proxies=self)

  def show_histogram_of_deltas(self,
        sites_cart,
        n_slots=5,
        f=None,
        prefix=""):
    return _show_histogram_of_deltas_impl(O=self,
      proxy_label="chiral volume",
      format_cutoffs="%8.3f",
      unit_cell=None,
      sites_cart=sites_cart, n_slots=n_slots, f=f, prefix=prefix)

  def show_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        f=None,
        prefix="",
        max_items=None):
    _show_sorted_impl(O=self,
        proxy_type=chirality,
        proxy_label="Chirality",
        item_label="chirality",
        by_value=by_value, unit_cell=None, sites_cart=sites_cart,
        site_labels=site_labels, f=f, prefix=prefix, max_items=max_items)

  def get_sorted (self,
        by_value,
        sites_cart,
        site_labels=None,
        max_items=None):
    return _get_sorted_impl(O=self,
        proxy_type=chirality,
        by_value=by_value, unit_cell=None, sites_cart=sites_cart,
        site_labels=site_labels, max_items=max_items,
        get_restraints_only=False)

class _shared_planarity_proxy(boost.python.injector, shared_planarity_proxy):

  def deltas_rms(O, sites_cart, unit_cell=None):
    if unit_cell is None:
      return planarity_deltas_rms(sites_cart=sites_cart, proxies=O)
    else:
      return planarity_deltas_rms(
        unit_cell=unit_cell, sites_cart=sites_cart, proxies=O)

  def residuals(O, sites_cart, unit_cell=None):
    if unit_cell is None:
      return planarity_residuals(sites_cart=sites_cart, proxies=O)
    else:
      return planarity_residuals(
        unit_cell=unit_cell, sites_cart=sites_cart, proxies=O)

  def get_sorted (O,
        by_value,
        sites_cart,
        site_labels=None,
        unit_cell=None,
        max_items=None):
    assert by_value in ["residual", "rms_deltas"]
    assert site_labels is None or len(site_labels) == sites_cart.size()
    if (O.size() == 0): return None, None
    if (max_items is not None and max_items <= 0): return None, None
    if (by_value == "residual"):
      if unit_cell is None:
        data_to_sort = O.residuals(sites_cart=sites_cart)
      else:
        data_to_sort = O.residuals(unit_cell=unit_cell,sites_cart=sites_cart)
    elif (by_value == "rms_deltas"):
      if unit_cell is None:
        data_to_sort = O.deltas_rms(sites_cart=sites_cart)
      else:
        data_to_sort = O.deltas_rms(unit_cell=unit_cell, sites_cart=sites_cart)
    else:
      raise AssertionError
    i_proxies_sorted = flex.sort_permutation(data=data_to_sort, reverse=True)
    if (max_items is not None):
      i_proxies_sorted = i_proxies_sorted[:max_items]
    sorted_table = []
    for i_proxy in i_proxies_sorted:
      proxy = O[i_proxy]
      len_max = 0
      labels = []
      for i_seq in proxy.i_seqs:
        if (site_labels is None): l = str(i_seq)
        else:                     l = site_labels[i_seq]
        len_max = max(len_max, len(l))
        labels.append(l)
      if unit_cell is None:
        restraint = planarity(sites_cart=sites_cart, proxy=proxy)
      else:
        restraint = planarity(unit_cell=unit_cell, sites_cart=sites_cart,
                              proxy=proxy)
      restraint_atoms = []
      for i, (i_seq,weight,delta,l) in enumerate(zip(proxy.i_seqs, proxy.weights,
                                      restraint.deltas(), labels)):
        sym_op = ""
        if proxy.sym_ops:
          rt_mx = proxy.sym_ops[i]
          if not rt_mx.is_unit_mx():
            sym_op = "%s %s" %(rdr_spacer, rt_mx.as_xyz())
        restraint_atoms.append((l, delta, weight_as_sigma(weight=weight),
          weight, sym_op))
      sorted_table.append((restraint_atoms, restraint.rms_deltas(),
        restraint.residual()))
    n_not_shown = O.size() - i_proxies_sorted.size()
    return sorted_table, n_not_shown

  # TODO: convert this to use get_sorted
  def show_sorted(O,
        by_value,
        sites_cart,
        site_labels=None,
        unit_cell=None,
        f=None,
        prefix="",
        max_items=None):
    assert by_value in ["residual", "rms_deltas"]
    assert site_labels is None or len(site_labels) == sites_cart.size()
    if (f is None): f = sys.stdout
    print >> f, "%sPlanarity restraints: %d" % (prefix, O.size())
    if (O.size() == 0): return
    if (max_items is not None and max_items <= 0): return
    if (by_value == "residual"):
      if unit_cell is None:
        data_to_sort = O.residuals(sites_cart=sites_cart)
      else:
        data_to_sort = O.residuals(unit_cell=unit_cell,sites_cart=sites_cart)
    elif (by_value == "rms_deltas"):
      if unit_cell is None:
        data_to_sort = O.deltas_rms(sites_cart=sites_cart)
      else:
        data_to_sort = O.deltas_rms(unit_cell=unit_cell, sites_cart=sites_cart)
    else:
      raise AssertionError
    i_proxies_sorted = flex.sort_permutation(data=data_to_sort, reverse=True)
    if (max_items is not None):
      i_proxies_sorted = i_proxies_sorted[:max_items]
    print >> f, "%sSorted by %s:" % (prefix, by_value)
    for i_proxy in i_proxies_sorted:
      proxy = O[i_proxy]
      len_max = 0
      ls = []
      for i_seq in proxy.i_seqs:
        if (site_labels is None): l = str(i_seq)
        else:                     l = site_labels[i_seq]
        len_max = max(len_max, len(l))
        ls.append(l)
      if unit_cell is None:
        restraint = planarity(sites_cart=sites_cart, proxy=proxy)
        sym_op_label = ""
      else:
        restraint = planarity(unit_cell=unit_cell, sites_cart=sites_cart,
                              proxy=proxy)
        sym_op_label = " sym.op."
      print >> f, \
        "%s      %s    delta    sigma   weight rms_deltas residual%s" % (
          prefix, " "*len_max, sym_op_label)
      s = "plane"
      rdr = None
      for i, (i_seq,weight,delta,l) in enumerate(zip(proxy.i_seqs, proxy.weights,
                                      restraint.deltas(), ls)):
        if (rdr is None):
          rdr = "   %6.2e %6.2e" % (
            restraint.rms_deltas(), restraint.residual())
          rdr_spacer = ""
        sym_op = ""
        if proxy.sym_ops:
          rt_mx = proxy.sym_ops[i]
          if not rt_mx.is_unit_mx():
            sym_op = "%s %s" %(rdr_spacer, rt_mx.as_xyz())
        print >> f, "%s%5s %s  %7.3f %6.2e %6.2e%s%s" % (
          prefix, s, l+" "*(len_max-len(l)),
          delta, weight_as_sigma(weight=weight), weight, rdr, sym_op)
        rdr = ""
        rdr_spacer = " "*20
        s = ""
    n_not_shown = O.size() - i_proxies_sorted.size()
    if (n_not_shown != 0):
      print >> f, prefix + "... (remaining %d not shown)" % n_not_shown

class _shared_bond_similarity_proxy(
  boost.python.injector, shared_bond_similarity_proxy):

  def deltas_rms(self, sites_cart, unit_cell=None):
    if unit_cell is None:
      return bond_similarity_deltas_rms(sites_cart=sites_cart, proxies=self)
    else:
      return bond_similarity_deltas_rms(
        unit_cell=unit_cell, sites_cart=sites_cart, proxies=self)

  def residuals(self, sites_cart, unit_cell=None):
    if unit_cell is None:
      return bond_similarity_residuals(sites_cart=sites_cart, proxies=self)
    else:
      return bond_similarity_residuals(
        unit_cell=unit_cell, sites_cart=sites_cart, proxies=self)

  def show_sorted(O,
      by_value,
      sites_cart,
      site_labels=None,
      unit_cell=None,
      f=None,
      prefix="",
      max_items=None):
    assert by_value in ["residual", "rms_deltas"]
    assert site_labels is None or len(site_labels) == sites_cart.size()
    if (f is None): f = sys.stdout
    print >> f, "%sBond similarity restraints: %d" % (prefix, O.size())
    if (O.size() == 0): return
    if (max_items is not None and max_items <= 0): return
    if (by_value == "residual"):
      if unit_cell is None:
        data_to_sort = O.residuals(sites_cart=sites_cart)
      else:
        data_to_sort = O.residuals(unit_cell=unit_cell,sites_cart=sites_cart)
    elif (by_value == "rms_deltas"):
      if unit_cell is None:
        data_to_sort = O.deltas_rms(sites_cart=sites_cart)
      else:
        data_to_sort = O.deltas_rms(unit_cell=unit_cell, sites_cart=sites_cart)
    else:
      raise AssertionError
    i_proxies_sorted = flex.sort_permutation(data=data_to_sort, reverse=True)
    if (max_items is not None):
      i_proxies_sorted = i_proxies_sorted[:max_items]
    print >> f, "%sSorted by %s:" % (prefix, by_value)
    for i_proxy in i_proxies_sorted:
      proxy = O[i_proxy]
      len_max = 0
      ls = []
      for pair in proxy.i_seqs:
        if (site_labels is None):
          l = "%s-%s" %(str(pair[0]), str(pair[1]))
        else:
          l = "%s-%s" %(site_labels[pair[0]], site_labels[pair[1]])
        len_max = max(len_max, len(l))
        ls.append(l)
      if unit_cell is None:
        restraint = bond_similarity(sites_cart=sites_cart, proxy=proxy)
        sym_op_label = ""
      else:
        restraint = bond_similarity(unit_cell=unit_cell, sites_cart=sites_cart,
                              proxy=proxy)
        sym_op_label = " sym.op."
      print >> f, \
        "%s     %s    delta    sigma   weight rms_deltas residual%s" % (
          prefix, " "*len_max, sym_op_label)
      s = "bond"
      rdr = None
      for i, (i_seq,weight,delta,l) in enumerate(zip(proxy.i_seqs, proxy.weights,
                                      restraint.deltas(), ls)):
        if (rdr is None):
          rdr = "   %6.2e %6.2e" % (
            restraint.rms_deltas(), restraint.residual())
          rdr_spacer = ""
        sym_op = ""
        if proxy.sym_ops:
          rt_mx = proxy.sym_ops[i]
          if not rt_mx.is_unit_mx():
            sym_op = "%s %s" %(rdr_spacer, rt_mx.as_xyz())
        print >> f, "%s%4s %s  %7.3f %6.2e %6.2e%s%s" % (
          prefix, s, l+" "*(len_max-len(l)),
          delta, weight_as_sigma(weight=weight), weight, rdr, sym_op)
        rdr = ""
        rdr_spacer = " "*20
        s = ""
    n_not_shown = O.size() - i_proxies_sorted.size()
    if (n_not_shown != 0):
      print >> f, prefix + "... (remaining %d not shown)" % n_not_shown

def _show_histogram_of_deltas_impl(O,
        proxy_label,
        format_cutoffs,
        unit_cell,
        sites_cart,
        n_slots,
        f,
        prefix):
    if (O.size() == 0): return
    if (f is None): f = sys.stdout
    print >> f, "%sHistogram of %s deviations from ideal:" % (
      prefix, proxy_label)
    if unit_cell is None:
      data = flex.abs(O.deltas(sites_cart=sites_cart))
    else:
      data = flex.abs(
        O.deltas(unit_cell=unit_cell, sites_cart=sites_cart))
    histogram = flex.histogram(
      data=data,
      n_slots=n_slots)
    fmt = "%%s  %s - %s: %%d" % (format_cutoffs, format_cutoffs)
    low_cutoff = histogram.data_min()
    for i,n in enumerate(histogram.slots()):
      high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
      print >> f, fmt % (prefix, low_cutoff, high_cutoff, n)
      low_cutoff = high_cutoff
    return histogram

def _get_sorted_impl(O,
        proxy_type,
        by_value,
        unit_cell,
        sites_cart,
        site_labels,
        max_items,
        get_restraints_only=True):
  assert by_value in ["residual", "delta"]
  assert site_labels is None or len(site_labels) == sites_cart.size()
  if (O.size() == 0): return None, None
  if (max_items is not None and max_items <= 0): return None, None
  if (by_value == "residual"):
    if unit_cell is None:
      data_to_sort = O.residuals(sites_cart=sites_cart)
    else:
      data_to_sort = O.residuals(unit_cell=unit_cell, sites_cart=sites_cart)
  elif (by_value == "delta"):
    if unit_cell is None:
      data_to_sort = flex.abs(O.deltas(sites_cart=sites_cart))
    else:
      data_to_sort = flex.abs(
        O.deltas(unit_cell=unit_cell, sites_cart=sites_cart))
  else:
    raise AssertionError
  i_proxies_sorted = flex.sort_permutation(data=data_to_sort, reverse=True)
  if (max_items is not None):
    i_proxies_sorted = i_proxies_sorted[:max_items]
  sorted_table = []
  for i_proxy in i_proxies_sorted:
    proxy = O[i_proxy]
    labels = []
    for n, i_seq in enumerate(proxy.i_seqs):
      if (site_labels is None): l = str(i_seq)
      else:                     l = site_labels[i_seq]
      if unit_cell and proxy.sym_ops:
        sym_op = proxy.sym_ops[n]
        if not sym_op.is_unit_mx():
          l += "  %s" %sym_op.as_xyz()
      labels.append(l)
    if unit_cell is None:
      restraint = proxy_type(
        sites_cart=sites_cart,
        proxy=proxy)
    else:
      restraint = proxy_type(
        unit_cell=unit_cell,
        sites_cart=sites_cart,
        proxy=proxy)
    if get_restraints_only :
      sorted_table.append((labels, restraint))
    else :
      restraint_info = restraint._get_sorted_item()
      sorted_table.append([labels] + restraint_info)
  n_not_shown = O.size() - i_proxies_sorted.size()
  return sorted_table, n_not_shown

def _show_sorted_impl(O,
        proxy_type,
        proxy_label,
        item_label,
        by_value,
        unit_cell,
        sites_cart,
        site_labels,
        f,
        prefix,
        max_items):
  if (f is None): f = sys.stdout
  sorted_table, n_not_shown = _get_sorted_impl(O,
        proxy_type=proxy_type,
        by_value=by_value,
        unit_cell=unit_cell,
        sites_cart=sites_cart,
        site_labels=site_labels,
        max_items=max_items,
        get_restraints_only=True)
  print >> f, "%s%s restraints: %d" % (prefix, proxy_label, O.size())
  if (O.size() == 0): return
  if (proxy_type is dihedral):
    n_harmonic = O.count_harmonic()
    n_sinusoidal = O.size() - n_harmonic
    print >> f, prefix+"  sinusoidal: %d" % n_sinusoidal
    print >> f, prefix+"    harmonic: %d" % n_harmonic
  if (max_items is not None and max_items <= 0): return
  item_label_blank = " " * len(item_label)
  print >> f, "%sSorted by %s:" % (prefix, by_value)
  for (labels, restraint) in sorted_table :
    s = item_label
    for l in labels :
      print >> f, "%s%s %s" % (prefix, s, l)
      s = item_label_blank
    restraint._show_sorted_item(f=f, prefix=prefix)
  if (n_not_shown != 0):
    print >> f, prefix + "... (remaining %d not shown)" % n_not_shown

class pair_proxies(object):

  def __init__(self,
        flags=None,
        bond_params_table=None,
        shell_asu_tables=None,
        model_indices=None,
        conformer_indices=None,
        sym_excl_indices=None,
        nonbonded_params=None,
        nonbonded_types=None,
        nonbonded_distance_cutoff_plus_buffer=None,
        min_cubicle_edge=5):
    self.bond_proxies = None
    self.nonbonded_proxies = None
    if (bond_params_table is not None
        and (flags is None or flags.bond)):
      if (shell_asu_tables is None):
        self.bond_proxies = bond_sorted_asu_proxies(
          bond_params_table=bond_params_table)
      else:
        assert len(shell_asu_tables) > 0
        self.bond_proxies = bond_sorted_asu_proxies(
          bond_params_table=bond_params_table,
          bond_asu_table=shell_asu_tables[0])
    if (nonbonded_types is not None
        and (flags is None or flags.nonbonded)):
      assert nonbonded_params is not None
      assert nonbonded_distance_cutoff_plus_buffer is not None
      assert shell_asu_tables is not None
      assert len(shell_asu_tables) > 0
      self.nonbonded_proxies = nonbonded_sorted_asu_proxies(
        model_indices=model_indices,
        conformer_indices=conformer_indices,
        sym_excl_indices=sym_excl_indices,
        nonbonded_params=nonbonded_params,
        nonbonded_types=nonbonded_types,
        nonbonded_distance_cutoff_plus_buffer=\
          nonbonded_distance_cutoff_plus_buffer,
        min_cubicle_edge=min_cubicle_edge,
        shell_asu_tables=shell_asu_tables)

class _motif(boost.python.injector, ext.motif):

  def show(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    print >> out, prefix+"geometry_restraints.motif {"
    print >> out, prefix+"  id = %s" % show_string(self.id)
    print >> out, prefix+"  description = %s" % show_string(self.description)
    for info in self.info:
      print >> out, prefix+"  info = %s" % show_string(info)
    for manipulation_id in self.manipulation_ids:
      print >> out, prefix+"  manipulation_id = %s" % (
        show_string(manipulation_id))
    self.show_atoms(out=out, prefix=prefix+"  ")
    self.show_bonds(out=out, prefix=prefix+"  ")
    self.show_angles(out=out, prefix=prefix+"  ")
    self.show_dihedrals(out=out, prefix=prefix+"  ")
    self.show_chiralities(out=out, prefix=prefix+"  ")
    self.show_planarities(out=out, prefix=prefix+"  ")
    print >> out, prefix+"}"

  def show_atoms(self, out=None, prefix=""):
    atoms = self.atoms_as_list()
    if (len(atoms) > 0):
      print >> out, prefix+"atom = " \
        "[name scattering_type nonbonded_type partial_charge]"
      for atom in atoms:
        print >> out, prefix+"atom = %s %s %s %.6g" % (
          show_string(atom.name),
          show_string(atom.scattering_type),
          show_string(atom.nonbonded_type),
          atom.partial_charge)

  def show_bonds(self, out=None, prefix=""):
    bonds = self.bonds_as_list()
    if (len(bonds) > 0):
      print >> out, prefix+"bond = " \
        "[atom_name*2 type distance_ideal weight id]"
      for bond in bonds:
        atom_names = bond.atom_names
        print >> out, prefix+"bond = %s %s %s %.6g %.6g %s" % (
          show_string(atom_names[0]),
          show_string(atom_names[1]),
          show_string(bond.type),
          bond.distance_ideal,
          bond.weight,
          show_string(bond.id))

  def show_angles(self, out=None, prefix=""):
    angles = self.angles_as_list()
    if (len(angles) > 0):
      print >> out, prefix+"angle = " \
        "[atom_name*3 angle_ideal weight id]"
      for angle in angles:
        atom_names = angle.atom_names
        print >> out, prefix+"angle = %s %s %s %.6g %.6g %s" % (
          show_string(atom_names[0]),
          show_string(atom_names[1]),
          show_string(atom_names[2]),
          angle.angle_ideal,
          angle.weight,
          show_string(angle.id))

  def show_dihedrals(self, out=None, prefix=""):
    dihedrals = self.dihedrals_as_list()
    if (len(dihedrals) > 0):
      print >> out, prefix+"dihedral = " \
        "[atom_name*4 angle_ideal weight periodicity id]"
      for dihedral in dihedrals:
        atom_names = dihedral.atom_names
        print >> out, prefix+"dihedral = %s %s %s %s %.6g %.6g %d %s" % (
          show_string(atom_names[0]),
          show_string(atom_names[1]),
          show_string(atom_names[2]),
          show_string(atom_names[3]),
          dihedral.angle_ideal,
          dihedral.weight,
          dihedral.periodicity,
          show_string(dihedral.id))

  def show_chiralities(self, out=None, prefix=""):
    chiralities = self.chiralities_as_list()
    if (len(chiralities) > 0):
      print >> out, prefix+"chirality = " \
        "[atom_name*4 volume_sign both_signs volume_ideal weight id]"
      for chirality in chiralities:
        atom_names = chirality.atom_names
        if (chirality.both_signs): both_signs = "True"
        else:                      both_signs = "False"
        print >> out, prefix+"chirality = %s %s %s %s %s %s %.6g %.6g %s" % (
          show_string(atom_names[0]),
          show_string(atom_names[1]),
          show_string(atom_names[2]),
          show_string(atom_names[3]),
          show_string(chirality.volume_sign),
          both_signs,
          chirality.volume_ideal,
          chirality.weight,
          show_string(chirality.id))

  def show_planarities(self, out=None, prefix=""):
    planarities = self.planarities_as_list()
    if (len(planarities) > 0):
      for planarity in planarities:
        print >> out, prefix+"planarity {"
        print >> out, prefix+"  id = %s" % show_string(planarity.id)
        assert planarity.weights.size() == planarity.atom_names.size()
        print >> out, prefix+"  atom = [name weight]"
        for an,w in zip(planarity.atom_names, planarity.weights):
          print >> out, prefix+"  atom = %s %.6g" % (show_string(an), w)
        print >> out, prefix+"}"

class _motif_alteration(boost.python.injector, ext.motif_alteration):

  def show(self, out=None, prefix="", previous_help=None):
    if (out is None): out = sys.stdout
    action = self.action
    operand = self.operand
    if (operand == "atom"):
      assert len(self.motif_ids) == 1
      atom = self.atom
      attr = "name scattering_type nonbonded_type partial_charge"
      if (action == "add"):
        help = prefix+"atom = add [motif_id %s]" % attr
        if (help != previous_help): print >> out, help
        print >> out, prefix+"atom = add %s %s %s %s %s" % (
          show_string(self.motif_ids[0]),
          show_string(atom.name),
          show_string(atom.scattering_type),
          show_string(atom.nonbonded_type),
          atom.partial_charge)
      elif (action == "change"):
        help = prefix+"atom = change [motif_id motif_atom_name \\\n" \
                    + prefix+"               %s]" % attr
        if (help != previous_help): print >> out, help
        print >> out, prefix+"atom = change %s %s \\" % (
          show_string(self.motif_ids[0]),
          show_string(self.motif_atom_name))
        if (not self.change_partial_charge()):
          partial_charge = "None"
        else:
          partial_charge = "%.6g" % atom.partial_charge
        print >> out, prefix+"              %s %s %s %s" % (
          show_string(atom.name),
          show_string(atom.scattering_type),
          show_string(atom.nonbonded_type),
          partial_charge)
      else:
        assert action == "delete"
        help = prefix+"atom = delete [motif_id motif_atom_name]"
        if (help != previous_help): print >> out, help
        print >> out, prefix+"atom = delete %s %s" % (
          show_string(self.motif_ids[0]),
          show_string(self.motif_atom_name))
    elif (operand == "bond"):
      assert len(self.motif_ids) == 2
      bond = self.bond
      help_lead = "bond = %s [(motif_id atom_name)*2" % action
      data_lead = "bond = %s %s %s %s %s" % (
        action,
        show_string(self.motif_ids[0]),
        show_string(bond.atom_names[0]),
        show_string(self.motif_ids[1]),
        show_string(bond.atom_names[1]))
      if (action == "delete"):
        help = prefix+"%s]" % help_lead
        if (help != previous_help): print >> out, help
        print >> out, prefix+"%s" % data_lead
      else:
        if (action == "change" and not self.change_distance_ideal()):
          distance_ideal = "None"
        else:
          distance_ideal = "%.6g" % bond.distance_ideal
        if (action == "change" and not self.change_weight()):
          weight = "None"
        else:
          weight = "%.6g" % bond.weight
        help = prefix+"%s type distance_ideal weight id]" % help_lead
        if (help != previous_help): print >> out, help
        print >> out, prefix+"%s %s %s %s" % (
          data_lead, distance_ideal, weight, show_string(bond.id))
    elif (operand == "angle"):
      assert len(self.motif_ids) == 3
      angle = self.angle
      help_lead = "angle = %s [(motif_id atom_name)*3" % action
      data_lead = "angle = %s %s %s %s %s %s %s" % (
        action,
        show_string(self.motif_ids[0]),
        show_string(angle.atom_names[0]),
        show_string(self.motif_ids[1]),
        show_string(angle.atom_names[1]),
        show_string(self.motif_ids[2]),
        show_string(angle.atom_names[2]))
      if (action == "delete"):
        help = prefix+"%s]" % help_lead
        if (help != previous_help): print >> out, help
        print >> out, prefix+"%s" % data_lead
      else:
        if (action == "change" and not self.change_angle_ideal()):
          angle_ideal = "None"
        else:
          angle_ideal = "%.6g" % angle.angle_ideal
        if (action == "change" and not self.change_weight()):
          weight = "None"
        else:
          weight = "%.6g" % angle.weight
        help = prefix+"%s type angle_ideal weight id]" % help_lead
        if (help != previous_help): print >> out, help
        print >> out, prefix+"%s %s %s %s" % (
          data_lead, angle_ideal, weight, show_string(angle.id))
    elif (operand == "dihedral"):
      assert len(self.motif_ids) == 4
      dihedral = self.dihedral
      help_lead = "dihedral = %s [(motif_id atom_name)*4" % action
      data_lead = "dihedral = %s %s %s %s %s %s %s %s %s" % (
        action,
        show_string(self.motif_ids[0]),
        show_string(dihedral.atom_names[0]),
        show_string(self.motif_ids[1]),
        show_string(dihedral.atom_names[1]),
        show_string(self.motif_ids[2]),
        show_string(dihedral.atom_names[2]),
        show_string(self.motif_ids[3]),
        show_string(dihedral.atom_names[3]))
      if (action == "delete"):
        help = prefix+"%s]" % help_lead
        if (help != previous_help): print >> out, help
        print >> out, prefix+"%s" % data_lead
      else:
        if (action == "change" and not self.change_angle_ideal()):
          angle_ideal = "None"
        else:
          angle_ideal = "%.6g" % dihedral.angle_ideal
        if (action == "change" and not self.change_weight()):
          weight = "None"
        else:
          weight = "%.6g" % dihedral.weight
        if (action == "change" and not self.change_periodicity()):
          periodicity = "None"
        else:
          periodicity = "%d" % dihedral.periodicity
        help = prefix+"%s angle_ideal weight periodicity id]" % help_lead
        if (help != previous_help): print >> out, help
        print >> out, prefix+"%s %s %s %s %s" % (
          data_lead, angle_ideal, weight, periodicity,
          show_string(dihedral.id))
    elif (operand == "chirality"):
      assert len(self.motif_ids) == 4
      chirality = self.chirality
      help_lead = "chirality = %s [(motif_id atom_name)*4" % action
      data_lead = "chirality = %s %s %s %s %s %s %s %s %s" % (
        action,
        show_string(self.motif_ids[0]),
        show_string(chirality.atom_names[0]),
        show_string(self.motif_ids[1]),
        show_string(chirality.atom_names[1]),
        show_string(self.motif_ids[2]),
        show_string(chirality.atom_names[2]),
        show_string(self.motif_ids[3]),
        show_string(chirality.atom_names[3]))
      if (action == "delete"):
        help = prefix+"%s]" % help_lead
        if (help != previous_help): print >> out, help
        print >> out, prefix+"%s" % data_lead
      else:
        if (action == "change" and not self.change_volume_ideal()):
          volume_ideal = "None"
        else:
          volume_ideal = "%.6g" % chirality.volume_ideal
        if (action == "change" and not self.change_weight()):
          weight = "None"
        else:
          weight = "%.6g" % chirality.weight
        help = prefix+"%s \\\n" % help_lead \
             + prefix+" "*(14+len(action)) \
             + "volume_sign volume_ideal weight id]"
        if (help != previous_help): print >> out, help
        print >> out, prefix+"%s \\\n%s%s%s %s %s %s" % (
          data_lead, prefix, " "*(13+len(action)),
          show_string(chirality.volume_sign),
          volume_ideal, weight, show_string(chirality.id))
    elif (operand == "planarity"):
      planarity = self.planarity
      print >> out, prefix+"planarity {"
      print >> out, prefix+"  action = %s" % action
      print >> out, prefix+"  motif_id = %s" % show_string(
        self.planarity_motif_id)
      print >> out, prefix+"  id = %s" % show_string(planarity.id)
      if (action == "add"):
        print >> out, prefix+"  atom = [motif_id name weight]"
        assert planarity.weights.size() == planarity.atom_names.size()
        assert self.motif_ids.size() == planarity.atom_names.size()
        for mi,an,w in zip(self.motif_ids,
                           planarity.atom_names,
                           planarity.weights):
          print >> out, prefix+"  atom = %s %s %.6g" % (
            show_string(mi), show_string(an), w)
      elif (action == "change"):
        assert planarity.weights.size() == planarity.atom_names.size()
        assert self.motif_ids.size() == planarity.atom_names.size()
        actions = self.planarity_atom_actions_as_list()
        assert len(actions) == planarity.atom_names.size()
        previous_help = None
        for ac,mi,an,w in zip(actions,
                              self.motif_ids,
                              planarity.atom_names,
                              planarity.weights):
          if (ac != "delete"):
            help = prefix+"  atom = %s [motif_id name weight]" % ac
            if (help != previous_help): print >> out, help
            print >> out, prefix+"  atom = %s %s %s %.6g" % (
              ac, show_string(mi), show_string(an), w)
          else:
            help = prefix+"  atom = delete [motif_id name]"
            if (help != previous_help): print >> out, help
            print >> out, prefix+"  atom = %s %s %s" % (
              ac, show_string(mi), show_string(an))
          previous_help = help
      print >> out, prefix+"}"
      help = None
    else:
      raise RuntimeError("Internal Error: unknown operand: %s" % operand)
    return help

class _motif_manipulation(boost.python.injector, ext.motif_manipulation):

  def show(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    print >> out, prefix+"geometry_restraints.motif_manipulation {"
    print >> out, prefix+"  id = %s" % show_string(self.id)
    print >> out, prefix+"  description = %s" % show_string(self.description)
    for info in self.info:
      print >> out, prefix+"  info = %s" % show_string(info)
    previous_help = None
    for alteration in self.alterations_as_list():
      previous_help = alteration.show(
        out=out, prefix=prefix+"  ", previous_help=previous_help)
    print >> out, prefix+"}"
