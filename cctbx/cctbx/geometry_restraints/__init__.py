import cctbx.crystal.direct_space_asu
from cctbx.array_family import flex
import scitbx.array_family.shared
from scitbx.python_utils.misc import adopt_init_args
from libtbx.test_utils import approx_equal

import boost.python
ext = boost.python.import_ext("cctbx_geometry_restraints_ext")
from cctbx_geometry_restraints_ext import *

import scitbx.stl.map
from stdlib import math
import sys

nonbonded_radius_table = scitbx.stl.map.stl_string_double

nonbonded_distance_table = scitbx.stl.map.stl_string_stl_map_stl_string_double
nonbonded_distance_dict = scitbx.stl.map.stl_string_double

def angle_delta_deg(angle_1, angle_2, periodicity=1):
  half_period = 180./max(1,periodicity)
  d = math.fmod(angle_2-angle_1, 2*half_period)
  if   (d < -half_period): d += 2*half_period
  elif (d >  half_period): d -= 2*half_period
  return d

class proxy_registry_add_result:

  def __init__(self, tabulated_proxy=None, is_new=False, is_conflicting=False):
    adopt_init_args(self, locals())

class proxy_registry_mixin:

  def discard_table(self):
    self.table = None

class bond_simple_proxy_registry(proxy_registry_mixin):

  def __init__(self, n_seq):
    self.table = [{} for i in xrange(n_seq)]
    self.proxies = shared_bond_simple_proxy()
    self.counts = flex.size_t()

  def process(self, proxy, tolerance=1.e-6):
    result = proxy_registry_add_result()
    proxy = proxy.sort_i_seqs()
    if (not self.table[proxy.i_seqs[0]].has_key(proxy.i_seqs[1])):
      self.table[proxy.i_seqs[0]][proxy.i_seqs[1]] = self.proxies.size()
      self.proxies.append(proxy)
      self.counts.append(1)
      result.tabulated_proxy = proxy
      result.is_new = True
    else:
      i_list = self.table[proxy.i_seqs[0]][proxy.i_seqs[1]]
      result.tabulated_proxy = self.proxies[i_list]
      if (   abs(result.tabulated_proxy.distance_ideal
                                - proxy.distance_ideal) > tolerance
          or abs(result.tabulated_proxy.weight - proxy.weight)
               > tolerance):
        result.is_conflicting = True
      else:
        self.counts[i_list] += 1
    return result

class angle_proxy_registry(proxy_registry_mixin):

  def __init__(self):
    self.table = {}
    self.proxies = shared_angle_proxy()
    self.counts = flex.size_t()

  def process(self, proxy, tolerance=1.e-6):
    result = proxy_registry_add_result()
    proxy = proxy.sort_i_seqs()
    tab_i_seq_1 = self.table.setdefault(proxy.i_seqs[1], {})
    i_seqs_0_2 = (proxy.i_seqs[0], proxy.i_seqs[2])
    if (not tab_i_seq_1.has_key(i_seqs_0_2)):
      tab_i_seq_1[i_seqs_0_2] = self.proxies.size()
      self.proxies.append(proxy)
      self.counts.append(1)
      result.tabulated_proxy = proxy
      result.is_new = True
    else:
      i_list = tab_i_seq_1[i_seqs_0_2]
      result.tabulated_proxy = self.proxies[i_list]
      if (   abs(angle_delta_deg(result.tabulated_proxy.angle_ideal,
                                 proxy.angle_ideal)) > tolerance
          or abs(result.tabulated_proxy.weight - proxy.weight)
               > tolerance):
        result.is_conflicting = True
      else:
        self.counts[i_list] += 1
    return result

class dihedral_proxy_registry(proxy_registry_mixin):

  def __init__(self):
    self.table = {}
    self.proxies = shared_dihedral_proxy()
    self.counts = flex.size_t()

  def process(self, proxy, tolerance=1.e-6):
    result = proxy_registry_add_result()
    proxy = proxy.sort_i_seqs()
    tab_i_seq_0 = self.table.setdefault(proxy.i_seqs[0], {})
    i_seqs_1_2_3 = (proxy.i_seqs[1], proxy.i_seqs[2], proxy.i_seqs[3])
    if (not tab_i_seq_0.has_key(i_seqs_1_2_3)):
      tab_i_seq_0[i_seqs_1_2_3] = self.proxies.size()
      self.proxies.append(proxy)
      self.counts.append(1)
      result.tabulated_proxy = proxy
      result.is_new = True
    else:
      i_list = tab_i_seq_0[i_seqs_1_2_3]
      result.tabulated_proxy = self.proxies[i_list]
      if (   abs(angle_delta_deg(result.tabulated_proxy.angle_ideal,
                                 proxy.angle_ideal,
                                 proxy.periodicity)) > tolerance
          or abs(result.tabulated_proxy.weight - proxy.weight)
               > tolerance
          or result.tabulated_proxy.periodicity != proxy.periodicity):
        result.is_conflicting = True
      else:
        self.counts[i_list] += 1
    return result

class chirality_proxy_registry(proxy_registry_mixin):

  def __init__(self):
    self.table = {}
    self.proxies = shared_chirality_proxy()
    self.counts = flex.size_t()

  def process(self, proxy, tolerance=1.e-6):
    result = proxy_registry_add_result()
    proxy = proxy.sort_i_seqs()
    tab_i_seq_0 = self.table.setdefault(proxy.i_seqs[0], {})
    i_seqs_1_2_3 = (proxy.i_seqs[1], proxy.i_seqs[2], proxy.i_seqs[3])
    if (not tab_i_seq_0.has_key(i_seqs_1_2_3)):
      tab_i_seq_0[i_seqs_1_2_3] = self.proxies.size()
      self.proxies.append(proxy)
      self.counts.append(1)
      result.tabulated_proxy = proxy
      result.is_new = True
    else:
      i_list = tab_i_seq_0[i_seqs_1_2_3]
      result.tabulated_proxy = self.proxies[i_list]
      if (   abs(result.tabulated_proxy.volume_ideal - proxy.volume_ideal)
               > tolerance
          or abs(result.tabulated_proxy.weight - proxy.weight)
               > tolerance):
        result.is_conflicting = True
      else:
        self.counts[i_list] += 1
    return result

class planarity_proxy_registry(proxy_registry_mixin):

  def __init__(self):
    self.table = {}
    self.proxies = shared_planarity_proxy()
    self.counts = flex.size_t()

  def process(self, proxy, tolerance=1.e-6):
    assert proxy.i_seqs.size() > 0
    result = proxy_registry_add_result()
    proxy = proxy.sort_i_seqs()
    tab_i_seq_0 = self.table.setdefault(proxy.i_seqs[0], {})
    i_seqs_1_up = tuple(proxy.i_seqs[1:])
    if (not tab_i_seq_0.has_key(i_seqs_1_up)):
      tab_i_seq_0[i_seqs_1_up] = self.proxies.size()
      self.proxies.append(proxy)
      self.counts.append(1)
      result.tabulated_proxy = proxy
      result.is_new = True
    else:
      i_list = tab_i_seq_0[i_seqs_1_up]
      result.tabulated_proxy = self.proxies[i_list]
      if (not approx_equal(result.tabulated_proxy.weights, proxy.weights,
                           eps=tolerance)):
        result.is_conflicting = True
      else:
        self.counts[i_list] += 1
    return result

class _bond_sorted_asu_proxies(boost.python.injector, bond_sorted_asu_proxies):

  def show_sorted_by_residual(self,
        sites_cart,
        labels=None,
        f=None,
        prefix="",
        max_lines=None):
    assert labels is None or len(labels) == sites_cart.size()
    if (self.n_total() == 0): return
    if (f is None): f = sys.stdout
    residuals = bond_residuals(
      sites_cart=sites_cart,
      sorted_asu_proxies=self)
    i_proxies_sorted = flex.sort_permutation(data=residuals, reverse=True)
    if (max_lines is not None and i_proxies_sorted.size() > max_lines+1):
      i_proxies_sorted = i_proxies_sorted[:max_lines]
    n_asu = 0
    n_simple = self.simple.size()
    if (labels is not None):
      max_label_lengths = [0,0]
      for i_proxy in i_proxies_sorted:
        if (i_proxy < n_simple):
          i_seqs = self.simple[i_proxy].i_seqs
        else:
          proxy = self.asu[i_proxy-n_simple]
          i_seqs = (proxy.i_seq, proxy.j_seq)
          n_asu += 1
        max_label_lengths = [max(m, len(labels[i_seq]))
          for m,i_seq in zip(max_label_lengths, i_seqs)]
      if (max(max_label_lengths) == 0):
        labels = None
      else:
        label_label_format = \
          " - ".join(["%%-%ds" % m for m in max_label_lengths]) + " "
        atom_i_atom_j_format = label_label_format.replace("-", "", 1)
    if (self.asu.size() == 0):
      asu_mappings = None
    else:
      asu_mappings = self.asu_mappings()
    show_legend = True
    atom_i_atom_j = ""
    label_label = ""
    for i_proxy in i_proxies_sorted:
      if (i_proxy < n_simple):
        proxy = self.simple[i_proxy]
        i_seq,j_seq = proxy.i_seqs
        rt_mx = None
        restraint = bond(
          sites_cart=sites_cart,
          proxy=proxy)
      else:
        proxy = self.asu[i_proxy-n_simple]
        i_seq,j_seq = proxy.i_seq,proxy.j_seq
        rt_mx = asu_mappings.get_rt_mx_i(pair=proxy).inverse().multiply(
                asu_mappings.get_rt_mx_j(pair=proxy))
        restraint = bond(
          sites_cart=sites_cart,
          asu_mappings=asu_mappings,
          proxy=proxy)
      if (show_legend):
        show_legend = False
        if (labels is not None):
          if (min(max_label_lengths) >= 6):
            atom_i_atom_j = atom_i_atom_j_format % ("atom i", "atom j")
          elif (min(max_label_lengths) >= 1):
            atom_i_atom_j = atom_i_atom_j_format % ("i", "j")
        print >> f, "%sBonded interactions sorted by residual:" % prefix
        print >> f, "%s%sideal  model  delta   weight residual" % (
          prefix, atom_i_atom_j),
        if (n_asu > 0):
          print >> f, "sym.op. j",
        print >> f
      if (labels is not None):
        label_label = label_label_format % (labels[i_seq], labels[j_seq])
      print >> f, "%s%s%5.3f %6.3f %6.3f %6.2e %6.2e" % (
        prefix, label_label,
        restraint.distance_ideal, restraint.distance_model, restraint.delta,
        restraint.weight, restraint.residual()),
      if (rt_mx is not None):
        print >> f, rt_mx,
      print >> f
    n_not_shown = residuals.size() - i_proxies_sorted.size()
    if (n_not_shown != 0):
      print >> f, prefix + "... (remaining %d not shown)" % n_not_shown
