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

  def __init__(self, tabulated_proxy=None, is_new=00000, is_conflicting=00000):
    adopt_init_args(self, locals())

class angle_proxy_registry:

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
      result.is_new = 0001
    else:
      i_list = tab_i_seq_1[i_seqs_0_2]
      result.tabulated_proxy = self.proxies[i_list]
      if (   abs(angle_delta_deg(result.tabulated_proxy.angle_ideal,
                                 proxy.angle_ideal)) > tolerance
          or abs(result.tabulated_proxy.weight - proxy.weight)
               > tolerance):
        result.is_conflicting = 0001
      else:
        self.counts[i_list] += 1
    return result

class dihedral_proxy_registry:

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
      result.is_new = 0001
    else:
      i_list = tab_i_seq_0[i_seqs_1_2_3]
      result.tabulated_proxy = self.proxies[i_list]
      if (   abs(angle_delta_deg(result.tabulated_proxy.angle_ideal,
                                 proxy.angle_ideal,
                                 proxy.periodicity)) > tolerance
          or abs(result.tabulated_proxy.weight - proxy.weight)
               > tolerance
          or result.tabulated_proxy.periodicity != proxy.periodicity):
        result.is_conflicting = 0001
      else:
        self.counts[i_list] += 1
    return result

class chirality_proxy_registry:

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
      result.is_new = 0001
    else:
      i_list = tab_i_seq_0[i_seqs_1_2_3]
      result.tabulated_proxy = self.proxies[i_list]
      if (   abs(result.tabulated_proxy.volume_ideal - proxy.volume_ideal)
               > tolerance
          or abs(result.tabulated_proxy.weight - proxy.weight)
               > tolerance):
        result.is_conflicting = 0001
      else:
        self.counts[i_list] += 1
    return result

class planarity_proxy_registry:

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
      result.is_new = 0001
    else:
      i_list = tab_i_seq_0[i_seqs_1_up]
      result.tabulated_proxy = self.proxies[i_list]
      if (not approx_equal(result.tabulated_proxy.weights, proxy.weights,
                           eps=tolerance)):
        result.is_conflicting = 0001
      else:
        self.counts[i_list] += 1
    return result
