import cctbx.crystal.direct_space_asu
from cctbx.array_family import flex
import scitbx.array_family.shared
from scitbx.python_utils.misc import adopt_init_args

import boost.python
ext = boost.python.import_ext("cctbx_geometry_restraints_ext")
from cctbx_geometry_restraints_ext import *

import scitbx.stl.map

nonbonded_radius_table = scitbx.stl.map.stl_string_double

nonbonded_distance_table = scitbx.stl.map.stl_string_stl_map_stl_string_double
nonbonded_distance_dict = scitbx.stl.map.stl_string_double

class proxy_registry_add_result:

  def __init__(self, tabulated_proxy=None, is_new=00000, is_conflicting=00000):
    adopt_init_args(self, locals())

class angle_proxy_registry:

  def __init__(self):
    self.table = {}
    self.proxies = shared_angle_proxy()
    self.counts = flex.size_t()

  def add(self, proxy, tolerance=1.e-6):
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
      if (   abs(result.tabulated_proxy.angle_ideal - proxy.angle_ideal)
               > tolerance
          or abs(result.tabulated_proxy.weight - proxy.weight)
               > tolerance):
        result.is_conflicting = 0001
      else:
        self.counts[i_list] += 1
    return result
