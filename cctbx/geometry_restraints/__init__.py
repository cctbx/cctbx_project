from __future__ import absolute_import, division, print_function
import cctbx.crystal.direct_space_asu # import dependency
from cctbx.array_family import flex
import scitbx.array_family.shared # import dependency
import cctbx.geometry # import dependency
from libtbx.test_utils import approx_equal
from libtbx.str_utils import show_string
from libtbx import group_args

import boost_adaptbx.boost.python as bp
from six.moves import range
from six.moves import zip
ext = bp.import_ext("cctbx_geometry_restraints_ext")
from cctbx_geometry_restraints_ext import *

import scitbx.stl.map
import math
import sys
import six

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

@bp.inject_into(bond_params_table)
class _():

  def lookup(self, i_seq, j_seq):
    if (i_seq > j_seq): i_seq, j_seq = j_seq, i_seq
    bond_params_by_j_seq = self[i_seq]
    if (j_seq not in bond_params_by_j_seq): return None
    return bond_params_by_j_seq[j_seq]

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

  def _append_proxy(self, source_info, proxy, process_result):
    self.proxies.append(proxy)
    self.source_labels.append(source_info.labels())
    self.source_n_expected_atoms.append(source_info.n_expected_atoms())
    process_result.tabulated_proxy = proxy
    process_result.is_new = True

  def _handle_conflict(self,
        source_info,
        proxy,
        i_list,
        process_result):
    # Almost never called.
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
  """
  self.table:
  [ {iseq1: Nproxy} ] , index in this array is iseq0
  """

  def __init__(self, n_seq, strict_conflict_handling):
    proxy_registry_base.__init__(self,
      proxies=shared_bond_simple_proxy(),
      strict_conflict_handling=strict_conflict_handling)
    self.n_seq = n_seq

  def expand_with_ncs(self, nrgl, masters_and_rest_iselection):
    # print("original proxies:", [p.i_seqs for p in self.proxies])
    # nrgl._show(brief=False)
    additional_proxies = []
    for i in range(len(self.proxies)):
      p = self.proxies[i]
      new_current_proxy_iseqs = (masters_and_rest_iselection[p.i_seqs[0]],
                                 masters_and_rest_iselection[p.i_seqs[1]])
      new_master_p = bond_simple_proxy(
          i_seqs = new_current_proxy_iseqs,
          distance_ideal=p.distance_ideal,
          weight=p.weight,
          slack=p.slack,
          limit=p.limit,
          top_out=p.top_out,
          origin_id=p.origin_id).sort_i_seqs()
      self.proxies[i] = new_master_p
      all_new_iseqs = nrgl.get_copy_iseqs(new_current_proxy_iseqs)
      # print('  i, p.i_seqs', i, p.i_seqs)
      # print('  all_new_iseqs', all_new_iseqs)
      if all_new_iseqs is not None:
        for new_iseqs in all_new_iseqs:
          new_proxy = bond_simple_proxy(
              i_seqs=new_iseqs,
              distance_ideal=p.distance_ideal,
              weight=p.weight,
              slack=p.slack,
              limit=p.limit,
              top_out=p.top_out,
              origin_id=p.origin_id).sort_i_seqs()
          # marking table
          self.table[new_iseqs[0]][new_iseqs[1]] = self.proxies.size()
          # ~ self._append_proxy
          additional_proxies.append(new_proxy)
          self.source_labels.append(self.source_labels[i])
          self.source_n_expected_atoms.append(self.source_n_expected_atoms[i])
    self.proxies.extend(additional_proxies)


  def initialize_table(self):
    proxy_registry_base.initialize_table(self)
    self.table = [{} for i in range(self.n_seq)]

  def is_proxy_set(self, i_seqs):
    return (i_seqs[1] in self.table[i_seqs[0]])

  def is_any_proxy_set(self, i_seqs, j_seqs):
    for i in i_seqs:
      for j in j_seqs:
        tmp = [i,j]
        tmp.sort()
        ps = self.is_proxy_set(tmp)
        if ps: return ps
    return False

  def process(self, source_info, proxy, tolerance=1.e-6, replace_in_place=False):
    result = proxy_registry_process_result()
    proxy = proxy.sort_i_seqs()
    if (proxy.i_seqs[1] not in self.table[proxy.i_seqs[0]]):
      self.table[proxy.i_seqs[0]][proxy.i_seqs[1]] = self.proxies.size()
      self._append_proxy(
        source_info=source_info,
        proxy=proxy,
        process_result=result)
    elif replace_in_place:
      i_list = self.table[proxy.i_seqs[0]][proxy.i_seqs[1]]
      tabulated_proxy = self.proxies[i_list]
      for attr in ['distance_ideal',
                   'limit',
                   'slack',
                   'top_out',
                   'weight']:
        setattr(tabulated_proxy, attr, getattr(proxy, attr))
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
    return result

class angle_proxy_registry(proxy_registry_base):
  #
  """
  self.table: nested dicts
    { iseq1:
      {
        (iseq0, iseq2) : Nproxy
      }
    }
  """

  def __init__(self, strict_conflict_handling):
    proxy_registry_base.__init__(self,
      proxies=shared_angle_proxy(),
      strict_conflict_handling=strict_conflict_handling)

  def expand_with_ncs(self, nrgl, masters_and_rest_iselection):
    additional_proxies = []
    # n_proxies = len(self.proxies)
    for i in range(len(self.proxies)):
      p = self.proxies[i]
      new_current_proxy_iseqs = [ masters_and_rest_iselection[iseq] for iseq in p.i_seqs]
      new_master_p = angle_proxy(
            i_seqs=new_current_proxy_iseqs,
            proxy=p).sort_i_seqs()
      self.proxies[i] = new_master_p
      all_new_iseqs = nrgl.get_copy_iseqs(new_current_proxy_iseqs)
      if all_new_iseqs is not None:
        for new_iseqs in all_new_iseqs:
          new_proxy = angle_proxy(
              i_seqs=new_iseqs,
              proxy=p).sort_i_seqs()
          # marking table
          tab_i_seq_1 = self.table.setdefault(new_proxy.i_seqs[1], {})
          tab_i_seq_1[(new_proxy.i_seqs[0], new_proxy.i_seqs[2])] = self.proxies.size()
          # ~ self._append_proxy
          additional_proxies.append(new_proxy)
          # self.proxies.append(new_proxy)
          self.source_labels.append(self.source_labels[i])
          self.source_n_expected_atoms.append(self.source_n_expected_atoms[i])
    self.proxies.extend(additional_proxies)

  def add_if_not_duplicated(self, proxy, tolerance=1.e-6):
    assert len(proxy.i_seqs) == 3
    proxy = proxy.sort_i_seqs()
    tab_i_seq_1 = self.table.setdefault(proxy.i_seqs[1], {})
    i_seqs_0_2 = (proxy.i_seqs[0], proxy.i_seqs[2])
    if (i_seqs_0_2 not in tab_i_seq_1):
      tab_i_seq_1[i_seqs_0_2] = self.proxies.size()
      self.proxies.append(proxy)
      return True
    return False

  def process(self, source_info, proxy, tolerance=1.e-6, replace_in_place=False):
    result = proxy_registry_process_result()
    proxy = proxy.sort_i_seqs()
    tab_i_seq_1 = self.table.setdefault(proxy.i_seqs[1], {})
    i_seqs_0_2 = (proxy.i_seqs[0], proxy.i_seqs[2])
    if (i_seqs_0_2 not in tab_i_seq_1):
      tab_i_seq_1[i_seqs_0_2] = self.proxies.size()
      self._append_proxy(
        source_info=source_info,
        proxy=proxy,
        process_result=result)
    elif replace_in_place:
      i_list = tab_i_seq_1[i_seqs_0_2]
      tabulated_proxy = self.proxies[i_list]
      for attr in ['angle_ideal', 'slack', 'weight']:
        setattr(tabulated_proxy, attr, getattr(proxy, attr))
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
  """
  self.table - similar to angle:
    { iseq0:
      {
        (iseq1, iseq2, iseq3) : Nproxy
      }
    }
  """

  def __init__(self, strict_conflict_handling):
    proxy_registry_base.__init__(self,
      proxies=shared_dihedral_proxy(),
      strict_conflict_handling=strict_conflict_handling)

  def expand_with_ncs(self, nrgl, masters_and_rest_iselection):
    additional_proxies = []
    for i in range(len(self.proxies)):
      p = self.proxies[i]
      new_current_proxy_iseqs = [ masters_and_rest_iselection[iseq] for iseq in p.i_seqs]
      new_master_p = dihedral_proxy(
            i_seqs=new_current_proxy_iseqs,
            proxy=p).sort_i_seqs()
      self.proxies[i] = new_master_p
      all_new_iseqs = nrgl.get_copy_iseqs(new_current_proxy_iseqs)
      if all_new_iseqs is not None:
        for new_iseqs in all_new_iseqs:
          new_proxy = dihedral_proxy(
              i_seqs=new_iseqs,
              proxy=p).sort_i_seqs()
          # marking table
          tab_i_seq_0 = self.table.setdefault(new_proxy.i_seqs[0], {})
          tab_i_seq_0[(new_proxy.i_seqs[1], new_proxy.i_seqs[2], new_proxy.i_seqs[3])] = self.proxies.size()
          # ~ self._append_proxy
          additional_proxies.append(new_proxy)
          # self.proxies.append(new_proxy)
          self.source_labels.append(self.source_labels[i])
          self.source_n_expected_atoms.append(self.source_n_expected_atoms[i])
    self.proxies.extend(additional_proxies)

  def add_if_not_duplicated(self, proxy, tolerance=1.e-6):
    assert len(proxy.i_seqs) == 4
    proxy = proxy.sort_i_seqs()
    tab_i_seq_0 = self.table.setdefault(proxy.i_seqs[0], {})
    i_seqs_1_2_3 = (proxy.i_seqs[1], proxy.i_seqs[2], proxy.i_seqs[3])
    if (i_seqs_1_2_3 not in tab_i_seq_0):
      tab_i_seq_0[i_seqs_1_2_3] = self.proxies.size()
      self.proxies.append(proxy)
      return True
    return False

  def process(self, source_info, proxy, tolerance=1.e-6, replace_in_place=False):
    result = proxy_registry_process_result()
    proxy = proxy.sort_i_seqs()
    tab_i_seq_0 = self.table.setdefault(proxy.i_seqs[0], {})
    i_seqs_1_2_3 = (proxy.i_seqs[1], proxy.i_seqs[2], proxy.i_seqs[3])
    if (i_seqs_1_2_3 not in tab_i_seq_0):
      tab_i_seq_0[i_seqs_1_2_3] = self.proxies.size()
      self._append_proxy(
        source_info=source_info,
        proxy=proxy,
        process_result=result)
    elif replace_in_place:
      i_list = tab_i_seq_0[i_seqs_1_2_3]
      tabulated_proxy = self.proxies[i_list]
      for attr in ['alt_angle_ideals',
                   'angle_ideal',
                   'limit',
                   'periodicity',
                   'slack',
                   'top_out',
                   'weight']:
        setattr(tabulated_proxy, attr, getattr(proxy, attr))
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
  """
  self.table - same as dihedral:
    { iseq0:
      {
        (iseq1, iseq2, iseq3) : Nproxy
      }
    }
  """

  def __init__(self, strict_conflict_handling):
    proxy_registry_base.__init__(self,
      proxies=shared_chirality_proxy(),
      strict_conflict_handling=strict_conflict_handling)

  def expand_with_ncs(self, nrgl, masters_and_rest_iselection):
    additional_proxies = []
    for i in range(len(self.proxies)):
      p = self.proxies[i]
      new_current_proxy_iseqs = [ masters_and_rest_iselection[iseq] for iseq in p.i_seqs]
      new_master_p = chirality_proxy(
            i_seqs=new_current_proxy_iseqs,
            proxy=p).sort_i_seqs()
      self.proxies[i] = new_master_p
      all_new_iseqs = nrgl.get_copy_iseqs(new_current_proxy_iseqs)
      if all_new_iseqs is not None:
        for new_iseqs in all_new_iseqs:
          new_proxy = chirality_proxy(
              i_seqs=new_iseqs,
              proxy=p).sort_i_seqs()
          # marking table
          tab_i_seq_0 = self.table.setdefault(new_proxy.i_seqs[0], {})
          tab_i_seq_0[(new_proxy.i_seqs[1], new_proxy.i_seqs[2], new_proxy.i_seqs[3])] = self.proxies.size()
          # ~ self._append_proxy
          additional_proxies.append(new_proxy)
          # self.proxies.append(new_proxy)
          self.source_labels.append(self.source_labels[i])
          self.source_n_expected_atoms.append(self.source_n_expected_atoms[i])
    self.proxies.extend(additional_proxies)

  def add_if_not_duplicated(self, proxy, tolerance=1.e-6):
    proxy = proxy.sort_i_seqs()
    tab_i_seq_0 = self.table.setdefault(proxy.i_seqs[0], {})
    i_seqs_1_2_3 = (proxy.i_seqs[1], proxy.i_seqs[2], proxy.i_seqs[3])
    if (i_seqs_1_2_3 not in tab_i_seq_0):
      tab_i_seq_0[i_seqs_1_2_3] = self.proxies.size()
      self.proxies.append(proxy)
      return True
    return False


  def process(self, source_info, proxy, tolerance=1.e-6):
    result = proxy_registry_process_result()
    proxy = proxy.sort_i_seqs()
    tab_i_seq_0 = self.table.setdefault(proxy.i_seqs[0], {})
    i_seqs_1_2_3 = (proxy.i_seqs[1], proxy.i_seqs[2], proxy.i_seqs[3])
    if (i_seqs_1_2_3 not in tab_i_seq_0):
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
    return result

class planarity_proxy_registry(proxy_registry_base):
  """
  self.table:

    self.table - similar to dihedral, chiralities, but undefined
      number of iseqs in the nested dictionary keys:
    { iseq0:
      {
        (iseq1, iseq2, iseq3, ... ) : Nproxy
      }
    }

  """

  def __init__(self, strict_conflict_handling):
    proxy_registry_base.__init__(self,
      proxies=shared_planarity_proxy(),
      strict_conflict_handling=strict_conflict_handling)

  def expand_with_ncs(self, nrgl, masters_and_rest_iselection):
    additional_proxies = []
    for i in range(len(self.proxies)):
      p = self.proxies[i]
      new_current_proxy_iseqs = [ masters_and_rest_iselection[iseq] for iseq in p.i_seqs]
      new_master_p = planarity_proxy(
            i_seqs=new_current_proxy_iseqs,
            proxy=p).sort_i_seqs()
      self.proxies[i] = new_master_p
      all_new_iseqs = nrgl.get_copy_iseqs(new_current_proxy_iseqs)
      if all_new_iseqs is not None:
        for new_iseqs in all_new_iseqs:
          new_proxy = planarity_proxy(
              i_seqs=new_iseqs,
              proxy=p).sort_i_seqs()
          # marking table
          tab_i_seq_0 = self.table.setdefault(new_proxy.i_seqs[0], {})
          tab_i_seq_0[(new_proxy.i_seqs[1:])] = self.proxies.size()
          # ~ self._append_proxy
          additional_proxies.append(new_proxy)
          # self.proxies.append(new_proxy)
          self.source_labels.append(self.source_labels[i])
          self.source_n_expected_atoms.append(self.source_n_expected_atoms[i])
    self.proxies.extend(additional_proxies)

  def add_if_not_duplicated(self, proxy, tolerance=1.e-6):
    assert proxy.i_seqs.size() > 2
    proxy = proxy.sort_i_seqs()
    tab_i_seq_0 = self.table.setdefault(proxy.i_seqs[0], {})
    i_seqs_1_up = tuple(proxy.i_seqs[1:])
    if (i_seqs_1_up not in tab_i_seq_0):
      tab_i_seq_0[i_seqs_1_up] = self.proxies.size()
      # saving proxy number in list
      self.proxies.append(proxy)
      return True
    return False

  def process(self, source_info, proxy, tolerance=1.e-6):
    assert proxy.i_seqs.size() > 0
    result = proxy_registry_process_result()
    proxy = proxy.sort_i_seqs()
    tab_i_seq_0 = self.table.setdefault(proxy.i_seqs[0], {})
    i_seqs_1_up = tuple(proxy.i_seqs[1:])
    if (i_seqs_1_up not in tab_i_seq_0):
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
    return result

class parallelity_proxy_registry(proxy_registry_base):
  """
  self.table:
  { ( (iseqs), (jseqs) ) : Nproxy }
  """

  def __init__(self, strict_conflict_handling):
    proxy_registry_base.__init__(self,
        proxies=shared_parallelity_proxy(),
        strict_conflict_handling=strict_conflict_handling)

  def expand_with_ncs(self, nrgl, masters_and_rest_iselection):
    additional_proxies = []
    for i in range(len(self.proxies)):
      p = self.proxies[i]
      new_current_proxy_iseqs = [ masters_and_rest_iselection[iseq] for iseq in p.i_seqs]
      new_current_proxy_jseqs = [ masters_and_rest_iselection[jseq] for jseq in p.j_seqs]
      new_master_p = parallelity_proxy(
            i_seqs=new_current_proxy_iseqs,
            j_seqs=new_current_proxy_jseqs,
            proxy=p).sort_i_seqs()
      self.proxies[i] = new_master_p
      # all_new_iseqs = nrgl.get_copy_iseqs(new_current_proxy_iseqs)
      all_new_iseqs = nrgl.get_copy_iseqs(new_current_proxy_iseqs)
      all_new_jseqs = nrgl.get_copy_iseqs(new_current_proxy_jseqs)
      if all_new_iseqs is not None and all_new_jseqs is not None:
        for new_iseqs, new_jseqs in zip(all_new_iseqs, all_new_jseqs):
          new_proxy = parallelity_proxy(
              i_seqs=new_iseqs,
              j_seqs=new_jseqs,
              proxy=p).sort_i_seqs()
          # marking table
          self.table[(tuple(new_proxy.i_seqs), tuple(new_proxy.j_seqs))] = self.proxies.size()
          # ~ self._append_proxy
          additional_proxies.append(new_proxy)
          # self.proxies.append(new_proxy)
          self.source_labels.append(self.source_labels[i])
          self.source_n_expected_atoms.append(self.source_n_expected_atoms[i])
    self.proxies.extend(additional_proxies)

  def add_if_not_duplicated(self, proxy, tolerance=1.e-6):
    assert proxy.i_seqs.size() > 2
    assert proxy.j_seqs.size() > 2
    proxy = proxy.sort_ij_seqs()
    tab_i_seqs = tuple(proxy.i_seqs)
    tab_j_seqs = tuple(proxy.j_seqs)
    if (self.table.get((tab_i_seqs, tab_j_seqs), -1) <0 and
        self.table.get((tab_j_seqs, tab_i_seqs), -1) <0):
      # saving proxy number in list
      self.table[(tab_i_seqs, tab_j_seqs)]= self.proxies.size()
      self.proxies.append(proxy)
      return True
    return False

  def process(self, source_info, proxy, tolerance=1.e-6):
    assert proxy.i_seqs.size() > 2
    assert proxy.j_seqs.size() > 2
    result = proxy_registry_process_result()
    proxy = proxy.sort_ij_seqs()
    # here we want to make sure that we don't have this proxy yet
    tab_i_seqs = tuple(proxy.i_seqs)
    tab_j_seqs = tuple(proxy.j_seqs)
    if (self.table.get((tab_i_seqs, tab_j_seqs), -1) <0 and
        self.table.get((tab_j_seqs, tab_i_seqs), -1) <0):
      # saving proxy number in list
      self.table[(tab_i_seqs, tab_j_seqs)]= self.proxies.size()
      self._append_proxy(
        source_info=source_info,
        proxy=proxy,
        process_result=result)
    else:
      # conflict, we have precisely the same proxy!!!
      i_list = max(self.table.get((tab_i_seqs, tab_j_seqs), -1),
                   self.table.get((tab_j_seqs, tab_i_seqs), -1))
      # number of duplicated proxy
      result.tabulated_proxy = self.proxies[i_list]
      if (not approx_equal(result.tabulated_proxy.weight, proxy.weight,
                           eps=tolerance)):
        self._handle_conflict(
          source_info=source_info,
          proxy=proxy,
          i_list=i_list,
          process_result=result)
    return result

@bp.inject_into(prolsq_repulsion_function)
class _():

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
                           max_items=None,
                           origin_id=None,
                           return_result = False):
  if return_result:
      result = group_args(group_args_type = 'Bond restraints',
        by_value = by_value,
        max_items = max_items,
        value_list = [],
      )

  if unit_cell is None:
    sorted_table, n_not_shown = self.get_sorted(
      by_value=by_value,
      sites_cart=sites_cart,
      site_labels=site_labels,
      max_items=max_items,
      origin_id=origin_id)
  else:
    sorted_table, n_not_shown = self.get_sorted(
      by_value=by_value,
      sites_cart=sites_cart,
      unit_cell=unit_cell,
      site_labels=site_labels,
      max_items=max_items,
      origin_id=origin_id)
  len_sorted_table = 0 if sorted_table is None else len(sorted_table)
  if n_not_shown is None: n_not_shown = 0
  print("%sBond restraints: %d" % (prefix, len_sorted_table+n_not_shown), file=f)
  if (f is None): f = sys.stdout
  if len_sorted_table+n_not_shown==0: return
  print("%sSorted by %s:" % (prefix, by_value), file=f)
  if sorted_table is not None:
    for restraint_info in sorted_table :
      (i_seq,j_seq,  # Always
        labels, distance_ideal, distance_model, slack, delta, sigma, weight,
        residual, sym_op_j, rt_mx) = restraint_info
      s = "bond"
      for label in labels :
        print("%s%4s %s" % (prefix, s, label), file=f)
        s = ""
      if (slack == 0):
        l = ""
        v = ""
      else:
        l = "  slack"
        v = " %6.3f" % slack
      print("%s  ideal  model%s  delta    sigma   weight residual%s" % (
        prefix, l, sym_op_j), file=f)
      print("%s  %5.3f %6.3f%s %6.3f %6.2e %6.2e %6.2e" % (
        prefix, distance_ideal, distance_model, v, delta,
        sigma, weight, residual), end='', file=f)
      if (rt_mx is not None):
        print(" " + str(rt_mx), end='', file=f)
      print(file=f)
      if return_result:
          value = group_args(
            group_args_type =
             'Bond distance:  ',
            labels = labels,
            delta = delta,
            sigma = sigma,
            residual = residual,
            ideal = distance_ideal,
            model = distance_model,)
          result.value_list.append(value)


  if (n_not_shown != 0):
    print(prefix + "... (remaining %d not shown)" % n_not_shown, file=f)

  if return_result:
    return result

@bp.inject_into(shared_bond_asu_proxy)
class _():

  def get_proxies_without_origin_id(self, origin_id):
    result = shared_bond_asu_proxy()
    for p in self:
      if p.origin_id != origin_id:
        result.append(p)
    return result

def resid_to_pymol(resid):
  resid = resid.strip()
  if resid.startswith('-'):
    resid = '\\' + resid
  return resid

@bp.inject_into(shared_bond_simple_proxy)
class _():

  def _generate_proxy_and_atom_labels(self, pdb_hierarchy):
    pdb_atoms = pdb_hierarchy.atoms()
    for proxy in self:
      i_seq, j_seq = proxy.i_seqs
      atom1 = pdb_atoms[i_seq].fetch_labels()
      atom2 = pdb_atoms[j_seq].fetch_labels()
      yield proxy, atom1, atom2

  def as_csv(self, pdb_hierarchy):
    # chain1,resname1,resid1,altloc1,name1,chain2,resname2,resid2,altloc2,name2
    # A,ALA,1," ",O,A,GLY,5," ",N
    result='chain1,resname1,resid1,name1,chain2,resname2,resid2,name2\n'
    for proxy, atom1, atom2 in self._generate_proxy_and_atom_labels(pdb_hierarchy):
      line = [ atom1.chain_id, atom1.resseq, atom1.altloc, atom1.name,
               atom2.chain_id, atom2.resseq, atom2.altloc, atom2.name,
               ]
      result += '%s\n' % ",".join(line)
    return result

  def as_pymol_dashes(self, pdb_hierarchy):
    # copied from mmtbx/geometry_restraints/hbond.py due to deprecation of
    # hbond.py
    result = ""
    for proxy, atom1, atom2 in self._generate_proxy_and_atom_labels(pdb_hierarchy):
      base_sele = """chain "%s" and resi %s and name %s and alt '%s'"""
      sele1 = base_sele % (atom1.chain_id.strip(), resid_to_pymol(atom1.resid()), atom1.name, atom1.altloc)
      sele2 = base_sele % (atom2.chain_id.strip(), resid_to_pymol(atom2.resid()), atom2.name, atom2.altloc)
      result += "dist %s, %s\n" % (sele1, sele2)
    return result

  def as_refmac_restraints(self, pdb_hierarchy):
    # copied from mmtbx/geometry_restraints/hbond.py due to deprecation of
    # hbond.py
    result = ""
    for proxy, atom1, atom2 in self._generate_proxy_and_atom_labels(pdb_hierarchy):
      sigma = 0.05
      if proxy.weight > 1e-5:
        sigma = 1.0/(proxy.weight**0.5)
      cmd = (("exte dist first chain %s residue %s atom %s " +
              "second chain %s residue %s atom %s value %.3f sigma %.2f") %
        (atom1.chain_id, atom1.resseq, atom1.name, atom2.chain_id,
         atom2.resseq, atom2.name, proxy.distance_ideal, sigma))
      result += "%s\n" % cmd
    return result

  def as_kinemage(self, pdb_hierarchy):
    # copied from mmtbx/geometry_restraints/hbond.py due to deprecation of
    # hbond.py. Outputs something, not tested.
    pdb_atoms = pdb_hierarchy.atoms()
    result = ""
    result += """\
@group {PHENIX H-bonds}
@subgroup {H-bond dots} dominant"""
    for proxy in self:
      i_seq, j_seq = proxy.i_seqs
      a = pdb_atoms[i_seq].xyz
      b = pdb_atoms[j_seq].xyz
      ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
      result += """@dotlist {Drawn dots} color= green"""
      for x in range(1, 12):
        fac = float(x) / 12
        vec = (a[0] + (ab[0]*fac), a[1] + (ab[1]*fac), a[2] + (ab[2]*fac))
        if (x == 1):
          result += "{drawn} %.4f %.4f %.4f" % vec
        else :
          result += "{''} %.4f %.4f %.4f" % vec
    return result

  # shared_bond_simple_proxy
  def get_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        unit_cell=None,
        max_items=None,
        origin_id=None):
    assert origin_id is None # not implemented
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
      info = (i_seq, j_seq, labels, restraint.distance_ideal, restraint.distance_model,
         restraint.slack, restraint.delta,
         weight_as_sigma(weight=restraint.weight), restraint.weight,
         restraint.residual(), sym_op_j, rt_mx)
      sorted_table.append(info)
    n_not_shown = data_to_sort.size() - i_proxies_sorted.size()
    return sorted_table, n_not_shown

  # shared_bond_simple_proxy
  def show_sorted(self,
                  by_value,
                  sites_cart,
                  site_labels=None,
                  unit_cell=None,
                  f=None,
                  prefix="",
                  max_items=None,
                  origin_id=None,
                  return_result = False):
    if f is None: f = sys.stdout
    # print >> f, "%sBond restraints: %d" % (prefix, self.size())
    return _bond_show_sorted_impl(self, by_value,
                           sites_cart=sites_cart,
                           site_labels=site_labels,
                           unit_cell=unit_cell,
                           f=f,
                           prefix=prefix,
                           max_items=max_items,
                           origin_id=origin_id,
                           return_result = return_result)

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

  def get_proxies_with_origin_id(self, origin_id=0): # this is the default
    result = []
    for p in self:
      if p.origin_id == origin_id:
        result.append(p)
    return result

  def get_proxies_without_origin_id(self, origin_id):
    result = shared_bond_simple_proxy()
    for p in self:
      if p.origin_id != origin_id:
        result.append(p)
    return result

@bp.inject_into(bond_sorted_asu_proxies)
class _():

  def show_histogram_of_model_distances(self,
        sites_cart,
        n_slots=5,
        cutoff_warn_small=0.5,
        cutoff_warn_large=5,
        cutoff_warn_extreme=20,
        f=None,
        prefix="",
        origin_id=None):
    if (self.n_total() == 0): return None
    if (f is None): f = sys.stdout
    print("%sHistogram of bond lengths:" % prefix, file=f)
    hdata = None
    if origin_id is None:
      hdata = bond_distances_model(
        sites_cart=sites_cart,
        sorted_asu_proxies=self)
    else:
      selected_simple = self.simple.proxy_select(origin_id = origin_id)
      selected_asu = self.asu.proxy_select(origin_id = origin_id)
      hdata_simple = bond_distances_model(
        sites_cart=sites_cart,
        proxies=selected_simple)
      sap = bond_sorted_asu_proxies(asu_mappings=self.asu_mappings())
      sap.process(selected_asu)
      hdata_asu = bond_distances_model(
        sites_cart=sites_cart,
        sorted_asu_proxies=sap)
      hdata_simple.extend(hdata_asu)
      hdata = hdata_simple
    histogram = flex.histogram(
      data=flex.double(hdata),
      n_slots=n_slots)
    low_cutoff = histogram.data_min()
    for i,n in enumerate(histogram.slots()):
      high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
      print("%s  %8.2f - %8.2f: %d" % (
        prefix, low_cutoff, high_cutoff, n), file=f)
      low_cutoff = high_cutoff
    if (cutoff_warn_small is not None
        and histogram.data_min() < cutoff_warn_small):
      print("%sWarning: very small bond lengths." % prefix, file=f)
    if (cutoff_warn_extreme is not None
        and histogram.data_max() > cutoff_warn_extreme):
      print("%sWarning: extremely large bond lengths." % prefix, file=f)
    elif (cutoff_warn_large is not None
          and histogram.data_max() > cutoff_warn_large):
      print("%sWarning: very large bond lengths." % prefix, file=f)
    return histogram

  def deltas(self, sites_cart, origin_id=None):
    if origin_id is None:
      return bond_deltas(sites_cart=sites_cart, sorted_asu_proxies=self)
    return bond_deltas(sites_cart=sites_cart, sorted_asu_proxies=self,
        origin_id=origin_id)

  def residuals(self, sites_cart, origin_id=None):
    if origin_id is None:
      return bond_residuals(sites_cart=sites_cart, sorted_asu_proxies=self)
    return bond_residuals(sites_cart=sites_cart, sorted_asu_proxies=self,
        origin_id=origin_id)

  def get_proxies_with_origin_id(self, origin_id=0): # this is the default
    result = []
    for p in self.simple:
      if p.origin_id == origin_id:
        result.append(p)
    for p in self.asu:
      if p.origin_id == origin_id:
        result.append(p)
    return result

  def show_histogram_of_deltas(self,
        sites_cart,
        n_slots=5,
        f=None,
        prefix="",
        origin_id=None):
    if (self.n_total() == 0): return
    if (f is None): f = sys.stdout
    print("%sHistogram of bond deltas:" % prefix, file=f)
    hdata = None
    if origin_id is None:
      hdata = bond_deltas(
        sites_cart=sites_cart,
        sorted_asu_proxies=self)
    else:
      sorted_table, n_not_shown = self.get_sorted(
                        by_value="delta",
                        sites_cart=sites_cart,
                        origin_id=origin_id)
      hd = [x[6] for x in sorted_table]
      hdata = flex.double(hd)
    histogram = flex.histogram(
      data=flex.abs(hdata),
      n_slots=n_slots)
    low_cutoff = histogram.data_min()
    for i,n in enumerate(histogram.slots()):
      high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
      print("%s  %8.3f - %8.3f: %d" % (
        prefix, low_cutoff, high_cutoff, n), file=f)
      low_cutoff = high_cutoff
    return histogram

  # bond_sorted_asu_proxies
  def get_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        max_items=None,
        origin_id=None):
    assert by_value in ["residual", "delta"]
    assert site_labels is None or len(site_labels) == sites_cart.size()
    if (self.n_total() == 0): return None, None
    if (max_items is not None and max_items <0): return None, None
    if (by_value == "residual"):
      data_to_sort = self.residuals(sites_cart=sites_cart)
    elif (by_value == "delta"):
      data_to_sort = flex.abs(self.deltas(sites_cart=sites_cart))
    else:
      raise AssertionError
    i_proxies_sorted = flex.sort_permutation(data=data_to_sort, reverse=True)
    total_proxies = len(i_proxies_sorted)
    correct_id_proxies = total_proxies
    if origin_id is not None:
      correct_id_proxies = [i.origin_id for i in self.simple].count(origin_id)
      correct_id_proxies += [i.origin_id for i in self.asu].count(origin_id)
    if max_items is None:
      max_items = correct_id_proxies
    if (self.asu.size() == 0):
      asu_mappings = None
    else:
      asu_mappings = self.asu_mappings()
    n_simple = self.simple.size()
    sorted_table = []
    n = 0
    n_outputted = 0
    n_excluded = 0
    while n < total_proxies and n_outputted < max_items:
      i_proxy = i_proxies_sorted[n]
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
        i_seq,j_seq = proxy.i_seq, proxy.j_seq
        rt_mx = asu_mappings.get_rt_mx_ji(pair=proxy)
        sym_op_j = " sym.op."
        restraint = bond(
          sites_cart=sites_cart,
          asu_mappings=asu_mappings,
          proxy=proxy)
      if origin_id is None or origin_id == proxy.origin_id:
        labels = []
        for i in [i_seq, j_seq]:
          if (site_labels is None): l = str(i)
          else:                     l = site_labels[i]
          labels.append(l)
        info=(i_seq, j_seq, labels, restraint.distance_ideal, restraint.distance_model,
           restraint.slack, restraint.delta,
           weight_as_sigma(weight=restraint.weight), restraint.weight,
           restraint.residual(), sym_op_j, rt_mx)
        sorted_table.append(info)
        n_outputted += 1
      else:
        n_excluded += 1
      n += 1
    n_not_shown = correct_id_proxies - n_outputted
    return sorted_table, n_not_shown

  def get_outliers(self, sites_cart, sigma_threshold, origin_id=None):
    result = []
    from cctbx.geometry_restraints.linking_class import linking_class
    origin_ids = linking_class()
    if origin_id is None: origin_id=origin_ids.get_origin_id('covalent geometry')
    vals = self.get_sorted(
        by_value="delta",
        sites_cart=sites_cart,
        origin_id=origin_id)[0]
    if(vals is None): return result
    for it in vals:
      i,j = it[0],it[1]
      delta = abs(it[6])
      sigma = it[7]
      if(delta > sigma*sigma_threshold):
        result.append([i,j])
    return result

  def get_filtered_deltas(self,
      sites_cart,
      origin_id=None):
    n_proxies = self.n_total()
    if (n_proxies == 0): return None
    if (self.asu.size() == 0):
      asu_mappings = None
    else:
      asu_mappings = self.asu_mappings()
    if origin_id is None:
      return self.deltas(sites_cart=sites_cart)
    else:
      result = flex.double()
      n_simple = self.simple.size()
      for i in range(n_simple):
        proxy = self.simple[i]
        if proxy.origin_id == origin_id:
          rt_mx = None
          sym_op_j = ""
          restraint = bond(
            sites_cart=sites_cart,
            proxy=proxy)
          result.append(restraint.delta)
      for i in range(n_simple, n_proxies):
        proxy = self.asu[i-n_simple]
        if proxy.origin_id == origin_id:
          rt_mx = asu_mappings.get_rt_mx_ji(pair=proxy)
          sym_op_j = " sym.op."
          restraint = bond(
            sites_cart=sites_cart,
            asu_mappings=asu_mappings,
            proxy=proxy)
          result.append(restraint.delta)
      return result if len(result) > 0 else None

  # bond_sorted_asu_proxies
  def show_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        f=None,
        prefix="",
        max_items=None,
        origin_id=None,
        return_result = False,):
    if f is None: f = sys.stdout
    # print >> f, "%sBond restraints: %d" % (prefix, self.n_total())
    return _bond_show_sorted_impl(self, by_value,
                          sites_cart=sites_cart,
                          site_labels=site_labels,
                          f=f,
                          prefix=prefix,
                          max_items=max_items,
                          origin_id=origin_id,
                          return_result = return_result)

@bp.inject_into(nonbonded_sorted_asu_proxies)
class _():

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
    print("%sHistogram of nonbonded interaction distances:" % prefix, file=f)
    histogram = flex.histogram(
      data=self.deltas(sites_cart=sites_cart),
      n_slots=n_slots)
    low_cutoff = histogram.data_min()
    for i,n in enumerate(histogram.slots()):
      high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
      print("%s  %8.2f - %8.2f: %d" % (
        prefix, low_cutoff, high_cutoff, n), file=f)
      low_cutoff = high_cutoff
    if (cutoff_warn_small is not None
        and histogram.data_min() < cutoff_warn_small):
      print("%sWarning: very small nonbonded interaction distances." % (
        prefix), file=f)
    return histogram

  def get_sorted_i_proxies(self,
                           by_value,
                           sites_cart,
                           site_labels=None,
                           max_items=None,
                           ):
    assert by_value in ["delta"]
    deltas = nonbonded_deltas(sites_cart=sites_cart, sorted_asu_proxies=self)
    if (deltas.size() == 0): return
    if (max_items is not None and max_items <= 0): return
    i_proxies_sorted = flex.sort_permutation(data=deltas)
    if (max_items is not None):
      i_proxies_sorted = i_proxies_sorted[:max_items]
    return i_proxies_sorted

  def sorted_value_proxies_generator(self,
                           by_value,
                           sites_cart,
                           cutoff = 100):
    assert by_value in ["delta"]
    deltas = nonbonded_deltas(sites_cart=sites_cart, sorted_asu_proxies=self)
    if (deltas.size() == 0): return
    i_proxies_sorted = flex.sort_permutation(data=deltas)
    n_proxies = deltas.size()
    n_simple = self.simple.size()
    asu_mappings = self.asu_mappings()
    i = 0
    while i < n_proxies and deltas[i_proxies_sorted[i]] < cutoff:
      i_proxy = i_proxies_sorted[i]
      if (i_proxy < n_simple):
        i_seq, j_seq = self.simple[i_proxy].i_seqs
        yield (
            i_seq,
            j_seq,
            deltas[i_proxy],
            None,
            "",
            self.simple[i_proxy],
            )
      else:
        i_seq, j_seq = self.asu[i_proxy-n_simple].i_seq, self.asu[i_proxy-n_simple].j_seq
        yield (
            i_seq,
            j_seq,
            deltas[i_proxy],
            asu_mappings.get_rt_mx_ji(pair=self.asu[i_proxy-n_simple]),
            " sym.op.",
            self.asu[i_proxy-n_simple])
      i += 1

  def get_sorted_proxies(self,
                         by_value,
                         sites_cart,
                         site_labels=None,
                         max_items=None,
                         ):
    assert by_value in ["delta"]
    i_proxies_sorted = self.get_sorted_i_proxies(by_value=by_value,
                                                 sites_cart=sites_cart,
                                                 site_labels=site_labels,
                                                 max_items=max_items,
      )
    if (self.asu.size() == 0):
      asu_mappings = None
    else:
      asu_mappings = self.asu_mappings()
    n_simple = self.simple.size()
    sorted_proxies=[]
    for i_proxy in i_proxies_sorted:
      if (i_proxy < n_simple):
        proxy = self.simple[i_proxy]
      else:
        proxy = self.asu[i_proxy-n_simple]
      sorted_proxies.append(proxy)
    return sorted_proxies

  def get_sorted(self,
                 by_value,
                 sites_cart,
                 site_labels=None,
                 max_items=None,
                 include_proxy=False,
                 ):
    assert by_value in ["delta"]
    deltas = nonbonded_deltas(sites_cart=sites_cart, sorted_asu_proxies=self)
    if (deltas.size() == 0): return [], 0
    if (max_items is not None and max_items <= 0): return [], deltas.size()
    i_proxies_sorted = flex.sort_permutation(data=deltas)
    if (max_items is not None):
      i_proxies_sorted = i_proxies_sorted[:max_items]
    if (self.asu.size() == 0):
      asu_mappings = None
    else:
      asu_mappings = self.asu_mappings()
    n_simple = self.simple.size()
    sorted_table=[]
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
      labels = []
      for i in [i_seq, j_seq]:
        if (site_labels is None): l = str(i)
        else:                     l = site_labels[i]
        labels.append(l)
      m, v = deltas[i_proxy], proxy.vdw_distance
      if include_proxy:
        sorted_table.append(
          (labels,
           i_seq,
           j_seq,
           deltas[i_proxy],
           proxy.vdw_distance,
           sym_op_j,
           rt_mx,
           proxy)
          )
      else:
        sorted_table.append(
          (labels,
           i_seq,
           j_seq,
           deltas[i_proxy],
           proxy.vdw_distance,
           sym_op_j,
           rt_mx)
          )
    n_not_shown = deltas.size() - i_proxies_sorted.size()
    return sorted_table, n_not_shown

  def get_symmetry_interacting_indices_unique(self, sites_cart):
    """
    See output of show_sorted below to understand what this function does.
    Find atom indices in the main copy that bump with the main copy via symmetry
    or periodicity. There may be several symmetry/periodicity related copies
    corresponding to given atom in main copy. The function returns a dictionary
    with keys being atom indices in main copy and values being the list of
    symmetry operations to generate interacting mates.
    """
    result = {}
    deltas = nonbonded_deltas(sites_cart=sites_cart, sorted_asu_proxies=self)
    if (deltas.size() == 0): return
    i_proxies_sorted = flex.sort_permutation(data=deltas)
    if (self.asu.size() == 0):
      asu_mappings = None
    else:
      asu_mappings = self.asu_mappings()
    n_simple = self.simple.size()
    for i_proxy in i_proxies_sorted:
      if (i_proxy >= n_simple):
        proxy = self.asu[i_proxy-n_simple]
        i_seq,j_seq = proxy.i_seq,proxy.j_seq
        rt_mx = asu_mappings.get_rt_mx_ji(pair=proxy)
        #if str(rt_mx)=="x,y,z": continue # why do I need this?
        result.setdefault(j_seq, []).append(rt_mx)
    for k,v in six.iteritems(result):
      result[k] = list(set(v))
    return result

  def show_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        f=None,
        prefix="",
        max_items=None,
        suppress_model_minus_vdw_greater_than=0.2,
        but_show_all_model_up_to=3.5,
        return_result = False):

    if return_result:
      result = group_args(group_args_type = 'Non-bonded restraints',
        by_value = by_value,
        max_items = max_items,
        value_list = [],
      )

    assert by_value in ["delta"]
    sorted_table, n_not_shown = self.get_sorted(
        by_value=by_value,
        sites_cart=sites_cart,
        site_labels=site_labels,
        max_items=max_items,
        include_proxy=False)
    if (f is None): f = sys.stdout
    print("%sNonbonded interactions: %d" % (prefix, len(sorted_table)+n_not_shown), file=f)
    if len(sorted_table) == 0: return
    print("%sSorted by model distance:" % prefix, file=f)

    for info in sorted_table:
      labels, i_seq, j_seq, delta, vdw_distance, sym_op_j, rt_mx = info
      def suppress():
        m, v = delta, vdw_distance
        if (suppress_model_minus_vdw_greater_than is None): return False
        if (m-v <= suppress_model_minus_vdw_greater_than): return False
        if (but_show_all_model_up_to is None): return True
        if (m <= but_show_all_model_up_to): return False
        return True
      if (suppress()): continue
      s = "nonbonded"
      for l in labels:
        print("%s%9s %s" % (prefix, s, l), file=f)
        s = ""
      print("%s   model   vdw%s" % (prefix, sym_op_j), file=f)
      print("%s  %6.3f %5.3f" % (prefix, delta, vdw_distance), end='', file=f)
      if (rt_mx is not None):
        print(" " + str(rt_mx), end='', file=f)
      print(file=f)
      if return_result:
          value = group_args(
            group_args_type =
             'Non-bonded distance:  ideal is vdw_distance, '+
                 'model is delta (actual)',
            labels = labels,
            delta = None,
            sigma = None,
            residual = None,
            ideal = vdw_distance,
            model = delta,)
          result.value_list.append(value)

    if (n_not_shown != 0):
      print(prefix + "... (remaining %d not shown)" % n_not_shown, file=f)

    if return_result:
      return result

@bp.inject_into(angle)
class _():

  def _show_sorted_item(O, f, prefix):
    print("%s    ideal   model   delta" \
      "    sigma   weight residual" % prefix, file=f)
    print("%s  %7.2f %7.2f %7.2f %6.2e %6.2e %6.2e" % (
      prefix,
      O.angle_ideal, O.angle_model, O.delta,
      weight_as_sigma(weight=O.weight), O.weight, O.residual()), file=f)

  def _get_sorted_item(O):
    return [O.angle_ideal, O.angle_model, O.delta,
            weight_as_sigma(weight=O.weight), O.weight, O.residual()]

@bp.inject_into(shared_angle_proxy)
class _():

  def as_pymol_dashes(self, pdb_hierarchy):
    # copied from mmtbx/geometry_restraints/hbond.py due to deprecation of
    # hbond.py
    result = ""
    pdb_atoms = pdb_hierarchy.atoms()
    i = 0
    for proxy in self:
      i_seq, j_seq, k_seq = proxy.i_seqs
      atom1 = pdb_atoms[i_seq].fetch_labels()
      atom2 = pdb_atoms[j_seq].fetch_labels()
      atom3 = pdb_atoms[k_seq].fetch_labels()
      base_sele = """chain "%s" and resi %s and name %s and alt '%s'"""
      sele1 = base_sele % (atom1.chain_id.strip(), resid_to_pymol(atom1.resid()), atom1.name, atom1.altloc)
      sele2 = base_sele % (atom2.chain_id.strip(), resid_to_pymol(atom2.resid()), atom2.name, atom2.altloc)
      sele3 = base_sele % (atom3.chain_id.strip(), resid_to_pymol(atom3.resid()), atom3.name, atom3.altloc)
      result += "angle a%d, %s, %s, %s\n" % (i, sele1, sele2, sele3)
      i += 1
    return result

  def deltas(self, sites_cart, unit_cell=None, origin_id=None):
    if unit_cell is None:
      if origin_id is None:
        return angle_deltas(sites_cart=sites_cart, proxies=self)
      return angle_deltas(sites_cart=sites_cart, proxies=self, origin_id=origin_id)
    else:
      if origin_id is None:
        return angle_deltas(
          unit_cell=unit_cell, sites_cart=sites_cart, proxies=self)
      return angle_deltas(
        unit_cell=unit_cell, sites_cart=sites_cart, proxies=self, origin_id=origin_id)

  def residuals(self, sites_cart, unit_cell=None):
    if unit_cell is None:
      return angle_residuals(sites_cart=sites_cart, proxies=self)
    else:
      return angle_residuals(
        unit_cell=unit_cell, sites_cart=sites_cart, proxies=self)

  def get_filtered_deltas(self,
      sites_cart,
      origin_id=None):
    # ANGLE
    n_proxies = self.size()
    if (n_proxies == 0): return None
    if origin_id is None:
      if unit_cell is None:
        result = flex.abs(O.deltas(sites_cart=sites_cart))
      else:
        result = flex.abs(
          O.deltas(unit_cell=unit_cell, sites_cart=sites_cart))
    else:
      result = flex.double()
      for i in range(n_proxies):
        proxy = self[i]
        if proxy.origin_id == origin_id:
          restraint = angle(
            sites_cart=sites_cart,
            proxy=proxy)
          result.append(restraint.delta)
    return result if len(result) > 0 else None

  def show_histogram_of_deltas(self,
        sites_cart,
        unit_cell=None,
        n_slots=5,
        proxy_label="bond angle",
        f=None,
        prefix="",
        origin_id=None):
    return _show_histogram_of_deltas_impl(O=self,
      proxy_label=proxy_label,
      format_cutoffs="%8.2f",
      unit_cell=unit_cell,
      sites_cart=sites_cart, n_slots=n_slots, f=f, prefix=prefix,
      origin_id=origin_id)

  def show_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        proxy_label="Bond angle",
        unit_cell=None,
        f=None,
        prefix="",
        max_items=None,
        origin_id=None,
        return_result = False):
    return _show_sorted_impl(O=self,
        proxy_type=angle,
        proxy_label=proxy_label,
        item_label="angle",
        by_value=by_value, unit_cell=unit_cell, sites_cart=sites_cart,
        site_labels=site_labels, f=f, prefix=prefix, max_items=max_items,
        origin_id=origin_id,
        return_result = return_result)

  def get_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        unit_cell=None,
        max_items=None,
        origin_id=None):
    return _get_sorted_impl(O=self,
        proxy_type=angle,
        by_value=by_value, unit_cell=unit_cell, sites_cart=sites_cart,
        site_labels=site_labels, max_items=max_items,
        get_restraints_only=False, origin_id=origin_id)

  def get_outliers(self, sites_cart, sigma_threshold, origin_id=0):
    result = []
    vals = self.get_sorted(by_value="delta",
                           origin_id=origin_id,
                           sites_cart=sites_cart)[0]
    if(vals is None): return result
    for it in vals:
      i,j,k = [int(i) for i in it[0]]
      delta = abs(it[3])
      sigma = it[4]
      if(delta > sigma*sigma_threshold):
        result.append([i,j,k])
    return result

@bp.inject_into(dihedral)
class _():

  def _show_sorted_item(O, f, prefix):
    print("%s    ideal   model   delta" \
      " %        s    sigma   weight residual" % (
        prefix, {False: "sinusoidal", True: " harmonic "}[O.periodicity <= 0]), file=f)
    angle_ideal = O.angle_model+O.delta
    if angle_ideal<-180: angle_ideal+=360
    print("%s  %7.2f %7.2f %7.2f %5d      %6.2e %6.2e %6.2e" % (
      prefix,
      angle_ideal, O.angle_model, O.delta, O.periodicity,
      weight_as_sigma(weight=O.weight), O.weight, O.residual()), file=f)

  def _get_sorted_item(O):
    return [O.angle_ideal, O.angle_model, O.delta, O.periodicity,
            weight_as_sigma(weight=O.weight), O.weight, O.residual()]

@bp.inject_into(shared_dihedral_proxy)
class _():

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
        proxy_label="dihedral angle",
        f=None,
        prefix=""):
    return _show_histogram_of_deltas_impl(O=self,
      proxy_label=proxy_label,
      format_cutoffs="%8.2f",
      unit_cell=unit_cell,
      sites_cart=sites_cart, n_slots=n_slots, f=f, prefix=prefix)

  def show_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        proxy_label="Dihedral angle",
        unit_cell=None,
        f=None,
        prefix="",
        max_items=None,
        origin_id=None,
        return_result = False):
    return _show_sorted_impl(O=self,
        proxy_type=dihedral,
        proxy_label=proxy_label,
        item_label="dihedral",
        by_value=by_value, unit_cell=unit_cell, sites_cart=sites_cart,
        site_labels=site_labels, f=f, prefix=prefix, max_items=max_items,
        origin_id=origin_id,
        return_result = return_result)

  def get_sorted(self,
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

  def get_outliers(self, sites_cart, sigma_threshold):
    result = []
    vals = self.get_sorted(by_value="delta", sites_cart=sites_cart)[0]
    if(vals is None): return result
    for it in vals:
      ind = [int(i) for i in it[0]]
      delta = abs(it[3])
      sigma = it[5]
      if(delta > sigma*sigma_threshold):
        result.append(ind)
    return result

@bp.inject_into(chirality)
class _():

  def _show_sorted_item(O, f, prefix):
    print("%s  both_signs  ideal   model" \
      "   delta    sigma   weight residual" % prefix, file=f)
    print("%s    %-5s   %7.2f %7.2f %7.2f %6.2e %6.2e %6.2e" % (
      prefix,
      str(O.both_signs), O.volume_ideal, O.volume_model, O.delta,
      weight_as_sigma(weight=O.weight), O.weight, O.residual()), file=f)

  def _get_sorted_item(O):
    return [str(O.both_signs), O.volume_ideal, O.volume_model, O.delta,
      weight_as_sigma(weight=O.weight), O.weight, O.residual()]

@bp.inject_into(shared_chirality_proxy)
class _():

  def deltas(self, sites_cart, unit_cell=None):
    if unit_cell is None:
      return chirality_deltas(sites_cart=sites_cart, proxies=self)
    else:
      return chirality_deltas(unit_cell=unit_cell, sites_cart=sites_cart, proxies=self)

  def residuals(self, sites_cart, unit_cell=None):
    if unit_cell is None:
      return chirality_residuals(sites_cart=sites_cart, proxies=self)
    else:
      return chirality_residuals(unit_cell=unit_cell, sites_cart=sites_cart, proxies=self)

  def show_histogram_of_deltas(self,
        sites_cart,
        n_slots=5,
        proxy_label="chiral volume",
        f=None,
        prefix=""):
    return _show_histogram_of_deltas_impl(O=self,
      proxy_label=proxy_label,
      format_cutoffs="%8.3f",
      unit_cell=None,
      sites_cart=sites_cart, n_slots=n_slots, f=f, prefix=prefix)

  def show_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        unit_cell=None,
        proxy_label="Chirality",
        f=None,
        prefix="",
        max_items=None,
        origin_id=None,
        return_result = False,
    ):

    return _show_sorted_impl(O=self,
        proxy_type=chirality,
        proxy_label=proxy_label,
        item_label="chirality",
        by_value=by_value, unit_cell=unit_cell, sites_cart=sites_cart,
        site_labels=site_labels, f=f, prefix=prefix, max_items=max_items,
        origin_id=origin_id,
        return_result = return_result)

  def get_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        max_items=None):
    return _get_sorted_impl(O=self,
        proxy_type=chirality,
        by_value=by_value, unit_cell=None, sites_cart=sites_cart,
        site_labels=site_labels, max_items=max_items,
        get_restraints_only=False)

  def get_outliers(self, sites_cart, sigma_threshold):
    result = []
    vals = self.get_sorted(by_value="delta", sites_cart=sites_cart)[0]
    if(vals is None): return result
    for it in vals:
      # it: [[iseqs], both_signs(bool), ideal, model, delta, sigma, weight, residual]
      ind = [int(i) for i in it[0]]
      delta = abs(it[4])
      sigma = it[5]
      if(delta > sigma*sigma_threshold):
        result.append(ind)
    return result

@bp.inject_into(shared_planarity_proxy)
class _():

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

  def get_sorted(O,
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
      labels = []
      for i_seq in proxy.i_seqs:
        if (site_labels is None): l = str(i_seq)
        else:                     l = site_labels[i_seq]
        labels.append(l)
      if unit_cell is None:
        restraint = planarity(sites_cart=sites_cart, proxy=proxy)
      else:
        restraint = planarity(unit_cell=unit_cell, sites_cart=sites_cart,
                              proxy=proxy)
      restraint_atoms = []
      for i,(i_seq,weight,delta,l) in enumerate(zip(proxy.i_seqs, proxy.weights,
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
        proxy_label=None, # not used yet
        unit_cell=None,
        f=None,
        prefix="",
        max_items=None,
        origin_id=None,
        return_result = False,
    ):
    if return_result:
      result = group_args(group_args_type = 'Planarity restraints',
        by_value = by_value,
        max_items = max_items,
        value_list = [],
      )



    assert by_value in ["residual", "rms_deltas"]
    assert site_labels is None or len(site_labels) == sites_cart.size()
    if (f is None): f = sys.stdout
    outl = ''
    if (O.size() == 0):
      print("%sPlanarity restraints: %d" % (prefix, O.size()), file=f)
      return
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
    n_planes=0
    for i_proxy in i_proxies_sorted:
      proxy = O[i_proxy]
      if origin_id is not None and proxy.origin_id!=origin_id: continue
      n_planes+=1
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
      outl +=  "%s      %s    delta    sigma   weight rms_deltas residual%s\n" % (
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
        outl += "%s%5s %s  %7.3f %6.2e %6.2e%s%s\n" % (
          prefix, s, l+" "*(len_max-len(l)),
          delta, weight_as_sigma(weight=weight), weight, rdr, sym_op)
        rdr = ""
        rdr_spacer = " "*20
        s = ""
        if return_result:
          sigma = weight_as_sigma(weight=weight)
          value = group_args(
            group_args_type =
             'Planarity restraint (energy for one atom, all atoms: %s)' %(
              str(ls)),
            labels = [l],
            delta = delta,
            sigma = sigma,
            residual = (delta/max(1.e-10,sigma))**2,
            ideal = None,
            model = None,)
          result.value_list.append(value)

    n_not_shown = O.size() - i_proxies_sorted.size()
    if (n_not_shown != 0):
      outl += prefix + "... (remaining %d not shown)\n" % n_not_shown
    if origin_id is None:
      n_planes = O.size()
    print("%sPlanarity restraints: %d" % (prefix, n_planes), file=f)
    if outl:
      print("%sSorted by %s:" % (prefix, by_value), file=f)
      print(outl[:-1], file=f)

    if return_result:
      return result

@bp.inject_into(parallelity)
class _():

  def _show_sorted_item(O, f, prefix):
    assert 0
    print("%s    ideal   model   delta" \
      "    sigma   weight residual" % prefix, file=f)

  def _get_sorted_item(self):
    return [self.residual(), self.weight, self.delta]

@bp.inject_into(shared_parallelity_proxy)
class _():

  def deltas(O, sites_cart, unit_cell=None):
    if unit_cell is None:
      return parallelity_deltas(sites_cart=sites_cart, proxies=O)
    else:
      return parallelity_deltas(
        unit_cell=unit_cell, sites_cart=sites_cart, proxies=O)

  def residuals(O, sites_cart, unit_cell=None):
    if unit_cell is None:
      return parallelity_residuals(sites_cart=sites_cart, proxies=O)
    else:
      return parallelity_residuals(
        unit_cell=unit_cell, sites_cart=sites_cart, proxies=O)

  def show_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        proxy_label=None, # not used yet
        f=None,
        prefix="",
        max_items=None,
        origin_id=None):
    if (f is None): f = sys.stdout
    sorted_table, n_not_shown = self.get_sorted(
          by_value=by_value,
          sites_cart=sites_cart,
          site_labels=site_labels,
          max_items=max_items,
          origin_id=origin_id)
    if sorted_table is None: return
    print("Parallelity restraints: %d" % (len(sorted_table)), file=f)
    if (self.size() == 0): return
    if len(sorted_table)==0: return
    if (max_items is not None and max_items <= 0): return
    print("%sSorted by %s:" % (prefix, by_value), file=f)
    for info in sorted_table:
      residual = info[1]
      delta_deg = delta = info[3]
      # delta_deg = math.degrees(math.acos(1-delta))
      weight = info[2]
      sigma = math.sqrt(1./weight)
      print("    plane 1                plane 2  "+\
          "              residual  delta(deg) sigma", file=f)
      r_info = "  %.2e %8.4f  %8.4f" % (residual, delta_deg, sigma)
      i_labels = info[0][0]
      j_labels = info[0][1]
      i = 0
      long_len = max(len(i_labels), len(j_labels))
      while i < long_len:
        print("    %s  %s%s" % \
            (i_labels[i] if i < len(i_labels) else " "*21,
             j_labels[i] if i < len(j_labels) else "",
             r_info), file=f)
        r_info = ""
        i += 1
      print(file=f)
    if (n_not_shown != 0):
      print(prefix + "... (remaining %d not shown)" % n_not_shown, file=f)

  def get_sorted(self,
        by_value,
        sites_cart,
        site_labels=None,
        max_items=None,
        origin_id=None):
    return _get_sorted_impl(O=self,
        proxy_type=parallelity,
        by_value=by_value, unit_cell=None, sites_cart=sites_cart,
        site_labels=site_labels, max_items=max_items,
        get_restraints_only=False,
        origin_id=origin_id)

@bp.inject_into(shared_bond_similarity_proxy)
class _():

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
      max_items=None,
      origin_id=None):
    assert by_value in ["residual", "rms_deltas"]
    assert site_labels is None or len(site_labels) == sites_cart.size()
    if (f is None): f = sys.stdout
    print("%sBond similarity restraints: %d" % (prefix, O.size()), file=f)
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
    print("%sSorted by %s:" % (prefix, by_value), file=f)
    for i_proxy in i_proxies_sorted:
      proxy = O[i_proxy]
      if origin_id is not None and proxy.origin_id!=origin_id: continue
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
      print("%s     %s    delta    sigma   weight rms_deltas residual%s" % (
          prefix, " "*len_max, sym_op_label), file=f)
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
        print("%s%4s %s  %7.3f %6.2e %6.2e%s%s" % (
          prefix, s, l+" "*(len_max-len(l)),
          delta, weight_as_sigma(weight=weight), weight, rdr, sym_op), file=f)
        rdr = ""
        rdr_spacer = " "*20
        s = ""
    n_not_shown = O.size() - i_proxies_sorted.size()
    if (n_not_shown != 0):
      print(prefix + "... (remaining %d not shown)" % n_not_shown, file=f)

def _show_histogram_of_deltas_impl(O,
        proxy_label,
        format_cutoffs,
        unit_cell,
        sites_cart,
        n_slots,
        f,
        prefix,
        origin_id=None):
    if (O.size() == 0): return
    if (f is None): f = sys.stdout
    print("%sHistogram of %s deviations from ideal:" % (
      prefix, proxy_label), file=f)
    selected_proxies = O
    if origin_id is not None:
      selected_proxies = O.proxy_select(origin_id=origin_id)
    if unit_cell is None:
      data = flex.abs(selected_proxies.deltas(sites_cart=sites_cart))
    else:
      data = flex.abs(
        selected_proxies.deltas(unit_cell=unit_cell, sites_cart=sites_cart))
    histogram = flex.histogram(
      data=data,
      n_slots=n_slots)
    fmt = "%%s  %s - %s: %%d" % (format_cutoffs, format_cutoffs)
    low_cutoff = histogram.data_min()
    for i,n in enumerate(histogram.slots()):
      high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
      print(fmt % (prefix, low_cutoff, high_cutoff, n), file=f)
      low_cutoff = high_cutoff
    return histogram

def _get_sorted_impl(O,
        proxy_type,
        by_value,
        unit_cell,
        sites_cart,
        site_labels,
        max_items,
        get_restraints_only=True,
        origin_id=None):
  assert by_value in ["residual", "delta"]
  assert site_labels is None or len(site_labels) == sites_cart.size()
  if (O.size() == 0): return None, None
  if (max_items is not None and max_items < 0): return None, None
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
  n_total_proxies = len(i_proxies_sorted)
  if max_items is None:
    max_items = n_total_proxies
  sorted_table = []
  n_added_proxies = 0
  correct_id_proxies = n_total_proxies
  if origin_id is not None:
    correct_id_proxies = [i.origin_id for i in O].count(origin_id)

  i = 0
  while i < n_total_proxies and n_added_proxies < max_items:
    proxy = O[i_proxies_sorted[i]]
    if (origin_id is None or
       (origin_id is not None and hasattr(proxy, "origin_id")
        and proxy.origin_id==origin_id)):
      labels = []
      labels_j = []
      for n, i_seq in enumerate(proxy.i_seqs):
        if (site_labels is None): l = str(i_seq)
        else:                     l = site_labels[i_seq]
        if unit_cell and proxy.sym_ops:
          sym_op = proxy.sym_ops[n]
          if not sym_op.is_unit_mx():
            l += "  %s" %sym_op.as_xyz()
        labels.append(l)
      if proxy_type==parallelity:
        for n, i_seq in enumerate(proxy.j_seqs):
          if (site_labels is None): l = str(i_seq)
          else:                     l = site_labels[i_seq]
          if unit_cell and proxy.sym_ops:
            sym_op = proxy.sym_ops[n]
            if not sym_op.is_unit_mx():
              l += "  %s" %sym_op.as_xyz()
          labels_j.append(l)
        labels = [labels,labels_j]
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
      n_added_proxies += 1
    i += 1
  n_not_shown = correct_id_proxies - n_added_proxies
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
        max_items,
        origin_id=None,
        return_result = False,
      ):
  if return_result:
    result = group_args(group_args_type = '%s restraints' %(proxy_label),
        by_value = by_value,
        max_items = max_items,
        value_list = [],)

  if (f is None): f = sys.stdout
  sorted_table, n_not_shown = _get_sorted_impl(O,
        proxy_type=proxy_type,
        by_value=by_value,
        unit_cell=unit_cell,
        sites_cart=sites_cart,
        site_labels=site_labels,
        max_items=max_items,
        get_restraints_only=True,
        origin_id=origin_id)
  len_sorted_table = 0 if sorted_table is None else len(sorted_table)
  if n_not_shown is None: n_not_shown = 0
  print("%s%s restraints: %d" % (prefix, proxy_label, len_sorted_table+n_not_shown), file=f)
  if (O.size() == 0): return
  if len_sorted_table+n_not_shown==0: return
  n_harmonic = 0
  n_sinusoidal = 0
  if (max_items is not None and max_items <= 0): return
  item_label_blank = " " * len(item_label)
  outl = six.StringIO()
  for (labels, restraint) in sorted_table :
    if (proxy_type is dihedral) and (origin_id==0 or origin_id is None):
      if restraint.periodicity<=0: n_harmonic+=1
      else: n_sinusoidal+=1
    s = item_label
    for l in labels :
      print("%s%s %s" % (prefix, s, l), file=outl)
      s = item_label_blank
    restraint._show_sorted_item(f=outl, prefix=prefix)

    if return_result:
          def value_from_restraint(restraint, key_ending = None):
            for x in dir(restraint):
              if not x.startswith("__") and x.endswith(key_ending):
                return getattr(restraint,x)

          #angle_ideal, angle_model, delta, variance**0.5, weight, residual
          delta = restraint.delta
          sigma = weight_as_sigma(weight = restraint.weight)
          ideal = value_from_restraint(restraint,'_ideal')
          model = value_from_restraint(restraint,'_model')
          residual = restraint.residual()
          if proxy_label == 'Dihedral angle':
            ideal = model + delta
            if ideal<-180: ideal+=360

          value = group_args(group_args_type = '%s result' %proxy_label,
            labels = labels,
            delta = delta,
            sigma = sigma,
            ideal = ideal,
            model = model,
            residual = residual)
          result.value_list.append(value)

  if (n_not_shown != 0):
    if (proxy_type is dihedral):
      n_harmonic = O.count_harmonic()
      n_sinusoidal = O.size() - n_harmonic
    print(prefix + "... (remaining %d not shown)" % n_not_shown, file=outl)
  #
  if (proxy_type is dihedral) and (origin_id==0 or origin_id is None):
    print(prefix+"  sinusoidal: %d" % n_sinusoidal, file=f)
    print(prefix+"    harmonic: %d" % n_harmonic, file=f)
  print("%sSorted by %s:" % (prefix, by_value), file=f)
  print(outl.getvalue()[:-1], file=f)

  if return_result:
    return result

class pair_proxies(object):

  def __init__(self,
        flags=None,
        bond_params_table=None,
        shell_asu_tables=None,
        model_indices=None,
        conformer_indices=None,
        sym_excl_indices=None,
        donor_acceptor_excl_groups=None,
        nonbonded_params=None,
        nonbonded_types=None,
        nonbonded_charges=None,
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
      if (nonbonded_charges is None):
        nonbonded_charges = flex.int(nonbonded_types.size(), 0)
      assert (nonbonded_types.size() == nonbonded_charges.size())
      self.nonbonded_proxies = nonbonded_sorted_asu_proxies(
        model_indices=model_indices,
        conformer_indices=conformer_indices,
        sym_excl_indices=sym_excl_indices,
        donor_acceptor_excl_groups=donor_acceptor_excl_groups,
        nonbonded_params=nonbonded_params,
        nonbonded_types=nonbonded_types,
        nonbonded_charges=nonbonded_charges,
        nonbonded_distance_cutoff_plus_buffer=\
          nonbonded_distance_cutoff_plus_buffer,
        min_cubicle_edge=min_cubicle_edge,
        shell_asu_tables=shell_asu_tables)

@bp.inject_into(ext.motif)
class _():

  def show(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    print(prefix+"geometry_restraints.motif {", file=out)
    print(prefix+"  id = %s" % show_string(self.id), file=out)
    print(prefix+"  description = %s" % show_string(self.description), file=out)
    for info in self.info:
      print(prefix+"  info = %s" % show_string(info), file=out)
    for manipulation_id in self.manipulation_ids:
      print(prefix+"  manipulation_id = %s" % (
        show_string(manipulation_id)), file=out)
    self.show_atoms(out=out, prefix=prefix+"  ")
    self.show_bonds(out=out, prefix=prefix+"  ")
    self.show_angles(out=out, prefix=prefix+"  ")
    self.show_dihedrals(out=out, prefix=prefix+"  ")
    self.show_chiralities(out=out, prefix=prefix+"  ")
    self.show_planarities(out=out, prefix=prefix+"  ")
    print(prefix+"}", file=out)

  def show_atoms(self, out=None, prefix=""):
    atoms = self.atoms_as_list()
    if (len(atoms) > 0):
      print(prefix+"atom = " \
        "[name scattering_type nonbonded_type partial_charge]", file=out)
      for atom in atoms:
        print(prefix+"atom = %s %s %s %.6g" % (
          show_string(atom.name),
          show_string(atom.scattering_type),
          show_string(atom.nonbonded_type),
          atom.partial_charge), file=out)

  def show_bonds(self, out=None, prefix=""):
    bonds = self.bonds_as_list()
    if (len(bonds) > 0):
      print(prefix+"bond = " \
        "[atom_name*2 type distance_ideal weight id]", file=out)
      for bond in bonds:
        atom_names = bond.atom_names
        print(prefix+"bond = %s %s %s %.6g %.6g %s" % (
          show_string(atom_names[0]),
          show_string(atom_names[1]),
          show_string(bond.type),
          bond.distance_ideal,
          bond.weight,
          show_string(bond.id)), file=out)

  def show_angles(self, out=None, prefix=""):
    angles = self.angles_as_list()
    if (len(angles) > 0):
      print(prefix+"angle = " \
        "[atom_name*3 angle_ideal weight id]", file=out)
      for angle in angles:
        atom_names = angle.atom_names
        print(prefix+"angle = %s %s %s %.6g %.6g %s" % (
          show_string(atom_names[0]),
          show_string(atom_names[1]),
          show_string(atom_names[2]),
          angle.angle_ideal,
          angle.weight,
          show_string(angle.id)), file=out)

  def show_dihedrals(self, out=None, prefix=""):
    dihedrals = self.dihedrals_as_list()
    if (len(dihedrals) > 0):
      print(prefix+"dihedral = " \
        "[atom_name*4 angle_ideal weight periodicity id]", file=out)
      for dihedral in dihedrals:
        atom_names = dihedral.atom_names
        print(prefix+"dihedral = %s %s %s %s %.6g %.6g %d %s" % (
          show_string(atom_names[0]),
          show_string(atom_names[1]),
          show_string(atom_names[2]),
          show_string(atom_names[3]),
          dihedral.angle_ideal,
          dihedral.weight,
          dihedral.periodicity,
          show_string(dihedral.id)), file=out)

  def show_chiralities(self, out=None, prefix=""):
    chiralities = self.chiralities_as_list()
    if (len(chiralities) > 0):
      print(prefix+"chirality = " \
        "[atom_name*4 volume_sign both_signs volume_ideal weight id]", file=out)
      for chirality in chiralities:
        atom_names = chirality.atom_names
        if (chirality.both_signs): both_signs = "True"
        else:                      both_signs = "False"
        print(prefix+"chirality = %s %s %s %s %s %s %.6g %.6g %s" % (
          show_string(atom_names[0]),
          show_string(atom_names[1]),
          show_string(atom_names[2]),
          show_string(atom_names[3]),
          show_string(chirality.volume_sign),
          both_signs,
          chirality.volume_ideal,
          chirality.weight,
          show_string(chirality.id)), file=out)

  def show_planarities(self, out=None, prefix=""):
    planarities = self.planarities_as_list()
    if (len(planarities) > 0):
      for planarity in planarities:
        print(prefix+"planarity {", file=out)
        print(prefix+"  id = %s" % show_string(planarity.id), file=out)
        assert planarity.weights.size() == planarity.atom_names.size()
        print(prefix+"  atom = [name weight]", file=out)
        for an,w in zip(planarity.atom_names, planarity.weights):
          print(prefix+"  atom = %s %.6g" % (show_string(an), w), file=out)
        print(prefix+"}", file=out)

@bp.inject_into(ext.motif_alteration)
class _():

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
        if (help != previous_help): print(help, file=out)
        print(prefix+"atom = add %s %s %s %s %s" % (
          show_string(self.motif_ids[0]),
          show_string(atom.name),
          show_string(atom.scattering_type),
          show_string(atom.nonbonded_type),
          atom.partial_charge), file=out)
      elif (action == "change"):
        help = prefix+"atom = change [motif_id motif_atom_name \\\n" \
                    + prefix+"               %s]" % attr
        if (help != previous_help): print(help, file=out)
        print(prefix+"atom = change %s %s \\" % (
          show_string(self.motif_ids[0]),
          show_string(self.motif_atom_name)), file=out)
        if (not self.change_partial_charge()):
          partial_charge = "None"
        else:
          partial_charge = "%.6g" % atom.partial_charge
        print(prefix+"              %s %s %s %s" % (
          show_string(atom.name),
          show_string(atom.scattering_type),
          show_string(atom.nonbonded_type),
          partial_charge), file=out)
      else:
        assert action == "delete"
        help = prefix+"atom = delete [motif_id motif_atom_name]"
        if (help != previous_help): print(help, file=out)
        print(prefix+"atom = delete %s %s" % (
          show_string(self.motif_ids[0]),
          show_string(self.motif_atom_name)), file=out)
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
        if (help != previous_help): print(help, file=out)
        print(prefix+"%s" % data_lead, file=out)
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
        if (help != previous_help): print(help, file=out)
        print(prefix+"%s %s %s %s" % (
          data_lead, distance_ideal, weight, show_string(bond.id)), file=out)
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
        if (help != previous_help): print(help, file=out)
        print(prefix+"%s" % data_lead, file=out)
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
        if (help != previous_help): print(help, file=out)
        print(prefix+"%s %s %s %s" % (
          data_lead, angle_ideal, weight, show_string(angle.id)), file=out)
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
        if (help != previous_help): print(help, file=out)
        print(prefix+"%s" % data_lead, file=out)
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
        if (help != previous_help): print(help, file=out)
        print(prefix+"%s %s %s %s %s" % (
          data_lead, angle_ideal, weight, periodicity,
          show_string(dihedral.id)), file=out)
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
        if (help != previous_help): print(help, file=out)
        print(prefix+"%s" % data_lead, file=out)
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
        if (help != previous_help): print(help, file=out)
        print(prefix+"%s \\\n%s%s%s %s %s %s" % (
          data_lead, prefix, " "*(13+len(action)),
          show_string(chirality.volume_sign),
          volume_ideal, weight, show_string(chirality.id)), file=out)
    elif (operand == "planarity"):
      planarity = self.planarity
      print(prefix+"planarity {", file=out)
      print(prefix+"  action = %s" % action, file=out)
      print(prefix+"  motif_id = %s" % show_string(
        self.planarity_motif_id), file=out)
      print(prefix+"  id = %s" % show_string(planarity.id), file=out)
      if (action == "add"):
        print(prefix+"  atom = [motif_id name weight]", file=out)
        assert planarity.weights.size() == planarity.atom_names.size()
        assert self.motif_ids.size() == planarity.atom_names.size()
        for mi,an,w in zip(self.motif_ids,
                           planarity.atom_names,
                           planarity.weights):
          print(prefix+"  atom = %s %s %.6g" % (
            show_string(mi), show_string(an), w), file=out)
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
            if (help != previous_help): print(help, file=out)
            print(prefix+"  atom = %s %s %s %.6g" % (
              ac, show_string(mi), show_string(an), w), file=out)
          else:
            help = prefix+"  atom = delete [motif_id name]"
            if (help != previous_help): print(help, file=out)
            print(prefix+"  atom = %s %s %s" % (
              ac, show_string(mi), show_string(an)), file=out)
          previous_help = help
      print(prefix+"}", file=out)
      help = None
    else:
      raise RuntimeError("Internal Error: unknown operand: %s" % operand)
    return help

@bp.inject_into(ext.motif_manipulation)
class _():

  def show(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    print(prefix+"geometry_restraints.motif_manipulation {", file=out)
    print(prefix+"  id = %s" % show_string(self.id), file=out)
    print(prefix+"  description = %s" % show_string(self.description), file=out)
    for info in self.info:
      print(prefix+"  info = %s" % show_string(info), file=out)
    previous_help = None
    for alteration in self.alterations_as_list():
      previous_help = alteration.show(
        out=out, prefix=prefix+"  ", previous_help=previous_help)
    print(prefix+"}", file=out)
