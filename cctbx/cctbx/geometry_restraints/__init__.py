import cctbx.crystal.direct_space_asu
from cctbx.array_family import flex
import scitbx.array_family.shared
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
        self.source_labels[i_list].append(source_info.labels())
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

class pair_labels_formatter:

  def __init__(self, sorted_proxies, i_proxies_sorted, labels):
    self.labels = labels
    n_simple = sorted_proxies.simple.size()
    n_asu = 0
    max_label_lengths = [0,0]
    for i_proxy in i_proxies_sorted:
      if (i_proxy < n_simple):
        i_seqs = sorted_proxies.simple[i_proxy].i_seqs
      else:
        proxy = sorted_proxies.asu[i_proxy-n_simple]
        i_seqs = (proxy.i_seq, proxy.j_seq)
        n_asu += 1
      if (labels is not None):
        max_label_lengths = [max(m, len(labels[i_seq]))
          for m,i_seq in zip(max_label_lengths, i_seqs)]
    if (max(max_label_lengths) == 0):
      self.labels = None
    else:
      self.label_label_format = \
        " - ".join(["%%-%ds" % m for m in max_label_lengths]) + " "
      self.atom_i_atom_j_format = self.label_label_format.replace("-", "", 1)
    self.max_label_lengths = max_label_lengths
    if (n_asu > 0):
      self.sym_op_j = " sym.op. j"
    else:
      self.sym_op_j = ""

  def atom_i_atom_j(self):
    if (min(self.max_label_lengths) >= 6):
      return self.atom_i_atom_j_format % ("atom i", "atom j")
    if (min(self.max_label_lengths) >= 1):
      return self.atom_i_atom_j_format % ("i", "j")
    return ""

  def label_label(self, i_seq, j_seq):
    if (self.labels is None): return ""
    return self.label_label_format % (self.labels[i_seq], self.labels[j_seq])

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

  def show_histogram_of_deltas(self,
        sites_cart,
        n_slots=5,
        f=None,
        prefix=""):
    if (self.n_total() == 0): return
    if (f is None): f = sys.stdout
    print >> f, "%sHistogram of bond deltas:" % prefix
    histogram = flex.histogram(
      data=flex.abs(bond_deltas(
        sites_cart=sites_cart,
        sorted_asu_proxies=self)),
      n_slots=n_slots)
    low_cutoff = histogram.data_min()
    for i,n in enumerate(histogram.slots()):
      high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
      print >> f, "%s  %8.3f - %8.3f: %d" % (
        prefix, low_cutoff, high_cutoff, n)
      low_cutoff = high_cutoff
    return histogram

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
    plf = pair_labels_formatter(
      sorted_proxies=self,
      i_proxies_sorted=i_proxies_sorted,
      labels=labels)
    if (self.asu.size() == 0):
      asu_mappings = None
    else:
      asu_mappings = self.asu_mappings()
    n_simple = self.simple.size()
    show_legend = True
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
        rt_mx = asu_mappings.get_rt_mx_ji(pair=proxy)
        restraint = bond(
          sites_cart=sites_cart,
          asu_mappings=asu_mappings,
          proxy=proxy)
      if (show_legend):
        show_legend = False
        print >> f, "%sBond restraints sorted by residual:" % prefix
        print >> f, "%s%sideal  model  delta   weight residual%s" % (
          prefix, plf.atom_i_atom_j(), plf.sym_op_j)
      print >> f, "%s%s%5.3f %6.3f %6.3f %6.2e %6.2e" % (
        prefix, plf.label_label(i_seq, j_seq),
        restraint.distance_ideal, restraint.distance_model, restraint.delta,
        restraint.weight, restraint.residual()),
      if (rt_mx is not None):
        print >> f, rt_mx,
      print >> f
    n_not_shown = residuals.size() - i_proxies_sorted.size()
    if (n_not_shown != 0):
      print >> f, prefix + "... (remaining %d not shown)" % n_not_shown

class _nonbonded_sorted_asu_proxies(boost.python.injector,
        nonbonded_sorted_asu_proxies):

  def deltas(self, sites_cart):
    return nonbonded_deltas(
      sites_cart=sites_cart,
      sorted_asu_proxies=self,
      function=prolsq_repulsion_function())

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

  def show_sorted_by_model_distance(self,
        sites_cart,
        labels=None,
        f=None,
        prefix="",
        max_lines=None):
    if (f is None): f = sys.stdout
    deltas = nonbonded_deltas(
      sites_cart=sites_cart,
      sorted_asu_proxies=self,
      function=prolsq_repulsion_function())
    i_proxies_sorted = flex.sort_permutation(data=deltas)
    if (max_lines is not None and i_proxies_sorted.size() > max_lines+1):
      i_proxies_sorted = i_proxies_sorted[:max_lines]
    plf = pair_labels_formatter(
      sorted_proxies=self,
      i_proxies_sorted=i_proxies_sorted,
      labels=labels)
    if (self.asu.size() == 0):
      asu_mappings = None
    else:
      asu_mappings = self.asu_mappings()
    show_legend = True
    n_simple = self.simple.size()
    for i_proxy in i_proxies_sorted:
      if (i_proxy < n_simple):
        proxy = self.simple[i_proxy]
        i_seq,j_seq = proxy.i_seqs
        rt_mx = None
      else:
        proxy = self.asu[i_proxy-n_simple]
        i_seq,j_seq = proxy.i_seq,proxy.j_seq
        rt_mx = asu_mappings.get_rt_mx_ji(pair=proxy)
      delta = deltas[i_proxy]
      if (show_legend):
        show_legend = False
        print >> f, "%sNonbonded interactions sorted by model distance:" % (
          prefix)
        print >> f, "%s%s model   vdw%s" % (
          prefix, plf.atom_i_atom_j(), plf.sym_op_j)
      print >> f, "%s%s%6.3f %5.3f" % (
        prefix, plf.label_label(i_seq, j_seq), delta, proxy.vdw_distance),
      if (rt_mx is not None):
        print >> f, rt_mx,
      print >> f
    n_not_shown = deltas.size() - i_proxies_sorted.size()
    if (n_not_shown != 0):
      print >> f, prefix + "... (remaining %d not shown)" % n_not_shown

class _shared_dihedral_proxy(boost.python.injector, shared_dihedral_proxy):

  def show_histogram_of_deltas(self,
        sites_cart,
        n_slots=5,
        f=None,
        prefix=""):
    if (self.size() == 0): return
    if (f is None): f = sys.stdout
    print >> f, "%sHistogram of dihedral angle deviations from ideal:" % prefix
    histogram = flex.histogram(
      data=flex.abs(dihedral_deltas(
        sites_cart=sites_cart,
        proxies=self)),
      n_slots=n_slots)
    low_cutoff = histogram.data_min()
    for i,n in enumerate(histogram.slots()):
      high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
      print >> f, "%s  %8.2f - %8.2f: %d" % (
        prefix, low_cutoff, high_cutoff, n)
      low_cutoff = high_cutoff
    return histogram

  def show_sorted_by_residual(self,
        sites_cart,
        labels=None,
        f=None,
        prefix="",
        max_lines=None):
    assert labels is None or len(labels) == sites_cart.size()
    if (self.size() == 0): return
    if (f is None): f = sys.stdout
    residuals = dihedral_residuals(
      sites_cart=sites_cart,
      proxies=self)
    i_proxies_sorted = flex.sort_permutation(data=residuals, reverse=True)
    if (max_lines is not None and i_proxies_sorted.size() > max_lines+1):
      i_proxies_sorted = i_proxies_sorted[:max_lines]
    print >> f, "%sDihedral angle restraints sorted by residual:" % prefix
    for i_proxy in i_proxies_sorted:
      proxy = self[i_proxy]
      restraint = dihedral(
        sites_cart=sites_cart,
        proxy=proxy)
      if (labels is not None):
        for i_seq in proxy.i_seqs:
          print >> f, "%s%s" % (prefix, labels[i_seq])
      print >> f, "%s    ideal   model   delta" \
        " periodicty    weight residual" % prefix
      print >> f, "%s  %7.2f %7.2f %7.2f %5d       %6.2e %6.2e" % (
        prefix,
        restraint.angle_ideal, restraint.angle_model, restraint.delta,
        restraint.periodicity, restraint.weight, restraint.residual())
    n_not_shown = self.size() - i_proxies_sorted.size()
    if (n_not_shown != 0):
      print >> f, prefix + "... (remaining %d not shown)" % n_not_shown

class pair_proxies(object):

  def __init__(self,
        flags=None,
        bond_params_table=None,
        shell_asu_tables=None,
        model_indices=None,
        conformer_indices=None,
        nonbonded_params=None,
        nonbonded_types=None,
        nonbonded_distance_cutoff_plus_buffer=None):
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
      unit_cell = shell_asu_tables[0].asu_mappings().unit_cell()
      if (  nonbonded_distance_cutoff_plus_buffer**2
          > unit_cell.shortest_vector_sq()):
        raise RuntimeError(
          "Nonbonded distance cutoff + buffer"
          " > shortest lattice translation vector: %.6g > %.6g" % (
            nonbonded_distance_cutoff_plus_buffer,
            unit_cell.shortest_vector_sq()**0.5))
      self.nonbonded_proxies = nonbonded_sorted_asu_proxies(
        model_indices=model_indices,
        conformer_indices=conformer_indices,
        nonbonded_params=nonbonded_params,
        nonbonded_types=nonbonded_types,
        nonbonded_distance_cutoff_plus_buffer=\
          nonbonded_distance_cutoff_plus_buffer,
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
