from iotbx.kriber import strudat
import iotbx.pdb
from iotbx.option_parser import iotbx_option_parser
from cctbx import crystal
from scitbx.python_utils.misc import adopt_init_args
import libtbx.itertbx
import math
import sys, os

pair_asu_table = crystal.pair_asu_table

class _pair_asu_table:

  def __init__(self, asu_mappings):
    self.asu_mappings = asu_mappings
    self.table = crystal.pair_asu_table_table(
      self.asu_mappings.mappings().size())

  def add_all_pairs(self, distance_cutoff, verbose=0):
    pair_generator = crystal.neighbors_fast_pair_generator(
      asu_mappings=self.asu_mappings,
      distance_cutoff=distance_cutoff,
      minimal=0001)
    for pair in pair_generator:
      rt_mx_i = self.asu_mappings.get_rt_mx_i(pair=pair)
      rt_mx_j = self.asu_mappings.get_rt_mx_j(pair=pair)
      self.add_pair(
        i_seq=pair.i_seq,
        j_seq=pair.j_seq,
        rt_mx_ji=rt_mx_i.inverse().multiply(rt_mx_j),
        verbose=verbose)
    return self

  def add_pair_sym_table(self, sym_table, verbose=0):
    for i_seq,sym_dict in enumerate(sym_table):
      for j_seq,rt_mx_list in sym_dict.items():
        for rt_mx_ji in rt_mx_list:
          self.add_pair(i_seq=i_seq, j_seq=j_seq, rt_mx_ji=rt_mx_ji)
    return self

  def add_pair(self, i_seq, j_seq, rt_mx_ji, verbose=0):
    is_new = self._process_pair(
      i_seq=i_seq,
      j_seq=j_seq,
      rt_mx_ji=rt_mx_ji,
      verbose=verbose)
    if (is_new and i_seq != j_seq):
      is_new = self._process_pair(
        i_seq=j_seq,
        j_seq=i_seq,
        rt_mx_ji=rt_mx_ji.inverse_cancel(),
        verbose=verbose)
      assert is_new
    return self

  def _process_pair(self, i_seq, j_seq, rt_mx_ji, verbose=0):
    rt_mx_i = self.asu_mappings.get_rt_mx(i_seq=i_seq, i_sym=0)
    rt_mx_j = rt_mx_i.multiply(rt_mx_ji)
    j_sym = self.asu_mappings.find_i_sym(i_seq=j_seq, rt_mx=rt_mx_j)
    assert j_sym >= 0
    j_sym_groups = self.table[i_seq].setdefault(j_seq, [])
    for j_sym_group in j_sym_groups:
      if (j_sym in j_sym_group):
        return 00000
    j_sym_groups.append([])
    j_syms = j_sym_groups[-1]
    if (0 or verbose):
      print "primary:     i_seq, j_seq, j_sym", i_seq, j_seq, j_sym
    for mi in self.asu_mappings.site_symmetry_table().get(i_seq).matrices():
      if (i_seq == j_seq):
        rt_mx_j_eq = rt_mx_i.multiply(rt_mx_ji.multiply(mi).inverse_cancel())
        j_sym_eq = self.asu_mappings.find_i_sym(i_seq=j_seq, rt_mx=rt_mx_j_eq)
        assert j_sym_eq >= 0
        if (not j_sym_eq in j_syms):
          j_syms.append(j_sym_eq)
          if (0 or verbose):
            print "    equiv a: i_seq, j_seq, j_sym", i_seq, j_seq, j_sym_eq
      rt_mx_j_eq = rt_mx_i.multiply(mi.multiply(rt_mx_ji))
      j_sym_eq = self.asu_mappings.find_i_sym(i_seq=j_seq, rt_mx=rt_mx_j_eq)
      assert j_sym_eq >= 0
      if (not j_sym_eq in j_syms):
        j_syms.append(j_sym_eq)
        if (0 or verbose):
          print "    equiv b: i_seq, j_seq, j_sym", i_seq, j_seq, j_sym_eq
    return 0001

  def __contains__(self, pair):
    j_sym_groups = self.table[pair.i_seq].get(pair.j_seq, None)
    if (j_sym_groups is not None):
      for j_sym_group in j_sym_groups:
        if (pair.j_sym in j_sym_group):
          return 0001
    return 00000

  def extract_pair_sym_table(self):
    sym_table = crystal.pair_sym_table(self.asu_mappings.mappings().size())
    for i_seq,j_seq_dict in enumerate(self.table):
      rt_mx_i_inv = self.asu_mappings.get_rt_mx(i_seq=i_seq, i_sym=0).inverse()
      for j_seq,j_sym_groups in j_seq_dict.items():
        if (j_seq < i_seq): continue
        for j_sym_group in j_sym_groups:
          j_sym = j_sym_group[0]
          rt_mx_j = self.asu_mappings.get_rt_mx(i_seq=j_seq, i_sym=j_sym)
          sym_table[i_seq].setdefault(j_seq).append(
            rt_mx_i_inv.multiply(rt_mx_j))
    return sym_table

def is_sym_equiv_interaction_simple(unit_cell,
                                    i_seq,
                                    site_frac_i,
                                    j_seq,
                                    site_frac_j,
                                    special_op_j,
                                    rt_mx_ji_1,
                                    rt_mx_ji_2):
  f = unit_cell.shortest_vector_sq()**.5*.1
  trial_shifts = [f*x for x in [math.sqrt(2),math.sqrt(3),math.sqrt(5)]]
  frac = unit_cell.fractionalize
  orth = unit_cell.orthogonalize
  dist = unit_cell.distance
  for shifts in [[0,0,0], trial_shifts]:
    site_j_mod = special_op_j * frac([x+s
      for x,s in zip(orth(site_frac_j),shifts)])
    if (shifts == [0,0,0] or j_seq != i_seq):
      site_i_mod = site_frac_i
    else:
      site_i_mod = site_j_mod
    d1 = dist(rt_mx_ji_1 * site_j_mod, site_i_mod)
    d2 = dist(rt_mx_ji_2 * site_j_mod, site_i_mod)
    if (shifts == [0,0,0]):
      if (abs(d1-d2) >= 1.e-3):
        return 00000
  return abs(d1-d2) < 1.e-3

def check_sym_equiv(structure, bond_asu_table, weak=00000):
  unit_cell = structure.unit_cell()
  asu_mappings = bond_asu_table.asu_mappings()
  sites_frac = structure.scatterers().extract_sites()
  for i_seq,records in enumerate(bond_asu_table.table()):
    rt_mx_i_inv = asu_mappings.get_rt_mx(i_seq, 0).inverse()
    for j_seq,j_sym_groups in records.items():
      i_group_rt_mx_jis = []
      for i_group,j_sym_group in enumerate(j_sym_groups):
        for j_sym in j_sym_group:
          rt_mx_ji = rt_mx_i_inv.multiply(asu_mappings.get_rt_mx(j_seq, j_sym))
          i_group_rt_mx_jis.append((i_group,rt_mx_ji))
      for gi,ri in i_group_rt_mx_jis:
        for gj,rj in i_group_rt_mx_jis:
          is_sym_equiv = is_sym_equiv_interaction_simple(
            unit_cell=unit_cell,
            i_seq=i_seq,
            site_frac_i=sites_frac[i_seq],
            j_seq=j_seq,
            site_frac_j=sites_frac[j_seq],
            special_op_j=asu_mappings.special_op(j_seq),
            rt_mx_ji_1=ri,
            rt_mx_ji_2=rj)
          if (is_sym_equiv):
            if (not weak): assert gi == gj
          else:
            assert gi != gj

def check_connectivities(bond_asu_table, connectivities, verbose=0):
  n_mismatches = 0
  for records,connectivity in zip(bond_asu_table.table(), connectivities):
    n = 0
    for j_seq,j_sym_groups in records.items():
      for j_sym_group in j_sym_groups:
        n += len(j_sym_group)
    if (0 or verbose):
      print "n, connectivity:", n, connectivity
    assert n == connectivity

def exercise(
      structure,
      distance_cutoff,
      connectivities=None,
      weak_check_sym_equiv=00000,
      verbose=0):
  if (0 or verbose):
    print "distance_cutoff:", distance_cutoff
  asu_mappings = structure.asu_mappings(buffer_thickness=distance_cutoff)
  for i_pass in xrange(2):
    if (i_pass == 0):
      bond_asu_table = pair_asu_table(asu_mappings=asu_mappings)
      bond_asu_table.add_all_pairs(
        distance_cutoff=distance_cutoff)
    else:
      bond_sym_table = bond_asu_table.extract_pair_sym_table()
      bond_asu_table = pair_asu_table(asu_mappings=asu_mappings)
      bond_asu_table.add_pair_sym_table(
        sym_table=bond_sym_table)
    if (connectivities is not None):
      check_connectivities(bond_asu_table, connectivities, verbose)
    check_sym_equiv(
      structure=structure,
      bond_asu_table=bond_asu_table,
      weak=weak_check_sym_equiv)

def run():
  command_line = (iotbx_option_parser(
    usage="python neighbors.py [options] studat_file [...]",
    description="Example: python neighbors.py strudat --tag=SOD")
    .enable_symmetry_comprehensive()
    .option(None, "--tag",
      action="store",
      type="string",
      dest="tag",
      help="tag as it appears in the strudat file")
    .option(None, "--distance_cutoff",
      action="store",
      type="float",
      dest="distance_cutoff",
      metavar="FLOAT")
    .option(None, "--verbose",
      action="store_true",
      dest="verbose",
      help="produce output")
  ).process()
  default_distance_cutoff=3.5
  file_names = command_line.args
  if (len(file_names) == 0):
    regression_misc = os.path.join(
      os.environ["LIBTBX_DIST_ROOT"], "regression", "misc")
    for file_name in ["strudat_zeolite_atlas", "strudat_special_bonds"]:
      path = os.path.join(regression_misc, file_name)
      if (os.path.isfile(path)):
        file_names.append(path)
  for file_name in file_names:
    try:
      strudat_entries = strudat.read_all_entries(open(file_name))
    except:
      strudat_entries = None
    if (strudat_entries is not None and len(strudat_entries.entries) > 0):
      for entry in strudat_entries.entries:
        if (    command_line.options.tag is not None
            and command_line.options.tag != entry.tag):
          continue
        if (0 or command_line.options.verbose):
          print "strudat tag:", entry.tag
        structure = entry.as_xray_structure()
        if (0 or command_line.options.verbose):
          structure.show_summary().show_scatterers()
        if (command_line.options.distance_cutoff is not None):
          distance_cutoff = command_line.options.distance_cutoff
        elif (entry.title.startswith("cutoff")):
          distance_cutoff = float(entry.title.split()[1])
        else:
          distance_cutoff = default_distance_cutoff
        weak_check_sym_equiv = (
          entry.reference.find("weak_check_sym_equiv") >= 0)
        connectivities = entry.connectivities(all_or_nothing=0001)
        exercise(
          structure=structure,
          distance_cutoff=distance_cutoff,
          connectivities=connectivities,
          weak_check_sym_equiv=weak_check_sym_equiv,
          verbose=command_line.options.verbose)
        if (0 or command_line.options.verbose):
          print
    else:
      try:
        structure = iotbx.pdb.as_xray_structure(
          file_name=file_name,
          crystal_symmetry=command_line.symmetry)
      except:
        raise RuntimeError("Coordinate file %s: unknown format." % file_name)
      else:
        if (0 or command_line.options.verbose):
          structure.show_summary().show_scatterers()
        exercise(
          structure=structure,
          distance_cutoff=distance_cutoff,
          verbose=command_line.options.verbose)
  print "OK"

if (__name__ == "__main__"):
  run()
