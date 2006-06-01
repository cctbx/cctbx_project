from iotbx.kriber import strudat
from cctbx import geometry_restraints
from cctbx import crystal
from cctbx.array_family import flex
import scitbx.math
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import libtbx.load_env
from cStringIO import StringIO
import math
import sys, os

def exercise_icosahedron(max_level=2, verbose=0):
  for level in xrange(0,max_level+1):
    if (0 or verbose):
      print "level:", level
    icosahedron = scitbx.math.icosahedron(level=level)
    try:
      distance_cutoff = icosahedron.next_neighbors_distance()*(1+1.e-3)
      estimated_distance_cutoff = False
    except RuntimeError, e:
      assert str(e) == "next_neighbors_distance not known."
      distance_cutoff = 0.4/(2**(level-1))
      estimated_distance_cutoff = True
    asu_mappings = crystal.direct_space_asu.non_crystallographic_asu_mappings(
      sites_cart=icosahedron.sites)
    pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
    pair_asu_table.add_all_pairs(distance_cutoff=distance_cutoff)
    if (0 or verbose):
      ps = pair_asu_table.show_distances(sites_cart=icosahedron.sites)
      print "level", level, "min", flex.min(ps.distances)
      print "     ", " ",   "max", flex.max(ps.distances)
      assert ps.pair_counts.all_eq(pair_asu_table.pair_counts())
      if (level == 0):
        for d in ps.distances:
          assert approx_equal(d, 1.0514622242382672)
    elif (level < 2):
      s = StringIO()
      ps = pair_asu_table.show_distances(sites_cart=icosahedron.sites, out=s)
      assert ps.pair_counts.all_eq(pair_asu_table.pair_counts())
      assert len(s.getvalue().splitlines()) == [72,320][level]
      del s
    if (level == 0):
      assert pair_asu_table.pair_counts().all_eq(5)
    else:
      assert pair_asu_table.pair_counts().all_eq(3)
    del pair_asu_table
    max_distance = crystal.neighbors_fast_pair_generator(
      asu_mappings=asu_mappings,
      distance_cutoff=distance_cutoff).max_distance_sq()**.5
    if (0 or verbose):
      print "max_distance:", max_distance
    if (not estimated_distance_cutoff):
      assert approx_equal(max_distance, icosahedron.next_neighbors_distance())
      assert approx_equal(max_distance/icosahedron.next_neighbors_distance(),1)

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
        return False
  return abs(d1-d2) < 1.e-3

def check_sym_equiv(structure, bond_asu_table, weak=False):
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

def exercise_incremental_pairs(
      structure,
      distance_cutoff,
      reference_pair_asu_table):
  ip = structure.incremental_pairs(distance_cutoff=distance_cutoff)
  for site_frac in structure.sites_frac():
    ip.process_site_frac(original_site=site_frac)
  assert ip.pair_asu_table().pair_counts().all_eq(
    reference_pair_asu_table.pair_counts())
  assert ip.pair_asu_table() == reference_pair_asu_table

def exercise_site_cluster_analysis(
      structure,
      distance_cutoff,
      reference_pair_asu_table):
  pat_selection = flex.size_t()
  pat_keep = []
  for i_seq,pair_asu_dict in enumerate(reference_pair_asu_table.table()):
    for j_seq,pair_asu_j_sym_groups in pair_asu_dict.items():
      if (j_seq == i_seq):
        for j_sym_group in pair_asu_j_sym_groups:
          assert 0 not in j_sym_group
        pat_keep.append(False)
        break
      if (j_seq < i_seq and pat_keep[j_seq]):
        pat_keep.append(False)
        break
    else:
      pat_keep.append(True)
      pat_selection.append(i_seq)
  assert reference_pair_asu_table.cluster_pivot_selection().all_eq(
    pat_selection)
  #
  sca = structure.site_cluster_analysis(distance_cutoff=distance_cutoff)
  sca_selection = flex.size_t()
  for i_seq,site_frac in enumerate(structure.sites_frac()):
    if (sca.process_site_frac(original_site=site_frac)):
      sca_selection.append(i_seq)
  assert sca_selection.all_eq(pat_selection)
  #
  sca = structure.site_cluster_analysis(distance_cutoff=distance_cutoff)
  sca_selection = sca.process_sites_frac(
    original_sites=structure.sites_frac(),
    site_symmetry_table=structure.site_symmetry_table())
  assert sca_selection.all_eq(pat_selection)
  #
  sca = structure.site_cluster_analysis(distance_cutoff=distance_cutoff)
  sca_selection = sca.process_sites_frac(
    original_sites=structure.sites_frac(),
    site_symmetry_table=structure.site_symmetry_table(),
    max_clusters=3)
  assert sca_selection.all_eq(pat_selection[:3])
  #
  sca = structure.site_cluster_analysis(distance_cutoff=distance_cutoff)
  sca_selection = sca.process_sites_frac(
    original_sites=structure.sites_frac())
  assert sca_selection.all_eq(pat_selection)
  #
  sca = structure.site_cluster_analysis(distance_cutoff=distance_cutoff)
  sca_selection = sca.process_sites_frac(
    original_sites=structure.sites_frac(),
    max_clusters=3)
  assert sca_selection.all_eq(pat_selection[:3])
  #
  sca = structure.site_cluster_analysis(distance_cutoff=distance_cutoff)
  sca_selection = sca.process_sites_cart(
    original_sites=structure.sites_cart(),
    site_symmetry_table=structure.site_symmetry_table())
  assert sca_selection.all_eq(pat_selection)
  #
  sca = structure.site_cluster_analysis(distance_cutoff=distance_cutoff)
  sca_selection = sca.process_sites_cart(
    original_sites=structure.sites_cart(),
    site_symmetry_table=structure.site_symmetry_table(),
    max_clusters=3)
  assert sca_selection.all_eq(pat_selection[:3])
  #
  sca = structure.site_cluster_analysis(distance_cutoff=distance_cutoff)
  sca_selection = sca.process_sites_cart(
    original_sites=structure.sites_cart())
  assert sca_selection.all_eq(pat_selection)
  #
  sca = structure.site_cluster_analysis(distance_cutoff=distance_cutoff)
  sca_selection = sca.process_sites_cart(
    original_sites=structure.sites_cart(),
    max_clusters=3)
  assert sca_selection.all_eq(pat_selection[:3])

def exercise(
      structure,
      distance_cutoff,
      connectivities=None,
      weak_check_sym_equiv=False,
      verbose=0):
  if (0 or verbose):
    print "distance_cutoff:", distance_cutoff
  asu_mappings = structure.asu_mappings(buffer_thickness=distance_cutoff)
  for i_pass in xrange(2):
    if (i_pass == 0):
      bond_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
      bond_asu_table.add_all_pairs(
        distance_cutoff=distance_cutoff)
      exercise_incremental_pairs(
        structure=structure,
        distance_cutoff=distance_cutoff,
        reference_pair_asu_table=bond_asu_table)
      exercise_site_cluster_analysis(
        structure=structure,
        distance_cutoff=distance_cutoff,
        reference_pair_asu_table=bond_asu_table)
    else:
      bond_sym_table = bond_asu_table.extract_pair_sym_table()
      bond_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
      bond_asu_table.add_pair_sym_table(
        sym_table=bond_sym_table)
    if (connectivities is not None):
      check_connectivities(bond_asu_table, connectivities, verbose)
    check_sym_equiv(
      structure=structure,
      bond_asu_table=bond_asu_table,
      weak=weak_check_sym_equiv)

def exercise_bond_sorted_asu_proxies(
      structure,
      distance_cutoff):
  asu_mappings = structure.asu_mappings(buffer_thickness=distance_cutoff)
  bond_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  bond_asu_table.add_all_pairs(distance_cutoff=distance_cutoff)
  bond_sym_table = bond_asu_table.extract_pair_sym_table()
  assert bond_sym_table.full_simple_connectivity().size() \
      == bond_sym_table.size()
  bond_params_table = geometry_restraints.bond_params_table(
    structure.scatterers().size())
  for i_seq,bond_sym_dict in enumerate(bond_sym_table):
    for j_seq in bond_sym_dict.keys():
      if (i_seq > j_seq):
        j_seq,i_seq = i_seq,j_seq
      bond_params_table[i_seq][j_seq] = geometry_restraints.bond_params(
        distance_ideal=3.1, weight=1)
  proxies_fast = geometry_restraints.bond_sorted_asu_proxies(
    bond_params_table=bond_params_table,
    bond_asu_table=bond_asu_table)
  pair_generator = crystal.neighbors_simple_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=distance_cutoff,
    minimal=False)
  proxies_slow = geometry_restraints.bond_sorted_asu_proxies(
    asu_mappings=asu_mappings)
  for pair in pair_generator:
    proxies_slow.process(geometry_restraints.bond_asu_proxy(
      pair=pair,
      distance_ideal=3.1,
      weight=1))
  assert proxies_slow.simple.size() == proxies_fast.simple.size()
  assert proxies_slow.asu.size() == proxies_fast.asu.size()
  ctrl = {}
  for proxy in proxies_slow.simple:
    assert not ctrl.has_key(proxy.i_seqs)
    ctrl[proxy.i_seqs] = 0
  for proxy in proxies_fast.simple:
    assert ctrl.has_key(proxy.i_seqs)
    ctrl[proxy.i_seqs] += 1
  assert ctrl.values() == [1]*len(ctrl)
  ctrl = {}
  for proxy in proxies_slow.asu:
    key = proxy.i_seq,proxy.j_seq,proxy.j_sym
    assert not ctrl.has_key(key)
    ctrl[key] = 0
  for proxy in proxies_fast.asu:
    key = proxy.i_seq,proxy.j_seq,proxy.j_sym
    assert ctrl.has_key(key)
    ctrl[key] += 1
  assert ctrl.values() == [1]*len(ctrl)

def exercise_all():
  verbose = "--verbose" in sys.argv[1:]
  exercise_icosahedron(verbose=verbose)
  default_distance_cutoff=3.5
  regression_misc = libtbx.env.find_in_repositories("regression/misc")
  if (regression_misc is None):
    print "Skipping exercise_all(): regression/misc not available"
    return
  file_names = []
  for file_name in ["strudat_zeolite_atlas", "strudat_special_bonds"]:
    path = os.path.join(regression_misc, file_name)
    if (not os.path.isfile(path)):
      print "Skipping %s test: input file not available" % file_name
    else:
      file_names.append(path)
  for file_name in file_names:
    strudat_entries = strudat.read_all_entries(open(file_name))
    for entry in strudat_entries.entries:
      if (0 or verbose):
        print "strudat tag:", entry.tag
      structure = entry.as_xray_structure()
      if (0 or verbose):
        structure.show_summary().show_scatterers()
      if (entry.title.startswith("cutoff")):
        distance_cutoff = float(entry.title.split()[1])
      else:
        distance_cutoff = default_distance_cutoff
      weak_check_sym_equiv = (
        entry.reference.find("weak_check_sym_equiv") >= 0)
      connectivities = entry.connectivities(all_or_nothing=True)
      if (1):
        exercise(
          structure=structure,
          distance_cutoff=distance_cutoff,
          connectivities=connectivities,
          weak_check_sym_equiv=weak_check_sym_equiv,
          verbose=verbose)
      if (0 or verbose):
        print
      if (file_name.endswith("strudat_zeolite_atlas")):
        exercise_bond_sorted_asu_proxies(
          structure=structure,
          distance_cutoff=distance_cutoff)

def run():
  exercise_all()
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
