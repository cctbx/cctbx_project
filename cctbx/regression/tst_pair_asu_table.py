from iotbx.kriber import strudat
from cctbx import geometry_restraints
from cctbx import crystal
from cctbx.array_family import flex
import scitbx.math
from scitbx import matrix
from libtbx.test_utils import approx_equal, show_diff
from libtbx.utils import format_cpu_times
import libtbx.load_env
from libtbx import dict_with_default_0
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
      ps = pair_asu_table.show_distances(sites_cart=icosahedron.sites) \
        .distances_info
      print "level", level, "min", flex.min(ps.distances)
      print "     ", " ",   "max", flex.max(ps.distances)
      assert ps.pair_counts.all_eq(pair_asu_table.pair_counts())
      if (level == 0):
        for d in ps.distances:
          assert approx_equal(d, 1.0514622242382672)
    elif (level < 2):
      s = StringIO()
      ps = pair_asu_table.show_distances(sites_cart=icosahedron.sites, out=s) \
        .distances_info
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
  assert reference_pair_asu_table.cluster_pivot_selection(
    max_clusters=3).all_eq(pat_selection[:3])
  #
  sca = structure.site_cluster_analysis(min_distance=distance_cutoff)
  sca_selection = flex.size_t()
  for i_seq,site_frac in enumerate(structure.sites_frac()):
    if (sca.process_site_frac(original_site=site_frac)):
      sca_selection.append(i_seq)
  assert sca_selection.all_eq(pat_selection)
  #
  sca = structure.site_cluster_analysis(min_distance=distance_cutoff)
  sca_selection = sca.process_sites_frac(
    original_sites=structure.sites_frac(),
    site_symmetry_table=structure.site_symmetry_table())
  assert sca_selection.all_eq(pat_selection)
  #
  sca = structure.site_cluster_analysis(min_distance=distance_cutoff)
  sca_selection = sca.process_sites_frac(
    original_sites=structure.sites_frac(),
    site_symmetry_table=structure.site_symmetry_table(),
    max_clusters=3)
  assert sca_selection.all_eq(pat_selection[:3])
  #
  sca = structure.site_cluster_analysis(min_distance=distance_cutoff)
  sca_selection = sca.process_sites_frac(
    original_sites=structure.sites_frac())
  assert sca_selection.all_eq(pat_selection)
  #
  sca = structure.site_cluster_analysis(min_distance=distance_cutoff)
  sca_selection = sca.process_sites_frac(
    original_sites=structure.sites_frac(),
    max_clusters=3)
  assert sca_selection.all_eq(pat_selection[:3])
  #
  sca = structure.site_cluster_analysis(min_distance=distance_cutoff)
  sca_selection = sca.process_sites_cart(
    original_sites=structure.sites_cart(),
    site_symmetry_table=structure.site_symmetry_table())
  assert sca_selection.all_eq(pat_selection)
  #
  sca = structure.site_cluster_analysis(min_distance=distance_cutoff)
  sca_selection = sca.process_sites_cart(
    original_sites=structure.sites_cart(),
    site_symmetry_table=structure.site_symmetry_table(),
    max_clusters=3)
  assert sca_selection.all_eq(pat_selection[:3])
  #
  sca = structure.site_cluster_analysis(min_distance=distance_cutoff)
  sca_selection = sca.process_sites_cart(
    original_sites=structure.sites_cart())
  assert sca_selection.all_eq(pat_selection)
  #
  sca = structure.site_cluster_analysis(min_distance=distance_cutoff)
  sca_selection = sca.process_sites_cart(
    original_sites=structure.sites_cart(),
    max_clusters=3)
  assert sca_selection.all_eq(pat_selection[:3])
  #
  sca = structure.site_cluster_analysis(
    min_distance=distance_cutoff,
    general_positions_only=True)
  sca_selection = sca.process_sites_frac(
    original_sites=structure.sites_frac(),
    site_symmetry_table=structure.site_symmetry_table())
  pat_selection = reference_pair_asu_table.cluster_pivot_selection(
    general_positions_only=True)
  assert sca_selection.all_eq(pat_selection)

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
      def exercise_symmetry_equivalent_pair_interactions():
        asu_mappings = bond_asu_table.asu_mappings()
        for i_seq, j_seq_dict in enumerate(bond_asu_table.table()):
          rt_mx_i = asu_mappings.get_rt_mx(i_seq, 0)
          rt_mx_i_inv = rt_mx_i.inverse()
          for j_seq,j_sym_group in j_seq_dict.items():
            scs = structure.scatterers()
            def get_coords(symops):
              result = []
              for s in symops:
                result.append(numstr(s * scs[j_seq].site))
              result.sort()
              return result
            prev_equiv_rt_mx_ji = None
            for j_syms in j_sym_group:
              equiv_rt_mx_ji = []
              for j_sym in j_syms:
                rt_mx_ji = rt_mx_i_inv.multiply(
                  asu_mappings.get_rt_mx(j_seq, j_sym))
                equiv_rt_mx_ji.append(rt_mx_ji)
              old_coords = get_coords(equiv_rt_mx_ji)
              all_sepi = set()
              for rt_mx_ji in equiv_rt_mx_ji:
                _ = asu_mappings.site_symmetry_table()
                sepi_obj = _.symmetry_equivalent_pair_interactions(
                  i_seq=i_seq, j_seq=j_seq, rt_mx_ji=rt_mx_ji)
                sepi = sepi_obj.get()
                new_coords = get_coords(sepi)
                assert new_coords == old_coords
                all_sepi.add(";".join([str(_) for _ in sepi]))
                for _ in equiv_rt_mx_ji:
                  assert sepi_obj.is_equivalent(rt_mx_ji=_)
                if (prev_equiv_rt_mx_ji is not None):
                  for _ in prev_equiv_rt_mx_ji:
                    assert not sepi_obj.is_equivalent(rt_mx_ji=_)
              assert len(all_sepi) == 1
              prev_equiv_rt_mx_ji = equiv_rt_mx_ji
      exercise_symmetry_equivalent_pair_interactions()
      def exercise_pair_sym_table_tidy_and_full_connectivity():
        def check_one_way(pst):
          for sym_pair in pst.iterator():
            i_seq, j_seq = sym_pair.i_seqs()
            assert i_seq <= j_seq
            assert len(pst[i_seq][j_seq]) > 0
            if (i_seq != j_seq):
              assert i_seq not in pst[j_seq]
        def check_two_way(pst):
          for sym_pair in pst.iterator():
            i_seq, j_seq = sym_pair.i_seqs()
            assert len(pst[i_seq][j_seq]) > 0
            assert len(pst[j_seq][i_seq]) > 0
        pst_extracted = bond_sym_table
        check_one_way(pst_extracted)
        sio_extracted = StringIO()
        structure.pair_sym_table_show(pst_extracted, out=sio_extracted)
        pst = pst_extracted.tidy(
          site_symmetry_table=structure.site_symmetry_table())
        check_one_way(pst)
        sio = StringIO()
        structure.pair_sym_table_show(pst, out=sio)
        assert not show_diff(sio.getvalue(), sio_extracted.getvalue())
        pst = pst_extracted.full_connectivity()
        check_two_way(pst)
        pst_full = pst_extracted.full_connectivity(
          site_symmetry_table=structure.site_symmetry_table())
        check_two_way(pst_full)
        sio = StringIO()
        structure.pair_sym_table_show(
          pst_full, is_full_connectivity=True, out=sio)
        assert sio.getvalue().find("sym. equiv.") < 0
        pst = pst_full.tidy(
          site_symmetry_table=structure.site_symmetry_table())
        check_one_way(pst)
        sio = StringIO()
        structure.pair_sym_table_show(pst, out=sio)
        assert not show_diff(sio.getvalue(), sio_extracted.getvalue())
        pst_full2 = pst_full.full_connectivity(
          site_symmetry_table=structure.site_symmetry_table())
        check_two_way(pst_full2)
        pst = pst_full2.tidy(
          site_symmetry_table=structure.site_symmetry_table())
        check_one_way(pst)
        sio = StringIO()
        structure.pair_sym_table_show(pst, out=sio)
        assert not show_diff(sio.getvalue(), sio_extracted.getvalue())
      exercise_pair_sym_table_tidy_and_full_connectivity()
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
  el = bond_sym_table.simple_edge_list()
  es = bond_sym_table.full_simple_connectivity()
  assert es.size() == bond_sym_table.size()
  for i,j in el:
    assert j in es[i]
    assert i in es[j]
  npis = bond_sym_table.number_of_pairs_involving_symmetry()
  assert len(list(bond_sym_table.iterator())) == len(el) + npis
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
  proxies_conservative = geometry_restraints.bond_sorted_asu_proxies(
    pair_asu_table=bond_asu_table)
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
  def compare_proxies(proxies_1, proxies_2):
    assert proxies_1.simple.size() == proxies_2.simple.size()
    assert proxies_1.asu.size() == proxies_2.asu.size()
    ctrl = {}
    for proxy in proxies_1.simple:
      assert not ctrl.has_key(proxy.i_seqs)
      ctrl[proxy.i_seqs] = 0
    for proxy in proxies_2.simple:
      assert ctrl.has_key(proxy.i_seqs)
      ctrl[proxy.i_seqs] += 1
    assert ctrl.values() == [1]*len(ctrl)
    ctrl = {}
    for proxy in proxies_1.asu:
      key = proxy.i_seq,proxy.j_seq,proxy.j_sym
      assert not ctrl.has_key(key)
      ctrl[key] = 0
    for proxy in proxies_2.asu:
      key = proxy.i_seq,proxy.j_seq,proxy.j_sym
      assert ctrl.has_key(key)
      ctrl[key] += 1
    assert ctrl.values() == [1]*len(ctrl)
  compare_proxies(proxies_1=proxies_fast, proxies_2=proxies_conservative)
  compare_proxies(proxies_1=proxies_fast, proxies_2=proxies_slow)
  sites_cart = structure.sites_cart()
  for proxy in proxies_conservative.simple:
    i,j = proxy.i_seqs
    assert approx_equal(
      abs(matrix.col(sites_cart[i]) - matrix.col(sites_cart[j])),
      proxy.distance_ideal)
    assert proxy.weight == 1
  distance = proxies_conservative.asu_mappings().unit_cell().distance
  get_rt_mx_ji = proxies_conservative.asu_mappings().get_rt_mx_ji
  sites_frac = structure.sites_frac()
  for proxy in proxies_conservative.asu:
    assert approx_equal(
      distance(
        sites_frac[proxy.i_seq],
        get_rt_mx_ji(pair=proxy) * sites_frac[proxy.j_seq]),
      proxy.distance_ideal)
    assert proxy.weight == 1

def py_pair_asu_table_angle_pair_asu_table(self):
  asu_mappings = self.asu_mappings()
  result = crystal.pair_asu_table(asu_mappings=asu_mappings)
  for i_seq,asu_dict in enumerate(self.table()):
    pair_list = []
    for j_seq,j_sym_groups in asu_dict.items():
      for i_group,j_sym_group in enumerate(j_sym_groups):
        for j_sym in j_sym_group:
          pair_list.append((j_seq,j_sym))
    for i_jj1 in xrange(0,len(pair_list)-1):
      jj1 = pair_list[i_jj1]
      rt_mx_jj1_inv = asu_mappings.get_rt_mx(*jj1).inverse()
      for i_jj2 in xrange(i_jj1+1,len(pair_list)):
        jj2 = pair_list[i_jj2]
        result.add_pair(
          i_seq=jj1[0],
          j_seq=jj2[0],
          rt_mx_ji=rt_mx_jj1_inv.multiply(asu_mappings.get_rt_mx(*jj2)))
  return result

def exercise_angle_pair_asu_table(
      structure,
      distance_cutoff,
      connectivities,
      reference_apatanl,
      reference_cppc):
  sg_asu_mappings = structure.asu_mappings(
    buffer_thickness=2*distance_cutoff)
  sg_pat = crystal.pair_asu_table(asu_mappings=sg_asu_mappings)
  sg_pat.add_all_pairs(
    distance_cutoff=distance_cutoff,
    min_cubicle_edge=0)
  # compare connectivities with reference
  assert list(sg_pat.pair_counts()) == connectivities
  #
  p1_structure = structure.expand_to_p1()
  p1_asu_mappings = p1_structure.asu_mappings(
    buffer_thickness=2*distance_cutoff)
  p1_pat = crystal.pair_asu_table(asu_mappings=p1_asu_mappings)
  p1_pat.add_all_pairs(
    distance_cutoff=distance_cutoff,
    min_cubicle_edge=0)
  sg_labels = structure.scatterers().extract_labels()
  p1_labels = p1_structure.scatterers().extract_labels()
  label_connect = dict(zip(sg_labels, sg_pat.pair_counts()))
  for l,c in zip(p1_labels, p1_pat.pair_counts()):
    # compare connectivities in original space group and in P1
    assert label_connect[l] == c
  #
  sg_apat_py = py_pair_asu_table_angle_pair_asu_table(self=sg_pat)
  sg_apat = sg_pat.angle_pair_asu_table()
  assert sg_apat.as_nested_lists() == sg_apat_py.as_nested_lists()
  sg_counts = {}
  for i_seq,pair_asu_dict in enumerate(sg_apat.table()):
    lbl_i = sg_labels[i_seq]
    for j_seq,pair_asu_j_sym_groups in pair_asu_dict.items():
      lbl_j = sg_labels[j_seq]
      for j_sym_group in pair_asu_j_sym_groups:
        sg_counts.setdefault(lbl_i, dict_with_default_0())[
                             lbl_j] += len(j_sym_group)
  p1_apat = p1_pat.angle_pair_asu_table()
  p1_counts = {}
  for i_seq,pair_asu_dict in enumerate(p1_apat.table()):
    lbl_i = p1_labels[i_seq]
    for j_seq,pair_asu_j_sym_groups in pair_asu_dict.items():
      lbl_j = p1_labels[j_seq]
      for j_sym_group in pair_asu_j_sym_groups:
        p1_counts.setdefault(lbl_i, dict_with_default_0())[
                             lbl_j] += len(j_sym_group)
  # self-consistency check
  multiplicities = {}
  for sc in structure.scatterers():
    multiplicities[sc.label] = sc.multiplicity()
  assert sorted(p1_counts.keys()) == sorted(sg_counts.keys())
  for lbl_i,sg_lc in sg_counts.items():
    p1_lc = p1_counts[lbl_i]
    assert sorted(p1_lc.keys()) == sorted(sg_lc.keys())
    for lbl_j,sg_c in sg_lc.items():
      p1_c = p1_lc[lbl_j]
      assert p1_c == sg_c * multiplicities[lbl_i]
  # compare with reference
  apatanl = str(sg_apat.as_nested_lists()).replace(" ","")
  if (reference_apatanl is not None):
    assert apatanl == reference_apatanl
  #
  counts = []
  for conserve_angles in [False, True]:
    proxies = structure.conservative_pair_proxies(
      bond_sym_table=sg_pat.extract_pair_sym_table(),
      conserve_angles=conserve_angles)
    counts.extend([proxies.bond.simple.size(), proxies.bond.asu.size()])
    if (not conserve_angles):
      assert proxies.angle is None
    else:
      counts.extend([proxies.angle.simple.size(), proxies.angle.asu.size()])
  cppc = ",".join([str(c) for c in counts])
  if (reference_cppc is not None):
    assert cppc == reference_cppc

def exercise_all():
  verbose = "--verbose" in sys.argv[1:]
  exercise_icosahedron(verbose=verbose)
  default_distance_cutoff = 3.5
  regression_misc = libtbx.env.find_in_repositories("phenix_regression/misc")
  if (regression_misc is None):
    print "Skipping exercise_all(): phenix_regression/misc not available"
    return
  def get_reference_dict(file_name):
    path = os.path.join(regression_misc, file_name)
    if (not os.path.isfile(path)):
      print "Skipping some tests: reference file not available:", path
      return None
    result = {}
    for line in open(path).read().splitlines():
      tag, data = line.split()
      assert not tag in result
      result[tag] = data
    return result
  reference_apatanl_dict = get_reference_dict(
    "angle_pair_asu_tables_as_nested_lists")
  reference_cppc_dict = get_reference_dict(
    "conservative_pair_proxies_counts")
  file_names = []
  for file_name in ["strudat_zeolite_atlas", "strudat_special_bonds"]:
    path = os.path.join(regression_misc, file_name)
    if (not os.path.isfile(path)):
      print "Skipping %s test: input file not available" % file_name
    else:
      file_names.append(path)
  for file_name in file_names:
    strudat_entries = strudat.read_all_entries(open(file_name))
    for i_entry,entry in enumerate(strudat_entries.entries):
      if (    file_name.endswith("strudat_zeolite_atlas")
          and not ("--full" in sys.argv[1:] or i_entry % 20 == 0)):
        continue
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
      if (reference_apatanl_dict is None):
        reference_apatanl = None
      else:
        assert entry.tag in reference_apatanl_dict
        reference_apatanl = reference_apatanl_dict[entry.tag]
      if (reference_cppc_dict is None):
        reference_cppc = None
      else:
        assert entry.tag in reference_cppc_dict
        reference_cppc = reference_cppc_dict[entry.tag]
      exercise_angle_pair_asu_table(
        structure=structure,
        distance_cutoff=distance_cutoff,
        connectivities=connectivities,
        reference_apatanl=reference_apatanl,
        reference_cppc=reference_cppc)

def run():
  exercise_all()
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
