from cctbx import restraints
from cctbx import crystal
from cctbx import sgtbx
from scitbx.python_utils.misc import adopt_init_args

class bond_record:

  def __init__(self, i_seqs, rt_mx, distance_ideal, weight):
    adopt_init_args(self, locals())

def from_bond_proxies(sorted_proxies):
  result = []
  rt_mx_ji = sgtbx.rt_mx(1,1)
  for proxy in sorted_proxies.proxies:
    result.append(bond_record(
      i_seqs=proxy.i_seqs,
      rt_mx=rt_mx_ji,
      distance_ideal=proxy.distance_ideal,
      weight=proxy.weight))
    result.append(bond_record(
      i_seqs=[proxy.i_seqs[1], proxy.i_seqs[0]],
      rt_mx=rt_mx_ji,
      distance_ideal=proxy.distance_ideal,
      weight=proxy.weight))
  asu_mappings = sorted_proxies.asu_mappings()
  for proxy in sorted_proxies.sym_proxies:
    p = proxy.pair
    rt_mx_i_inverse = asu_mappings.get_rt_mx(i_seq=p.i_seq,i_sym=0).inverse()
    rt_mx_j = asu_mappings.get_rt_mx(i_seq=p.j_seq, i_sym=p.j_sym)
    rt_mx_ji = rt_mx_i_inverse.multiply(rt_mx_j)
    result.append(bond_record(
      i_seqs=[p.i_seq, p.j_seq],
      rt_mx=rt_mx_ji,
      distance_ideal=proxy.distance_ideal,
      weight=proxy.weight))
  return result

def as_bond_proxies(sorted_proxies, records):
  asu_mappings = sorted_proxies.asu_mappings()
  for record in records:
    rt_mx_i = asu_mappings.get_rt_mx(record.i_seqs[0], 0)
    i_sym = asu_mappings.find_i_sym(
      i_seq=record.i_seqs[1],
      rt_mx=rt_mx_i.multiply(record.rt_mx))
    assert i_sym >= 0
    i_seq, j_seq, j_sym = record.i_seqs[0], record.i_seqs[1], i_sym
    if (i_seq < j_seq or j_sym != 0):
      pair = asu_mappings.make_pair(i_seq=i_seq, j_seq=j_seq, j_sym=j_sym)
      proxy = restraints.bond_sym_proxy(
        pair=pair,
        distance_ideal=record.distance_ideal,
        weight=record.weight)
    sorted_proxies.process(proxy=proxy)

def check_bond_proxies(structure, sites_cart, sorted_proxies_a,
                                              sorted_proxies_b):
  scatterers = structure.scatterers()
  for proxy_pair in zip(sorted_proxies_a.proxies, sorted_proxies_b.proxies):
    for proxy in proxy_pair:
      print "check bond proxy:",
      for i_seq in proxy.i_seqs:
        print "%s(%d)" % (scatterers[i_seq].label, i_seq),
      print "w=%.6g" % proxy.weight,
      r = restraints.bond(
            sites_cart=sites_cart,
            proxy=proxy)
      print "ideal,model=%.6g %.6g %.6g" % (
        proxy.distance_ideal, r.distance_model, r.residual())
    print
  assert len(sorted_proxies_a.sym_proxies) \
      == len(sorted_proxies_b.sym_proxies)
  for i_proxy in xrange(len(sorted_proxies_a.sym_proxies)):
    for sorted_proxies in (sorted_proxies_a, sorted_proxies_b):
      proxy = sorted_proxies.sym_proxies[i_proxy]
      print "check bond proxy:",
      pair = proxy.pair
      for i_seq in [pair.i_seq, pair.j_seq]:
        print "%s(%d)" % (scatterers[i_seq].label, i_seq),
      print "j_sym=%d" % pair.j_sym,
      print "w=%.6g" % proxy.weight,
      r = restraints.bond(
            sites_cart=sites_cart,
            asu_mappings=sorted_proxies.asu_mappings(),
            proxy=proxy)
      print "ideal,model=%.6g %.6g %.6g" % (
        proxy.distance_ideal, r.distance_model, r.residual())
    print

class recompute_interaction_proxies:

  def __init__(self, structure,
                     original_asu_mappings,
                     bond_records,
                     sites_cart,
                     nonbonded_cutoff=None,
                     nonbonded_cutoff_tolerance=1.e-4):
    structure.scatterers().set_sites(
      structure.unit_cell().fractionalization_matrix() * sites_cart)
    nonbonded_cutoff_plus = nonbonded_cutoff * (1+nonbonded_cutoff_tolerance)
    asu_mappings = crystal.direct_space_asu.asu_mappings(
      space_group=original_asu_mappings.space_group(),
      asu=original_asu_mappings.asu(),
      buffer_thickness=nonbonded_cutoff_plus,
      min_distance_sym_equiv=original_asu_mappings.min_distance_sym_equiv())
    asu_mappings.process_sites_cart(original_sites=sites_cart)
    self.bond_sorted_proxies = restraints.bond_sorted_proxies(
      asu_mappings=asu_mappings)
    assert original_asu_mappings.mappings().size() \
                 == asu_mappings.mappings().size()
    as_bond_proxies(
      sorted_proxies=self.bond_sorted_proxies,
      records=bond_records)
    repulsion_n_shells = 3
    from cctbx.crystal import distance_ls
    term_table, bond_registries = distance_ls.coordination_sequences_sorted(
      structure=structure,
      sorted_proxies=self.bond_sorted_proxies,
      n_shells=repulsion_n_shells,
      heterogeneous_pairs_only=0001)
    pair_generator = crystal.neighbors_fast_pair_generator(
      asu_mappings=asu_mappings,
      distance_cutoff=nonbonded_cutoff)
    self.repulsion_proxies = distance_ls.start_repulsion_proxies(
      bond_registries=bond_registries,
      pair_generator=pair_generator,
      vdw_radius=-1,
      vdw_radius_1_4=-2)
    distance_ls.set_vdw_radii(
      structure=structure,
      repulsion_proxies=self.repulsion_proxies)
    print "repulsion_proxies n_total:", \
          self.repulsion_proxies.sorted_proxies.n_total()
