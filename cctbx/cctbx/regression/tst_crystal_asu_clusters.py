from cctbx import crystal
import cctbx.crystal.direct_space_asu
from cctbx.array_family import flex

def exercise_non_crystallographic():
  sites_cart = flex.vec3_double([
    (7.767, 5.853, 7.671),
    (6.935, 5.032, 8.622),
    (5.918, 4.176, 8.140),
    (7.161, 5.107, 10.012),
    (5.126, 3.395, 9.038),
    (6.382, 4.336, 10.930),
    (5.360, 3.476, 10.439),
    (7.956, 7.811, 6.133),
    (8.506, 7.237, 5.169),
    (8.143, 9.010, 6.428),
    (6.253, 5.840, 5.439),
    (5.364, 7.253, 5.745),
    (5.875, 6.461, 6.183),
    (5.216, 5.927, 6.782),
    (7.000, 7.000, 7.000)])
  for i_trial in xrange(10):
    if (i_trial > 0):
      sites_cart = sites_cart.select(
        flex.random_permutation(size=sites_cart.size()))
    asu_mappings = crystal.direct_space_asu.non_crystallographic_asu_mappings(
      sites_cart=sites_cart)
    for distance_cutoff, expected_cluster_sizes in [
          (1.0, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
          (1.1, [4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
          (1.5, [6, 5, 3, 1]),
          (1.6, [15])]:
      pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
      pair_asu_table.add_all_pairs(distance_cutoff=distance_cutoff)
      clusters = crystal.asu_clusters(pair_asu_table).sort().clusters
      assert [cluster.size() for cluster in clusters] == expected_cluster_sizes

def exercise_crystallographic():
  crystal_symmetry = crystal.symmetry(
    unit_cell=(10, 10, 10, 90, 90, 90),
    space_group_symbol="P 1 1 2")
  sites_frac = flex.vec3_double([
    (0.1, 0.1, 0.0),
    (0.9, 0.1, 0.0)])
  for distance_cutoff in [1,2]:
    asu_mappings = crystal_symmetry.asu_mappings(
      buffer_thickness=distance_cutoff).process_sites_frac(
        original_sites=sites_frac)
    pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
    pair_asu_table.add_all_pairs(distance_cutoff=distance_cutoff)
    for strictly_in_asu in [True, False]:
      clusters = crystal.asu_clusters(
        pair_asu_table=pair_asu_table,
        strictly_in_asu=strictly_in_asu).sort().clusters
      cluster_sizes = [cluster.size() for cluster in clusters]
      if (distance_cutoff == 1 or strictly_in_asu):
        assert cluster_sizes == [1, 1]
      else:
        assert cluster_sizes == [2]
  sites_frac = flex.vec3_double([
    (0.1, 0.1, 0.0),
    (0.2, 0.2, 0.0),
    (0.1, 0.3, 0.0),
    (0.9, 0.1, 0.0),
    (0.8, 0.2, 0.0)])
  for i_trial in xrange(10):
    if (i_trial > 0):
      sites_frac = sites_frac.select(
        flex.random_permutation(size=sites_frac.size()))
    for distance_cutoff in [1.5,2]:
      asu_mappings = crystal_symmetry.asu_mappings(
        buffer_thickness=distance_cutoff).process_sites_frac(
          original_sites=sites_frac)
      pair_asu_table = crystal.pair_asu_table(
        asu_mappings=asu_mappings).add_all_pairs(
          distance_cutoff=distance_cutoff)
      for strictly_in_asu in [True, False]:
        clusters = crystal.asu_clusters(
          pair_asu_table=pair_asu_table,
          strictly_in_asu=strictly_in_asu).sort().clusters
        cluster_sizes = [cluster.size() for cluster in clusters]
        if (distance_cutoff == 1.5 or strictly_in_asu):
          assert cluster_sizes == [3, 2]
        else:
          assert cluster_sizes == [5]

def exercise():
  exercise_non_crystallographic()
  exercise_crystallographic()
  print "OK"

if (__name__ == "__main__"):
  exercise()
