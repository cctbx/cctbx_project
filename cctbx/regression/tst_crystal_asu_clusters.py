from __future__ import absolute_import, division, print_function
from cctbx import crystal
import cctbx.crystal.direct_space_asu
from cctbx.array_family import flex
from six.moves import range

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
  have_one_unsorted = False
  for i_trial in range(10):
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
      clusters = crystal.asu_clusters(
        pair_asu_table=pair_asu_table).sort_index_groups_by_size()
      assert [cluster.size() for cluster in clusters.index_groups] \
          == expected_cluster_sizes
      if (distance_cutoff == 1.5):
        for cluster in clusters.index_groups:
          sorted_indices = list(cluster)
          sorted_indices.sort()
          if (list(cluster) != sorted_indices):
            have_one_unsorted = True
        clusters.sort_indices_in_each_group()
        for cluster in clusters.index_groups:
          sorted_indices = list(cluster)
          sorted_indices.sort()
          assert list(cluster) == sorted_indices
  assert have_one_unsorted

def exercise_crystallographic():
  crystal_symmetry = crystal.symmetry(
    unit_cell=(10, 10, 10, 90, 90, 90),
    space_group_symbol="P 1 1 2")
  sites_frac = flex.vec3_double([
    (0.1, 0.1, 0.0),
    (0.9, 0.1, 0.0)])
  for distance_cutoff in [1,2]:
    pair_asu_table = \
      crystal_symmetry.special_position_settings().pair_asu_table(
        distance_cutoff=distance_cutoff,
        sites_frac=sites_frac)
    for strictly_in_asu in [True, False]:
      cluster = crystal.asu_clusters(
        pair_asu_table=pair_asu_table,
        strictly_in_asu=strictly_in_asu).sort_index_groups_by_size()
      cluster_sizes = [cluster.size() for cluster in cluster.index_groups]
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
  for i_trial in range(10):
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
        cluster = crystal.asu_clusters(
          pair_asu_table=pair_asu_table,
          strictly_in_asu=strictly_in_asu).sort_index_groups_by_size()
        cluster_sizes = [cluster.size() for cluster in cluster.index_groups]
        if (distance_cutoff == 1.5 or strictly_in_asu):
          assert cluster_sizes == [3, 2]
        else:
          assert cluster_sizes == [5]

def exercise():
  exercise_non_crystallographic()
  exercise_crystallographic()
  print("OK")

if (__name__ == "__main__"):
  exercise()
