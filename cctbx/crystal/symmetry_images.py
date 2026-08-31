"""Symmetry images lying within a radius of selected sites."""
from __future__ import absolute_import, division, print_function
from six.moves import range
import collections

# The two tables `images_near_seeds` reads.  Neither records which sites it
# was built from; one paired with a table built from other sites can give
# wrong images.
tables = collections.namedtuple(
  "tables", ["pair_sym_table", "site_symmetry_table"])

# A site counts as lying on a symmetry element only when its symmetry mate is
# within this distance.  Coordinates written to three decimals put a site on
# an element within 0.002 A of its mate.  A site taken to lie on an element is
# moved onto it, its own image then coincides with it and is dropped, and its
# neighbours' images are measured from the moved position.  The move is not
# bounded by this distance, since the site symmetry is re-derived from the
# moved site until it settles: up to 0.023 A where the site symmetry has order
# 48, and 0.005 A where it has order 2.
coincident_distance = 1.e-2

def tables_within_radius(sites_cart, crystal_symmetry, radius):
  """The tables `images_near_seeds` reads, covering every contact within
  *radius*.

  Neither depends on the seed selection.

  Parameters
  ----------
  sites_cart : scitbx.array_family.flex.vec3_double
      Cartesian coordinates of the atoms.
  crystal_symmetry : cctbx.crystal.symmetry
      Cell and space group.
  radius : float
      Distance threshold in Angstrom.

  Returns
  -------
  tables
      ``pair_sym_table``, each pair listed from both ends with every
      symmetry-equivalent interaction, and the ``site_symmetry_table`` it was
      built with.
  """
  # An empty crystal_symmetry is truthy, so test the members.
  assert crystal_symmetry is not None, "no crystal_symmetry"
  assert crystal_symmetry.unit_cell() is not None, "no unit cell"
  assert crystal_symmetry.space_group() is not None, "no space group"
  assert 0 <= radius < float("inf"), f"radius is {radius}"
  special_position_settings = crystal_symmetry.special_position_settings(
    min_distance_sym_equiv=coincident_distance)
  site_symmetry_table = special_position_settings.site_symmetry_table(
    sites_cart=sites_cart)
  # add_all_pairs searches to distance_cutoff*(1+1e-6); the buffer matches.
  pair_asu_table = special_position_settings.pair_asu_table(
    distance_cutoff               = radius,
    sites_cart                    = sites_cart,
    site_symmetry_table           = site_symmetry_table,
    asu_mappings_buffer_thickness = radius * (1 + 1.e-6))
  # skip_j_seq_less_than_i_seq=False lists each pair from both ends.
  # all_interactions_from_inside_asu=True keeps every operator of a group of
  # symmetry-equivalent interactions; the default keeps one of each group.
  return tables(
    pair_sym_table = pair_asu_table.extract_pair_sym_table(
      skip_j_seq_less_than_i_seq       = False,
      all_interactions_from_inside_asu = True),
    site_symmetry_table = site_symmetry_table)

def images_near_seeds(symmetry_tables, seed_indices):
  """Symmetry images in *symmetry_tables* that lie near a seed.

  The radius is the one the tables were built with.

  Parameters
  ----------
  symmetry_tables : tables
      As returned by `tables_within_radius`.
  seed_indices : iterable of int
      Indices into the sites the tables were built from.

  Returns
  -------
  dict
      ``{j_seq: [rt_mx_ji, ...]}``.  Applying ``rt_mx_ji`` to the fractional
      coordinates of site ``j_seq`` places that image near a seed.  Operators
      carry denominators (1, 12), so equal operators hash equally.

  Notes
  -----
  The result identifies images, not pairs: which seed an image was near is not
  recorded.

  An image is reported once however many operators produce it, and an image
  that coincides with the site it came from is not reported.

  Both hold exactly for a site on a symmetry element.  For a site moved onto
  one, the images the move brings together are merged, at most twice the move
  apart, and an image that close to the radius can be reported or left out.
  """
  matrices_by_j_seq = {}
  by_j_seq = {}
  pair_sym_table = symmetry_tables.pair_sym_table
  site_symmetry_table = symmetry_tables.site_symmetry_table
  assert site_symmetry_table.indices().size() == pair_sym_table.size(), (
    "the two tables cover different numbers of sites")
  for i_seq in set(seed_indices):
    assert i_seq >= 0, f"seed index is {i_seq}"
    for j_seq, rt_mx_ji_list in pair_sym_table[i_seq].items():
      by_key = None
      # stl_vector_rt_mx has no __iter__, so a for-in loop falls back to the
      # __getitem__ protocol and costs one C++ exception per vector.
      for i_op in range(len(rt_mx_ji_list)):
        rt_mx_ji = rt_mx_ji_list[i_op]
        # An operator of site j's own symmetry is recorded in the table as the
        # identity, so this drops those too.
        if (rt_mx_ji.is_unit_mx()): continue
        if (by_key is None):
          site_symmetry_matrices = matrices_by_j_seq.get(j_seq)
          if (site_symmetry_matrices is None):
            site_symmetry_matrices = site_symmetry_table.get(j_seq).matrices()
            matrices_by_j_seq[j_seq] = site_symmetry_matrices
          one_matrix = len(site_symmetry_matrices) == 1
          by_key = by_j_seq.get(j_seq)
          if (by_key is None):
            by_key = by_j_seq[j_seq] = {}
        if (one_matrix):
          key = rt_mx_ji.as_xyz()
          if (key not in by_key):
            by_key[key] = rt_mx_ji.new_denominators(1, 12)
        else:
          coset = {}
          for matrix in site_symmetry_matrices:
            image = rt_mx_ji.multiply(matrix)
            coset[image.as_xyz()] = image
          key = min(coset)
          if (key not in by_key):
            # The canonical member of the coset, not whichever was seen
            # first, so a union of seeds gives the union of the results.
            by_key[key] = coset[key].new_denominators(1, 12)
  result = {}
  for j_seq, by_key in by_j_seq.items():
    result[j_seq] = list(by_key.values())
  return result

def images_within_radius(sites_cart, crystal_symmetry, radius, seed_indices):
  """Symmetry images within *radius* of the seed atoms.

  `tables_within_radius` followed by `images_near_seeds`, which document the
  arguments and the result.
  """
  return images_near_seeds(
    tables_within_radius(sites_cart, crystal_symmetry, radius), seed_indices)
