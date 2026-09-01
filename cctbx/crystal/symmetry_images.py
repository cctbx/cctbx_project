"""Symmetry images lying within a radius of selected sites."""
from __future__ import absolute_import, division, print_function
import collections

# Both tables must be built from the same sites; neither records which, so a
# mismatched pair gives wrong images.
tables = collections.namedtuple(
  "tables", ["pair_sym_table", "site_symmetry_table"])

# A site whose symmetry mate is within this distance is taken to lie on a
# symmetry element and is moved onto it.  Coordinates written to three
# decimals leave such a site within 0.002 A of its mate.
coincident_distance = 1.e-2

def tables_within_radius(sites_cart, crystal_symmetry, radius):
  """Pair and site symmetry tables covering every contact within *radius*.

  Neither depends on the seed selection.

  Parameters
  ----------
  sites_cart : flex.vec3_double
  crystal_symmetry : cctbx.crystal.symmetry
  radius : float
      In Angstrom.

  Returns
  -------
  tables
      Each pair listed from both ends, with every symmetry-equivalent
      interaction.
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
  # add_all_pairs searches to distance_cutoff*(1+is_inside_epsilon()), which
  # defaults to 1e-6; the buffer matches.
  pair_asu_table = special_position_settings.pair_asu_table(
    distance_cutoff               = radius,
    sites_cart                    = sites_cart,
    site_symmetry_table           = site_symmetry_table,
    asu_mappings_buffer_thickness = radius * (1 + 1.e-6))
  # all_interactions_from_inside_asu keeps every operator of a group of
  # symmetry-equivalent interactions.
  return tables(
    pair_sym_table = pair_asu_table.extract_pair_sym_table(
      skip_j_seq_less_than_i_seq       = False,
      all_interactions_from_inside_asu = True),
    site_symmetry_table = site_symmetry_table)

def images_near_seeds(symmetry_tables, seed_indices):
  """Symmetry images near the sites *seed_indices* selects.

  Parameters
  ----------
  symmetry_tables : tables
  seed_indices : iterable of int
      Indices into the sites the tables were built from.

  Returns
  -------
  dict
      ``{j_seq: [rt_mx_ji, ...]}``.  Applying *rt_mx_ji* to the fractional
      coordinates of site *j_seq* places that image within the tables'
      radius of a seed; which seed is not recorded.  Operators carry
      denominators (1, 12), so equal operators hash equally.

  Notes
  -----
  One operator per image; an image coincident with its site is dropped.  For
  a moved site the merge is approximate, so an image near the radius can be
  reported or left out.
  """
  matrices_by_j_seq = {}
  by_j_seq = {}
  pair_sym_table = symmetry_tables.pair_sym_table
  site_symmetry_table = symmetry_tables.site_symmetry_table
  assert site_symmetry_table.indices().size() == pair_sym_table.size(), (
    "the two tables cover different numbers of sites")
  for i_seq in sorted(set(seed_indices)):
    assert 0 <= i_seq < pair_sym_table.size(), f"seed index is {i_seq}"
    for j_seq, rt_mx_ji_list in pair_sym_table[i_seq].items():
      by_key = None
      # stl_vector_rt_mx has no __iter__, so a for-in loop falls back to the
      # __getitem__ protocol and costs one C++ exception per vector.
      for i_op in range(len(rt_mx_ji_list)):
        rt_mx_ji = rt_mx_ji_list[i_op]
        # Site j's own symmetry operators are recorded as the identity, so
        # this drops them as well.
        if (rt_mx_ji.is_unit_mx()): continue
        if (by_key is None):
          # matrices() builds a fresh tuple on each call, so cache it; None
          # records a general position.
          if (j_seq in matrices_by_j_seq):
            site_symmetry_matrices = matrices_by_j_seq[j_seq]
          elif (site_symmetry_table.is_special_position(j_seq)):
            site_symmetry_matrices = site_symmetry_table.get(j_seq).matrices()
            matrices_by_j_seq[j_seq] = site_symmetry_matrices
          else:
            site_symmetry_matrices = matrices_by_j_seq[j_seq] = None
          by_key = by_j_seq.get(j_seq)
          if (by_key is None):
            by_key = by_j_seq[j_seq] = {}
        if (site_symmetry_matrices is not None):
          # The canonical coset member, so a union of seeds gives the union
          # of the results.  min() orders by rt_mx.__lt__, the ordering of
          # space_group.make_tidy().
          rt_mx_ji = min([rt_mx_ji.multiply(matrix)
                          for matrix in site_symmetry_matrices])
        # Equal operators hash equally only on equal denominators, so the
        # normalised operator is the key.
        by_key[rt_mx_ji.new_denominators(1, 12)] = None
  result = {}
  for j_seq, by_key in by_j_seq.items():
    result[j_seq] = list(by_key)
  return result

def images_within_radius(sites_cart, crystal_symmetry, radius, seed_indices):
  """Symmetry images within *radius* of the seeds."""
  return images_near_seeds(
    tables_within_radius(sites_cart, crystal_symmetry, radius), seed_indices)
