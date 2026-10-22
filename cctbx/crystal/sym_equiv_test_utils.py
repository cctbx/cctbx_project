"""Fixtures and brute-force oracles for tests of symmetry-equivalent sites."""
from cctbx import crystal
from cctbx import sgtbx
from scitbx import matrix
from scitbx.array_family import flex

# Fixtures whose first site lies on a symmetry element, with the sites to
# select.
special_position_fixtures = [
  ("C 1 2 1", (6, 7, 28, 90, 90, 90),                 # two-fold
   [(0.0, 0.2, 0.0), (0.13, 0.31, 0.07), (0.21, 0.44, 0.19)], [0, 1]),
  ("F m m 2", (9.673, 12.574, 16.444, 90, 90, 90),    # mirror, F centred
   [(0.0, 1 / 6., 0.0), (0.113, 0.207, 0.331), (0.717, 0.311, 0.613)], [1, 2]),
  ("R 3 :H", (12.1, 12.1, 16.7, 90, 90, 120),         # three-fold
   [(0.0, 0.0, 0.11), (0.213, 0.331, 0.07), (0.41, 0.62, 0.28)], [1, 2]),
  ("P 43 21 2", (12.1, 12.1, 16.7, 90, 90, 90),       # two-fold at (x, x, 0)
   [(0.19, 0.19, 0.0), (0.113, 0.407, 0.231), (0.617, 0.311, 0.513)], [1, 2]),
  ("P -1", (7, 8, 9, 90, 90, 90),                     # inversion centre
   [(0.0, 0.0, 0.0), (0.21, 0.33, 0.14), (0.44, 0.19, 0.37)], [1, 2]),
  ("P m -3 m", (9, 9, 9, 90, 90, 90),                 # site symmetry m-3m
   [(0.0, 0.0, 0.0), (0.213, 0.331, 0.07), (0.41, 0.28, 0.62)], [1, 2]),
]

def general_position_structure(space_group_info=None, n_sites=12):
  """A structure with all sites in general positions and fixed coordinates."""
  if (space_group_info is None):
    space_group_info = sgtbx.space_group_info("P 21 21 21")
  cs = crystal.symmetry(
    unit_cell=space_group_info.any_compatible_unit_cell(volume=3000.),
    space_group_info=space_group_info)
  sites_frac = flex.vec3_double()
  for i in range(n_sites):
    sites_frac.append((
      (0.0417 + 0.0713 * i) % 0.83,
      (0.1103 + 0.0891 * i) % 0.79,
      (0.0629 + 0.1097 * i) % 0.87))
  return cs, cs.unit_cell().orthogonalize(sites_frac)

def brute_force_sym_equiv(sites_cart, crystal_symmetry, radius,
                          primary_selection, n_shifts=2):
  """Equivalents within *radius* of a selected site, from all operators over
  a block of lattice translations, without asu_mappings.

  Compares operators, so it is valid only for sites in general positions.
  """
  unit_cell = crystal_symmetry.unit_cell()
  sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
  selected_fracs = [sites_frac[i_seq] for i_seq in primary_selection]
  identity_as_xyz = sgtbx.rt_mx().as_xyz()
  result = {}
  for rotation_op in crystal_symmetry.space_group().all_ops():
    for i in range(-n_shifts, n_shifts + 1):
      for j in range(-n_shifts, n_shifts + 1):
        for k in range(-n_shifts, n_shifts + 1):
          op = rotation_op + (i, j, k)
          as_xyz = op.as_xyz()
          if (as_xyz == identity_as_xyz): continue
          for j_seq in range(sites_frac.size()):
            equiv_frac = op * sites_frac[j_seq]
            for selected_frac in selected_fracs:
              if (unit_cell.distance(equiv_frac, selected_frac) <= radius):
                if (j_seq not in result):
                  result[j_seq] = set()
                result[j_seq].add(as_xyz)
                break
  return result

def unique_positions(positions, tolerance=1.e-4):
  """Positions, with duplicates closer than *tolerance* dropped.

  Compares against every position kept: two positions can share a coordinate
  to within rounding, so sorting does not bring duplicates together.
  """
  result = []
  for position in positions:
    is_duplicate = False
    for kept in result:
      if ((matrix.col(position) - matrix.col(kept)).length() < tolerance):
        is_duplicate = True
        break
    if (not is_duplicate):
      result.append(position)
  return result

def brute_force_sym_equiv_positions(sites_cart, crystal_symmetry, radius,
                                primary_selection, n_shifts=2):
  """Equivalent positions within *radius* of a selected site, from all
  operators over a block of lattice translations, without asu_mappings.

  Positions rather than operators, so this holds for sites on special
  positions too.  Equivalents coincident with their site are dropped.
  """
  unit_cell = crystal_symmetry.unit_cell()
  sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
  selected_fracs = [sites_frac[i_seq] for i_seq in primary_selection]
  result = {}
  for rotation_op in crystal_symmetry.space_group().all_ops():
    for i in range(-n_shifts, n_shifts + 1):
      for j in range(-n_shifts, n_shifts + 1):
        for k in range(-n_shifts, n_shifts + 1):
          op = rotation_op + (i, j, k)
          for j_seq in range(sites_frac.size()):
            equiv_frac = op * sites_frac[j_seq]
            if (unit_cell.distance(equiv_frac, sites_frac[j_seq]) <= 1.e-6):
              continue
            for selected_frac in selected_fracs:
              if (unit_cell.distance(equiv_frac, selected_frac) <= radius):
                if (j_seq not in result):
                  result[j_seq] = []
                result[j_seq].append(unit_cell.orthogonalize(equiv_frac))
                break
  for j_seq in result:
    result[j_seq] = unique_positions(result[j_seq])
  return result

def assert_same_positions(expected, obtained, info, tolerance=1.e-4):
  """The two position tables list the same equivalents, in any order.

  Two positions can share a coordinate to within rounding, so their order
  after sorting is not the same on every platform; each position is matched
  to one of its own instead.
  """
  assert sorted(obtained) == sorted(expected), (
    info, sorted(expected), sorted(obtained))
  for j_seq, positions in expected.items():
    unmatched = list(obtained[j_seq])
    assert len(unmatched) == len(positions), (
      info, j_seq, len(positions), len(unmatched))
    for position in positions:
      i_match = None
      for i_candidate, candidate in enumerate(unmatched):
        if ((matrix.col(position) - matrix.col(candidate)).length()
            < tolerance):
          i_match = i_candidate
          break
      if (i_match is None):
        raise AssertionError((info, j_seq, position, obtained[j_seq]))
      del unmatched[i_match]
