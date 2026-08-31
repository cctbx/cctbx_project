from __future__ import absolute_import, division, print_function
from cctbx import crystal
from cctbx import sgtbx
from cctbx.crystal import symmetry_images
from scitbx.array_family import flex
from six.moves import range

# Space groups exercised against the brute-force enumeration, each with a cell
# its lattice allows.
space_groups_and_cells = [
  ("P 1",        (12.1, 14.3, 16.7, 90, 90,  90)),
  ("P -1",       (12.1, 14.3, 16.7, 90, 90,  90)),
  ("P 21",       (12.1, 14.3, 16.7, 90, 90,  90)),
  ("P 1 21/c 1", (12.1, 14.3, 16.7, 90, 105, 90)),
  ("C 1 2 1",    (12.1, 14.3, 16.7, 90, 90,  90)),
  ("P 21 21 21", (12.1, 14.3, 16.7, 90, 90,  90)),
  ("I 41",       (12.1, 12.1, 16.7, 90, 90,  90)),
  ("P 43 21 2",  (12.1, 12.1, 16.7, 90, 90,  90)),
  ("R 3 :H",     (12.1, 12.1, 16.7, 90, 90, 120)),
  ("P 63",       (12.1, 12.1, 16.7, 90, 90, 120)),
]

def general_position_structure(
      space_group_symbol="P 21 21 21",
      unit_cell=(12.1, 14.3, 16.7, 90, 90, 90),
      n_sites=12):
  """A structure with all sites in general positions and fixed coordinates."""
  cs = crystal.symmetry(unit_cell=unit_cell,
                        space_group_symbol=space_group_symbol)
  sites_frac = flex.vec3_double()
  for i in range(n_sites):
    sites_frac.append((
      (0.0417 + 0.0713 * i) % 0.83,
      (0.1103 + 0.0891 * i) % 0.79,
      (0.0629 + 0.1097 * i) % 0.87))
  return cs, cs.unit_cell().orthogonalize(sites_frac)

def brute_force_images(sites_cart, crystal_symmetry, radius, seed_indices,
                       n_shifts=2):
  """Images within *radius* of a seed, from all operators over a block of
  lattice translations.

  Uses no asu_mappings, so it does not repeat the library path.  Valid only
  for sites in general positions: an operator in a site's own symmetry is
  recorded in the table as the identity.
  """
  unit_cell = crystal_symmetry.unit_cell()
  sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
  seed_fracs = [sites_frac[i_seq] for i_seq in seed_indices]
  identity_as_xyz = sgtbx.rt_mx().as_xyz()
  result = {}
  for rotation_op in crystal_symmetry.space_group().all_ops():
    for i in range(-n_shifts, n_shifts + 1):
      for j in range(-n_shifts, n_shifts + 1):
        for k in range(-n_shifts, n_shifts + 1):
          op = sgtbx.rt_mx(
            rotation_op.r(),
            rotation_op.t().plus(sgtbx.tr_vec((i, j, k), 1)))
          as_xyz = op.as_xyz()
          if (as_xyz == identity_as_xyz): continue
          for j_seq in range(sites_frac.size()):
            image_frac = op * sites_frac[j_seq]
            for seed_frac in seed_fracs:
              if (unit_cell.distance(image_frac, seed_frac) <= radius):
                result.setdefault(j_seq, set()).add(as_xyz)
                break
  return result

def as_xyz_sets(images):
  """{j_seq: set of operator strings}."""
  result = {}
  for j_seq, ops in images.items():
    result[j_seq] = set([op.as_xyz() for op in ops])
  return result

def image_positions(unit_cell, sites_frac, images):
  """{j_seq: set of Cartesian image positions, rounded to 4 decimals}."""
  result = {}
  for j_seq, ops in images.items():
    result[j_seq] = set([
      tuple(round(x, 4) + 0.0
            for x in unit_cell.orthogonalize(op * sites_frac[j_seq]))
      for op in ops])
  return result

def brute_force_image_positions(sites_cart, crystal_symmetry, radius,
                                seed_indices, n_shifts=2):
  """Image positions within *radius* of a seed, from all operators over a
  block of lattice translations.

  Reports positions rather than operators, so unlike `brute_force_images` it
  holds for sites on special positions, where several operators give one
  image.  Images coincident with the site are dropped.
  """
  unit_cell = crystal_symmetry.unit_cell()
  sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
  seed_fracs = [sites_frac[i_seq] for i_seq in seed_indices]
  result = {}
  for rotation_op in crystal_symmetry.space_group().all_ops():
    for i in range(-n_shifts, n_shifts + 1):
      for j in range(-n_shifts, n_shifts + 1):
        for k in range(-n_shifts, n_shifts + 1):
          op = sgtbx.rt_mx(
            rotation_op.r(),
            rotation_op.t().plus(sgtbx.tr_vec((i, j, k), 1)))
          for j_seq in range(sites_frac.size()):
            image_frac = op * sites_frac[j_seq]
            if (unit_cell.distance(image_frac, sites_frac[j_seq]) <= 1.e-6):
              continue
            for seed_frac in seed_fracs:
              if (unit_cell.distance(image_frac, seed_frac) <= radius):
                result.setdefault(j_seq, set()).add(
                  tuple(round(x, 4) + 0.0
                        for x in unit_cell.orthogonalize(image_frac)))
                break
  return result

def differences(expected, obtained):
  """Per j_seq, what is missing from and extra in *obtained*."""
  result = {}
  for j_seq in set(expected) | set(obtained):
    missing = expected.get(j_seq, set()) - obtained.get(j_seq, set())
    extra = obtained.get(j_seq, set()) - expected.get(j_seq, set())
    if (missing or extra):
      result[j_seq] = (sorted(missing), sorted(extra))
  return result

def exercise_matches_brute_force():
  """Every image a direct search finds is found, and no others.

  The larger radii reach far enough for most of these groups to give an atom
  images of itself, of which a table extracted without
  all_interactions_from_inside_asu keeps one per symmetry-equivalent group.
  """
  seed_indices = [0, 3, 7]
  for space_group_symbol, unit_cell in space_groups_and_cells:
    cs, sites_cart = general_position_structure(
      space_group_symbol=space_group_symbol, unit_cell=unit_cell)
    n_images = 0
    for radius in [1.0, 2.0, 5.0, 8.0]:
      obtained = as_xyz_sets(symmetry_images.images_within_radius(
        sites_cart       = sites_cart,
        crystal_symmetry = cs,
        radius           = radius,
        seed_indices     = seed_indices))
      expected = brute_force_images(
        sites_cart       = sites_cart,
        crystal_symmetry = cs,
        radius           = radius,
        seed_indices     = seed_indices)
      assert obtained == expected, (
        space_group_symbol, radius, differences(expected, obtained))
      n_images += len(expected)
    # The smallest radii can find nothing; the sweep as a whole must not.
    assert n_images > 0, space_group_symbol

def exercise_brute_force_block_is_large_enough():
  """The reference enumeration covers enough lattice translations."""
  seed_indices = [0, 3, 7]
  for space_group_symbol, unit_cell in space_groups_and_cells:
    cs, sites_cart = general_position_structure(
      space_group_symbol=space_group_symbol, unit_cell=unit_cell)
    for radius in [1.0, 2.0, 5.0, 8.0]:
      near = brute_force_images(
        sites_cart, cs, radius, seed_indices, n_shifts=2)
      far = brute_force_images(
        sites_cart, cs, radius, seed_indices, n_shifts=3)
      assert near == far, (space_group_symbol, radius)

def exercise_seeding_is_a_subset():
  """Seeding on fewer atoms returns part of the whole-structure result."""
  cs, sites_cart = general_position_structure()
  radius = 5.0
  everything = as_xyz_sets(symmetry_images.images_within_radius(
    sites_cart, cs, radius, seed_indices=list(range(sites_cart.size()))))
  subset = as_xyz_sets(symmetry_images.images_within_radius(
    sites_cart, cs, radius, seed_indices=[0, 1]))
  assert len(subset) > 0
  assert len(subset) < len(everything)
  for j_seq, ops in subset.items():
    assert j_seq in everything
    assert ops.issubset(everything[j_seq])

def exercise_both_pair_directions():
  """Images are reported whether the seed is the first or second atom of the
  underlying pair.
  """
  cs, sites_cart = general_position_structure()
  radius = 5.0
  symmetry_tables = symmetry_images.tables_within_radius(
    sites_cart, cs, radius)
  n_checked = 0
  for i_seq in range(sites_cart.size()):
    forward = symmetry_images.images_near_seeds(symmetry_tables, [i_seq])
    for j_seq, ops in forward.items():
      if (j_seq == i_seq): continue
      backward = symmetry_images.images_near_seeds(
        symmetry_tables, [j_seq])
      assert i_seq in backward, (i_seq, j_seq)
      backward_as_xyz = set([op.as_xyz() for op in backward[i_seq]])
      for op in ops:
        inverse_as_xyz = op.inverse().new_denominators(1, 12).as_xyz()
        assert inverse_as_xyz in backward_as_xyz, (
          i_seq, j_seq, op.as_xyz(), sorted(backward_as_xyz))
        n_checked += 1
  assert n_checked > 0

def exercise_distances_hold():
  """A returned operator places its atom within radius of a seed.

  No site here is moved onto a symmetry element.
  """
  cs, sites_cart = general_position_structure()
  radius = 5.0
  seed_indices = [0, 5]
  sites_frac = cs.unit_cell().fractionalize(sites_cart=sites_cart)
  images = symmetry_images.images_within_radius(
    sites_cart, cs, radius, seed_indices=seed_indices)
  assert len(images) > 0
  for j_seq, ops in images.items():
    for op in ops:
      image_frac = op * sites_frac[j_seq]
      distances = [cs.unit_cell().distance(image_frac, sites_frac[i_seq])
                   for i_seq in seed_indices]
      assert min(distances) <= radius * (1 + 1.e-6), (j_seq, op.as_xyz())

def exercise_identity_excluded():
  """The identity is never returned."""
  cs, sites_cart = general_position_structure()
  images = symmetry_images.images_within_radius(
    sites_cart, cs, radius=5.0, seed_indices=list(range(sites_cart.size())))
  assert len(images) > 0
  identity_as_xyz = sgtbx.rt_mx().as_xyz()
  for j_seq, ops in images.items():
    for op in ops:
      assert op.as_xyz() != identity_as_xyz

def special_position_structure(unit_cell=(6, 7, 28, 90, 90, 90)):
  """A structure whose first site lies on the two-fold of C 1 2 1."""
  cs = crystal.symmetry(unit_cell=unit_cell, space_group_symbol="C 1 2 1")
  sites_frac = flex.vec3_double([
    (0.0,   0.2,  0.0),
    (0.13,  0.31, 0.07),
    (0.21,  0.44, 0.19)])
  return cs, sites_frac, cs.unit_cell().orthogonalize(sites_frac)

def exercise_special_position():
  """An atom on a rotation axis yields no image on top of itself."""
  cs, sites_frac, sites_cart = special_position_structure()
  images = symmetry_images.images_within_radius(
    sites_cart, cs, radius=5.0, seed_indices=[0])
  ops = images.get(0, [])
  assert len(ops) > 0
  for op in ops:
    image_frac = op * sites_frac[0]
    distance = cs.unit_cell().distance(image_frac, sites_frac[0])
    assert distance > 1.e-6, op.as_xyz()

def exercise_self_images_of_a_seed():
  """An atom's own symmetry images are reported.

  A table extracted without all_interactions_from_inside_asu keeps one image
  per group of symmetry-equivalent interactions and misses some of these.
  Site 1 is in a general position.
  """
  cs, sites_frac, sites_cart = special_position_structure()
  for seed in [0, 1]:
    images = symmetry_images.images_within_radius(
      sites_cart, cs, radius=5.0, seed_indices=[seed])
    obtained = image_positions(cs.unit_cell(), sites_frac, images)
    expected = brute_force_image_positions(sites_cart, cs, 5.0, [seed])
    # Four images for seed 0 and five for seed 1; an extraction that
    # keeps one per group returns 1 and 3.
    assert len(expected[seed]) > 1, seed
    assert obtained.get(seed) == expected[seed], (
      seed, differences(expected, obtained))

def exercise_prebuilt_table():
  """A separately built table gives the same answer as building it inline."""
  cs, sites_cart = general_position_structure()
  radius = 5.0
  seed_indices = [0, 4, 9]
  direct = symmetry_images.images_within_radius(
    sites_cart, cs, radius, seed_indices=seed_indices)
  staged = symmetry_images.images_near_seeds(
    symmetry_images.tables_within_radius(sites_cart, cs, radius), seed_indices)
  assert len(direct) > 0
  assert as_xyz_sets(direct) == as_xyz_sets(staged)

def exercise_operators_are_hash_stable():
  """Returned operators survive a set round trip.

  sgtbx.rt_mx has an inconsistent __hash__ and __eq__: operators that compare
  equal but carry different translation denominators hash differently, so an
  un-normalised operator can enter a set twice.
  """
  cs, sites_cart = general_position_structure()
  images = symmetry_images.images_within_radius(
    sites_cart, cs, radius=5.0, seed_indices=[0, 2, 6])
  assert len(images) > 0
  for j_seq, ops in images.items():
    assert len(set([op.as_xyz() for op in ops])) == len(ops)
    ops_set = set(ops)
    assert len(ops_set) == len(ops), j_seq
    for op in ops:
      assert op in ops_set
      assert sgtbx.rt_mx(op.as_xyz()).new_denominators(1, 12) in ops_set

def near_special_position_structure():
  """A structure whose first site sits 0.2 A off the two-fold of C 1 2 1.

  Its symmetry mate is 0.4 A away, far enough that the site is not taken to
  lie on the axis.
  """
  cs = crystal.symmetry(unit_cell=(20, 24, 28, 90, 90, 90),
                        space_group_symbol="C 1 2 1")
  sites_frac = flex.vec3_double([
    (0.01, 0.2,  0.0),
    (0.13, 0.31, 0.07)])
  return cs, sites_frac, cs.unit_cell().orthogonalize(sites_frac)

def exercise_a_site_near_an_element_is_general():
  """A site near a symmetry element keeps every image of its own.

  Taking it to lie on the element merges it with its mate, and the image is
  then dropped as coincident.
  """
  cs, sites_frac, sites_cart = near_special_position_structure()
  symmetry_tables = symmetry_images.tables_within_radius(
    sites_cart, cs, 5.0)
  special = list(
    symmetry_tables.site_symmetry_table.special_position_indices())
  assert special == [], special
  images = symmetry_images.images_within_radius(
    sites_cart, cs, 5.0, seed_indices=[0])
  obtained = image_positions(cs.unit_cell(), sites_frac, images)
  expected = brute_force_image_positions(sites_cart, cs, 5.0, [0])
  assert len(expected) > 0
  assert obtained == expected, differences(expected, obtained)

def exercise_a_neighbour_near_an_element_keeps_its_images():
  """An image of a site near an element is reported when it is in range.

  Taking the site to lie on the element loses the image.
  """
  cs = crystal.symmetry(unit_cell=(9.579, 10.153, 11.64, 90, 90, 90),
                        space_group_symbol="P -1")
  # Site 1 lies 0.1775 A off the inversion centre at (1/2, 0, 0).  Under
  # x,y+1,z it is 4.99 A from the seed, under -x+1,-y+1,-z it is 5.34 A.
  sites_frac = flex.vec3_double([
    (0.672825, 0.536167, 0.113296),
    (0.506371, -0.016354, 0.001298),
    (0.493881, 0.352158, 0.718093)])
  sites_cart = cs.unit_cell().orthogonalize(sites_frac)
  images = symmetry_images.images_within_radius(
    sites_cart, cs, 5.0, seed_indices=[0])
  obtained = image_positions(cs.unit_cell(), sites_frac, images)
  expected = brute_force_image_positions(sites_cart, cs, 5.0, [0])
  assert len(expected) > 0
  assert obtained == expected, differences(expected, obtained)

def exercise_result_is_lists_of_operators():
  """The keys are plain ints and the values are lists of operators."""
  cs, sites_cart = general_position_structure()
  images = symmetry_images.images_within_radius(
    sites_cart, cs, 5.0, seed_indices=[0, 3, 7])
  assert isinstance(images, dict)
  assert len(images) > 0
  for j_seq, ops in images.items():
    assert isinstance(j_seq, int), type(j_seq)
    assert isinstance(ops, list), (j_seq, type(ops))
    assert len(ops) > 0
    for op in ops:
      assert isinstance(op, sgtbx.rt_mx)

def exercise_operators_from_mixed_denominators():
  """Equal operators arriving on different denominators collapse to one.

  No library path produces this, so the table is built by hand.
  """
  table = crystal.pair_sym_table(2)
  ops = sgtbx.stl_vector_rt_mx()
  rt_mx = sgtbx.rt_mx("-x,y+1/2,-z")
  ops.append(rt_mx.new_denominators(1, 12))
  ops.append(rt_mx.new_denominators(1, 24))
  ops.append(sgtbx.rt_mx("x+1,y,z"))
  table[0][1] = ops
  cs = crystal.symmetry(unit_cell=(30, 30, 30, 90, 90, 90),
                        space_group_symbol="P 1")
  site_symmetry_table = cs.special_position_settings().site_symmetry_table(
    sites_cart=flex.vec3_double([(1, 2, 3), (4, 5, 6)]))
  images = symmetry_images.images_near_seeds(
    symmetry_images.tables(table, site_symmetry_table), [0])
  assert sorted(images) == [1], sorted(images)
  obtained = images[1]
  assert len(obtained) == 2, [op.as_xyz() for op in obtained]
  assert len(set(obtained)) == len(obtained)
  assert set([op.t().den() for op in obtained]) == set([12])

def exercise_seed_indices_contract():
  """Seeds are taken as a set; no seeds gives no images.

  A negative index is rejected rather than counted from the end.
  """
  cs, sites_cart = general_position_structure()
  radius = 5.0
  once = as_xyz_sets(symmetry_images.images_within_radius(
    sites_cart, cs, radius, seed_indices=[0, 3]))
  repeated = as_xyz_sets(symmetry_images.images_within_radius(
    sites_cart, cs, radius, seed_indices=[3, 0, 3, 0]))
  assert len(once) > 0
  assert once == repeated
  assert symmetry_images.images_within_radius(
    sites_cart, cs, radius, seed_indices=[]) == {}
  try:
    symmetry_images.images_within_radius(
      sites_cart, cs, radius, seed_indices=[-1])
  except AssertionError:
    pass
  else:
    raise AssertionError("accepted a negative seed index")

def exercise_every_seed_contributes():
  """Seeding on many atoms is the union of seeding on each of them.

  Two seeds can reach one image of an atom on a special position by different
  operators, so which operator is reported must not depend on the seeds.
  """
  cs, sites_cart = general_position_structure()
  cases = [(cs, sites_cart, list(range(sites_cart.size())))]
  for space_group_symbol, unit_cell, fracs, seeds in special_position_fixtures:
    cs = crystal.symmetry(unit_cell=unit_cell,
                          space_group_symbol=space_group_symbol)
    cases.append((cs, cs.unit_cell().orthogonalize(flex.vec3_double(fracs)),
                  seeds))
  radius = 5.0
  for cs, sites_cart, seed_indices in cases:
    together = as_xyz_sets(symmetry_images.images_within_radius(
      sites_cart, cs, radius, seed_indices=seed_indices))
    apart = {}
    for i_seq in seed_indices:
      for j_seq, ops in as_xyz_sets(symmetry_images.images_within_radius(
          sites_cart, cs, radius, seed_indices=[i_seq])).items():
        apart.setdefault(j_seq, set()).update(ops)
    assert len(together) > 0
    assert together == apart, differences(apart, together)

def exercise_requires_crystal_symmetry():
  """None and an empty crystal.symmetry are rejected."""
  cs, sites_cart = general_position_structure()
  placeholder = crystal.symmetry()
  assert placeholder.unit_cell() is None
  assert placeholder.space_group_info() is None
  for bad in [None, placeholder]:
    try:
      symmetry_images.images_within_radius(sites_cart, bad, 5.0, [0])
    except AssertionError:
      pass
    else:
      raise AssertionError(f"accepted crystal_symmetry {bad!r}")

def exercise_one_operator_per_image():
  """Each image of an atom is reported once.

  Two seeds can reach one image of a third atom under operators that differ
  by that atom's site symmetry.
  """
  cs = crystal.symmetry(unit_cell=(9.673, 12.574, 16.444, 90, 90, 90),
                        space_group_symbol="F m m 2")
  # Site 0 lies on the mirror at x = 0; the seeds are general.
  sites_frac = flex.vec3_double([
    (0.0,   1 / 6., 0.0),
    (0.113, 0.207,  0.331),
    (0.717, 0.311,  0.613)])
  sites_cart = cs.unit_cell().orthogonalize(sites_frac)
  images = symmetry_images.images_within_radius(
    sites_cart, cs, radius=5.0, seed_indices=[1, 2])
  assert len(images) > 0
  for j_seq, ops in images.items():
    positions = set()
    for op in ops:
      image_frac = op * sites_frac[j_seq]
      positions.add(tuple(round(x, 6) for x in image_frac))
    assert len(positions) == len(ops), (
      j_seq, sorted([op.as_xyz() for op in ops]))

# Fixtures whose first site lies on a symmetry element, with the seeds to use.
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
]

def near_special_position_seed_structure():
  """A structure whose first site sits 0.15 A off an inversion centre.

  The other sites are in general positions.
  """
  space_group_info = sgtbx.space_group_info("P 1 21/c 1")
  cs = crystal.symmetry(
    unit_cell=space_group_info.any_compatible_unit_cell(volume=700.),
    space_group_info=space_group_info)
  # The images of sites 1 and 2 fall just inside the radius used below.
  sites_frac = flex.vec3_double([
    (0.00632,  0.005007, 0.012803),
    (0.845806, 0.150766, 0.661288),
    (0.698725, 0.188033, 0.675395),
    (0.87933,  0.832146, 0.150299)])
  return cs, sites_frac, cs.unit_cell().orthogonalize(sites_frac)

def exercise_near_special_seed_keeps_its_neighbours():
  """A seed lying near a symmetry element still finds the images around it.

  At the library default tolerance the seed is taken to lie on the element,
  and images inside the radius are then lost.
  """
  cs, sites_frac, sites_cart = near_special_position_seed_structure()
  radius = 4.0
  images = symmetry_images.images_within_radius(
    sites_cart, cs, radius, seed_indices=[0])
  obtained = image_positions(cs.unit_cell(), sites_frac, images)
  expected = brute_force_image_positions(sites_cart, cs, radius, [0])
  assert len(expected) > 1
  assert obtained == expected, differences(expected, obtained)
  # Nothing is reported beyond the radius either.
  for j_seq, ops in images.items():
    for op in ops:
      distance = cs.unit_cell().distance(
        op * sites_frac[j_seq], sites_frac[0])
      assert distance <= radius * (1 + 1.e-6), (j_seq, op.as_xyz(), distance)

def exercise_a_radius_that_is_not_a_distance_is_rejected():
  """A radius that is not a length is rejected."""
  cs, sites_cart = general_position_structure()
  for radius in [-1.0, float("inf"), float("nan")]:
    try:
      symmetry_images.tables_within_radius(sites_cart, cs, radius)
    except AssertionError:
      pass
    else:
      raise AssertionError(f"accepted a radius of {radius}")

def exercise_tables_must_cover_the_same_sites():
  """Two tables built for different numbers of sites are rejected."""
  cs, sites_cart = general_position_structure()
  symmetry_tables = symmetry_images.tables_within_radius(sites_cart, cs, 5.0)
  fewer = symmetry_images.tables_within_radius(
    sites_cart[:4], cs, 5.0).site_symmetry_table
  try:
    symmetry_images.images_near_seeds(
      symmetry_images.tables(symmetry_tables.pair_sym_table, fewer), [0])
  except AssertionError:
    pass
  else:
    raise AssertionError("accepted tables covering different sites")

def exercise_a_site_written_on_an_element_is_coincident():
  """A site on a symmetry element is recognised at the precision it is written.

  Three decimals leave a site on the element within 0.002 A of its mate.  The
  0.03 case puts the site 0.042 A off the axis and 0.085 A from its mate, so
  it is not taken to lie on the axis.
  """
  cs = crystal.symmetry(unit_cell=(20, 24, 28, 90, 90, 90),
                        space_group_symbol="C 1 2 1")
  # The two-fold runs along (0, y, 0).
  for offset, on_the_element in [(0.0005, True), (0.03, False)]:
    sites_cart = flex.vec3_double([
      (offset, 4.8, offset), (2.6, 7.44, 1.96)])
    symmetry_tables = symmetry_images.tables_within_radius(
      sites_cart, cs, 5.0)
    special = list(
      symmetry_tables.site_symmetry_table.special_position_indices())
    assert special == ([0] if on_the_element else []), (offset, special)

def exercise_special_positions_are_complete():
  """Every image of an atom on a symmetry element is reported, exactly once.

  `brute_force_images` compares operators and so cannot be used here; this
  compares positions.
  """
  n_checked = 0
  for space_group_symbol, unit_cell, fracs, seeds in special_position_fixtures:
    cs = crystal.symmetry(unit_cell=unit_cell,
                          space_group_symbol=space_group_symbol)
    sites_frac = flex.vec3_double(fracs)
    sites_cart = cs.unit_cell().orthogonalize(sites_frac)
    for radius in [4.0, 6.0]:
      images = symmetry_images.images_within_radius(
        sites_cart, cs, radius, seed_indices=seeds)
      obtained = image_positions(cs.unit_cell(), sites_frac, images)
      expected = brute_force_image_positions(sites_cart, cs, radius, seeds)
      assert obtained == expected, (
        space_group_symbol, radius, differences(expected, obtained))
      assert len(expected) > 0, (space_group_symbol, radius)
      for j_seq, ops in images.items():
        assert len(ops) == len(obtained[j_seq]), (
          space_group_symbol, radius, j_seq)
        n_checked += len(ops)
  assert n_checked > 0

def exercise_tolerance_admits_a_three_decimal_worst_case():
  """A site written on an element at three decimals is recognised.

  Rounding moves a site up to 0.0005 A along each of three axes, so across an
  inversion centre its mate is 2*sqrt(3)*0.0005 = 0.0017 A away.  The
  tolerance has to clear that, not just the 0.0014 A of a two-axis offset.
  """
  cs = crystal.symmetry(unit_cell=(21, 23, 25, 90, 90, 90),
                        space_group_symbol="P -1")
  sites_cart = flex.vec3_double([(0.0005, 0.0005, 0.0005), (7.3, 8.1, 9.7)])
  sites_frac = cs.unit_cell().fractionalize(sites_cart)
  mate = cs.unit_cell().distance(
    sites_frac[0], tuple(-x for x in sites_frac[0]))
  assert abs(mate - 0.0017320) < 1.e-6, mate
  symmetry_tables = symmetry_images.tables_within_radius(
    sites_cart, cs, 5.0)
  special = list(
    symmetry_tables.site_symmetry_table.special_position_indices())
  assert special == [0], special

def exercise_a_contact_exactly_at_the_radius_is_reported():
  """A contact lying exactly at the radius is reported.

  A buffer below the radius makes the library raise.
  """
  n_checked = 0
  # The cell edge is the radius, so a lattice translation lands on it.
  for edge in [4.0, 6.0, 7.0, 8.0]:
    cs = crystal.symmetry(unit_cell=(edge, 7.0, 28.0, 90, 90, 90),
                          space_group_symbol="C 1 2 1")
    sites_frac = flex.vec3_double([(0.0, 0.2, 0.0), (0.31, 0.44, 0.19)])
    sites_cart = cs.unit_cell().orthogonalize(sites_frac)
    images = symmetry_images.images_within_radius(
      sites_cart, cs, edge, seed_indices=[0])
    obtained = image_positions(cs.unit_cell(), sites_frac, images)
    expected = brute_force_image_positions(sites_cart, cs, edge, [0])
    assert obtained == expected, (edge, differences(expected, obtained))
    at_the_edge = [
      op for op in images[0]
      if abs(cs.unit_cell().distance(op * sites_frac[0], sites_frac[0])
             - edge) < 1.e-9]
    assert len(at_the_edge) > 0, edge
    n_checked += len(at_the_edge)
  assert n_checked > 0

def exercise_seed_order_does_not_change_the_result():
  """The result does not depend on the order the seeds are given in."""
  cs, sites_cart = general_position_structure()
  forward = symmetry_images.images_within_radius(
    sites_cart, cs, 5.0, seed_indices=[0, 3, 7])
  reverse = symmetry_images.images_within_radius(
    sites_cart, cs, 5.0, seed_indices=[7, 3, 0])
  assert len(forward) > 0
  assert list(forward) == list(reverse)
  for j_seq in forward:
    assert ([op.as_xyz() for op in forward[j_seq]]
            == [op.as_xyz() for op in reverse[j_seq]]), j_seq

def run():
  exercise_matches_brute_force()
  exercise_brute_force_block_is_large_enough()
  exercise_seeding_is_a_subset()
  exercise_both_pair_directions()
  exercise_distances_hold()
  exercise_identity_excluded()
  exercise_special_position()
  exercise_self_images_of_a_seed()
  exercise_prebuilt_table()
  exercise_operators_are_hash_stable()
  exercise_operators_from_mixed_denominators()
  exercise_seed_indices_contract()
  exercise_requires_crystal_symmetry()
  exercise_one_operator_per_image()
  exercise_a_radius_that_is_not_a_distance_is_rejected()
  exercise_tables_must_cover_the_same_sites()
  exercise_a_site_written_on_an_element_is_coincident()
  exercise_tolerance_admits_a_three_decimal_worst_case()
  exercise_a_contact_exactly_at_the_radius_is_reported()
  exercise_seed_order_does_not_change_the_result()
  exercise_special_positions_are_complete()
  exercise_near_special_seed_keeps_its_neighbours()
  exercise_a_site_near_an_element_is_general()
  exercise_a_neighbour_near_an_element_keeps_its_images()
  exercise_result_is_lists_of_operators()
  exercise_every_seed_contributes()
  print("OK")

if (__name__ == "__main__"):
  run()
