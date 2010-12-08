from cctbx.eltbx.distance_based_connectivity import expected_bond_lengths_by_element_pair
from cctbx.eltbx.van_der_waals_radii import vdw
from scitbx.stl import map
import cctbx.crystal
import cctbx.uctbx
from scitbx.array_family import flex, shared

expected_bond_lengths = map.stl_string_double()
vdw_radii = map.stl_string_double()
for (e1, e2), length in expected_bond_lengths_by_element_pair.iteritems() :
  expected_bond_lengths[e1+e2] = length
for k, v in vdw.table.items() :
  vdw_radii[k.upper()] = v

# XXX severe duplication: cctbx/eltbx/distance_based_connectivity.py
def build_bond_list (
      sites_cart,
      elements,
      conformer_indices=None,
      search_max_distance=None,
      tolerance_factor_expected_bond_length=1.3,
      fallback_expected_bond_length=2.0,
      fallback_search_max_distance=3.0) :
  assert sites_cart.size() == elements.size()
  assert (isinstance(tolerance_factor_expected_bond_length, float) or
          isinstance(tolerance_factor_expected_bond_length, int))
  assert (isinstance(fallback_expected_bond_length, float) or
          isinstance(fallback_expected_bond_length, int))
  assert (isinstance(fallback_search_max_distance, float) or
          isinstance(fallback_search_max_distance, int))
  if (conformer_indices is None) :
    conformer_indices = flex.size_t(sites_cart.size(), 0)
  stripped_elements = elements.strip().upper()
  if (search_max_distance is None):
    search_max_distance = 2 * max([vdw_radii.get(e, 0.0) for e in elements])
    if (search_max_distance == 0.0):
      search_max_distance = fallback_search_max_distance
    else:
      search_max_distance *= tolerance_factor_expected_bond_length
  else :
    assert isinstance(search_max_distance, float)
  box_symmetry = cctbx.crystal.symmetry(
    unit_cell=cctbx.uctbx.non_crystallographic_unit_cell(
      sites_cart=sites_cart,
      buffer_layer=search_max_distance*(1+1.e-4)), # uncritical tolerance
    space_group_symbol="P1").special_position_settings()
  asu_mappings = box_symmetry.asu_mappings(
    buffer_thickness=search_max_distance,
    sites_cart=sites_cart)
  pair_generator = cctbx.crystal.neighbors_fast_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=search_max_distance,
    minimal=True)
  result = shared.stl_set_unsigned(sites_cart.size())
  bonds = pair_generator.distance_based_connectivity(
    elements=stripped_elements,
    conformer_indices=conformer_indices,
    expected_bond_lengths=expected_bond_lengths,
    vdw_radii=vdw_radii,
    fallback_expected_bond_length=fallback_expected_bond_length,
    tolerance_factor_expected_bond_length=tolerance_factor_expected_bond_length)
  return bonds
#  pair_generator.restart()
#  for pair in pair_generator :
#    pair_elems = tuple(sorted(
#      [stripped_elements[i] for i in [pair.i_seq, pair.j_seq]]))
#    elem_key = "%s%s" % (pair_elems[0], pair_elems[1])
#    ebl = expected_bond_lengths.get(elem_key, None)
#    if (ebl == 0.0):
#      continue
#    if (ebl is None):
#      ebl = max([vdw_table.get(e, 0.0) for e in pair_elems])
#      if (ebl == 0.0):
#        ebl = fallback_expected_bond_length
#        if (ebl is None):
#          continue
#    cutoff_sq = (ebl * tolerance_factor_expected_bond_length)**2
#    if (pair.dist_sq > cutoff_sq):
#      continue
#    result[pair.i_seq].append(pair.j_seq)
#    result[pair.j_seq].append(pair.i_seq)
#  return result

def exercise () :
  # caffeine
  sites_cart = flex.vec3_double([
    (-2.986, 0.015, 1.643),
    (-1.545, 0.015, 1.643),
    (-0.733, 0.015, 2.801),
    (0.592, 0.015, 2.395),
    (0.618, 0.015, 1.034),
    (1.758, 0.015, 0.102),
    (3.092, -0.06, 0.694),
    (1.525, 0.015, -1.360),
    (2.489, -0.024, -2.139),
    (0.158, 0.015, -1.888),
    (-0.025, 0.024, -3.330),
    (-0.986, 0.015, -0.959),
    (-2.155, 0.008, -1.408),
    (-0.733, 0.015, 0.565),
    (-3.346, 0.016, 2.662),
    (-3.347, 0.896, 1.133),
    (-3.347, -0.868, 1.136),
    (-1.083, 0.02, 3.822),
    (3.184, -0.975, 1.26),
    (3.245, 0.785, 1.348),
    (3.835, -0.047, -0.09),
    (0.508, 0.861, -3.756),
    (-1.076, 0.113, -3.560),
    (0.358, -0.896, -3.748)
  ])
  elements = flex.std_string([
    ' C', ' N', ' C', ' N', ' C', ' N', ' C', ' C', ' O', ' N', ' C', ' C',
    ' O', ' C', ' H', ' H', ' H', ' H', ' H', ' H', ' H', ' H', ' H', ' H'])
  bonds = build_bond_list(
    sites_cart=sites_cart,
    elements=elements)
  assert bonds.size() == sites_cart.size()
  #print list(bonds[0])
  assert list(bonds[0]) == [1, 14, 15, 16]
  print "OK"

if __name__ == "__main__" :
  exercise()
