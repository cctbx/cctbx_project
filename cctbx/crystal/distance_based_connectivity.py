from cctbx.eltbx.distance_based_connectivity import \
  expected_bond_lengths_by_element_pair
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

def build_simple_two_way_bond_sets (
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
    search_max_distance = 2 * max([vdw_radii.get(e, 0.0)
      for e in stripped_elements])
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
  return pair_generator.distance_based_simple_two_way_bond_sets(
    elements=stripped_elements,
    conformer_indices=conformer_indices,
    expected_bond_lengths=expected_bond_lengths,
    vdw_radii=vdw_radii,
    fallback_expected_bond_length=fallback_expected_bond_length,
    tolerance_factor_expected_bond_length=tolerance_factor_expected_bond_length)
