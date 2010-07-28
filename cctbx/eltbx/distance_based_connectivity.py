expected_bond_lengths_by_element_pair = {
  # based on elbow/chemistry/BondLengths.py rev. 42
  # max of averages, rounded to one decimal
('H', 'H'): 0.0,
('AL', 'F'): 1.8,
('AS', 'C'): 2.0,
('AS', 'O'): 1.7,
('B', 'C'): 1.6,
('B', 'O'): 1.5,
('BR', 'C'): 2.0,
('C', 'C'): 1.5,
('C', 'CL'): 1.8,
('C', 'F'): 1.3,
('C', 'H'): 1.1,
('C', 'HG'): 2.3,
('C', 'N'): 1.4,
('C', 'O'): 1.4,
('C', 'P'): 1.7,
('C', 'S'): 1.7,
('C', 'SE'): 1.9,
('CO', 'N'): 2.0,
('CU', 'N'): 2.1,
('CU', 'O'): 1.8,
('F', 'O'): 1.8,
('FE', 'FE'): 2.2,
('FE', 'N'): 2.0,
('FE', 'O'): 2.0,
('FE', 'S'): 2.2,
('H', 'N'): 1.0,
('H', 'O'): 1.0,
('H', 'S'): 1.0,
('HG', 'O'): 2.3,
('MG', 'N'): 2.0,
('MG', 'O'): 2.2,
('N', 'N'): 1.3,
('N', 'NI'): 2.1,
('N', 'O'): 1.4,
('N', 'P'): 1.6,
('N', 'RU'): 1.8,
('N', 'S'): 1.6,
('O', 'O'): 1.4,
('O', 'P'): 1.6,
('O', 'S'): 1.5,
('O', 'U'): 1.8,
('O', 'V'): 2.0,
('O', 'W'): 2.0,
('P', 'S'): 1.7,
('S', 'S'): 2.0}

def build_edge_list(
      sites_cart,
      elements,
      search_max_distance=None,
      tolerance_factor_expected_bond_length=1.3,
      fallback_expected_bond_length=2,
      fallback_search_max_distance=3):
  result = []
  if (sites_cart.size() == 0):
    return result
  import cctbx.crystal
  from cctbx.eltbx.van_der_waals_radii import vdw
  vdw_table = dict([(k.upper(),v) for k,v in vdw.table.items()])
  elements_strip_upper = [e.strip().upper() for e in elements]
  if (search_max_distance is None):
    search_max_distance = max([vdw_table.get(e, 0.0)
      for e in elements_strip_upper])
    if (search_max_distance == 0.0):
      search_max_distance = fallback_search_max_distance
    else:
      search_max_distance *= tolerance_factor_expected_bond_length
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
  for pair in pair_generator:
    pair_elems = tuple(sorted([elements_strip_upper[i]
      for i in [pair.i_seq, pair.j_seq]]))
    ebl = expected_bond_lengths_by_element_pair.get(pair_elems)
    if (ebl == 0.0):
      continue
    if (ebl is None):
      ebl = max([vdw_table.get(e, 0.0) for e in pair_elems])
      if (ebl == 0.0):
        ebl = fallback_expected_bond_length
        if (ebl is None):
          continue
    cutoff_sq = (ebl * tolerance_factor_expected_bond_length)**2
    if (pair.dist_sq > cutoff_sq):
      continue
    result.append(tuple(sorted((pair.i_seq, pair.j_seq))))
  return result
