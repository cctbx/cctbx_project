from __future__ import absolute_import, division, print_function
import libtbx.load_env
import iotbx.map_manager
import mmtbx.model
from iotbx.map_model_manager import map_model_manager
from scitbx.array_family import flex

# get_any_map_manager() and its callers with half maps only.
# Managers are used fresh: routines that reference map_manager() (and so
#  create the averaged full map) are not called before the call under test.

data_ccp4 = libtbx.env.under_dist('iotbx',
  'regression/data/non_zero_origin_map.ccp4')

expected_shift_cart = (-74.70333, -72.30750, -73.74375)
expected_site = (-66.70333, -62.30750, -61.74375)

def get_half_maps(mm):
  h1 = mm.deep_copy()
  h2 = mm.deep_copy()
  md = h2.map_data().as_1d()  # shares data with h2
  md.set_selected(flex.bool(md.size(), True), md * 2. + 1.)
  assert not same_map_values(h1, h2)  # h2 differs from h1
  return h1, h2

def new_model(mm):
  m = mmtbx.model.manager.from_sites_cart(atom_name = ' CA ',
    resname = 'ALA', chain_id = 'A', b_iso = 30., occ = 1.,
    scatterer = 'C', sites_cart = flex.vec3_double(((8, 10, 12),)),
    crystal_symmetry = mm.crystal_symmetry())
  m.set_shift_cart((0, 0, 0))
  return m

def close(a, b, eps = 1.e-3):
  return max([abs(x - y) for x, y in zip(a, b)]) < eps

def same_map_values(a, b, eps = 1.e-6):
  d = a.map_data().as_1d() - b.map_data().as_1d()
  return flex.max(flex.abs(d)) < eps

def map_dict_identity(mmm):
  return [(k, id(v)) for k, v in mmm.map_dict().items()]

def exercise():
  print("Data:", data_ccp4)
  mm = iotbx.map_manager.map_manager(data_ccp4)
  failures = []

  def half_maps_only():
    h1, h2 = get_half_maps(mm)
    return map_model_manager(map_manager_1 = h1, map_manager_2 = h2)

  # (1) Main map present: returns the main map_manager object
  for label, mmm in (
      ("main map", map_model_manager(map_manager = mm.deep_copy())),
      ("main map and half maps", map_model_manager(
         map_manager = mm.deep_copy(),
         map_manager_1 = get_half_maps(mm)[0],
         map_manager_2 = get_half_maps(mm)[1]))):
    main = mmm.map_dict().get('map_manager')
    if main is None or mmm.get_any_map_manager() is not main:
      failures.append("B1(1) %s: get_any_map_manager() is not the main map" %(
        label))

  # (2) Half maps only, map_manager() never referenced
  mmm = half_maps_only()
  before = map_dict_identity(mmm)
  existing = list(mmm.map_dict().values())
  any_mm = mmm.get_any_map_manager()
  if any_mm is None:
    failures.append("B1(2) half maps only: get_any_map_manager() returned None")
  elif not [v for v in existing if v is any_mm]:
    failures.append("B1(2) half maps only: returned map is not an existing map")
  elif not (any_mm is mmm.map_dict().get('map_manager_1') or
            any_mm is mmm.map_dict().get('map_manager_2')):
    failures.append("B1(2) half maps only: returned map is not a half map")
  if mmm.map_dict().get('map_manager') is not None:
    failures.append("B1(2) half maps only: an averaged main map was created")
  if map_dict_identity(mmm) != before:
    failures.append("B1(2) half maps only: map_dict changed")

  # (3) Empty manager, and model only (dummy map from the model)
  mmm = map_model_manager()
  try:
    any_mm = mmm.get_any_map_manager()
    if any_mm is not None:
      failures.append("B1(3) empty manager: did not return None")
  except Exception as e:
    failures.append("B1(3) empty manager: raised %s: %s" %(
      type(e).__name__, e))
  mmm = map_model_manager(model = new_model(mm))
  any_mm = mmm.get_any_map_manager()
  if any_mm is None or not any_mm.is_dummy_map_manager() or \
      mmm.map_dict().get('map_manager') is not any_mm:
    failures.append("B1(3) model only: no dummy map made from the model")

  # (4) shift_any_model_to_match on a fresh half-maps-only manager,
  #   compared with a manager where map_manager() was referenced first
  ref = half_maps_only()
  ref.map_manager()
  ref_model = new_model(mm)
  ref.shift_any_model_to_match(ref_model)
  ref_site = tuple(ref_model.get_sites_cart()[0])
  ref_shift_cart = tuple(ref_model.shift_cart())
  print("Reference (map_manager() referenced): site", ref_site,
     "shift_cart", ref_shift_cart)
  if not (close(ref_site, expected_site) and
          close(ref_shift_cart, expected_shift_cart)):
    failures.append("B1(4) reference: unexpected site %s shift_cart %s" %(
      str(ref_site), str(ref_shift_cart)))
  mmm = half_maps_only()
  model = new_model(mm)
  mmm.shift_any_model_to_match(model)
  site = tuple(model.get_sites_cart()[0])
  shift_cart = tuple(model.shift_cart())
  print("Fresh half maps only: site", site, "shift_cart", shift_cart)
  if close(site, (8, 10, 12)):
    failures.append(
      "B1(4) shift_any_model_to_match half maps only: model not moved")
  elif not (close(site, ref_site) and close(shift_cart, ref_shift_cart)):
    failures.append("B1(4) shift_any_model_to_match half maps only: " +
      "site %s shift_cart %s differ from reference" %(
      str(site), str(shift_cart)))

  # (5) get_map_model_manager_with_selected on a fresh half-maps-only manager
  mmm = half_maps_only()
  before = map_dict_identity(mmm)
  h1 = mmm.map_dict()['map_manager_1']
  h2 = mmm.map_dict()['map_manager_2']
  average = half_maps_only().map_manager()  # what map_manager() would give
  try:
    w = mmm.get_map_model_manager_with_selected(
      map_id_list = ['map_manager_1', 'map_manager_2'])
  except Exception as e:
    w = None
    failures.append("B1(5) get_map_model_manager_with_selected half maps " +
      "only: raised %s: %s" %(type(e).__name__, e))
  if map_dict_identity(mmm) != before:
    failures.append("B1(5) source map_dict changed")
  if w is not None:
    for map_id, h in (('map_manager_1', h1), ('map_manager_2', h2)):
      sel = w.map_dict().get(map_id)
      if sel is None or not same_map_values(sel, h):
        failures.append("B1(5) selected manager lacks %s" %(map_id))
    main = w.map_dict().get('map_manager')
    if main is not None and (main is h1 or main is h2):
      failures.append("B1(5) selected manager presents a half map as its " +
        "full map")
    main = w.map_manager()
    if main is None or not same_map_values(main, average):
      failures.append("B1(5) selected manager full map is not the average " +
        "of the half maps")

  if failures:
    for f in failures:
      print("FAIL:", f)
    raise AssertionError("B1 FAILURES: %s" %(len(failures)))
  print("All B1 checks passed")

if __name__ == "__main__":
  exercise()
  print("OK")
