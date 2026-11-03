from __future__ import absolute_import, division, print_function
import sys
import threading
import libtbx.load_env
import iotbx.pdb
import iotbx.map_manager
import mmtbx.model
from cctbx import crystal
from scitbx.array_family import flex
from six.moves import cStringIO as StringIO
from iotbx.map_model_manager import map_model_manager, \
  _match_map_model_ncs, is_same_shift_cart

# Small defects in map_model_manager (M1-M6). Each function can be run
#  alone: libtbx.python tst_map_model_manager_cleanup.py m1 m3 ...

data_ccp4 = libtbx.env.under_dist('iotbx',
  'regression/data/non_zero_origin_map.ccp4')

def zero_origin_map():
  cs = crystal.symmetry((30, 40, 50, 90, 90, 90), 1)
  return iotbx.map_manager.map_manager(
    map_data = flex.double(flex.grid((30, 40, 50)), 1.0),
    unit_cell_grid = (30, 40, 50), unit_cell_crystal_symmetry = cs,
    wrapping = True)

def new_model(crystal_symmetry, site = (8, 10, 12)):
  return mmtbx.model.manager.from_sites_cart(atom_name = ' CA ',
    resname = 'ALA', chain_id = 'A', b_iso = 30., occ = 1.,
    scatterer = 'C', sites_cart = flex.vec3_double((site,)),
    crystal_symmetry = crystal_symmetry)

def model_without_symmetry():
  lines = ["ATOM      1  CA  ALA A   1       8.000  10.000  12.000"
    "  1.00 30.00           C", "END"]
  return mmtbx.model.manager(
    model_input = iotbx.pdb.input(source_info = None, lines = lines))

def m1():
  '''_generate_new_model_id returns temp_0, then temp_1 when temp_0 exists'''
  mmm = map_model_manager(map_manager = zero_origin_map())
  assert mmm._generate_new_model_id() == 'temp_0'
  mmm.add_model_by_id(new_model(mmm.crystal_symmetry()), 'temp_0')
  # Hang guard on every platform: call it in a daemon thread and wait 10 s
  result = []
  def call():
    result.append(mmm._generate_new_model_id())
  thread = threading.Thread(target = call)
  thread.daemon = True
  thread.start()
  thread.join(10)
  if thread.is_alive():
    raise AssertionError("M1 hang guard: _generate_new_model_id() did not "
      "return within 10 s with 'temp_0' present")
  assert len(result) == 1, "M1: _generate_new_model_id() raised an exception"
  new_id = result[0]
  assert new_id == 'temp_1', new_id

def m2():
  '''read_ncs_file with a one-operator ncs_spec file'''
  from mmtbx.ncs.ncs import ncs, is_identity
  unit = ncs()
  unit.set_unit_ncs()
  f = StringIO()
  unit.format_all_for_group_specification(out = f)
  file_name = 'tst_map_model_manager_cleanup_unit.ncs_spec'
  with open(file_name, 'w') as out:
    out.write(f.getvalue())
  match = _match_map_model_ncs()
  match.add_map_manager(zero_origin_map())
  assert match.map_manager().ncs_object() is None
  match.read_ncs_file(file_name)
  ncs_object = match.map_manager().ncs_object()
  assert ncs_object is not None
  assert ncs_object.max_operators() == 1, ncs_object.max_operators()
  group = ncs_object.ncs_groups()[0]
  assert is_identity(group.rota_matrices()[0], group.translations_orth()[0])

def m3():
  '''add_crystal_symmetry_to_model_if_necessary without map_manager'''
  mm = iotbx.map_manager.map_manager(data_ccp4)
  mmm = map_model_manager(map_manager = mm)
  model = model_without_symmetry()
  assert not model.crystal_symmetry()
  mmm.add_crystal_symmetry_to_model_if_necessary(model)
  map_manager = mmm.map_manager()
  assert model.crystal_symmetry().is_similar_symmetry(
    map_manager.crystal_symmetry())
  assert model.unit_cell_crystal_symmetry().is_similar_symmetry(
    map_manager.unit_cell_crystal_symmetry())

def m4():
  '''deep_copy and customized_copy keep the tolerances'''
  mmm = map_model_manager(map_manager = zero_origin_map(),
    absolute_angle_tolerance = 5, absolute_length_tolerance = 5,
    shift_tol = 0.5)
  for label, copy in (('deep_copy', mmm.deep_copy()),
      ('customized_copy', mmm.customized_copy())):
    values = (copy._absolute_angle_tolerance,
      copy._absolute_length_tolerance, copy._shift_tol)
    assert values == (5, 5, 0.5), "M4 %s tolerances %s, expected %s" %(
      label, str(values), str((5, 5, 0.5)))

def manager_with_incompatible_maps():
  mmm = map_model_manager(map_manager = zero_origin_map())
  other = mmm.map_manager().deep_copy()
  other.set_original_origin_and_gridding(original_origin = (1, 2, 3))
  assert not other.is_similar(mmm.map_manager())
  mmm._map_dict['incompatible_map'] = other # not checked on purpose
  return mmm

def m5a():
  '''check_consistency names the incompatible map managers'''
  first = "Consistency check failure in map_model_manager 'map_model_manager'"
  ids = ("'incompatible_map'", "'map_manager'")
  mmm = manager_with_incompatible_maps()
  try:
    mmm.check_consistency(stop_on_errors = True)
  except AssertionError as e:
    text = str(e)
  else:
    raise AssertionError("M5a: check_consistency did not raise")
  assert text.startswith(first), text
  for id in ids:
    assert text.find(id) > -1, "M5a: AssertionError text lacks map id %s: %s" %(
      id, text)
  mmm = manager_with_incompatible_maps()
  log = StringIO()
  mmm.set_log(log)
  assert mmm.check_consistency(stop_on_errors = False) is False
  for id in ids:
    assert log.getvalue().find(id) > -1, \
      "M5a: printed log lacks map id %s: %s" %(id, log.getvalue())

def m5b_empty():
  '''check_consistency of an empty manager'''
  assert map_model_manager().check_consistency() is True

def models_only_manager(crystal_symmetries):
  # No real map; add_model_by_id may make a placeholder map from 'model'
  mmm = map_model_manager()
  for i, cs in enumerate(crystal_symmetries):
    mmm.add_model_by_id(new_model(cs), 'model' if i == 0 else 'model_%d' %(i))
  assert len(mmm.models()) == len(crystal_symmetries)
  for mm in mmm.map_managers():
    assert mm.is_dummy_map_manager()
  return mmm

def m5b_model_only():
  '''check_consistency of a manager with only a model is True, no map made'''
  cs = crystal.symmetry((30, 40, 50, 90, 90, 90), 1)
  for mmm in (map_model_manager(model = new_model(cs)),
      models_only_manager([cs])):
    map_ids = mmm.map_id_list()
    assert map_ids == [], map_ids
    assert mmm.check_consistency(print_errors = False) is True
    assert mmm.map_id_list() == map_ids, (map_ids, mmm.map_id_list())

def m5b_models_agree():
  '''check_consistency of two models with the same symmetry is True'''
  cs = crystal.symmetry((30, 40, 50, 90, 90, 90), 1)
  mmm = models_only_manager([cs, cs])
  map_ids = mmm.map_id_list()
  assert mmm.check_consistency(print_errors = False) is True
  assert mmm.map_id_list() == map_ids, (map_ids, mmm.map_id_list())

def m5b_models_disagree():
  '''check_consistency of two models with different symmetry fails'''
  first = "Consistency check failure in map_model_manager"
  css = [crystal.symmetry((30, 40, 50, 90, 90, 90), 1),
    crystal.symmetry((35, 40, 50, 90, 90, 90), 1)]
  mmm = models_only_manager(css)
  assert mmm.check_consistency(stop_on_errors = False,
    print_errors = False) is False
  mmm = models_only_manager(css)
  try:
    mmm.check_consistency(stop_on_errors = True, print_errors = False)
  except AssertionError as e:
    text = str(e)
  else:
    raise AssertionError("M5b: check_consistency did not raise")
  assert text.startswith(first), text

def m5b_model_then_add():
  '''models added after 'model' (a placeholder map exists): models compared'''
  cs = crystal.symmetry((30, 40, 50, 90, 90, 90), 1)
  other_cs = crystal.symmetry((35, 40, 50, 90, 90, 90), 1)
  for how in ('constructor', 'add_model_by_id'):
    for b_cs, expected in ((cs, True), (other_cs, False)):
      if how == 'constructor':
        mmm = map_model_manager(model = new_model(cs))
      else:
        mmm = map_model_manager()
        mmm.add_model_by_id(new_model(cs), 'model')
      mmm.add_model_by_id(new_model(b_cs), 'b')
      assert len(mmm.models()) == 2
      assert len(mmm.map_managers()) == 1 and \
        mmm.map_managers()[0].is_dummy_map_manager(), (how, mmm.map_id_list())
      map_ids = mmm.map_id_list()
      if expected:
        assert mmm.check_consistency(print_errors = False) is True
      else:
        assert mmm.check_consistency(stop_on_errors = False,
          print_errors = False) is False, (how, str(b_cs))
      assert mmm.map_id_list() == map_ids, (map_ids, mmm.map_id_list())

def m5b_shifted_models():
  '''models with symmetry: equal shift_cart -> True, different -> False'''
  cs = crystal.symmetry((30, 40, 50, 90, 90, 90), 1)
  for shift_b, expected in (((1, 2, 3), True), ((1, 2, 4), False)):
    mmm = models_only_manager([cs, cs])
    mmm.get_model_by_id('model').set_shift_cart((1, 2, 3))
    mmm.get_model_by_id('model_1').set_shift_cart(shift_b)
    for m in mmm.models():
      assert m.shift_cart() is not None and \
        m.unit_cell_crystal_symmetry() is not None
    map_ids = mmm.map_id_list()
    if expected:
      assert mmm.check_consistency(print_errors = False) is True
    else:
      assert mmm.check_consistency(stop_on_errors = False,
        print_errors = False) is False, shift_b
    assert mmm.map_id_list() == map_ids, (map_ids, mmm.map_id_list())

def m5b_model_without_symmetry_unchanged():
  '''models without crystal symmetry: the old path is used and returns False
  (as before); the printed line lets the base and candidate be compared'''
  cs = crystal.symmetry((30, 40, 50, 90, 90, 90), 1)
  for case in ('alone', 'with_other'):
    mmm = map_model_manager()
    if case == 'with_other':  # added first: no placeholder made yet
      mmm.add_model_by_id(new_model(cs), 'b')
    mmm.add_model_by_id(model_without_symmetry(), 'model')
    assert mmm.model().crystal_symmetry() is None
    try:
      result = "returned %s" %(str(mmm.check_consistency(
        stop_on_errors = False, print_errors = False)))
    except Exception as e:
      result = "%s: %s" %(type(e).__name__, str(e).splitlines()[0])
    print("M5B_NO_SYMMETRY_RESULT %s: %s" %(case, result))
    assert result == "returned False", (case, result)

def m5b_with_map_unchanged():
  '''check_consistency with a map and a model: True, and False on mismatch'''
  mmm = map_model_manager(map_manager = zero_origin_map(),
    model = new_model(crystal.symmetry((30, 40, 50, 90, 90, 90), 1)))
  assert mmm.check_consistency(print_errors = False) is True
  mmm = map_model_manager(map_manager = zero_origin_map())
  mmm._model_dict['model'] = new_model(  # not checked on purpose
    crystal.symmetry((35, 40, 50, 90, 90, 90), 1))
  assert mmm.check_consistency(stop_on_errors = False,
    print_errors = False) is False

def m6():
  '''is_same_shift_cart is equivalent to the comparison it replaces'''
  def old_differs(a, b):
    return max([abs(x - y) for x, y in zip(a, b)]) > 0.001
  pairs = []
  for base in ((0, 0, 0), (10, 20, 30), (-74.70333, -72.3075, -73.74375)):
    for d in (0, 0.0005, 0.001, 0.0011, 0.002, -0.0005, -0.001, -0.0011,
        -0.002):
      for axis in range(3):
        other = list(base)
        other[axis] += d
        pairs.append((base, tuple(other)))
    pairs.append((base, tuple([x + 0.0005 for x in base])))
    pairs.append((base, tuple([x - 0.002 for x in base])))
    pairs.append((base, (base[0] + 0.0005, base[1] - 0.0011, base[2])))
  for a, b in pairs:
    assert (not is_same_shift_cart(a, b)) == old_differs(a, b), (a, b)
  # add_model_by_id shifts a model whose shift_cart differs by 0.002 only
  for d, expect_shift in ((0.002, True), (0.0005, False)):
    mmm = map_model_manager(map_manager = zero_origin_map())
    assert mmm.map_manager().shift_cart() == (0, 0, 0)
    model = new_model(mmm.crystal_symmetry())
    model.set_shift_cart((d, 0, 0))
    x_before = model.get_sites_cart()[0][0]
    mmm.add_model_by_id(model, 'm')
    dx = mmm.get_model_by_id('m').get_sites_cart()[0][0] - x_before
    shift_x = mmm.get_model_by_id('m').shift_cart()[0]
    if expect_shift:
      assert abs(dx + d) < 1.e-6 and abs(shift_x) < 1.e-6, (d, dx, shift_x)
    else:
      assert abs(dx) < 1.e-6 and abs(shift_x - d) < 1.e-6, (d, dx, shift_x)

all_tests = ['m1', 'm2', 'm3', 'm4', 'm5a', 'm5b_empty', 'm5b_model_only',
  'm5b_models_agree', 'm5b_models_disagree', 'm5b_model_then_add',
  'm5b_shifted_models', 'm5b_model_without_symmetry_unchanged',
  'm5b_with_map_unchanged', 'm6']

def run(names):
  for name in names:
    assert name in all_tests, name
    print("RUN %s" %(name))
    sys.stdout.flush()
    globals()[name]()
    print("OK %s" %(name))
  print("OK")

if __name__ == "__main__":
  run(sys.argv[1:] or all_tests)
