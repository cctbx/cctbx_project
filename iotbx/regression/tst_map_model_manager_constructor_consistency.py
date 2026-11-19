from __future__ import absolute_import, division, print_function
# Tests that the map_model_manager constructor treats extra models like the
#  main model (crystal_symmetry, unit_cell_crystal_symmetry, shift_cart, written
#  output, half maps only) and gives a supplied ncs_object to every map.
#  Extra models are also placed once when the map is already shifted.
import sys
import libtbx.load_env
from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx.array_family import flex
from iotbx.map_manager import map_manager
from iotbx.map_model_manager import map_model_manager
from mmtbx.model import manager as model_manager
from mmtbx.ncs.ncs import ncs
import iotbx.cif
import iotbx.pdb

# Sites of the input model and the sites the main model has after the
#  constructor shifts the origin of non_zero_origin_map.ccp4 (origin 100,100,100)
input_sites = ((8, 10, 12), (14, 15, 16))
expected_shifted_sites = ((-66.70333, -62.30750, -61.74375),
                          (-60.70333, -57.30750, -57.74375))

def get_map_file_name():
  return libtbx.env.under_dist(
      module_name = "iotbx",
      path = "regression/data/non_zero_origin_map.ccp4")

def make_model(crystal_symmetry):
  return model_manager.from_sites_cart(
         atom_name = ' CA ',
         resname = 'ALA',
         chain_id = 'A',
         b_iso = 30.,
         occ = 1.,
         scatterer = 'C',
         sites_cart = flex.vec3_double(input_sites),
         crystal_symmetry = crystal_symmetry)

def cryst1_lines(model):
  return [x for x in model.model_as_pdb().splitlines()
     if x.startswith("CRYST1")]

def pdb_sites(model):
  inp = iotbx.pdb.input(source_info = None,
     lines = flex.split_lines(model.model_as_pdb()))
  return inp.atoms().extract_xyz()

def cif_block(model):
  cif_model = iotbx.cif.reader(input_string = model.model_as_mmcif()).model()
  return list(cif_model.values())[0]

def cif_cell(model):
  block = cif_block(model)
  return [block["_cell.%s" %(key)] for key in ("length_a", "length_b",
     "length_c", "angle_alpha", "angle_beta", "angle_gamma")]

def cif_sites(model):
  block = cif_block(model)
  x, y, z = [flex.double([float(v) for v in block["_atom_site.Cartn_%s" %(k)]])
     for k in ("x", "y", "z")]
  return flex.vec3_double(x, y, z)

def sites_equal(a, b, eps):
  return a.size() == b.size() and approx_equal(
     list(a), list(b), eps = eps, out = None)

class failure_list(object):
  def __init__(self):
    self.failures = []
  def check(self, ok, text):
    print("%s: %s" %("OK  " if ok else "FAIL", text))
    if not ok:
      self.failures.append(text)

def check_models_match(mmm, f, label):
  '''
  Criterion items (1) and (2): 'extra' has the same symmetries and shift_cart
  as 'model', the manager is consistent, and neither model was moved.
  '''
  model = mmm.get_model_by_id('model')
  extra = mmm.get_model_by_id('extra')
  for name, a, b in (
      ("crystal_symmetry",
         model.crystal_symmetry(), extra.crystal_symmetry()),
      ("unit_cell_crystal_symmetry",
         model.unit_cell_crystal_symmetry(),
         extra.unit_cell_crystal_symmetry()),):
    f.check(a.is_similar_symmetry(b,
       absolute_angle_tolerance = mmm._absolute_angle_tolerance,
       absolute_length_tolerance = mmm._absolute_length_tolerance),
       "A1(1) %s: extra %s equals model %s (%s vs %s)" %(
       label, name, name, b.unit_cell().parameters(),
       a.unit_cell().parameters()))
  f.check(approx_equal(extra.shift_cart(), model.shift_cart(),
       eps = mmm._shift_tol, out = None),
       "A1(1) %s: extra shift_cart equals model shift_cart (%s vs %s)" %(
       label, str(extra.shift_cart()), str(model.shift_cart())))
  f.check(mmm.check_consistency(stop_on_errors = False, print_errors = True),
       "A1(1) %s: check_consistency() is True" %(label))
  for m_id in ('model', 'extra'):
    sites = mmm.get_model_by_id(m_id).get_sites_cart()
    f.check(sites_equal(sites, flex.vec3_double(expected_shifted_sites),
       eps = 1.e-3),
       "A1(2) %s: sites of '%s' are the sites the main model gets (%s)" %(
       label, m_id, str(sites[0])))

def tst_extra_models_match_main_model():
  ''' Criterion A1: extra models are treated like the main model '''
  file_name = get_map_file_name()
  f = failure_list()
  mm = map_manager(file_name)
  box_cs = mm.crystal_symmetry()
  full_cs = mm.unit_cell_crystal_symmetry()
  assert not box_cs.is_similar_symmetry(full_cs)  # the map is boxed

  # Models with box crystal_symmetry (the failing case) and full cell control
  for label, cs in (('box cell', box_cs), ('full cell', full_cs)):
    mmm = map_model_manager(map_manager = map_manager(file_name),
      model = make_model(cs),
      extra_model_list = [make_model(cs)],
      extra_model_id_list = ['extra'])
    check_models_match(mmm, f, label)
    if label != 'box cell':
      continue
    model = mmm.get_model_by_id('model')
    extra = mmm.get_model_by_id('extra')

    # (3) Written cell metadata
    f.check(cryst1_lines(extra) == cryst1_lines(model),
       "A1(3) PDB CRYST1 of extra equals model (%s vs %s)" %(
        cryst1_lines(extra), cryst1_lines(model)))
    f.check(cif_cell(extra) == cif_cell(model),
       "A1(3) mmCIF cell of extra equals model (%s vs %s)" %(
        cif_cell(extra), cif_cell(model)))

    # (4) remove origin shift and unit_cell_crystal_symmetry
    sites_before = dict([(m_id, mmm.get_model_by_id(m_id).get_sites_cart(
       ).deep_copy()) for m_id in ('model', 'extra')])
    mmm.remove_origin_shift_and_unit_cell_crystal_symmetry()
    for m_id in ('model', 'extra'):
      m = mmm.get_model_by_id(m_id)
      f.check(approx_equal(m.shift_cart(), (0, 0, 0), eps = 1.e-3,
         out = None),
         "A1(4) shift_cart of '%s' is (0,0,0) after remove_origin_shift (%s)"
         %(m_id, str(m.shift_cart())))
      f.check(sites_equal(m.get_sites_cart(), sites_before[m_id],
         eps = 1.e-5),
         "A1(4) in-memory sites of '%s' unchanged by remove_origin_shift" %(
         m_id))
    model_sites = model.get_sites_cart()
    for fmt, get_sites in (('PDB', pdb_sites), ('mmCIF', cif_sites)):
      written_model = get_sites(model)
      written_extra = get_sites(extra)
      f.check(sites_equal(written_extra, written_model, eps = 1.e-3),
         "A1(4) %s coordinates written for extra equal model (%s vs %s)" %(
         fmt, str(written_extra[0]), str(written_model[0])))
      f.check(sites_equal(written_extra, model_sites, eps = 1.e-3),
         "A1(4) %s coordinates written for extra equal in-memory sites" %(
         fmt))

  # (5) Half maps only, then (1) and (2)
  try:
    mmm = map_model_manager(map_manager_1 = map_manager(file_name),
      map_manager_2 = map_manager(file_name),
      model = make_model(box_cs),
      extra_model_list = [make_model(box_cs)],
      extra_model_id_list = ['extra'])
  except AttributeError as e:
    mmm = None
    f.check(False,
       "A1(5) half maps only: constructor raised AttributeError: %s" %(e))
  if mmm is not None:
    f.check(True, "A1(5) half maps only: constructor ran without error")
    check_models_match(mmm, f, 'half maps only')

  if f.failures:
    raise AssertionError("A1 FAILURES:\n  " + "\n  ".join(f.failures))
  print("A1 extra models: OK")

def make_two_fold_ncs(center = (85., 86., 85.)):
  ''' NCS object with identity and a two-fold about z through center '''
  c = matrix.col(center)
  r_identity = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))
  r_two_fold = matrix.sqr((-1, 0, 0, 0, -1, 0, 0, 0, 1))
  ncs_object = ncs()
  ncs_object.import_ncs_group(
    ncs_rota_matr = [r_identity, r_two_fold],
    center_orth = [c, c],
    trans_orth = [matrix.col((0, 0, 0)), c - r_two_fold * c])
  return ncs_object

def ncs_operators(ncs_object):
  ''' rotations, translations and centers of all groups as flat list '''
  values = []
  for g in ncs_object.ncs_groups():
    for r, t, c in zip(g.rota_matrices(), g.translations_orth(),
        g.centers()):
      values += list(r) + list(t) + list(c)
  return values

def check_all_maps_have_ncs(mmm, map_ids, supplied_ncs, f, label):
  main_ncs = mmm.map_manager().ncs_object()
  # Supplied NCS moved by the coordinate shift of the maps and models
  expected_ncs = supplied_ncs.coordinate_offset(
     coordinate_offset = mmm.map_manager().shift_cart())
  f.check(main_ncs is not None, "A2 %s: main map has NCS" %(label))
  if main_ncs is None:
    return
  # The main map NCS is the supplied one shifted with the map
  f.check(approx_equal(ncs_operators(main_ncs), ncs_operators(expected_ncs),
     eps = 1.e-3, out = None),
     "A2 %s: main map NCS operators are the supplied ones shifted with map"
      %(label))
  for map_id in map_ids:
    m = mmm.get_map_manager_by_id(map_id)
    n = m.ncs_object()
    f.check(n is not None, "A2 %s: map '%s' has an NCS object" %(
       label, map_id))
    if n is None:
      continue
    f.check(n.max_operators() == 2 and approx_equal(n.shift_cart(),
       m.shift_cart(), eps = 1.e-3, out = None),
       "A2 %s: NCS of map '%s' has 2 operators and the map's shift_cart" %(
       label, map_id))
    f.check(approx_equal(ncs_operators(n), ncs_operators(main_ncs),
       eps = 1.e-3, out = None),
       "A2 %s: NCS operators of map '%s' equal those of the main map" %(
       label, map_id))

def tst_supplied_ncs_object_reaches_all_maps():
  ''' Criterion A2: a supplied ncs_object is given to every map '''
  file_name = get_map_file_name()
  f = failure_list()
  mm = map_manager(file_name)
  supplied_ncs = make_two_fold_ncs()

  def get_mask():
    mask_mm = map_manager(file_name)
    mask_mm.set_is_mask(True)
    return mask_mm

  # Reference model coordinates without an ncs_object
  reference = map_model_manager(map_manager = map_manager(file_name),
    model = make_model(mm.crystal_symmetry()))
  reference_sites = reference.model().get_sites_cart()

  # Full map, half maps, an extra map and a mask
  mmm = map_model_manager(map_manager = map_manager(file_name),
    map_manager_1 = map_manager(file_name),
    map_manager_2 = map_manager(file_name),
    extra_map_manager_list = [map_manager(file_name), get_mask()],
    extra_map_manager_id_list = ['extra_1', 'map_manager_mask'],
    model = make_model(mm.crystal_symmetry()),
    ncs_object = supplied_ncs.deep_copy())
  check_all_maps_have_ncs(mmm, ['map_manager', 'map_manager_1',
     'map_manager_2', 'extra_1', 'map_manager_mask'],
     supplied_ncs, f, 'full map, half maps, extra maps')
  f.check(sites_equal(mmm.model().get_sites_cart(), reference_sites,
     eps = 1.e-5),
     "A2 supplying the NCS object does not change model coordinates")

  # A map that comes with its own (different) NCS object gets the supplied one
  mm2_with_ncs = map_manager(file_name)
  mm2_with_ncs.set_ncs_object(make_two_fold_ncs(center = (80., 80., 80.)))
  try:
    mmm = map_model_manager(map_manager = map_manager(file_name),
      map_manager_1 = map_manager(file_name),
      map_manager_2 = mm2_with_ncs,
      model = make_model(mm.crystal_symmetry()),
      ncs_object = supplied_ncs.deep_copy())
  except AssertionError as e:
    mmm = None
    f.check(False, "A2 map with its own NCS: constructor raised "+
      "AssertionError: %s" %(" ".join(str(e).split())))
  if mmm is not None:
    check_all_maps_have_ncs(mmm, ['map_manager', 'map_manager_1',
       'map_manager_2'], supplied_ncs, f, 'map with its own NCS')

  # Half maps only
  mmm = map_model_manager(
    map_manager_1 = map_manager(file_name),
    map_manager_2 = map_manager(file_name),
    extra_map_manager_list = [map_manager(file_name)],
    extra_map_manager_id_list = ['extra_1'],
    model = make_model(mm.crystal_symmetry()),
    ncs_object = supplied_ncs.deep_copy())
  check_all_maps_have_ncs(mmm, ['map_manager_1', 'map_manager_2',
     'extra_1', 'map_manager'], supplied_ncs, f, 'half maps only')
  f.check(sites_equal(mmm.model().get_sites_cart(), reference_sites,
     eps = 1.e-5),
     "A2 half maps only: supplying the NCS object does not change model "+
     "coordinates")

  if f.failures:
    raise AssertionError("A2 FAILURES:\n  " + "\n  ".join(f.failures))
  print("A2 supplied NCS object: OK")

# shift_cart of non_zero_origin_map.ccp4 after its origin is shifted
expected_shift_cart = (-74.70333, -72.30750, -73.74375)

def check_shifted_map_case(mmm, f, label, sites_before = None):
  '''
  Criterion C1: 'extra' has the sites, shift_cart and symmetries of 'model',
  the sites are those the main model gets, and the manager is consistent.
  Optionally the sites are unchanged from sites_before (no coordinate moves).
  '''
  model = mmm.get_model_by_id('model')
  extra = mmm.get_model_by_id('extra')
  for name, a, b in (
      ("crystal_symmetry",
         model.crystal_symmetry(), extra.crystal_symmetry()),
      ("unit_cell_crystal_symmetry",
         model.unit_cell_crystal_symmetry(),
         extra.unit_cell_crystal_symmetry()),):
    f.check(a.is_similar_symmetry(b,
       absolute_angle_tolerance = mmm._absolute_angle_tolerance,
       absolute_length_tolerance = mmm._absolute_length_tolerance),
       "C1 %s: extra %s equals model %s" %(label, name, name))
  for m_id in ('model', 'extra'):
    m = mmm.get_model_by_id(m_id)
    f.check(approx_equal(m.shift_cart(), expected_shift_cart, eps = 1.e-3,
       out = None),
       "C1 %s: shift_cart of '%s' is the map's shift_cart (%s)" %(
       label, m_id, str(m.shift_cart())))
    sites = m.get_sites_cart()
    f.check(sites_equal(sites, flex.vec3_double(expected_shifted_sites),
       eps = 1.e-3),
       "C1 %s: sites of '%s' are the sites the main model gets (%s)" %(
       label, m_id, str(sites[0])))
    if sites_before is not None:
      f.check(sites_equal(sites, sites_before[m_id], eps = 1.e-5),
         "C1 %s: sites of '%s' were not moved" %(label, m_id))
  f.check(sites_equal(extra.get_sites_cart(), model.get_sites_cart(),
     eps = 1.e-3), "C1 %s: sites of extra equal sites of model" %(label))
  f.check(mmm.check_consistency(stop_on_errors = False, print_errors = True),
       "C1 %s: check_consistency() is True" %(label))

def tst_extra_models_with_shifted_map():
  '''
  Criterion C1: extra models are placed once, also when the map passed to
  the constructor is already shifted
  '''
  file_name = get_map_file_name()
  f = failure_list()

  def get_models(mm, in_shifted_frame = False):
    cs = mm.unit_cell_crystal_symmetry()
    models = {'model': make_model(cs), 'extra': make_model(cs)}
    if in_shifted_frame:
      for m in models.values():
        mm.shift_model_to_match_map(m)
    return models

  def construct(mm, models, label):
    try:
      return map_model_manager(map_manager = mm, model = models['model'],
        extra_model_list = [models['extra']],
        extra_model_id_list = ['extra'])
    except AssertionError as e:
      extra = models['extra']
      message = " ".join(str(e).split())
      f.check(False, "C1 %s: constructor raised AssertionError%s "
       "(extra model afterwards: sites %s shift_cart %s)" %(
       label, (": " + message) if message else "",
       str(extra.get_sites_cart()[0]), str(extra.shift_cart())))
      return None

  # (1) Map shifted before the constructor, models in the original frame
  label = '(1) pre-shifted map'
  mm = map_manager(file_name)
  mm.shift_origin()
  mmm = construct(mm, get_models(mm), label)
  if mmm is not None:
    check_shifted_map_case(mmm, f, label)

  # (2) Map shifted before the constructor, models already in shifted frame
  label = '(2) pre-shifted map, models in shifted frame'
  mm = map_manager(file_name)
  mm.shift_origin()
  models = get_models(mm, in_shifted_frame = True)
  sites_before = dict([(m_id, m.get_sites_cart().deep_copy())
     for m_id, m in models.items()])
  for m_id in ('model', 'extra'):
    f.check(sites_equal(sites_before[m_id],
      flex.vec3_double(expected_shifted_sites), eps = 1.e-3),
      "C1 %s: set-up: '%s' is in the shifted frame" %(label, m_id))
  mmm = construct(mm, models, label)
  if mmm is not None:
    check_shifted_map_case(mmm, f, label, sites_before = sites_before)

  # (3) The same map object used for a second manager (the first manager
  #  shifts the map in place)
  label = '(3) same map object reused'
  mm = map_manager(file_name)
  first = map_model_manager(map_manager = mm,
    model = make_model(mm.unit_cell_crystal_symmetry()))
  f.check(first.map_manager() is mm and approx_equal(mm.shift_cart(),
    expected_shift_cart, eps = 1.e-3, out = None),
    "C1 %s: set-up: the first manager shifted the map in place" %(label))
  mmm = construct(mm, get_models(mm), label)
  if mmm is not None:
    check_shifted_map_case(mmm, f, label)

  # (4) Control: map not shifted before the constructor
  label = '(4) unshifted map'
  mm = map_manager(file_name)
  mmm = construct(mm, get_models(mm), label)
  if mmm is not None:
    check_shifted_map_case(mmm, f, label)

  if f.failures:
    raise AssertionError("C1 FAILURES:\n  " + "\n  ".join(f.failures))
  print("C1 extra models with shifted map: OK")

if __name__ == "__main__":
  # Optional argument: models, ncs or shifted to run only that test
  args = sys.argv[1:]
  if not args or 'models' in args:
    tst_extra_models_match_main_model()
  if not args or 'ncs' in args:
    tst_supplied_ncs_object_reaches_all_maps()
  if not args or 'shifted' in args:
    tst_extra_models_with_shifted_map()
  print("OK")
