from __future__ import absolute_import, division, print_function
import os, sys
import libtbx.load_env
from iotbx.data_manager import DataManager
from iotbx.map_model_manager import map_model_manager

def test_01(method = 'model_sharpen',
   expected_results=None):

  # Source data

  data_dir = os.path.dirname(os.path.abspath(__file__))
  data_ccp4 = os.path.join(data_dir, 'data',
                          'non_zero_origin_map.ccp4')
  data_pdb = os.path.join(data_dir, 'data',
                          'non_zero_origin_model.pdb')
  data_ncs_spec = os.path.join(data_dir, 'data',
                          'non_zero_origin_ncs_spec.ncs_spec')

  # Read in data

  dm = DataManager(['ncs_spec','model', 'real_map', 'phil'])
  dm.set_overwrite(True)

  map_file=data_ccp4
  dm.process_real_map_file(map_file)
  mm = dm.get_real_map(map_file)

  model_file=data_pdb
  dm.process_model_file(model_file)
  model = dm.get_model(model_file)

  ncs_file=data_ncs_spec
  dm.process_ncs_spec_file(ncs_file)
  ncs = dm.get_ncs_spec(ncs_file)

  mmm=map_model_manager(
    model = model,
    map_manager_1 = mm.deep_copy(),
    map_manager_2 = mm.deep_copy(),
    ncs_object = ncs,
    wrapping = False)
  mmm.add_map_manager_by_id(
     map_id='external_map',map_manager=mmm.map_manager().deep_copy())
  if method == 'local_resolution_map':
    mmm.generate_map(map_id = 'map_manager_2')
    dd = mmm.resolution()
  else:
    mmm.set_resolution(3)
  mmm.set_log(sys.stdout)

  dc = mmm.deep_copy()

  sharpen_method = getattr(mmm,method)

  # sharpen by method (can be model_sharpen, half_map_sharpen or
  #     external_sharpen and also local_resolution_map)

  if method == 'local_resolution_map':
    mm = sharpen_method()
    x = mm.map_data().as_1d().min_max_mean()
    from libtbx.test_utils import approx_equal
    assert approx_equal((x.min,x.max,x.mean),
      (1.9835316316768663, 2.3146054110530336, 2.1600110833032353),
       eps = 0.01)
    mm = mmm.local_resolution_map(map_id_1 = None, map_id_2 = None,
     model_id = 'model', map_id = 'map_manager_1')
    xx = mm.map_data().as_1d().min_max_mean()
    from libtbx.test_utils import approx_equal
    assert approx_equal((x.mean, xx.mean),
      (2.1600110833032353, 2.16),
       eps = 0.01)

  else: # usual
    def clear_sharpened_map():
      # Remove the sharpened map left by any previous call, so that each
      #   check below can see only the map written by the call before it
      if 'map_manager_scaled' in mmm.map_id_list():
        mmm.remove_map_manager_by_id('map_manager_scaled')

    def sharpened_map_model_cc():
      # Sharpening writes its result to 'map_manager_scaled' and leaves
      #   'map_manager' unchanged, so check the sharpened map
      assert 'map_manager_scaled' in mmm.map_id_list()
      return mmm.map_model_cc(map_id = 'map_manager_scaled')

    clear_sharpened_map()
    sharpen_method(anisotropic_sharpen = False, n_bins=10)
    assert sharpened_map_model_cc() > 0.9
    clear_sharpened_map()
    sharpen_method(anisotropic_sharpen = False, n_bins=10,
       local_sharpen = True)
    assert sharpened_map_model_cc() > 0.9
    clear_sharpened_map()
    sharpen_method(anisotropic_sharpen = True, n_bins=10)
    assert sharpened_map_model_cc() > 0.9
    clear_sharpened_map()
    sharpen_method(anisotropic_sharpen = True, n_bins=10,
       local_sharpen = True, n_boxes = 1)
    assert sharpened_map_model_cc() > 0.9


def test_02():
  """Summary of anisotropic scaling when anisotropy is not available"""

  from libtbx.test_utils import approx_equal

  # Cheap map_model_manager: no data files, no chem_data
  mmm = map_model_manager()
  mmm.set_log(None)
  mmm.generate_map(d_min = 3)
  mmm.add_map_manager_by_id(map_id = 'previous_map',
     map_manager = mmm.map_manager().deep_copy())

  # Supply the anisotropy of each map ourselves, keyed on map_id so that
  #  the two calls made by _get_aniso_before_and_after are distinguishable
  aniso_by_map_id = {}
  def get_aniso_of_map(d_min = None, map_id = None):
    return aniso_by_map_id[map_id]
  mmm._get_aniso_of_map = get_aniso_of_map

  prev_b_cart = (10., 20., 30., 1., 2., 3.)
  new_b_cart = (4., 5., 6., 0., 1., 2.)

  # Both available: b_sharpen is previous minus new, element by element
  aniso_by_map_id['previous_map'] = prev_b_cart
  aniso_by_map_id['map_manager'] = new_b_cart
  info = mmm._get_aniso_before_and_after(d_min = 3,
     map_id = 'map_manager', previous_map_id = 'previous_map')
  assert approx_equal(tuple(info.b_sharpen), (6., 15., 24., 1., 1., 1.))
  assert info.text.find('B-cart') > -1
  assert info.text.find('B-sharpen') > -1
  assert info.text.find('not available') < 0

  # Not available for one map or for both: return normally, b_sharpen is
  #  None and the text says the information is not available
  for prev_value, new_value in [
      (None, new_b_cart),
      (prev_b_cart, None),
      (None, None)]:
    aniso_by_map_id['previous_map'] = prev_value
    aniso_by_map_id['map_manager'] = new_value
    info = mmm._get_aniso_before_and_after(d_min = 3,
       map_id = 'map_manager', previous_map_id = 'previous_map')
    assert info.b_sharpen is None
    assert info.text.find('not available') > -1


def test_03():
  """Removing anisotropy when the overall anisotropy is not available"""

  from libtbx.test_utils import approx_equal
  from scitbx.array_family import flex
  from six.moves import StringIO

  # Cheap map_model_manager: no data files, no chem_data.  Two non-mask maps
  #  so that removal from all maps means something
  mmm = map_model_manager()
  mmm.set_log(None)
  mmm.generate_map(d_min = 3)
  mmm.add_map_manager_by_id(map_id = 'previous_map',
     map_manager = mmm.map_manager().deep_copy())

  # Supply the anisotropy of the map ourselves.  None is what
  #  _get_aniso_of_map returns when the anisotropic scaling fit failed
  aniso_b_cart_to_return = [None]
  def get_aniso_of_map(d_min = None, map_id = None):
    return aniso_b_cart_to_return[0]
  mmm._get_aniso_of_map = get_aniso_of_map

  def get_all_map_data():
    map_data_dict = {}
    for map_id in mmm.map_id_list():
      map_data_dict[map_id] = mmm.get_map_manager_by_id(map_id
         ).map_data().deep_copy()
    return map_data_dict

  def assert_all_maps_unchanged(map_data_dict):
    # Any anisotropy correction moves map values by an amount comparable to
    #  the variation in the map itself, so the largest change anywhere in a
    #  map, in units of that map's own standard deviation, is zero only if
    #  the map was left alone
    assert len(map_data_dict) > 1
    for map_id in map_data_dict.keys():
      original = map_data_dict[map_id].as_1d()
      current = mmm.get_map_manager_by_id(map_id).map_data().as_1d()
      assert current.size() == original.size()
      sd = original.sample_standard_deviation()
      assert sd > 0
      biggest_change = flex.max(flex.abs(current - original))
      assert approx_equal(biggest_change/sd, 0., eps = 1.e-6)

  # Not available, remove from all maps: no map is touched, the log says
  #  why, and the log does not claim that anything was removed
  aniso_b_cart_to_return[0] = None
  map_data_dict = get_all_map_data()
  f = StringIO()
  mmm.set_log(f)
  result = mmm.remove_anisotropy(d_min = 3, b_iso = 30,
     map_id = 'map_manager', remove_from_all_maps = True)
  mmm.set_log(None)
  assert result is None
  assert_all_maps_unchanged(map_data_dict)
  assert f.getvalue().find('Unable to determine overall anisotropy') > -1
  assert f.getvalue().find('Removed anisotropy from map') < 0

  # Not available, returning a map instead: nothing to return and again
  #  no map is touched
  map_data_dict = get_all_map_data()
  result = mmm.remove_anisotropy(d_min = 3, b_iso = 30,
     map_id = 'map_manager', remove_from_all_maps = False)
  assert result is None
  assert_all_maps_unchanged(map_data_dict)

  # Available, returning a map: a map_manager on the same grid comes back
  aniso_b_cart_to_return[0] = (10., 20., 30., 1., 2., 3.)
  result = mmm.remove_anisotropy(d_min = 3, b_iso = 30,
     map_id = 'map_manager', remove_from_all_maps = False)
  assert result is not None
  assert result.map_data().all() == mmm.map_manager().map_data().all()

def test_04():
  """A supplied aniso_b_cart is removed whether b_iso is given, None or 0"""

  from libtbx.test_utils import approx_equal
  from libtbx.utils import Sorry
  from scitbx.array_family import flex
  from six.moves import StringIO
  from cctbx.maptbx import segment_and_split_map
  from cctbx.maptbx.segment_and_split_map import map_coeffs_as_fp_phi
  from cctbx.maptbx.refine_sharpening import analyze_aniso_object

  d_min = 3
  supplied_b_cart = (10., 20., 30., 1., 2., 3.)

  # Two different maps: the generated map (the reference, map_id
  #  'map_manager') and a blurred copy, so that the two maps give different
  #  estimates of b_iso
  def get_mmm():
    mmm = map_model_manager()
    mmm.set_log(None)
    mmm.generate_map(d_min = d_min)
    blurred_map_coeffs = mmm.map_manager().map_as_fourier_coefficients(
       d_min = d_min).apply_debye_waller_factors(b_iso = 60)
    mmm.add_map_manager_by_id(map_id = 'blurred_map',
       map_manager = mmm.map_manager().fourier_coefficients_as_map_manager(
       blurred_map_coeffs))
    return mmm

  def get_b_iso_of_map(mm):
    f_array,phases=map_coeffs_as_fp_phi(mm.map_as_fourier_coefficients(
      d_min = d_min))
    return segment_and_split_map.get_b_iso(f_array, d_min = d_min)

  # The expected result, calculated directly: remove supplied_b_cart and
  #  leave b_iso behind
  def get_expected_map_data(mmm, mm, b_iso):
    map_coeffs = mm.map_as_fourier_coefficients(d_min = d_min)
    f_array,phases=map_coeffs_as_fp_phi(map_coeffs)
    analyze_aniso = analyze_aniso_object()
    analyze_aniso.b_cart = supplied_b_cart
    analyze_aniso.b_cart_aniso_removed = [ -b_iso, -b_iso, -b_iso, 0, 0, 0]
    scaled_f_array = analyze_aniso.apply_aniso_correction(f_array=f_array)
    return mmm.map_manager().fourier_coefficients_as_map_manager(
      scaled_f_array.phase_transfer(phase_source=phases, deg=True)
      ).map_data()

  def assert_same_map_data(current, expected):
    sd = expected.as_1d().sample_standard_deviation()
    assert sd > 0
    biggest_change = flex.max(flex.abs(current.as_1d() - expected.as_1d()))
    assert approx_equal(biggest_change/sd, 0., eps = 1.e-4)

  mmm = get_mmm()
  estimated_b_iso = get_b_iso_of_map(mmm.map_manager())
  blurred_b_iso = get_b_iso_of_map(mmm.get_map_manager_by_id('blurred_map'))
  # A b_iso estimated from each map separately would not match
  assert abs(estimated_b_iso - blurred_b_iso) > 10
  # Nor would 0 treated as missing, or a missing b_iso treated as 30
  assert abs(estimated_b_iso) > 10
  assert abs(estimated_b_iso - 30) > 10

  for b_iso, b_iso_used in ((30, 30), (None, estimated_b_iso), (0, 0)):

    # Remove from all maps: both maps get supplied_b_cart and the same b_iso
    mmm = get_mmm()
    expected = {}
    for map_id in ('map_manager', 'blurred_map'):
      expected[map_id] = get_expected_map_data(mmm,
         mmm.get_map_manager_by_id(map_id), b_iso_used)
    f = StringIO()
    mmm.set_log(f)
    result = mmm.remove_anisotropy(d_min = d_min,
       aniso_b_cart = supplied_b_cart, b_iso = b_iso,
       map_id = 'map_manager', remove_from_all_maps = True)
    mmm.set_log(None)
    assert tuple(result) == supplied_b_cart
    for map_id in ('map_manager', 'blurred_map'):
      assert_same_map_data(mmm.get_map_manager_by_id(map_id).map_data(),
        expected[map_id])
    if b_iso is None:
      assert f.getvalue().find(
        'b_iso not supplied; estimated from the reference map') > -1

    # Returning a map: the reference map with supplied_b_cart removed
    mmm = get_mmm()
    result = mmm.remove_anisotropy(d_min = d_min,
       aniso_b_cart = supplied_b_cart, b_iso = b_iso,
       map_id = 'map_manager', remove_from_all_maps = False)
    assert_same_map_data(result.map_data(),
      get_expected_map_data(mmm, mmm.map_manager(), b_iso_used))

  # b_iso cannot be estimated: Sorry is raised and no map is touched
  original_get_b_iso = segment_and_split_map.get_b_iso
  def failed_get_b_iso(miller_array, d_min = None,
      return_aniso_scale_and_b = False, d_max = 100000.):
    return 0., None
  segment_and_split_map.get_b_iso = failed_get_b_iso
  try:
    for remove_from_all_maps in (True, False):
      mmm = get_mmm()
      map_data_dict = {}
      for map_id in mmm.map_id_list():
        map_data_dict[map_id] = mmm.get_map_manager_by_id(map_id
           ).map_data().deep_copy()
      try:
        mmm.remove_anisotropy(d_min = d_min,
          aniso_b_cart = supplied_b_cart, b_iso = None,
          map_id = 'map_manager', remove_from_all_maps = remove_from_all_maps)
      except Sorry as e:
        assert str(e).find('Unable to estimate b_iso') > -1
      else:
        raise AssertionError("Expected Sorry when b_iso cannot be estimated")
      for map_id in mmm.map_id_list():
        assert mmm.get_map_manager_by_id(map_id).map_data().as_1d(
          ).all_eq(map_data_dict[map_id].as_1d())
  finally:
    segment_and_split_map.get_b_iso = original_get_b_iso


# ----------------------------------------------------------------------------

if (__name__ == '__main__'):
  test_02()
  test_03()
  test_04()
  if libtbx.env.find_in_repositories(relative_path='chem_data') is not None:
    test_01(method = 'model_sharpen')
  else:
    print('Skip test_01, chem_data not available')
  print ("OK")

