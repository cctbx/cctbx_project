from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from six.moves import cStringIO as StringIO
from iotbx.map_model_manager import map_model_manager

def sharpened_map_data(mmm, model_id):
  mmm.model_sharpen(model_id = model_id, n_bins = 10)
  return mmm.get_map_manager_by_id('map_manager_scaled').map_data()

def exercise():
  '''
  model_sharpen(model_id = 'other') boxes, masks and finds k_sol and b_sol
  around the model 'other', so it gives the same map as the same model
  stored as 'model'.  Checked with a different model stored as 'model'
  (residues 219-221), and with no model 'model'.  Sharpening leaves the
  models it was given unchanged.  (Fix 398a7df82b.)
  '''
  start_mmm = map_model_manager()
  start_mmm.generate_map(d_min = 4, wrapping = False)
  start_mmm.set_log(StringIO())
  start_mmm.set_resolution(4)
  full_model = start_mmm.model().deep_copy()

  reference = sharpened_map_data(start_mmm.deep_copy(), 'model')
  tolerance = 1.e-6 * reference.standard_deviation_of_the_sample()

  # A fragment stored as 'model' must not be used
  mmm = start_mmm.deep_copy()
  mmm.add_model_by_id(model = full_model.deep_copy(), model_id = 'other')
  mmm.add_model_by_id(model = full_model.apply_selection_string(
     'resseq 219:221'), model_id = 'model', overwrite = True)
  assert mmm.get_model_by_id('model').get_number_of_atoms() > 0
  fragment_sites = mmm.get_model_by_id('model').get_sites_cart().deep_copy()
  with_fragment = sharpened_map_data(mmm, 'other')
  d = flex.max(flex.abs(with_fragment - reference))
  assert d < tolerance, d
  # The models given to sharpening are unchanged
  for model, sites in ((mmm.get_model_by_id('other'),
      full_model.get_sites_cart()),
      (mmm.get_model_by_id('model'), fragment_sites)):
    assert model.get_sites_cart().rms_difference(sites) == 0

  # No model 'model' at all
  mmm = start_mmm.deep_copy()
  mmm.add_model_by_id(model = full_model.deep_copy(), model_id = 'other')
  mmm.remove_model_by_id(model_id = 'model')
  without_model = sharpened_map_data(mmm, 'other')
  d = flex.max(flex.abs(without_model - reference))
  assert d < tolerance, d

if __name__ == "__main__":
  exercise()
  print("OK")
