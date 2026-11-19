from __future__ import absolute_import, division, print_function
import os, shutil, tempfile
from scitbx.array_family import flex
from six.moves import cStringIO as StringIO
from iotbx.map_model_manager import map_model_manager

def get_mmm(d_min = 3):
  mmm = map_model_manager()
  mmm.generate_map(d_min = d_min, wrapping = False)
  mmm.set_log(StringIO())
  mmm.set_resolution(d_min)
  return mmm

def sharpened_map_data(mmm, model_id, local_sharpen, k_sol, b_sol):
  if local_sharpen:
    mmm.model_sharpen(model_id = model_id, n_bins = 10, local_sharpen = True,
      n_boxes = 8, nproc = 1, k_sol = k_sol, b_sol = b_sol)
  else:
    mmm.model_sharpen(model_id = model_id, n_bins = 10,
      k_sol = k_sol, b_sol = b_sol)
  return mmm.get_map_manager_by_id('map_manager_scaled').map_data()

def exercise(local_sharpen, k_sol = None, b_sol = None):
  '''
  model_sharpen(model_id = 'other') boxes, masks and finds k_sol and b_sol
  around the model 'other', so it gives the same map as the same model
  stored as 'model', with global or local sharpening, and with k_sol and
  b_sol supplied (masking only).  Checked with a different model stored as
  'model' (residues 219-221), and with no model 'model'.  Sharpening leaves
  the models it was given unchanged.
  '''
  reference_mmm = get_mmm()
  full_model = reference_mmm.model().deep_copy()
  reference = sharpened_map_data(reference_mmm, 'model', local_sharpen,
    k_sol, b_sol)

  # A fragment stored as 'model' must not be used
  mmm = get_mmm()
  mmm.add_model_by_id(model = full_model.deep_copy(), model_id = 'other')
  mmm.add_model_by_id(model = full_model.apply_selection_string(
     'resseq 219:221'), model_id = 'model', overwrite = True)
  assert mmm.get_model_by_id('model').get_number_of_atoms() > 0
  fragment_sites = mmm.get_model_by_id('model').get_sites_cart().deep_copy()
  with_fragment = sharpened_map_data(mmm, 'other', local_sharpen,
    k_sol, b_sol)
  d = flex.max(flex.abs(with_fragment - reference))
  print("local_sharpen %s, k_sol %s, fragment stored as 'model': max "
    "difference %.3g (map rms %.3g)" %(local_sharpen, k_sol, d,
    reference.standard_deviation_of_the_sample()))
  assert d < 1.e-6 * reference.standard_deviation_of_the_sample(), d
  # The models given to sharpening are unchanged
  for model, sites in ((mmm.get_model_by_id('other'),
      full_model.get_sites_cart()),
      (mmm.get_model_by_id('model'), fragment_sites)):
    assert model.get_sites_cart().rms_difference(sites) == 0

  # No model 'model' at all
  mmm = get_mmm()
  mmm.add_model_by_id(model = full_model.deep_copy(), model_id = 'other')
  mmm.remove_model_by_id(model_id = 'model')
  without_model = sharpened_map_data(mmm, 'other', local_sharpen,
    k_sol, b_sol)
  d = flex.max(flex.abs(without_model - reference))
  print("local_sharpen %s, k_sol %s, no model 'model': max difference %.3g" %(
    local_sharpen, k_sol, d))
  assert d < 1.e-6 * reference.standard_deviation_of_the_sample(), d

def run_in_temp_dir(f):
  '''Run f() in a new temporary directory (library code may write files
     to the current directory), then remove it'''
  cwd = os.getcwd()
  tmp_dir = tempfile.mkdtemp()
  try:
    os.chdir(tmp_dir)
    return f()
  finally:
    os.chdir(cwd)
    shutil.rmtree(tmp_dir, ignore_errors = True)

if __name__ == "__main__":
  run_in_temp_dir(lambda: exercise(local_sharpen = False))
  run_in_temp_dir(lambda: exercise(local_sharpen = True))
  run_in_temp_dir(lambda: exercise(local_sharpen = False,
    k_sol = 0.35, b_sol = 46.))
  print("OK")
