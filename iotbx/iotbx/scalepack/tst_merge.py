from iotbx.scalepack import merge
from iotbx.regression.utils import random_f_calc
from cctbx import miller
from cctbx.array_family import flex
from cctbx.development import debug_utils
from libtbx.test_utils import approx_equal
import sys

def recycle(miller_array):
  merge.write(file_name="tmp.sca", miller_array=miller_array)
  read_back_file = merge.reader(file_handle=open("tmp.sca"))
  read_back_arrays = read_back_file.as_miller_arrays()
  assert len(read_back_arrays) == 1
  read_back_array = read_back_arrays[0]
  read_back_input_indexing = read_back_array.adopt_set(miller_array)
  if (miller_array.is_xray_amplitude_array()):
    read_back_input_indexing = read_back_input_indexing.f_sq_as_f()
  regression = flex.linear_regression(
    miller_array.data(),
    read_back_input_indexing.data())
  assert approx_equal(regression.slope(), 1, eps=1.e-3)
  assert abs(regression.y_intercept()) < 1
  regression = flex.linear_regression(
    miller_array.sigmas(),
    read_back_input_indexing.sigmas())
  if (miller_array.is_xray_intensity_array()):
    assert approx_equal(regression.slope(), 1, eps=1.e-3)
  else:
    assert approx_equal(regression.slope(), 1, eps=1.e-1)
  assert abs(regression.y_intercept()) < 1

def exercise(space_group_info, n_scatterers=8, d_min=2.5,
             anomalous_flag=False, verbose=0):
  f_calc = random_f_calc(
    space_group_info=space_group_info,
    n_scatterers=n_scatterers,
    d_min=d_min,
    anomalous_flag=anomalous_flag,
    verbose=verbose)
  if (f_calc is None): return
  data = flex.norm(f_calc.data())
  scale_factor = 9999999/flex.max(data)
  data = data * scale_factor
  f_calc = miller.array(
    miller_set=f_calc,
    data=data,
    sigmas=data/10).set_observation_type_xray_intensity()
  r = flex.random_double(size=data.size())
  p = flex.sort_permutation(r)
  f_calc = f_calc.apply_selection(p)
  recycle(miller_array=f_calc)
  recycle(miller_array=f_calc.f_sq_as_f())

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True):
    exercise(
      space_group_info,
      anomalous_flag=anomalous_flag,
      verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
