from __future__ import absolute_import, division, print_function
from iotbx.scalepack import merge
from iotbx.regression.utils import random_f_calc
from cctbx import miller
from cctbx.array_family import flex
from cctbx.development import debug_utils
from libtbx.test_utils import approx_equal
import sys

def recycle(miller_array):
  merge.write(file_name="tmp.sca", miller_array=miller_array)
  with open("tmp.sca") as f:
    read_back_file = merge.reader(file_handle=f)
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
  scale_factor = 9999998/flex.max(data)
  data = data * scale_factor + 1
  f_calc = miller.array(
    miller_set=f_calc,
    data=data,
    sigmas=data/10).set_observation_type_xray_intensity()
  f_calc = f_calc.select(flex.random_permutation(size=data.size()))
  recycle(miller_array=f_calc)
  recycle(miller_array=f_calc.f_sq_as_f())

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True):
    exercise(
      space_group_info,
      anomalous_flag=anomalous_flag,
      verbose=flags.Verbose)

def exercise_overloads():
  with open("overloads.sca", "w") as f:
    f.write("""\
    1
 -987
    50.000    50.000    80.000    90.000    90.000   120.000 p3121
  19   2   3 ******* 47482.6 ******* 26861.1
  19   2   2 16333.4 17143.9 38472.0 22574.8
  19   2   1 49448.5 24728.9 28427.6 18873.9
  19   2   0 34296.1 24479.1 25846.4 22660.4
  19   2  -1 27513.6 23318.8 30341.0 19273.1""")
  with open("overloads.sca") as f:
    file_in = merge.reader(file_handle=f)
  arrays = file_in.as_miller_arrays()
  assert (len(arrays[0].indices()) == 8)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  exercise_overloads()

if (__name__ == "__main__"):
  run()
