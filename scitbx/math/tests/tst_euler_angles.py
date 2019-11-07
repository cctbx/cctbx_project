from __future__ import absolute_import, division, print_function
from scitbx.math import euler_angles as euler
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import random
from six.moves import range

def exercise_core(angles_in):
  for euler_matrix, euler_angles in [
        (euler.xyz_matrix, euler.xyz_angles),
        (euler.yzx_matrix, euler.yzx_angles),
        (euler.zyz_matrix, euler.zyz_angles)]:
    m = euler_matrix(*angles_in)
    angles_out = euler_angles(m)
    m2 = euler_matrix(*angles_out)
    assert approx_equal(m2, m)

def exercise():
  random.seed(0)
  for a1_in in range(0, 400, 15):
    for a2_in in range(0, 400, 15):
      for a3_in in range(0, 400, 15):
        exercise_core((a1_in, a2_in, a3_in))
  for i_trial in range(1000):
    exercise_core([random.random()*360-180 for i in [0,1,2]])
  print(format_cpu_times())

if (__name__ == "__main__"):
  exercise()
