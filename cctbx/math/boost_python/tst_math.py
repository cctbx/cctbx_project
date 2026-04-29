from __future__ import absolute_import, division, print_function
from cctbx import math_module
from libtbx.test_utils import approx_equal
import math
from six.moves import range
import random

def exercise():
  n = 2**12
  cos_sin_table = math_module.ext.cos_sin_table(n)
  assert cos_sin_table.n_points() == n
  for i in range(n):
    u = i / float(n)
    x = u * 2 * math.pi
    assert approx_equal(cos_sin_table.get(u), complex(math.cos(x),math.sin(x)))
  print("OK")

def exercise2():
  n = 2**12
  cos_sin_table = math_module.ext.cos_sin_table(n, True)
  assert cos_sin_table.n_points() == n
  for i in range(n):
    x = random.uniform(0., 6.)
    u = x / (2 * math.pi)
    assert approx_equal(cos_sin_table.get(u), complex(math.cos(x),math.sin(x)))
  print("OK")

if (__name__ == "__main__"):
  exercise()
  exercise2()
