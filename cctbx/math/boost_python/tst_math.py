from cctbx import math_module
from scitbx.test_utils import approx_equal
import math

def exercise():
  n = 2**12
  cos_sin_table = math_module.ext.cos_sin_table(n)
  assert cos_sin_table.n_points() == n
  for i in xrange(n):
    u = i / float(n)
    x = u * 2 * math.pi
    assert approx_equal(cos_sin_table.get(u), complex(math.cos(x),math.sin(x)))

if (__name__ == "__main__"):
  exercise()
