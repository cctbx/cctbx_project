from __future__ import division, absolute_import
from __future__ import print_function
import boost.python
ext = boost.python.import_ext("boost_adaptbx_boost_thread_test_ext")

from libtbx.test_utils import approx_equal
import math

def sanity_test():
  (pi, e) = ext.test_boost_thread()
  assert approx_equal(pi, math.pi, eps=5), pi
  assert approx_equal(e, math.exp(1), eps=5), e

def run():
  sanity_test()
  print('OK')

if __name__ == '__main__':
  run()
