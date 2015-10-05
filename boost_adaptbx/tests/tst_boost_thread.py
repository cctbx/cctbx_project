from __future__ import division
import boost.python
ext = boost.python.import_ext("boost_adaptbx_boost_thread_test_ext")

from libtbx.test_utils import approx_equal
from stdlib import math

def sanity_test():
  (pi, e) = ext.test_boost_thread()
  assert approx_equal(pi, math.pi, eps=5), pi
  assert approx_equal(e, math.exp(1), eps=5), e

def run():
  sanity_test()
  print 'OK'

if __name__ == '__main__':
  run()
