from libtbx.utils import user_plus_sys_time
from libtbx.test_utils import approx_equal

from iotbx.shelx import hklf

def run(filename):
  print filename
  print "C++ reader"
  fast_timer = user_plus_sys_time()
  fast = hklf.reader(filename=filename)
  t_fast = fast_timer.elapsed()
  print "%.2f s" %t_fast

  print "Python reader"
  slow_timer = user_plus_sys_time()
  slow = hklf.python_reader(filename=filename)
  t_slow = slow_timer.elapsed()
  print "%.2f s" %t_slow

  print "%s reflections" %fast.indices().size()
  print "Speed-up : %.2f" %(t_slow/t_fast)
  print "-----------------"

  fast_array = fast.as_miller_arrays()[0]
  slow_array = slow.as_miller_arrays()[0]
  assert fast_array.indices().all_eq(slow_array.indices())
  assert fast_array.data().all_approx_equal(slow_array.data())
  assert fast_array.sigmas().all_approx_equal(slow_array.sigmas())

if __name__ == '__main__':
  import sys
  for f in sys.argv[1:]:
    run(f)
  print "OK"
  #run(sys.argv[1])
