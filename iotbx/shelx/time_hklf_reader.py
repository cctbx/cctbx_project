from libtbx.utils import user_plus_sys_time
import os

from iotbx.shelx import hklf

def run(filename):
  print filename
  filename = os.path.expanduser(filename)

  print "Spirit reader"
  fast_timer = user_plus_sys_time()
  fast = hklf.fast_reader(filename=filename)
  t_fast = fast_timer.elapsed()
  print "%.2f s" %t_fast

  print "C++ simple reader"
  simple_timer = user_plus_sys_time()
  simple = hklf.simple_reader(filename=filename)
  t_simple = simple_timer.elapsed()
  print "%.2f s" %t_simple

  print "Python reader"
  slow_timer = user_plus_sys_time()
  slow = hklf.python_reader(filename=filename)
  t_slow = slow_timer.elapsed()
  print "%.2f s" %t_slow

  print "%s reflections" %fast.indices().size()
  print "Speed-up : ",
  try:
    print "%.2f" %(t_slow/t_fast)
  except ZeroDivisionError:
    print "> 1/(timer resolution)"
  print "-----------------"

  fast_array = fast.as_miller_arrays()[0]
  simple_array = simple.as_miller_arrays()[0]
  slow_array = slow.as_miller_arrays()[0]
  for a in (fast_array, simple_array):
    assert a.indices().all_eq(slow_array.indices())
    assert a.data().all_approx_equal(slow_array.data())
    assert a.sigmas().all_approx_equal(slow_array.sigmas())

if __name__ == '__main__':
  import sys
  for f in sys.argv[1:]:
    run(f)
  print "OK"
