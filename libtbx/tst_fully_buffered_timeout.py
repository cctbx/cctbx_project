from __future__ import division
from libtbx.easy_run import fully_buffered
from time import time
from libtbx.test_utils import approx_equal

def test_command(cmd, to, expected_time):
  t0 = time()
  fully_buffered(command=cmd, timeout=to)
  t1 = time()
  assert approx_equal(t1-t0, expected_time, 0.1)

def exercise():
  test_command("sleep 1", 2, 1)
  test_command("sleep 1000", 2, 2)
  test_command("sleep 5", 4, 4)
  test_command("sleep 5", 10, 5)
  print "OK"

if __name__ == "__main__" :
  exercise()
