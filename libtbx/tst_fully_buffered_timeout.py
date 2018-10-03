from __future__ import division
from libtbx.easy_run import fully_buffered
from time import time
from libtbx.test_utils import approx_equal
import sys

def test_command(cmd, to, expected_time):
  t0 = time()
  fb = fully_buffered(command=cmd, timeout=to)
  t1 = time()
  if to > expected_time:
    assert fb.return_code == 0
  else:
    assert fb.return_code == -15
  assert approx_equal(t1-t0, expected_time, 0.5)

def exercise():
  test_command("sleep 1", 2, 1)
  test_command("sleep 1000", 2, 2)
  test_command("sleep 5", 3, 3)
  test_command("sleep 5", 10, 5)
  print "OK"

if __name__ == "__main__":
  if sys.platform != 'win32':
    exercise()
