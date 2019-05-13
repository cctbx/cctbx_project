from __future__ import division
from __future__ import print_function
from mmtbx.utils import run_reduce_with_timeout
from time import time
from mmtbx.regression import model_1yjp
import sys

def exercise(prefix="tst_reduce_timeout_1"):
  fn = "%s.pdb" % prefix
  with open(fn, 'w') as f:
    f.write(model_1yjp)
  t0 = time()
  rr = run_reduce_with_timeout(
      stdin_lines=None,
      file_name=fn,
      parameters="-oh -his -flip -keep -allalt -pen9999",
      override_auto_timeout_with=None)
  t1 = time()
  assert rr.return_code == 0
  stdout_lines_1 = rr.stdout_lines
  rr = run_reduce_with_timeout(
      stdin_lines=model_1yjp,
      file_name=None,
      parameters="-oh -his -flip -keep -allalt -pen9999 -",
      override_auto_timeout_with=None)
  t1 = time()
  stdout_lines_2 = rr.stdout_lines
  assert rr.return_code == 0
  assert len(stdout_lines_1) > 0
  assert len(stdout_lines_1) == len(stdout_lines_2)
  for l1, l2 in zip(stdout_lines_1, stdout_lines_2):
    assert l1 == l2
  assert 'ATOM      0  HB3 GLN A   5       3.844   2.489   7.623  1.00 11.96           H   new' in stdout_lines_1
  assert 'ATOM      0  HB3 GLN A   5       3.844   2.489   7.623  1.00 11.96           H   new' in stdout_lines_2
  t0 = time()
  rr = run_reduce_with_timeout(
      stdin_lines=None,
      file_name=fn,
      parameters="-oh -his -flip -keep -allalt -pen9999",
      override_auto_timeout_with=0.01)
  t1 = time()
  assert t1-t0 < 0.1
  assert rr.return_code == -15

  t0 = time()
  rr = run_reduce_with_timeout(
      stdin_lines=model_1yjp,
      file_name=None,
      parameters="-oh -his -flip -keep -allalt -pen9999 -",
      override_auto_timeout_with=0.01)
  t1 = time()
  assert t1-t0 < 0.1
  assert rr.return_code == -15

  print("OK")

if __name__ == "__main__" :
  if sys.platform != 'win32':
    exercise()
