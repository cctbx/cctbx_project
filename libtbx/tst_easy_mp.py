from libtbx import unpicklable
from libtbx.test_utils import Exception_expected
from libtbx import Auto
from cStringIO import StringIO
import sys

class potentially_large(unpicklable):

  def __init__(self, size):
    self.array = range(3, size+3)

  def __call__(self, i):
    return self.array[i]

def eval_parallel(
      data,
      func_wrapper=Auto,
      log=None,
      exercise_out_of_range=False,
      exercise_fail=False):
  size = len(data.array)
  args = range(size)
  if (exercise_out_of_range):
    args.append(size)
  from libtbx import easy_mp
  if (exercise_fail):
    mp_results = easy_mp.pool_map(func=data, args=args)
  else:
    mp_results = easy_mp.pool_map(
      fixed_func=data, args=args, func_wrapper=func_wrapper, log=log,
      buffer_stdout_stderr=exercise_out_of_range)
  if (not exercise_out_of_range):
    assert mp_results == range(3, size+3)
  else:
    assert mp_results[:size] == zip([""]*size, range(3, size+3))
    assert mp_results[size][0].startswith("CAUGHT EXCEPTION:")
    assert mp_results[size][0].find("IndexError: ") > 0
    assert mp_results[size][1] is None

def exercise(exercise_fail):
  from libtbx import easy_mp
  mp_problem = easy_mp.detect_problem()
  if (mp_problem is not None):
    print "Skipping tst_easy_mp.py: %s" % mp_problem
    return
  data = potentially_large(size=1000)
  eval_parallel(data)
  assert len(easy_mp.fixed_func_registry) == 0
  eval_parallel(data, func_wrapper=None)
  assert len(easy_mp.fixed_func_registry) == 1
  eval_parallel(data, func_wrapper=None, log=sys.stdout)
  assert len(easy_mp.fixed_func_registry) == 2
  sio = StringIO()
  eval_parallel(data, func_wrapper=None, log=sio)
  assert len(easy_mp.fixed_func_registry) == 3
  lines = sio.getvalue().splitlines()
  assert len(lines) == 2
  assert lines[0].startswith("multiprocessing pool size: ")
  assert lines[1].startswith("wall clock time: ")
  eval_parallel(data, exercise_out_of_range=True)
  if (exercise_fail):
    eval_parallel(data, exercise_fail=True)
    raise Exception_expected
  del data
  assert len(easy_mp.fixed_func_registry) == 0

def run(args):
  assert args in [[], ["--fail"]]
  exercise_fail = (len(args) != 0)
  exercise(exercise_fail)
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
