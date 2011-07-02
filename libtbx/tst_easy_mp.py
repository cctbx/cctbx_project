from libtbx import unpicklable
from libtbx.test_utils import Exception_expected
import sys

class potentially_large(unpicklable):

  def __init__(self, size):
    self.array = range(3, size+3)

  def __call__(self, i):
    return self.array[i]

def eval_parallel(data, exercise_fail=False):
  size = len(data.array)
  import multiprocessing
  n_proc = min(size, multiprocessing.cpu_count())
  if (exercise_fail):
    mp_pool = multiprocessing.Pool(processes=n_proc)
    mp_results = mp_pool.map(data, range(size))
  else:
    from libtbx import easy_mp
    mp_pool = easy_mp.Pool(processes=n_proc, fixed_func=data)
    mp_results = mp_pool.map_fixed_func(range(size))
  assert mp_results == range(3, size+3)

def exercise(exercise_fail):
  import os
  if (os.name == "nt"):
    print "Skipping tst_easy_mp.py: Windows is not a supported platform."
    return
  vers_info = sys.version_info[:2]
  if (vers_info < (2,6)):
    print "Skipping tst_easy_mp.py: Python 2.6 or higher is required."
    return
  data = potentially_large(size=1000)
  eval_parallel(data)
  from libtbx import easy_mp
  assert len(easy_mp.fixed_func_registry) == 1
  eval_parallel(data)
  assert len(easy_mp.fixed_func_registry) == 2
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
