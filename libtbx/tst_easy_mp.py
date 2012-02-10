from libtbx import unpicklable
from libtbx.test_utils import Exception_expected
from libtbx import Auto
from cStringIO import StringIO
import sys

def exercise_func_wrapper_sub_directories():
  from libtbx.easy_mp import func_wrapper_sub_directories as f
  w = f("")
  assert w.sub_name_format == "%03d"
  w = f("x")
  assert w.sub_name_format == "x%03d"
  w = f("%05d")
  assert w.sub_name_format == "%05d"
  w = f("#")
  assert w.sub_name_format == "%01d"
  w = f("##")
  assert w.sub_name_format == "%02d"
  w = f("y####")
  assert w.sub_name_format == "y%04d"
  w = f("#z###")
  assert w.sub_name_format == "#z###%03d"

class potentially_large(unpicklable):

  def __init__(self, size):
    self.array = range(3, size+3)

  def __call__(self, i):
    return self.array[i]

def eval_parallel(
      data,
      func_wrapper="simple",
      index_args=True,
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
    if (func_wrapper == "simple" and exercise_out_of_range):
      func_wrapper = "buffer_stdout_stderr"
    mp_results = easy_mp.pool_map(
      fixed_func=data,
      args=args,
      func_wrapper=func_wrapper,
      index_args=index_args,
      log=log)
  if (not exercise_out_of_range):
    assert mp_results == range(3, size+3)
  else:
    assert mp_results[:size] == zip([""]*size, range(3, size+3))
    assert mp_results[size][0].startswith("CAUGHT EXCEPTION:")
    assert mp_results[size][0].find("IndexError: ") > 0
    assert mp_results[size][1] is None

def exercise(exercise_fail):
  exercise_func_wrapper_sub_directories()
  from libtbx import easy_mp
  mp_problem = easy_mp.detect_problem()
  if (mp_problem is not None):
    print "Skipping tst_easy_mp.py: %s" % mp_problem
    return
  data = potentially_large(size=1000)
  eval_parallel(data)
  assert len(easy_mp.fixed_func_registry) == 0
  eval_parallel(data, func_wrapper=None, index_args=False)
  assert len(easy_mp.fixed_func_registry) == 1
  eval_parallel(data, func_wrapper=None, index_args=False, log=sys.stdout)
  assert len(easy_mp.fixed_func_registry) == 2
  sio = StringIO()
  eval_parallel(data, func_wrapper=None, index_args=False, log=sio)
  assert len(easy_mp.fixed_func_registry) == 3
  lines = sio.getvalue().splitlines()
  assert len(lines) == 2
  assert lines[0].startswith("multiprocessing pool size: ")
  assert lines[1].startswith("wall clock time: ")
  eval_parallel(data, exercise_out_of_range=True)
  if (exercise_fail):
    eval_parallel(data, exercise_fail=True)
    raise Exception_expected
  results = easy_mp.pool_map(fixed_func=data, args=range(1000), processes=Auto)
  del data
  assert len(easy_mp.fixed_func_registry) == 0
  #
  from libtbx.clear_paths import \
    remove_or_rename_files_and_directories_if_possible as clear
  from libtbx import only_element
  import os
  op = os.path
  def fixed_func(arg):
    print "hello world", arg
    return 10*arg
  def go():
    return easy_mp.pool_map(
      fixed_func=fixed_func,
      func_wrapper="sub_directories",
      args=[1,2])
  clear(paths=["mp000", "mp001"])
  results = go()
  assert results == [(None, 10), (None, 20)]
  for i in [1,2]:
    only_element(open("mp%03d/log" % (i-1)).read().splitlines()) \
      == "hello world %d" % i
  results = go()
  assert results == [
    ('sub-directory exists already: "mp000"', None),
    ('sub-directory exists already: "mp001"', None)]
  clear(paths=["mp001"])
  results = go()
  assert results == [
    ('sub-directory exists already: "mp000"', None),
    (None, 20)]
  clear(paths=["mp000", "mp001"])
  results = easy_mp.pool_map(
    fixed_func=fixed_func,
    func_wrapper=easy_mp.func_wrapper_sub_directories(makedirs_mode=0),
    args=[1,2])
  assert results == [
    ('cannot chdir to sub-directory: "mp000"', None),
    ('cannot chdir to sub-directory: "mp001"', None)]
  clear(paths=["mp000", "mp001"])
  clear(paths=["bf000", "bf001"])
  def bad_func(arg):
    raise RuntimeError(str(arg))
  results = easy_mp.pool_map(
    fixed_func=bad_func,
    func_wrapper="sub_directories:bf",
    args=[1,2])
  assert results == [
    ('CAUGHT EXCEPTION: "bf000/log"', None),
    ('CAUGHT EXCEPTION: "bf001/log"', None)]
  for i in [1,2]:
    assert open("bf%03d/log" % (i-1)).read().splitlines()[-1] \
      == "RuntimeError: %d" % i

def run(args):
  assert args in [[], ["--fail"]]
  exercise_fail = (len(args) != 0)
  exercise(exercise_fail)
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
