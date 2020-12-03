from __future__ import absolute_import, division, print_function
from libtbx import unpicklable
from libtbx.test_utils import Exception_expected
from libtbx import Auto
from six.moves import cStringIO as StringIO
import sys
from six.moves import range, zip

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
    self.array = list(range(3, size+3))

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
  args = list(range(size))
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
    assert mp_results == list(range(3, size+3))
  else:
    assert mp_results[:size] == list(zip([""]*size, range(3, size+3)))
    assert mp_results[size][0].startswith("CAUGHT EXCEPTION:")
    assert mp_results[size][0].find("IndexError: ") > 0
    assert mp_results[size][1] is None

def exercise(exercise_fail):
  exercise_func_wrapper_sub_directories()
  from libtbx import easy_mp
  mp_problem = easy_mp.detect_problem()
  if (mp_problem is not None):
    print("Skipping tst_easy_mp.py: %s" % mp_problem)
    return
  check_if_stacktrace_is_propagated_properly(method='threading', nproc=2)
  check_if_stacktrace_is_propagated_properly(method='multiprocessing', nproc=2)
  check_if_stacktrace_is_propagated_properly(method='threading', nproc=1)
  check_if_stacktrace_is_propagated_properly(method='multiprocessing', nproc=1)
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
  results = easy_mp.pool_map(fixed_func=data, args=list(range(1000)), processes=Auto)
  del data
  assert len(easy_mp.fixed_func_registry) == 0
  #
  from libtbx.clear_paths import \
    remove_or_rename_files_and_directories_if_possible as clear
  import os
  op = os.path
  def fixed_func(arg):
    print("hello world", arg)
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
    with open("mp%03d/log" % (i-1), 'r') as fh:
      assert fh.read().splitlines() == ["hello world %d" % i]
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
    with open("bf%03d/log" % (i-1)) as f:
      assert f.read().splitlines()[-1] == "RuntimeError: %d" % i
  out = StringIO()
  def simple_func(arg):
    from math import log
    x = float(arg)
    y = 0
    for i in range(2, 1000):
      y += log(x) / log(float(i))
    return y
  def cb(result):
    out.write("%.3f\n" % result)
  result = easy_mp.pool_map(
    func=simple_func,
    args=[1,2,3,4,5,6,7,8],
    call_back_for_serial_run=cb,
    processes=1)
  assert (out.getvalue() == """\
0.000
122.891
194.777
245.782
285.344
317.668
344.998
368.673
""")

def _may_divide_by_zero(divideby):
  '''
  A helper function for the check_if_stacktrace_is_propagated_properly test.
  Must be on module level for parallelization on Windows machines.
  '''
  return 7 / divideby

def check_if_stacktrace_is_propagated_properly(method, nproc):
  exception_seen = False
  from libtbx.easy_mp import parallel_map
  import traceback

  try:
    results = parallel_map(
      func=_may_divide_by_zero,
      iterable=[2,1,0],
      method=method,
      processes=nproc,
      preserve_exception_message=True)
  except ZeroDivisionError as e:
    exception_seen = True
    exc_type, exc_value, exc_traceback = sys.exc_info()
    assert "division by zero" in str(exc_value), "Exception value mismatch: '%s'" % exc_value

    stack_contains_fail_function = False
    # Two options: Either the original stack is available directly
    for (filename, line, function, text) in traceback.extract_tb(exc_traceback):
      if function == _may_divide_by_zero.__name__:
        stack_contains_fail_function = True
    # or it should be preserved in the string representation of the exception
    from libtbx.scheduling import stacktrace
    ex, st = stacktrace.exc_info()
    if ex is not None and _may_divide_by_zero.__name__ in "".join( st ):
      stack_contains_fail_function = True
    if not stack_contains_fail_function:
      print("Thrown exception: %s:" % str(e))
      traceback.print_tb(exc_traceback)
      print("")
      assert stack_contains_fail_function, "Stacktrace lost"
  except Exception as e:
    print("Exception type mismatch, expected ZeroDivisionError")
    raise
  assert exception_seen, "Expected exception not thrown"

def run(args):
  assert args in [[], ["--fail"]]
  exercise_fail = (len(args) != 0)
  exercise(exercise_fail)
  print("OK")

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
