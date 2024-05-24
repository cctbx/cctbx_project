from __future__ import absolute_import, division, print_function
import libtbx.test_utils.parallel
from libtbx.utils import Sorry, Usage
import libtbx.phil
import random
import os
import sys

master_phil = libtbx.phil.parse("""
directory = None
  .type = path
  .multiple = True
module = None
  .type = str
  .multiple = True
script = None
  .type = path
  .multiple = True
nproc = 1
  .type=  int
shuffle = False
  .type = bool
quiet = False
  .type = bool
verbosity = 1
  .type = int
stderr = False
  .type = bool
skip_missing = False
  .type = bool
run_in_tmp_dir = False
  .type = bool
max_time = 180
  .type = float(value_min=0)
  .help = "Print warning and timing for all tests that take longer"
          "than max_time (in seconds) to run."
slow_tests = False
  .type = bool
  .help = "If True, also run any tests marked as slow, if any"
""")

def run(args,
   return_list_of_tests=None,
   python_keyword_text="",
   max_tests=None,
   start_test=None,
   tests_to_skip=None,
   tests_to_run=None,
   expected_failures_from_phenix_regression=[],
   unstables_from_phenix_regression = [],
   supplied_list_of_tests = None):

  if (len(args) == 0):
    raise Usage("""libtbx.run_tests_parallel [module=NAME] [directory=path]""")
  user_phil = []
  for arg in args :
    if os.path.isdir(arg):
      user_phil.append(libtbx.phil.parse("directory=%s" % arg))
    else :
      try :
        arg_phil = libtbx.phil.parse(arg)
      except RuntimeError :
        raise Sorry("Unrecognized argument '%s'" % arg)
      else :
        user_phil.append(arg_phil)

  params = master_phil.fetch(sources=user_phil).extract()
  if params.run_in_tmp_dir:
    from libtbx.test_utils import open_tmp_directory
    run_dir = open_tmp_directory()
    print('Running tests in %s' % run_dir)
    os.chdir(run_dir)
  elif return_list_of_tests:
    pass # don't need to check anything
  else:
    cwd = os.getcwd()
    cwd_files = os.listdir(cwd)
    if cwd_files and cwd_files != ["default.profraw"]:
      raise Sorry("Please run this program in an empty directory.")
  if (len(params.directory) == 0) and (len(params.module) == 0):
    raise Sorry("Please specify modules and/or directories to test.")
  all_tests = []
  expected_failure_list = []
  expected_unstable_list = []
  parallel_list = []
  if (not return_list_of_tests) and (supplied_list_of_tests is None): # (this
       #fails with return_list_of_tests)
    all_tests.extend(libtbx.test_utils.parallel.make_commands(params.script,
      python_keyword_text=python_keyword_text))
  for dir_name in params.directory :
    if os.path.split(dir_name)[-1].find("cctbx_project")>-1:
      print('DANGER '*10)
      print('Using the directory option in cctbx_project can be very time consuming')
      print('DANGER '*10)
    dir_tests = libtbx.test_utils.parallel.find_tests(dir_name)
    all_tests.extend(libtbx.test_utils.parallel.make_commands(dir_tests,
      python_keyword_text=python_keyword_text))
  for module_name in params.module :
    module_tests = libtbx.test_utils.parallel.get_module_tests(module_name,
       slow_tests = params.slow_tests,
       python_keyword_text=python_keyword_text,
        skip_missing = params.skip_missing)
    fail_tests = libtbx.test_utils.parallel.\
      get_module_expected_test_failures(module_name,
        skip_missing = params.skip_missing)
    unstable_tests = libtbx.test_utils.\
      parallel.get_module_expected_unstable_tests(module_name,
        skip_missing = params.skip_missing)
    parallel_tests = libtbx.test_utils.parallel.\
      get_module_parallel_tests(module_name,
        skip_missing = params.skip_missing)
    all_tests.extend(module_tests)
    all_tests.extend(fail_tests)
    all_tests.extend(unstable_tests)
    expected_failure_list.extend(fail_tests)
    expected_unstable_list.extend(unstable_tests)
    parallel_list.extend(parallel_tests)

  # add expected failures from phenix regression
  for ef in expected_failures_from_phenix_regression:
    for t in all_tests:
      if t.find(ef) > -1:
        expected_failure_list.append(t)

  # add unstables from phenix regression
  for u in unstables_from_phenix_regression:
    for t in all_tests:
      if t.find(u) > -1:
        expected_unstable_list.append(t)
  # Run all above to get expected failures etc

  if (supplied_list_of_tests is not None):
    all_tests = supplied_list_of_tests  # just use supplied tests

  # run only specified tests:
  if tests_to_run:
      new_tests=[]
      for t in all_tests:
        keep=False
        for tts in tests_to_run:
          if t.find(tts)>-1:
            keep=True
        if keep:
          print ("Keeping the test %s" %(t))
          new_tests.append(t)
        else:
          pass
      all_tests=new_tests

  # remove any specified tests:
  if tests_to_skip:
      new_tests=[]
      for t in all_tests:
        ok=True
        for tts in tests_to_skip:
          if t.find(tts)>-1:
            ok=False
        if ok:
          new_tests.append(t)
        else:
          print ("Skipping the test %s" %(t))
      all_tests=new_tests


  # check that test lists are unique
  seen = set()
  duplicates = set()
  for t in all_tests:
      if t in seen:
        duplicates.add(t)
      else:
        seen.add(t)
  assert len(duplicates) == 0, "Duplicate tests found.\n%s" % list(duplicates)
  if start_test:
    all_tests=all_tests[start_test:]
    print ("Starting with test # %s " %(start_test))
  if max_tests:
    all_tests=all_tests[:max_tests]
    print("Running only %s tests" %(max_tests))

  if return_list_of_tests:
    return all_tests
  if (len(all_tests) == 0):
    if (supplied_list_of_tests is not None):
      raise Sorry("No tests to run")
    else: # usual
      raise Sorry("No test scripts found in %s." % params.directory)
  if (params.shuffle):
    random.shuffle(all_tests)
  if (params.quiet):
    params.verbosity = 0
  with open("run_tests_parallel_zlog", "w") as log:
    result = libtbx.test_utils.parallel.run_command_list(
      cmd_list=all_tests,
      expected_failure_list=expected_failure_list,
      expected_unstable_list=expected_unstable_list,
      parallel_list=parallel_list,
      nprocs=params.nproc,
      log=log,
      verbosity=params.verbosity,
      max_time=params.max_time)
  print("\nSee run_tests_parallel_zlog for full output.\n")
  if (result.failure > 0):
    print("")
    print("*" * 80)
    print("ERROR: %d TEST FAILURES.  PLEASE FIX BEFORE COMMITTING CODE." % \
      result.failure)
    print("*" * 80)
    print("")
  return result.failure

if (__name__ == "__main__"):
  if (run(sys.argv[1:]) > 0):
    sys.exit(1)
