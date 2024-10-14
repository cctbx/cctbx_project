from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from libtbx import test_utils
import libtbx.load_env
from libtbx.utils import import_python_object, multi_out
from libtbx import group_args, Auto
from multiprocessing import Pool, cpu_count
from six.moves import cStringIO as StringIO
import traceback
import time
import os
import re
import sys
import re
import codecs

try:
  from enum import IntEnum
except ImportError:
  IntEnum = object

class Status(IntEnum):
  OK = 0
  WARNING = 1
  FAIL = 2
  EXPECTED_FAIL = 3
  EXPECTED_UNSTABLE = 4

QUIET = 0
DEFAULT_VERBOSITY = 1
EXTRA_VERBOSE = 2 # for nightly builds

def get_test_list(module_name, test_type='tst_list', valgrind=False,
         python_keyword_text=None, skip_missing = False):
  dist_path = libtbx.env.dist_path(module_name)
  if dist_path is None:
    if skip_missing:
      return []
    else:
      raise AssertionError("'%s' is not a valid CCTBX module." % module_name)
  # XXX don't check for file name, because import conventions differ among
  # derived modules - dist_path can be either the sys.path entry or the actual
  # module contents.  If the file is missing the import below will fail, which
  # is okay for testing purposes.
  try:
    tst_list = list(import_python_object(
      import_path="%s.run_tests.%s" % (module_name, test_type),
      error_prefix="",
      target_must_be="",
      where_str="").object)
  except AttributeError:
    tst_list = []
  build_path = libtbx.env.under_build(module_name)
  assert (build_path is not None) and (dist_path is not None)
  commands = []
  co = group_args(
    verbose=False,
    quick=True,
    valgrind=valgrind,
    python_keyword_text=python_keyword_text)
  for cmd in test_utils.iter_tests_cmd(
      co=co,
      build_dir=build_path,
      dist_dir=dist_path,
      tst_list=tst_list):
    commands.append(cmd)
  return commands

def get_module_tests(module_name, valgrind=False, slow_tests=False,
      python_keyword_text=None, skip_missing = False):
  tst_list = get_test_list(module_name, valgrind=valgrind,
     python_keyword_text=python_keyword_text, skip_missing = skip_missing)
  if slow_tests:
    tst_list.extend(get_test_list(module_name, test_type='tst_list_slow',
                    valgrind=valgrind,
                    python_keyword_text=python_keyword_text,
                    skip_missing = skip_missing))
  return tst_list

def get_module_expected_test_failures(module_name, skip_missing = False):
  return get_test_list(module_name, test_type='tst_list_expected_failures',
                    skip_missing = skip_missing)

def get_module_expected_unstable_tests(module_name, skip_missing = False):
  return get_test_list(module_name, test_type='tst_list_expected_unstable',
                    skip_missing = skip_missing)

def get_module_parallel_tests(module_name, skip_missing = False):
  return get_test_list(module_name, test_type='tst_list_parallel',
                     skip_missing = skip_missing)

def find_tests(dir_name):
  if not os.path.isdir(dir_name):
    raise RuntimeError("'%s' is not a directory." % dir_name)
  all_tests = []
  for root, dirnames, filenames in os.walk(dir_name):
    for file_name in filenames :
      base, ext = os.path.splitext(file_name)
      if (file_name.startswith("tst")) and (ext in [".py", ".csh", ".sh"]):
        all_tests.append(os.path.join(root, file_name))
  return all_tests

def run_command(command,
                verbosity=DEFAULT_VERBOSITY):
  try:
    t0=time.time()
    cmd_result = easy_run.fully_buffered(
      command=command,
      )
    cmd_result.wall_time = time.time()-t0
    cmd_result.error_lines = phenix_separate_output(cmd_result)
    return cmd_result
  except KeyboardInterrupt:
    traceback_str = "\n".join(traceback.format_tb(sys.exc_info()[2]))
    sys.stdout.write(traceback_str)
    return None

def phenix_separate_output(cmd_result):
  """Separates warning and error lines from a commands output streams.

  This will only give results if the phenix_regression package exists,
  and returns an empty list otherwise.
  """
  try:
    from phenix_regression.command_line import find_errors_and_warnings
  except ImportError:
    return []
  out = cmd_result.stdout_lines
  err = cmd_result.stderr_lines
  bad_lines = []
  for i, line in enumerate(out+err):
    if find_errors_and_warnings.is_error_or_warning(line, line.lower()):
      bad_lines.append(line)
  return bad_lines

def reconstruct_test_name(command):
  if hasattr(command, 'test_class') and hasattr(command, 'test_name'):
    return command.test_class, command.test_name

  pattern = '^[^"]*"([^"]*)"([^"]*)$'
  m = re.search(pattern, command)
  # command =  libtbx.python "/file.py" 90000
  #            \-- ignore --/ \--m1--/ \-m2-/
  if m:
    file = m.group(1)
    parameter = m.group(2).lstrip()
    filtered_parameter = re.sub(r"[^A-Za-z0-9\-_]", '', parameter)
    if filtered_parameter == '':
      test_name = file
    else:
      test_name = "%s.%s" % (file, filtered_parameter)

    pattern2 = r'^/?(.*?/((modules|build)/(cctbx\_project/|xia2/Test/)?))?(.*)/([^/&]*?)(\.py)?$'
    m2 = re.search(pattern2, file)
    if m2:
#     print "M (%s) (%s) (%s) (%s) (%s) (%s) (%s)" % (m2.group(1,2,3,4,5,6,7))
      test_name = m2.group(6).replace('/', '.')
      test_class = m2.group(5).replace('/', '.')
      is_python_dot_identifier = r"^([^\d\W]\w*)\.([^\d\W]\w*)\Z"
      if re.search(is_python_dot_identifier, parameter):
        # if parameter consists of two joined python identifiers, use it as test name
        test_class = "%s.%s" % (test_class, test_name)
        test_name = parameter
      elif filtered_parameter != '':
        # otherwise append sanitized (filtered) parameter to test name
        # so that each test has again a unique name
        test_name = "%s.%s" % (test_name, filtered_parameter)
      return (test_class, test_name)
    return (file, test_name)
  return (command, command)

def write_JUnit_XML(results, output_filename="output.xml"):
  """Write a JUnit XML test report to a file if junit_xml is available."""
  try:
    import junit_xml
  except ImportError:
    return

  test_cases = []
  for result in results:
    test_name = reconstruct_test_name(result.command)
    tc = junit_xml.TestCase(classname=test_name[0],
                  name=test_name[1],
                  elapsed_sec=result.wall_time,
                  stdout='\n'.join(result.stdout_lines),
                  stderr='\n'.join(result.stderr_lines))
    if result.return_code == 0:
      # Identify skipped tests
      output = '\n'.join(result.stdout_lines + result.stderr_lines)
      if re.search('skip', output, re.IGNORECASE):
        # find first line including word 'skip' and use it as message
        skipline = re.search('^((.*)skip(.*))$', output, re.IGNORECASE | re.MULTILINE).group(1)
        tc.add_skipped_info(skipline)
    elif result.alert_status == Status.EXPECTED_FAIL:
      tc.add_skipped_info("Expected test failure")
    elif result.alert_status == Status.EXPECTED_UNSTABLE:
      tc.add_skipped_info("Expected test instability")
    else:
      # Test failed. Extract error message and stack trace if possible
      error_message = 'exit code %d' % result.return_code
      error_output = '\n'.join(result.stderr_lines)
      if result.stderr_lines:
        error_message = result.stderr_lines[-1]
        if len(result.stderr_lines) > 20:
          error_output = '\n'.join(result.stderr_lines[-20:])
      tc.add_failure_info(message=error_message, output=error_output)
    test_cases.append(tc)
  ts = junit_xml.TestSuite("libtbx.run_tests_parallel", test_cases=test_cases)
  with codecs.open(output_filename, "w", encoding="utf-8") as f:
    ts.to_file(f, [ts], prettyprint=True, encoding="utf-8")

class run_command_list(object):
  def __init__(self,
                cmd_list,
                expected_failure_list=None,
                expected_unstable_list=None,
                parallel_list=None,
                nprocs=1,
                out=sys.stdout,
                log=None,
                verbosity=DEFAULT_VERBOSITY,
                max_time=180):
    if (log is None) : log = null_out()
    self.out = multi_out()
    self.log = log
    self.out.register("stdout", out)
    self.out.register("log", log)
    self.verbosity = verbosity
    self.quiet = (verbosity == 0)
    self.results = []
    self.pool = None

    self.expected_failure_list = expected_failure_list
    if self.expected_failure_list is None:
      self.expected_failure_list = []
    self.expected_unstable_list = expected_unstable_list
    if self.expected_unstable_list is None:
      self.expected_unstable_list = []
    self.parallel_list = parallel_list
    if self.parallel_list is None:
      self.parallel_list = []

    # Filter cmd list for duplicates.
    self.cmd_list = []
    for cmd in cmd_list :
      if (not cmd in self.cmd_list):
        self.cmd_list.append(cmd)
      else :
        print("Test %s repeated, skipping"%cmd, file=self.out)

    # Set number of processors.
    if (nprocs is Auto):
      nprocs = cpu_count()
    if len(self.parallel_list) == 0:
      nprocs = min(nprocs, len(self.cmd_list))

    # Starting summary.
    if (self.verbosity > 0):
      print("Running %d tests on %s processors:"%
        (len(self.cmd_list) + len(self.parallel_list), nprocs), file=self.out)
      for cmd in self.cmd_list:
        print("  %s"%cmd, file=self.out)
      for cmd in self.parallel_list:
        print("  %s [Parallel]"%cmd, file=self.out)
      print("", file=self.out)
      self.out.flush()

    # Either run tests in parallel or run parallel tests, but
    # can't run parallel tests in parallel (cctbx#95)
    os.environ['OPENBLAS_NUM_THREADS'] = "1"

    t_start = time.time()

    if nprocs > 1:
      # Run the tests with multiprocessing pool.
      self.pool = Pool(processes=nprocs)
      for command in self.cmd_list:
        self.pool.apply_async(
          run_command,
          [command, verbosity],
          callback=self.save_result)
      try:
        self.pool.close()
      except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating", file=self.out)
        self.pool.terminate()
      finally:
        try:
          self.pool.join()
        except KeyboardInterrupt:
          pass
    else:
      # Run tests serially.
      self.run_serial(self.cmd_list)

    # run parallel tests with multiple processors per test
    os.environ['OMP_NUM_THREADS'] = str(nprocs)
    self.run_serial(self.parallel_list)
    os.environ['OMP_NUM_THREADS'] = "1"

    # Print ending summary.
    t_end = time.time()
    print("="*80, file=self.out)
    print("", file=self.out)
    print("Tests finished. Elapsed time: %.2fs" %(t_end-t_start), file=self.out)
    print("", file=self.out)

    # Process results for errors and warnings.
    extra_stderr = len([result for result in self.results if result.stderr_lines])
    longjobs = [result for result in self.results if result.wall_time > max_time]
    warnings = [result for result in self.results if result.alert_status == Status.WARNING]
    failures = [result for result in self.results if result.alert_status == Status.FAIL]
    expected_failures = [result for result in self.results if result.alert_status == Status.EXPECTED_FAIL]
    expected_unstable = [result for result in self.results if result.alert_status == Status.EXPECTED_UNSTABLE]
    self.finished = len(self.results)
    self.failure = len(failures)
    self.warning = len(warnings)
    self.expected_failures = len(expected_failures)
    self.expected_unstable = len(expected_unstable)

    # Try writing the XML result file
    write_JUnit_XML(self.results, "output.xml")

    # Run time distribution.
    if (libtbx.env.has_module("scitbx")):
      try:
        from scitbx.array_family import flex
        print("Distribution of test runtimes:", file=self.out)
        hist = flex.histogram(flex.double([result.wall_time for result in self.results]), n_slots=10)
        hist.show(f=self.out, prefix="  ", format_cutoffs="%.1fs")
        print("", file=self.out)
      except Exception:
        # Failing on histogram rendering is a poor way to end a test run
        print("(Failed to render histogram of test results)", file=self.out)

    # Long job warning.
    if longjobs:
      print("", file=self.out)
      print("Warning: the following jobs took at least %d seconds:"%max_time, file=self.out)
      for result in sorted(longjobs, key=lambda result:result.wall_time):
        print("  %s: %.1fs"%(result.command, result.wall_time), file=self.out)
    else:
      # Just print 5 worst offenders to encourage developers to check them out
      print("", file=self.out)
      print("Warning: the following are 5 longest jobs:", file=self.out)
      for result in sorted(self.results, key=lambda result:-result.wall_time)[:5]:
        print("  %s: %.1fs"%(result.command, result.wall_time), file=self.out)
    print("Please try to reduce overall runtime - consider splitting up these tests.", file=self.out)
    print("", file=self.out)


    # Failures.
    if failures:
      print("", file=self.out)
      print("Error: the following jobs returned non-zero exit codes or suspicious stderr output:", file=self.out)
      print("", file=self.out)
      for result in warnings:
        self.display_result(result, alert=Status.WARNING, out=self.out, log_return=self.out, log_stderr=self.out)
      for result in failures:
        self.display_result(result, alert=Status.FAIL, out=self.out, log_return=self.out, log_stderr=self.out)
      print("", file=self.out)
      print("Please verify these tests manually.", file=self.out)
      print("", file=self.out)

    # Summary
    print("Summary:", file=self.out)
    print("  Tests run                    :",self.finished, file=self.out)
    print("  Failures                     :",self.failure, file=self.out)
    print("  Warnings (possible failures) :",self.warning, file=self.out)
    print("  Known Failures (% 3d)         :" % len(self.expected_failure_list),
      self.expected_failures, file=self.out)
    print("  Known Unstable (% 3d)         :" % len(self.expected_unstable_list),
      self.expected_unstable, file=self.out)
    print("  Stderr output (discouraged)  :",extra_stderr, file=self.out)
    if (self.finished != len(self.parallel_list) + len(self.cmd_list)):
      print("*" * 80, file=self.out)
      print("  WARNING: NOT ALL TESTS FINISHED!", file=self.out)
      print("*" * 80, file=self.out)

  def determine_result_status(self, result):
    alert = Status.OK
    # Note: error_lines is populated when package phenix_regression is configured
    if result.error_lines:
      if self.verbosity == EXTRA_VERBOSE:
        print("ERROR BEGIN "*10, file=self.out)
        print(result.error_lines, file=self.out)
        print('-'*80, file=self.out)
        print(result.stderr_lines, file=self.out)
        print("ERROR -END- "*10, file=self.out)
      alert = Status.WARNING
    if result.return_code != 0:
      if self.verbosity == EXTRA_VERBOSE:
        print("RETURN CODE BEGIN "*5, file=self.out)
        print(result.return_code, file=self.out)
        print("RETURN CODE -END- "*5, file=self.out)
      alert = Status.FAIL
      if result.command in self.expected_failure_list:
        alert = Status.EXPECTED_FAIL
      elif result.command in self.expected_unstable_list:
        alert = Status.EXPECTED_UNSTABLE
    return alert

  def run_serial(self, command_list):
    for command in command_list:
      rc = run_command(command, verbosity=self.verbosity)
      if self.save_result(rc) == False:
        break

  def save_result(self, result):
    if result is None:
      print("ERROR: job returned None, assuming CTRL+C pressed", file=self.out)
      if self.pool: self.pool.terminate()
      return False
    result.alert_status = self.determine_result_status(result)
    # If we got an 'OK' status but have stderr output, we can ignore this if
    # "OK" was the last line (e.g. python's unittest does this)
    if result.alert_status == Status.OK and result.stderr_lines:
      test_ok = (result.stderr_lines[-1].strip().startswith("OK"))
      if test_ok:
        result.stderr_lines = []
    self.results.append(result)
    kw = {}
    kw['out'] = self.out
    kw['log_return'] = self.log
    kw['log_stderr'] = True
    kw['log_stdout'] = self.log
    if self.quiet:
      kw['out'] = self.log
      kw['log_stderr'] = False
    elif self.verbosity == EXTRA_VERBOSE:
      kw['log_return'] = self.out
      kw['log_stderr'] = True
      kw['log_stdout'] = self.out
    # For any "not good" result, print out some more details
    elif not result.alert_status == Status.OK:
      kw['log_return'] = self.out
      kw['log_stderr'] = True
    self.display_result(
      result,
      alert=result.alert_status,
      **kw
    )

  def display_result(self, result, alert, out=None, log_return=True, log_stderr=True, log_stdout=False):
    status = {Status.OK: 'OK', Status.WARNING: 'WARNING', Status.FAIL: 'FAIL',
              Status.EXPECTED_FAIL: 'EXPECTED FAIL',
              Status.EXPECTED_UNSTABLE: 'EXPECTED UNSTABLE'}
    if out:
      print("%s [%s] %.1fs"%(result.command, status[alert], result.wall_time), file=out)
      out.flush()
    if log_return:
      print("  Time: %5.2f"%result.wall_time, file=log_return)
      print("  Return code: %s"%result.return_code, file=log_return)
      print("  OKs:", len([x for x in result.stdout_lines if 'OK' in x]), file=log_return)
      log_return.flush()
    if log_stdout and (len(result.stdout_lines) > 0):
      print("  Standard out:", file=log_stdout)
      print("    "+"\n    ".join(result.stdout_lines), file=log_stdout)
      log_stdout.flush()
    if log_stderr and (len(result.stderr_lines) > 0):
      print("  Standard error:", file=sys.stderr)
      print("    "+"\n    ".join(result.stderr_lines), file=sys.stderr)
      sys.stderr.flush()

def make_commands(files,python_keyword_text=""):
  commands = []
  non_executable = []
  unrecognized = []
  for file_name in files :
    if file_name.endswith('.py'):
      if python_keyword_text:
        cmd = 'libtbx.python %s "%s"'%(python_keyword_text,file_name)
      else:
        cmd = 'libtbx.python "%s"'%(file_name) # usual
    elif file_name.endswith('.sh'):
      # interpreter = 'libtbx.bash'
      cmd = file_name
    elif file_name.endswith('.csh'):
      # interpreter = 'libtbx.csh'
      cmd = file_name
    else:
      unrecognized.append(file_name)
      continue
    # cmd = '%s "%s"'%(interpreter, file_name)
    if cmd not in commands:
      commands.append(cmd)
  if (len(unrecognized) > 0):
    raise RuntimeError("""\
The following files could not be recognized as programs:

  %s

Please use the extensions .py, .sh, or .csh for all tests.  Shell scripts will
also need to be made executable.""" % "\n  ".join(unrecognized))
  if (len(non_executable) > 0):
    raise RuntimeError("""\
The following shell scripts are not executable:

  %s

Please enable execution of these to continue.""" % "\n  ".join(non_executable))
  return commands

def exercise():
  log=open("test_utils_parallel_zlog", "wb")
  f0 = open("test_parallel.csh", "wb")
  f0.write("""\
#!/bin/csh

echo "hello world"

exit""")
  f0.close()
  easy_run.call("chmod 755 test_parallel.csh")
  f1 = open("fails.py", "wb")
  f1.write("""\
print "this will crash"
assert 0""")
  f1.close()
  f2 = open("works.py", "wb")
  f2.write("""\
print "hello world"
print "OK"
""")
  f2.close()
  out = StringIO()
  run_command_list([
    "libtbx.python fails.py",
    "libtbx.python works.py",
    "csh test_parallel.csh",
    "./test_parallel.csh",
    ],
    nprocs=1,
    log=log,
    out=out)
  log.close()
  assert ("ERROR - return code 1" in out.getvalue())
  assert ("  Tests run                    : 4" in out.getvalue())
  assert ("  Failures                     : 1" in out.getvalue())
  assert ("  Warnings (possible failures) : 0" in out.getvalue())

if (__name__ == "__main__"):
  exercise()
  print("OK")
