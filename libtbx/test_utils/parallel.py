
from __future__ import division
from libtbx import easy_run
from libtbx import test_utils
import libtbx.load_env
from libtbx.utils import import_python_object, Sorry, multi_out
from libtbx import group_args, Auto
from multiprocessing import Pool, cpu_count
from cStringIO import StringIO
import traceback
import time
import os
import sys

QUIET = 0
DEFAULT_VERBOSITY = 1
EXTRA_VERBOSE = 2 # for nightly builds

def get_module_tests (module_name, valgrind=False) :
  dist_path = libtbx.env.dist_path(module_name)
  if (dist_path is None) :
    raise Sorry("'%s' is not a valid CCTBX module." % module_name)
  # XXX don't check for file name, because import conventions differ among
  # derived modules - dist_path can be either the sys.path entry or the actual
  # module contents.  If the file is missing the import below will fail, which
  # is okay for testing purposes.
  tst_list = import_python_object(
    import_path="%s.run_tests.tst_list" % module_name,
    error_prefix="",
    target_must_be="",
    where_str="").object
  assert (isinstance(tst_list, tuple) or isinstance(tst_list, list))
  build_path = libtbx.env.under_build(module_name)
  assert (build_path is not None) and (dist_path is not None)
  commands = []
  co = group_args(
    verbose=False,
    quick=True,
    valgrind=valgrind)
  for cmd in test_utils.iter_tests_cmd(
      co=co,
      build_dir=build_path,
      dist_dir=dist_path,
      tst_list=tst_list) :
    commands.append(cmd)
  return commands

def find_tests (dir_name) :
  if (not os.path.isdir(dir_name)) :
    raise RuntimeError("'%s' is not a directory." % dir_name)
  all_tests = []
  for root, dirnames, filenames in os.walk(dir_name):
    for file_name in filenames :
      base, ext = os.path.splitext(file_name)
      if (file_name.startswith("tst")) and (ext in [".py", ".csh", ".sh"]) :
        all_tests.append(os.path.join(root, file_name))
  return all_tests

def run_command(command,
                verbosity=DEFAULT_VERBOSITY,
                out=None):
  if (out is None) :
    out = sys.stdout
  try :
    t0=time.time()
    cmd_result = easy_run.fully_buffered(
      command=command,
      #join_stdout_stderr=join_stdout_stderr,
      )
    if 0:
      test_utils._check_command_output(
        lines=cmd_result.stdout_lines,
        show_command_if_error=1, #show_command_if_error,
        sorry_expected=0, #sorry_expected,
        )
    cmd_result.wall_time = time.time()-t0
    cmd_result.error_lines = evaluate_output(cmd_result)
    return cmd_result
  except KeyboardInterrupt :
    traceback_str = "\n".join(traceback.format_tb(sys.exc_info()[2]))
    sys.stdout.write(traceback_str)
    return None

def evaluate_output (cmd_result) :
  have_regression = libtbx.env.has_module("phenix_regression")
  bad_lines = []
  if (have_regression) :
    from phenix_regression.command_line import find_errors_and_warnings
    out = cmd_result.stdout_lines
    err = cmd_result.stderr_lines
    for i, line in enumerate(out+err) :
      if (find_errors_and_warnings.is_error_or_warning(line, line.lower())) :
        bad_lines.append(line)
  return bad_lines

def reconstruct_test_name (command) :
  import re
  if "-m pytest" in command:
    m = re.search('-m pytest ([^:]*)::(.*)$', command)
    command = 'libtbx.python "%s" %s' % (m.group(1), m.group(2))
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

    pattern2  = '^/?(.*?/((modules|build)/(cctbx\_project/|xia2/Test/)?))?(.*)/([^/&]*?)(\.py)?$'
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

class run_command_list (object) :
  def __init__ (self,
                cmd_list,
                nprocs=1,
                out=sys.stdout,
                log=None,
                verbosity=DEFAULT_VERBOSITY,
                max_time=180) :
    if (log is None) : log = null_out()
    self.out = multi_out()
    self.log = log
    self.out.register("stdout", out)
    self.out.register("log", log)
    self.verbosity = verbosity
    self.quiet = (verbosity == 0)
    self.results = []

    # Filter cmd list for duplicates.
    self.cmd_list = []
    for cmd in cmd_list :
      if (not cmd in self.cmd_list) :
        self.cmd_list.append(cmd)
      else :
        print >> self.out, "Test %s repeated, skipping"%cmd

    # Set number of processors.
    if (nprocs is Auto) :
      nprocs = cpu_count()
    nprocs = min(nprocs, len(self.cmd_list))

    # Starting summary.
    if (self.verbosity > 0) :
      print >> self.out, "Running %d tests on %s processors:"%(len(self.cmd_list), nprocs)
      for cmd in self.cmd_list:
        print >> self.out, "  %s"%cmd
      print >> self.out, ""

    t_start = time.time()
    if nprocs > 1:
      # Run the tests with multiprocessing pool.
      pool = Pool(processes=nprocs)
      for command in self.cmd_list:
        pool.apply_async(
          run_command,
          [command, verbosity, out],
          callback=self.save_result)
      try:
        pool.close()
      except KeyboardInterrupt:
        print >> self.out, "Caught KeyboardInterrupt, terminating"
        pool.terminate()
      finally:
        pool.join()
    else:
      # Run tests serially.
      for command in self.cmd_list:
        rc = run_command(command, verbosity=verbosity, out=out)
        self.save_result(rc)

    # Print ending summary.
    t_end = time.time()
    print >> self.out, "="*80
    print >> self.out, ""
    print >> self.out, "Tests finished. Elapsed time: %.2fs" %(t_end-t_start)
    print >> self.out, ""
    test_cases = []
    # Process results for errors and warnings.
    extra_stderr = len([result for result in self.results if result.stderr_lines])
    longjobs = [result for result in self.results if result.wall_time > max_time]
    warnings = [result for result in self.results if self.check_alert(result) == 1]
    failures = [result for result in self.results if self.check_alert(result) == 2]
    self.finished = len(self.results)
    self.failure = len(failures)
    self.warning = len(warnings)

    # Output JUnit XML if possible
    try:
      from junit_xml import TestSuite, TestCase
      import re
      def decode_string(string):
        try:
          return string.encode('ascii', 'xmlcharrefreplace')
        except Exception: # intentional
          return unicode(string, errors='ignore').encode('ascii', 'xmlcharrefreplace')
      for result in self.results:
        test_name = reconstruct_test_name(result.command)
        plain_stdout = map(decode_string, result.stdout_lines)
        plain_stderr = map(decode_string, result.stderr_lines)
        output = '\n'.join(plain_stdout + plain_stderr)
        tc = TestCase(classname=test_name[0],
                      name=test_name[1],
                      elapsed_sec=result.wall_time,
                      stdout='\n'.join(plain_stdout),
                      stderr='\n'.join(plain_stderr))
        if result.return_code == 0:
          # Identify skipped tests
          if re.search('skip', output, re.IGNORECASE):
            # find first line including word 'skip' and use it as message
            skipline = re.search('^((.*)skip(.*))$', output, re.IGNORECASE | re.MULTILINE).group(1)
            tc.add_skipped_info(skipline)
        else:
          # Test failed. Extract error message and stack trace if possible
          error_message = 'exit code %d' % result.return_code
          error_output = '\n'.join(plain_stderr)
          if plain_stderr:
            error_message = plain_stderr[-1]
            if len(plain_stderr) > 20:
              error_output = '\n'.join(plain_stderr[-20:])
          tc.add_failure_info(message=error_message, output=error_output)
        test_cases.append(tc)
      ts = TestSuite("libtbx.run_tests_parallel", test_cases=test_cases)
      with open('output.xml', 'wb') as f:
        print >> f, TestSuite.to_xml_string([ts], prettyprint=True)
    except ImportError, e:
      pass

    # Run time distribution.
    if (libtbx.env.has_module("scitbx")) :
      from scitbx.array_family import flex
      print >> self.out, "Distribution of test runtimes:"
      hist = flex.histogram(flex.double([result.wall_time for result in self.results]), n_slots=10)
      hist.show(f=self.out, prefix="  ", format_cutoffs="%.1fs")
      print >> self.out, ""

    # Long job warning.
    if longjobs:
      print >> self.out, ""
      print >> self.out, "Warning: the following jobs took at least %d seconds:"%max_time
      for result in sorted(longjobs, key=lambda result:result.wall_time):
        print >> self.out, "  %s: %.1fs"%(result.command, result.wall_time)
      print >> self.out, "Please try to reduce overall runtime - consider splitting up these tests."

    # Failures.
    if failures:
      print >> self.out, ""
      print >> self.out, "Error: the following jobs returned non-zero exit codes or suspicious stderr output:"
      print >> self.out, ""
      for result in warnings:
        self.display_result(result, alert=1, out=self.out, log_return=self.out, log_stderr=self.out)
      for result in failures:
        self.display_result(result, alert=2, out=self.out, log_return=self.out, log_stderr=self.out)
      print >> self.out, ""
      print >> self.out, "Please verify these tests manually."
      print >> self.out, ""

    # Summary
    print >> self.out, "Summary:"
    print >> self.out, "  Tests run                    :",self.finished
    print >> self.out, "  Failures                     :",self.failure
    print >> self.out, "  Warnings (possible failures) :",self.warning
    print >> self.out, "  Stderr output (discouraged)  :",extra_stderr
    if (self.finished != len(self.cmd_list)) :
      print >> self.out, "*" * 80
      print >> self.out, "  WARNING: NOT ALL TESTS FINISHED!"
      print >> self.out, "*" * 80

  def check_alert(self, result):
    alert = 0
    if result.error_lines:
      if self.verbosity == EXTRA_VERBOSE:
        print >> self.out, "ERROR BEGIN "*10
        print >> self.out, result.error_lines
        print >> self.out, '-'*80
        print >> self.out, result.stderr_lines
        print >> self.out, "ERROR -END- "*10
      alert = 1
    if result.return_code != 0:
      if self.verbosity == EXTRA_VERBOSE:
        print >> self.out, "RETURN CODE BEGIN "*5
        print >> self.out, result.return_code
        print >> self.out, "RETURN CODE -END- "*5
      alert = 2
    return alert

  def save_result (self, result) :
    if (result is None ):
      print >> self.out, "ERROR: job returned None"
      return
    alert = self.check_alert(result)
    # Clear stderr when python unittest framework is used and test passes
    if alert == 0 and result.stderr_lines:
      test_ok = (result.stderr_lines[-1].strip() == 'OK')
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
    elif alert:
      kw['log_return'] = self.out
      kw['log_stderr'] = True
    self.display_result(
      result,
      alert=alert,
      **kw
    )

  def display_result (self, result, alert, out=None, log_return=True, log_stderr=True, log_stdout=False) :
    status = ['OK', 'WARNING', 'FAIL']
    if out:
      print >> out, "%s [%s]"%(result.command, status[alert])
      out.flush()
    if log_return:
      print >> log_return, "  Time: %5.2f"%result.wall_time
      print >> log_return, "  Return code: %s"%result.return_code
      print >> log_return, "  OKs:", len(filter(lambda x:'OK' in x, result.stdout_lines))
      log_return.flush()
    if log_stdout and (len(result.stdout_lines) > 0) :
      print >> log_stdout, "  Standard out:"
      print >> log_stdout, "    "+"\n    ".join(result.stdout_lines)
      log_stdout.flush()
    if log_stderr and (len(result.stderr_lines) > 0) :
      print >> sys.stderr, "  Standard error:"
      print >> sys.stderr, "    "+"\n    ".join(result.stderr_lines)
      sys.stderr.flush()

def make_commands (files) :
  commands = []
  non_executable = []
  unrecognized = []
  for file_name in files :
    if file_name.endswith('.py'):
      cmd = 'libtbx.python "%s"'%(file_name)
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
  if (len(unrecognized) > 0) :
    raise RuntimeError("""\
The following files could not be recognized as programs:

  %s

Please use the extensions .py, .sh, or .csh for all tests.  Shell scripts will
also need to be made executable.""" % "\n  ".join(unrecognized))
  if (len(non_executable) > 0) :
    raise RuntimeError("""\
The following shell scripts are not executable:

  %s

Please enable execution of these to continue.""" % "\n  ".join(non_executable))
  return commands

def exercise () :
  log=file("test_utils_parallel_zlog", "wb")
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

if (__name__ == "__main__") :
  exercise()
  print "OK"
