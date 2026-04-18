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
import errno
import shutil
import stat

try:
  from enum import IntEnum
except ImportError:
  IntEnum = object

def _chmod_then_retry(func, path, exc_info):
  """Recover from Windows read-only removal failures inside ``rmtree``.

  Designed to be passed as the ``onerror`` callback to
  :func:`shutil.rmtree`. When ``func`` is a removal function
  (``os.unlink``, ``os.remove``, or ``os.rmdir``) and the underlying
  exception is an ``OSError`` with ``errno`` equal to ``EACCES`` or
  ``EPERM``, the target path is ``chmod``-ed writable and the
  removal is retried. Any other error is ignored — the caller
  treats cleanup as best-effort.

  Parameters
  ----------
  func : callable
      The filesystem function that raised (one of ``os.unlink``,
      ``os.remove``, ``os.rmdir``).
  path : str
      The path the function was attempting to operate on.
  exc_info : tuple
      The exception triple as produced by ``sys.exc_info()``.

  Notes
  -----
  On POSIX systems a file's permissions do not prevent unlinking
  when its parent directory is writable, so this callback is
  effectively dead code there. On Windows, files copied out of a
  read-only source tree carry the read-only attribute; removing
  them requires clearing the attribute first.

  Never raises.
  """
  excvalue = exc_info[1]
  if (func in (os.unlink, os.remove, os.rmdir)
      and isinstance(excvalue, OSError)
      and excvalue.errno in (errno.EACCES, errno.EPERM)):
    try:
      os.chmod(path, stat.S_IWUSR | stat.S_IWRITE)
      func(path)
    except OSError:
      pass

def _reset_cwd_for_retry(cwd):
  """Empty a test's working directory in place before a retry.

  Removes every entry inside ``cwd`` while leaving the directory
  itself intact. Subdirectories are recursively removed via
  :func:`shutil.rmtree` with :func:`_chmod_then_retry` as the
  ``onerror`` callback; plain files are removed with an
  ``os.chmod`` + retry fallback for the Windows read-only case.

  Parameters
  ----------
  cwd : str or None
      Absolute path of the test's per-run working directory (as
      produced by ``_build_unique_dir_mapping``). ``None``, the
      empty string, and non-existent paths are silently tolerated.

  Notes
  -----
  Symlinks are unlinked rather than followed, so removal never
  escapes ``cwd``. The directory itself is not removed or
  recreated — this avoids a Windows race in which a freshly
  created directory is not yet visible to a subprocess spawned
  immediately afterward.

  Best-effort. Never raises.
  """
  if not cwd or not os.path.isdir(cwd):
    return
  for entry in os.listdir(cwd):
    path = os.path.join(cwd, entry)
    try:
      if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path, onerror=_chmod_then_retry)
      else:
        try:
          os.remove(path)
        except OSError:
          try:
            os.chmod(path, stat.S_IWUSR | stat.S_IWRITE)
            os.remove(path)
          except OSError:
            pass
    except OSError:
      pass

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
                verbosity=DEFAULT_VERBOSITY,
                cwd=None,
                max_retries=0,
                skip_retry=False):
  """Run a single test command, optionally retrying on failure.

  The command is executed via :func:`libtbx.easy_run.fully_buffered`.
  When the first attempt exits with a non-zero return code, up to
  ``max_retries`` additional attempts are performed with exponential
  backoff (1s, 2s, 4s, ..., capped at 60s) between them. If ``cwd``
  is supplied, it is emptied via :func:`_reset_cwd_for_retry` before
  each retry so the next attempt starts from a clean state.

  Parameters
  ----------
  command : str
      The shell command to execute (typically
      ``libtbx.python "path/to/tst_xxx.py"``).
  verbosity : int, optional
      Verbosity level. Accepted for API compatibility with older
      callers; not referenced by the retry path.
  cwd : str or None, optional
      Absolute path of the working directory in which to run the
      command. When set, the directory is emptied between retries.
      When ``None`` (the default), the command inherits the
      caller's cwd and no cleanup is performed.
  max_retries : int, optional
      Number of additional attempts to make after an initial
      failure. Total attempts = ``1 + max_retries``. Default ``0``
      (no retries, preserves legacy behavior).
  skip_retry : bool, optional
      If ``True``, perform exactly one attempt regardless of
      ``max_retries``. Used by callers for tests in the
      expected-failure list.

  Returns
  -------
  libtbx.easy_run.fully_buffered or None
      The result of the final attempt, with additional attributes
      ``wall_time`` (float, seconds of the final attempt),
      ``error_lines`` (list of suspicious stdout/stderr lines as
      classified by ``phenix_regression``), ``attempt`` (int,
      1-indexed index of the final attempt), and ``attempts_total``
      (int, number of attempts allowed for this command). Output
      from earlier attempts is discarded. Returns ``None`` only on
      KeyboardInterrupt.

  Raises
  ------
  ValueError
      If ``max_retries`` is negative.

  Notes
  -----
  Only the final attempt's stdout/stderr and timing are returned;
  output from earlier attempts is discarded. Callers that need to
  retain per-attempt output must capture it themselves.
  """
  if max_retries < 0:
    raise ValueError(
      "max_retries must be >= 0, got %d" % max_retries)
  try:
    attempts = 1 if skip_retry else 1 + max_retries
    last_result = None
    for attempt in range(1, attempts + 1):
      if attempt > 1:
        time.sleep(min(2 ** (attempt - 2), 60))
        if cwd is not None:
          _reset_cwd_for_retry(cwd)
      t0 = time.time()
      cmd_result = easy_run.fully_buffered(
        command=command,
        cwd=cwd,
        )
      cmd_result.wall_time = time.time() - t0
      cmd_result.error_lines = phenix_separate_output(cmd_result)
      cmd_result.attempt = attempt
      cmd_result.attempts_total = attempts
      last_result = cmd_result
      if cmd_result.return_code == 0:
        break
    return last_result
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
                max_time=180,
                cwd_map=None,
                max_retries=0):
    """Run a list of test commands with optional parallelism and retry.

    Parameters
    ----------
    cmd_list : list of str
        Shell commands to execute, typically produced by
        :func:`make_commands`. Duplicates are silently dropped.
    expected_failure_list : list of str, optional
        Commands known to fail; their non-zero exit codes are
        reclassified as ``Status.EXPECTED_FAIL``. Never retried.
    expected_unstable_list : list of str, optional
        Commands known to flicker; their non-zero exit codes are
        reclassified as ``Status.EXPECTED_UNSTABLE``. Retried
        normally.
    parallel_list : list of str, optional
        Commands that run one-at-a-time but each with multiple
        processors (``OMP_NUM_THREADS`` set to ``nprocs``).
    nprocs : int or libtbx.Auto, optional
        Worker count for the multiprocessing pool. ``Auto`` means
        :func:`multiprocessing.cpu_count`. Default ``1``.
    out : file, optional
        Stream for per-test progress output. Default ``sys.stdout``.
    log : file or None, optional
        Secondary log stream. Default ``None`` (drops output).
    verbosity : int, optional
        ``0`` = quiet, ``1`` = default, ``2`` = extra verbose.
    max_time : float, optional
        Tests taking longer than ``max_time`` seconds are flagged
        as "long jobs" in the summary.
    cwd_map : dict, optional
        ``{command: cwd}`` mapping used when tests are dispatched
        with per-test working directories (the
        ``run_in_unique_dirs`` mode).
    max_retries : int, optional
        Additional attempts for commands that exit with a non-zero
        return code. ``0`` (the default) preserves legacy
        behavior. Commands in ``expected_failure_list`` are never
        retried regardless of this value.
    """
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
    self.cwd_map = cwd_map if cwd_map is not None else {}
    self.max_retries = max_retries

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
        skip_retry = command in self.expected_failure_list
        self.pool.apply_async(
          run_command,
          [command, verbosity, self.cwd_map.get(command),
           self.max_retries, skip_retry],
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
    retry_results = [r for r in self.results if getattr(r, 'attempt', 1) > 1]
    if retry_results:
      extra_attempts = sum(r.attempt - 1 for r in retry_results)
      print("  Retries used                 : %d attempts across %d tests"
        % (extra_attempts, len(retry_results)), file=self.out)
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
    """Run a list of commands serially, forwarding retry configuration.

    Parameters
    ----------
    command_list : list of str
        Commands to dispatch. Each is run through
        :func:`run_command` with this instance's ``max_retries``
        and with ``skip_retry=True`` iff the command appears in
        ``expected_failure_list``.
    """
    for command in command_list:
      skip_retry = command in self.expected_failure_list
      rc = run_command(
        command,
        verbosity=self.verbosity,
        cwd=self.cwd_map.get(command),
        max_retries=self.max_retries,
        skip_retry=skip_retry)
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

  def display_result(self, result, alert, out=None, log_return=True,
                     log_stderr=True, log_stdout=False):
    """Print one test's outcome, with a retry annotation when relevant.

    Appends ``(passed on attempt N of M)`` or ``(failed after M
    attempts)`` to the status line whenever ``result.attempt`` is
    greater than 1.
    """
    status = {Status.OK: 'OK', Status.WARNING: 'WARNING', Status.FAIL: 'FAIL',
              Status.EXPECTED_FAIL: 'EXPECTED FAIL',
              Status.EXPECTED_UNSTABLE: 'EXPECTED UNSTABLE'}
    retry_note = ""
    attempt = getattr(result, 'attempt', 1)
    attempts_total = getattr(result, 'attempts_total', 1)
    if attempt > 1:
      if result.return_code == 0:
        retry_note = "  (passed on attempt %d of %d)" % (
          attempt, attempts_total)
      else:
        retry_note = "  (failed after %d attempts)" % attempt
    if out:
      print("%s [%s] %.1fs%s" % (
        result.command, status[alert], result.wall_time, retry_note),
        file=out)
      out.flush()
    if log_return:
      print("  Time: %5.2f"%result.wall_time, file=log_return)
      print("  Return code: %s"%result.return_code, file=log_return)
      print("  OKs:", len([x for x in result.stdout_lines if 'OK' in x]),
        file=log_return)
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
