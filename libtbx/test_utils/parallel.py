from libtbx import easy_run
from libtbx import test_utils
import libtbx.load_env
from libtbx.utils import import_python_object, Sorry, multi_out
from libtbx import group_args
from multiprocessing import Pool
from cStringIO import StringIO
import traceback
import time
import os
import sys

max_time = 200 # XXX the lower the better

def get_module_tests (module_name, valgrind=False) :
  dist_path = libtbx.env.dist_path(module_name)
  if (dist_path is None) :
    raise Sorry("'%s' is not a valid CCTBX module." % module_name)
  elif (not os.path.isfile(os.path.join(dist_path, "run_tests.py"))) :
    raise Sorry("%s/run_tests.py does not exist." % module_name)
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
                verbose=False,
                out=None):
  try :
    t0=time.time()
    sys.stdout.flush()
    sys.stderr.flush()
    cmd_result = easy_run.fully_buffered(
      command=command,
      #join_stdout_stderr=join_stdout_stderr,
      )
    #if (len(cmd_result.stderr_lines) != 0):
    #  if verbose :
    #    print '!'*80
    #    print "command"
    #    print command
    #    if (cmd_result.return_code != 0) :
    #      print "ERROR - return code %d" % cmd_result.return_code
    #    print "stderr"
    #    #print "\n".join(cmd_result.stdout_lines)
    #    print "\n".join(cmd_result.stderr_lines)
    #    #cmd_result.raise_if_errors()
    #    print '!'*80
    if 0:
      test_utils._check_command_output(
        lines=cmd_result.stdout_lines,
        show_command_if_error=1, #show_command_if_error,
        sorry_expected=0, #sorry_expected,
        )
    sys.stdout.flush()
    sys.stderr.flush()
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

class run_command_list (object) :
  def __init__ (self,
                cmd_list,
                nprocs=1,
                out=sys.stdout,
                log=None,
                quiet=False) :
    if (log is None) : log = null_out()
    self.out = multi_out()
    self.log = log
    self.out.register("stdout", out)
    self.out.register("log", log)
    self.quiet = quiet
    self.cmd_list = []
    for cmd in cmd_list :
      if (not cmd in self.cmd_list) :
        self.cmd_list.append(cmd)
      else :
        print >> self.out, "  test %s repeated, skipping" % cmd
    nprocs = min(nprocs, len(self.cmd_list))
    print >> self.out, "\n  Starting command list"
    print >> self.out, "    NProcs :",nprocs
    print >> self.out, "    Cmds   :",len(self.cmd_list)
    t_start = time.time()
    if nprocs>1:
      pool = Pool(processes=nprocs)
    self.results = []
    for command in self.cmd_list:
      if nprocs>1:
        pool.apply_async(
          run_command,
          [command, (not quiet), out],
          callback=self.save_result)
      else:
        rc = run_command(command, verbose=(not quiet), out=out)
        self.save_result(rc)
    if nprocs>1:
      try :
        try :
          pool.close()
        except KeyboardInterrupt :
          print >> self.out, "Caught KeyboardInterrupt, terminating"
          pool.terminate()
      finally :
        pool.join()
      print >> self.out, '\nProcesses have joined : %d\n' % len(self.results)
    t_end = time.time()
    print >> self.out, ""
    print >> self.out, "Elapsed time: %.2fs" %(t_end-t_start)
    print >> self.out, ""
    finished = 0
    warning = 0
    extra_stderr = 0
    failure = 0
    failures = []
    long_jobs = []
    long_runtimes = []
    runtimes = []
    for result in self.results :
      finished += 1
      runtimes.append(result.wall_time)
      if (result.return_code != 0) :
        failure += 1
        failures.append(result)
      else :
        if (len(result.error_lines) != 0) :
          warning += 1
          failures.append(result)
        if (len(result.stderr_lines) != 0):
          extra_stderr += 1
      if (result.wall_time > max_time) :
        long_jobs.append(result.command)
        long_runtimes.append(result.wall_time)
    if (libtbx.env.has_module("scitbx")) :
      from scitbx.array_family import flex
      print >> self.out, "Distribution of test runtimes:"
      hist = flex.histogram(flex.double(runtimes), n_slots=20)
      hist.show(f=self.out, prefix="  ", format_cutoffs="%.1fs")
      print >> self.out, ""
    if (len(long_jobs) > 0) :
      print >> self.out, ""
      print >> self.out, "WARNING: the following jobs took at least %d seconds:" % \
        max_time
      jobs_and_timings = list(zip(long_jobs, long_runtimes))
      jobs_and_timings.sort(lambda x,y: cmp(x[1], y[1]))
      for cmd, runtime in jobs_and_timings :
        print >> self.out, "  " + cmd + " : %.1fs" % runtime
      print >> self.out, "Please try to reduce overall runtime - consider splitting up these tests."
    if (len(failures) > 0) :
      print >> self.out, ""
      print >> self.out, "ERROR: the following jobs returned non-zero exit codes or suspicious stderr output:"
      for result in failures :
        print >> self.out, ""
        print >> self.out, result.command + "(exit code %d):" % result.return_code
        for line in result.stderr_lines :
          print >> self.out, "  " + line
        for line in result.error_lines :
          print >> self.out, "  " + line
        print >> self.out, ""
      print >> self.out, "Please verify these tests manually."
      print >> self.out, ""
    print >> self.out, "Summary:"
    print >> self.out, "  Tests run                    :",finished
    print >> self.out, "  Failures                     :",failure
    print >> self.out, "  Warnings (possible failures) :",warning
    print >> self.out, "  Stderr output (discouraged)  :",extra_stderr
    if (finished != len(self.cmd_list)) :
      print >> self.out, "*" * 80
      print >> self.out, "  WARNING: NOT ALL TESTS FINISHED!"
      print >> self.out, "*" * 80


  def save_result (self, result) :
    if (result is None ):
      print >> self.out, "ERROR: job returned None"
      return
    self.results.append(result)
    self.display_result(result)

  def display_result (self, result) :
    alert = False
    if (result.return_code != 0) or (len(result.error_lines) > 0) :
      alert = True
    if (not self.quiet) and (len(result.error_lines) > 0) :
      alert = True
    if (alert) :
      out = self.out
    else :
      out = self.log
    print >> out, '_'*80
    print >> out, '\ncommand : "%s"' % result.command
    print >> out, 'return_code : %s' % result.return_code
    oks = 0
    for line in result.stdout_lines:
      if line.find("OK")>-1:
        oks += 1
    print >> out, 'Found %d "OK"' % oks
    if (len(result.stderr_lines) != 0):
      print >> out, 'stderr-'*10
      print >> out, "\n".join(result.stderr_lines)
    if (len(result.error_lines) > 0) :
      print >> out, "possible errors:"
      print >> out, "\n".join(result.error_lines)
    print >> out, "time : %5.2fs" % result.wall_time
    if result.wall_time>60:
      print >> out, '!'*78
      print >> out, "!!","WARNING "*9,"!!"
      print >> out, "!!  %-71s !!" %"TEST TAKES MORE THAN A MINUTE"
      print >> out, "!!","WARNING "*9,"!!"
      print >> out, '!'*78

def make_commands (files) :
  commands = []
  non_executable = []
  unrecognized = []
  for file_name in files :
    if (file_name.endswith(".py")) :
      cmd = "libtbx.python \"%s\"" % file_name
      if (not cmd in commands) :
        commands.append(cmd)
    elif (file_name.endswith(".sh")) or (file_name.endswith(".csh")) :
      if (not os.access(file_name, os.X_OK)) :
        non_executable.append(file_name)
      elif (not file_name in commands) :
        commands.append(file_name)
    else :
      unrecognized.append(file_name)
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
