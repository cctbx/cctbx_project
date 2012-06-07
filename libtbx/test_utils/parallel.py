from libtbx import easy_run
from libtbx import test_utils
import libtbx.load_env
from libtbx.utils import import_python_object, Sorry, multi_out
from libtbx import group_args
from multiprocessing import Pool
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
                ):
  t0=time.time()
  sys.stdout.flush()
  sys.stderr.flush()
  cmd_result = easy_run.fully_buffered(
    command=command,
    #join_stdout_stderr=join_stdout_stderr,
    )
  if (len(cmd_result.stderr_lines) != 0):
    if verbose:
      print '!'*80
      print "command"
      print command
      if (cmd_result.return_code != 0) :
        print "ERROR - return code %d" % cmd_result.return_code
      print "stderr"
      #print "\n".join(cmd_result.stdout_lines)
      print "\n".join(cmd_result.stderr_lines)
      #cmd_result.raise_if_errors()
      print '!'*80
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

def evaluate_output (cmd_result) :
  have_regression = libtbx.env.has_module("phenix_regression")
  bad_lines = []
  if (have_regression) :
    from phenix_regression.command_line import find_errors_and_warnings
    out = cmd_result.stdout_lines
    err = cmd_result.stderr_lines
    for i, line in enumerate(out) :
      if (find_errors_and_warnings.is_error_or_warning(line, line.lower())) :
        bad_lines.append(line)
  return bad_lines

def display_result(result, log=sys.stdout):
  print >> log, '_'*80
  print >> log, '\ncommand : "%s"' % result.command
  print >> log, 'return_code : %s' % result.return_code
  print >> log, 'stdout-'*10
  stdout_lines = "\n".join(result.stdout_lines)
  print >> log, stdout_lines
  oks = 0
  for line in result.stdout_lines:
    if line.find("OK")>-1:
      oks += 1
  print >> log, 'Found %d "OK"' % oks
  if (len(result.stderr_lines) != 0):
    print >> log, 'stderr-'*10
    print >> log, "\n".join(result.stderr_lines)
  if (len(result.error_lines) > 0) :
    print >> log, "possible errors:"
    print >> log, "\n".join(result.error_lines)
  print >> log, "time : %5.2fs" % result.wall_time
  if result.wall_time>60:
    print >> log, '!'*78
    print >> log, "!!","WARNING "*9,"!!"
    print >> log, "!!  %-71s !!" %"TEST TAKES MORE THAN A MINUTE"
    print >> log, "!!","WARNING "*9,"!!"
    print >> log, '!'*78

def run_command_list(cmd_list,
                     nprocs=1,
                     out=sys.stdout,
                     log=None,
                     verbose=0):
  if (log is None) : log = null_out()
  out_ = multi_out()
  out_.register("stdout", out)
  out_.register("log", log)
  unfiltered = list(cmd_list)
  cmd_list = []
  for cmd in unfiltered :
    if (not cmd in cmd_list) :
      cmd_list.append(cmd)
    else :
      print >> out_, "  test %s repeated, skipping" % cmd
  nprocs = min(nprocs, len(cmd_list))
  print >> out_, "\n  Starting command list"
  print >> out_, "    NProcs :",nprocs
  print >> out_, "    Cmds   :",len(cmd_list)
  t_start = time.time()
  if nprocs>1:
    pool = Pool(processes=nprocs)

  results = []
  def save_result (result) :
    results.append(result)
    if (verbose > 0) :
      display_result(result, log=out_)
    else :
      display_result(result, log=log)
  for command in cmd_list:
    if nprocs>1:
      pool.apply_async(
        run_command,
        [command, True],
        callback=save_result)
    else:
      rc = run_command(command, verbose=True)
      display_result(rc, log=log)
      results.append(rc)
  if nprocs>1:
    try :
      try :
        pool.close()
      except KeyboardInterrupt :
        print >> out_, "Caught KeyboardInterrupt, terminating"
        pool.terminate()
    finally :
      pool.join()
    print >> out_, '\nProcesses have joined : %d\n' % len(results)
  t_end = time.time()
  print >> out_, ""
  print >> out_, "Elapsed time: %.2fs" %(t_end-t_start)
  print >> out_, ""
  finished = 0
  warning = 0
  extra_stderr = 0
  failure = 0
  failures = []
  long_jobs = []
  long_runtimes = []
  runtimes = []
  for result in results :
    finished += 1
    runtimes.append(result.wall_time)
    if (result.return_code != 0) :
      failure += 1
      failures.append(result.command)
    else :
      if (len(result.error_lines) != 0) :
        warning += 1
        failures.append(result.command)
      if (len(result.stderr_lines) != 0):
        extra_stderr += 1
    if (result.wall_time > max_time) :
      long_jobs.append(result.command)
      long_runtimes.append(result.wall_time)
  print >> out_, "Summary:"
  print >> out_, "  Tests run                    :",finished
  print >> out_, "  Failures                     :",failure
  print >> out_, "  Warnings (possible failures) :",warning
  print >> out_, "  Stderr output (discouraged)  :",extra_stderr
  if (libtbx.env.has_module("scitbx")) :
    from scitbx.array_family import flex
    print >> out_, "Distribution of test runtimes:"
    hist = flex.histogram(flex.double(runtimes), n_slots=20)
    hist.show(f=out_, prefix="  ", format_cutoffs="%.1fs")
    print >> out_, ""
  if (len(long_jobs) > 0) :
    print >> out_, ""
    print >> out_, "WARNING: the following jobs took at least %d seconds:" % \
      max_time
    jobs_and_timings = list(zip(long_jobs, long_runtimes))
    jobs_and_timings.sort(lambda x,y: cmp(x[1], y[1]))
    for cmd, runtime in jobs_and_timings :
      print >> out_, "  " + cmd + " : %.1fs" % runtime
    print >> out_, "Please try to reduce overall runtime - consider splitting up these tests."
  if (len(failures) > 0) :
    print >> out_, ""
    print >> out_, "ERROR: the following jobs returned non-zero exit codes or suspicious stderr output:"
    for command in failures :
      print >> out_, "  " + command
    print >> out_, "Please verify these tests manually."
    print >> out_, ""

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

if __name__=="__main__":
  cwd = os.path.join(os.environ["PHENIX"],
                     "cctbx_project",
                     "libtbx",
                     "test_utils",
                     )

  log=file("zlog", "wb")
  run_command_list([
    "libtbx.python %s" % os.path.join(cwd, "fails.py"),
    "libtbx.python %s" % os.path.join(cwd, "works.py"),
    "csh %s" % os.path.join(cwd, "test.csh"),
    "%s" % os.path.join(cwd, "test.csh"),
    ],
    nprocs=2,
    log=log,
    )
  log.close()
  #run(sys.argv[1])
