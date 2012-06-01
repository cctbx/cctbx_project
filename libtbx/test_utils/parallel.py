from libtbx import easy_run
from libtbx import test_utils
import libtbx.load_env
from libtbx.utils import import_python_object, Sorry
from multiprocessing import Pool
import time
import os
import sys


def get_module_tests (module_name) :
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
  for cmd in test_utils.iter_tests_cmd(
      co=None,
      build_dir=build_path,
      dist_dir=dist_path,
      tst_list=tst_list) :
    commands.append(cmd)
  return commands

def find_tests (dir_name) :
  assert os.path.isdir(dir_name)
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
  return cmd_result

def display_result(result, log=sys.stdout):
  #print dir(result)
  print >> log, '_'*80
  print >> log, '\ncommand : "%s"' % result.command
  print >> log, 'return_code : %s' % result.return_code
  print >> log, 'stdout-'*10
  print >> log, "\n".join(result.stdout_lines)
  if (len(result.stderr_lines) != 0):
    print >> log, 'stderr-'*10
    print >> log, "\n".join(result.stderr_lines)
  print >> log, "time : %5.2fs" % result.wall_time
  if result.wall_time>60:
    print >> log, '!'*78
    print >> log, "!!","WARNING "*9,"!!"
    print >> log, "!!  %-71s !!" %"TEST TAKES MORE THAN A MINUTE"
    print >> log, "!!","WARNING "*9,"!!"
    print >> log, '!'*78

def run_command_list(cmd_list,
                     nprocs=1,
                     log=sys.stdout,
                     ):
  nprocs = min(nprocs, len(cmd_list))
  print "\n  Starting command list"
  print "    NProcs :",nprocs
  print "    Cmds   :",len(cmd_list)
  if nprocs>1:
    pool = Pool(processes=nprocs)

  results = []
  def save_result (result) :
    results.append(result)
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
    pool.close()
    pool.join()
    print '\nProcesses have joined : %d\n' % len(results)
  finished = 0
  warning = 0
  failure = 0
  long_jobs = []
  for result in results :
    finished += 1
    if (result.return_code != 0) :
      failure += 1
    elif (len(result.stderr_lines) != 0):
      warning += 1
    if (result.wall_time > 60) :
      long_jobs.append(result.command)
  print 'Done with output'
  print "  NProcs    :",nprocs
  print "  Tests run :",finished
  print "  Failures  :",failure
  print "  Warnings  :",warning
  if (len(long_jobs) > 0) :
    print ""
    print "WARNING: the following jobs took at least 60 seconds each:"
    for cmd in long_jobs :
      print "  " + cmd
    print "Please try to reduce overall runtime - consider splitting up these tests."
  print >> log, '\n\nDone with output'
  print >> log, "  NProcs    :",nprocs
  print >> log, "  Tests run :",finished
  print >> log, "  Failures  :",failure
  print >> log, "  Warnings  :",warning

def make_commands (files) :
  commands = []
  non_executable = []
  unrecognized = []
  for file_name in files :
    if (file_name.endswith(".py")) :
      commands.append("libtbx.python %s" % file_name)
    elif (file_name.endswith(".sh")) or (file_name.endswith(".csh")) :
      if (not os.access(file_name, os.X_OK)) :
        non_executable.append(file_name)
      else :
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
