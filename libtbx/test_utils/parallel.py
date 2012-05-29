import os, sys
import time
from libtbx import easy_run
from libtbx import test_utils
from multiprocessing import Pool

results=[]
def process_callback(arg):
  results.append(arg)

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
  print "\n  Starting command list"
  print "    NProcs :",nprocs
  print "    Cmds   :",len(cmd_list)
  if nprocs>1:
    pool = Pool(processes=nprocs)

  success = 0
  failure = 0
  for command in cmd_list:
    if nprocs>1:
      rc = pool.apply_async(
        run_command,
        [command, True],
        callback=process_callback,
        )
    else:
      rc = run_command(command, verbose=True)
      display_result(rc, log=log)
      success+=1
      if (len(result.stderr_lines) != 0):
        failure+=1

  if nprocs>1:
    pool.close()
    pool.join()
    print '\nProcesses have joined : %d\n' % len(results)
    for result in results:
      display_result(result, log=log)
      success+=1
      if (len(result.stderr_lines) != 0):
        failure+=1
  print 'Done with output'
  print "  NProcs    :",nprocs
  print "  Tests run :",success
  print "  Failures  :",failure
  print >> log, '\n\nDone with output'
  print >> log, "  NProcs    :",nprocs
  print >> log, "  Tests run :",success
  print >> log, "  Failures  :",failure

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
