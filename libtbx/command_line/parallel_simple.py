import time
import sys, os
op = os.path

def show_traceback(file):
  import traceback
  print >> file
  traceback.print_exc(file=file)
  print >> file

def fmt_time(t0):
  return "JOB wall clock time: %.2f s" % (time.time() - t0)

def run_one_cmd(cmd_info):
  sys.stdout.flush()
  t0 = time.time()
  from libtbx import easy_run
  if (cmd_info.index == 0):
    try:
      easy_run.call(command=cmd_info.cmd)
    except: # intentional
      show_traceback(file=sys.stdout)
    print fmt_time(t0)
    print
  else:
    sio = None
    try:
      buffers = easy_run.fully_buffered(
        command=cmd_info.cmd,
        join_stdout_stderr=True,
        stdout_splitlines=False)
    except: # intentional
      from cStringIO import StringIO
      sio = StringIO()
      show_traceback(file=sio)
    else:
      f = open(cmd_info.log, "w")
      f.write(buffers.stdout_buffer)
      if (sio is not None):
        f.write(sio.getvalue())
      print >> f, fmt_time(t0)
      del f
  sys.stdout.flush()

def run_in_dir(cmd_info):
  d = op.dirname(cmd_info.log)
  os.mkdir(d)
  os.chdir(d)
  from libtbx.command_line import printenv
  printenv.show(out=open("os_environ_at_start", "w"))
  from libtbx.easy_run import subprocess
  log = open("log", "w")
  t0 = time.time()
  try:
    subprocess.Popen(
      args=cmd_info.cmd,
      shell=True,
      bufsize=-1,
      stdout=log,
      stderr=log,
      universal_newlines=True).wait()
  except: # intentional
    show_traceback(file=log)
  print >> log, fmt_time(t0)
  sys.stdout.flush()

def run(args):
  if (len(args) == 0): args = ["--help"]
  import libtbx.load_env
  from libtbx.option_parser import option_parser
  command_line = (option_parser(
    usage="%s [options] file_listing_commands [...]"
      % libtbx.env.dispatcher_name)
    .option(None, "--dirs",
      action="store",
      type="str",
      help="create a sub-directory for each run.")
    .option("-j", "--jobs",
      action="store",
      type="int",
      help="Maximum number of parallel jobs (default: all CPUs).")
  ).process(args=args)
  co = command_line.options
  #
  import libtbx.introspection
  import libtbx.utils
  from libtbx.utils import Sorry
  from libtbx.str_utils import show_string
  import multiprocessing
  #
  cmd_infos = []
  for file_listing_commands in command_line.args:
    for line in open(op.expandvars(file_listing_commands)).read().splitlines():
      ll = line.lstrip()
      if (ll.startswith("#")): continue
      if (ll.startswith("set ")): continue
      if (ll.startswith("setenv ")):
        flds = ll.split(None, 2)
        if (len(flds) == 2): flds.append("")
        os.environ[flds[1]] = flds[2]
        continue
      index = len(cmd_infos)
      cmd_infos.append(libtbx.group_args(index=index, cmd=line, log=None))
  n_proc = min(len(cmd_infos), libtbx.introspection.number_of_processors())
  if (co.jobs is not None):
    n_proc = max(1, min(co.jobs, n_proc))
  print "Number of processors:", n_proc
  print
  sys.stdout.flush()
  show_times = libtbx.utils.show_times(time_start="now")
  def show_log(cmd_info):
    need_log = (cmd_info.index != 0 or co.dirs is not None)
    if (cmd_info.log is None):
      if (need_log):
        print "MISSING: output of command with index %03d:" % cmd_info.index
        print "  %s" % cmd_info.cmd
    elif (not op.isfile(cmd_info.log)):
      if (need_log):
        print "MISSING:", cmd_info.log
    else:
      lines = open(cmd_info.log).read().splitlines()
      if (len(lines) > 10):
        print "@BEGIN"
      if (len(lines) != 0):
        print "\n".join(lines)
      if (len(lines) > 10):
        print "@END"
  def show_logs():
    for cmd_info in cmd_infos:
      if (cmd_info.index == 0 and co.dirs is None): continue
      print "command:", cmd_info.cmd
      show_log(cmd_info=cmd_info)
  if (co.dirs is None):
    for cmd_info in cmd_infos:
      cmd_info.log = "log%03d" % cmd_info.index
      libtbx.utils.remove_files(cmd_info.log)
    if (n_proc < 2):
      for cmd_info in cmd_infos:
        print "command:", cmd_info.cmd
        run_one_cmd(cmd_info=cmd_info)
        show_log(cmd_info=cmd_info)
    else:
      mp_pool = multiprocessing.Pool(processes=n_proc)
      mp_pool.map(run_one_cmd, cmd_infos)
      show_logs()
  else:
    is_clean = True
    for cmd_info in cmd_infos:
      d = "%s%03d" % (co.dirs, cmd_info.index)
      if (op.exists(d)):
        print >> sys.stderr, "exists already: %s" % show_string(d)
        is_clean = False
      cmd_info.log = op.join(d, "log")
    if (not is_clean):
      raise Sorry(
        "Please remove the existing directories or files,"
        " or use a different --dirs assignment.")
    mp_pool = multiprocessing.Pool(processes=n_proc)
    mp_pool.map(run_in_dir, cmd_infos)
    show_logs()
  print
  show_times()
  print
  sys.stdout.flush()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
