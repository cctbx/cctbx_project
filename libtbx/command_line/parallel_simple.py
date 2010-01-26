def run_one_cmd(cmd_info):
  import sys
  sys.stdout.flush()
  import time
  t0 = time.time()
  def fmt_time():
    return "wall clock time: %.2f s" % (time.time() - t0)
  from libtbx import easy_run
  if (cmd_info.index == 0):
    easy_run.call(command=cmd_info.cmd)
    print fmt_time()
    print
  else:
    buffers = easy_run.fully_buffered(
      command=cmd_info.cmd,
      join_stdout_stderr=True,
      stdout_splitlines=False)
    f = open(cmd_info.log, "w")
    f.write(buffers.stdout_buffer)
    print >> f, fmt_time()
    del f
  sys.stdout.flush()

def run(args):
  import libtbx.utils
  if (len(args) not in [1,2]):
    import libtbx.load_env
    raise libtbx.utils.Usage(
      "%s file_listing_commands [#processors]" %
        libtbx.env.dispatcher_name)
  file_listing_commands = args[0]
  if (len(args) == 1):
    max_proc = None
  else:
    max_proc = int(args[1])
    assert max_proc > 0
  import libtbx
  import os
  op = os.path
  cmd_infos = []
  for line in open(op.expandvars(file_listing_commands)).read().splitlines():
    ll = line.lstrip()
    if (ll.startswith("#")): continue
    if (ll.startswith("set ")): continue
    index = len(cmd_infos)
    cmd_infos.append(libtbx.group_args(
      index=index,
      log="log%03d"%index,
      cmd=line))
  for cmd_info in cmd_infos:
    libtbx.utils.remove_files(cmd_info.log)
  import libtbx.introspection
  n_proc = min(len(cmd_infos), libtbx.introspection.number_of_processors())
  if (max_proc is not None):
    n_proc = min(max_proc, n_proc)
  print "Number of processors:", n_proc
  print
  import sys
  sys.stdout.flush()
  show_times = libtbx.utils.show_times(time_start="now")
  def show_log(cmd_info):
    if (not op.isfile(cmd_info.log)):
      if (cmd_info.index != 0):
        print "MISSING:", cmd_info.log
    else:
      lines = open(cmd_info.log).read().splitlines()
      if (len(lines) > 10):
        print "@BEGIN"
      if (len(lines) != 0):
        print "\n".join(lines)
      if (len(lines) > 10):
        print "@END"
  if (n_proc < 2):
    for cmd_info in cmd_infos:
      print "command:", cmd_info.cmd
      run_one_cmd(cmd_info=cmd_info)
      show_log(cmd_info=cmd_info)
  else:
    import multiprocessing
    mp_pool = multiprocessing.Pool(processes=n_proc)
    mp_pool.map(run_one_cmd, cmd_infos)
    for cmd_info in cmd_infos:
      if (cmd_info.index == 0): continue
      print "command:", cmd_info.cmd
      show_log(cmd_info=cmd_info)
  print
  show_times()
  print
  sys.stdout.flush()

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
