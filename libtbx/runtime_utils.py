
# detached process management (used in Phenix GUI)
# This is intended to imitate running a process using libtbx.thread_utils,
# but split across truly independent processes (potentially on different
# systems).  I use something like this to start a server:
#
#   easy_run.call("libtbx.start_process run.pkl &")

# FIXME this duplicates code in libtbx.thread_utils; it is also not especially
# well tested except to the extent it is used daily in the Phenix GUI

from __future__ import absolute_import, division, print_function
from six.moves import cStringIO as StringIO
from libtbx.utils import Sorry, Abort, multi_out, host_and_user
from libtbx import easy_pickle
from libtbx import adopt_init_args, group_args
import libtbx.load_env
import libtbx.phil
import traceback
import signal
import stat
import time
import os
import sys
from six.moves import range

process_master_phil = libtbx.phil.parse("""
run_file = None
  .type = path
prefix = None
  .type = str
output_dir = None
  .type = path
tmp_dir = None
  .type = path
debug = False
  .type = bool
timeout = 200
  .type = int
buffer_stdout = False
  .type = bool
fsync = True
  .type = bool
""")

class simple_target(object):
  def __init__(self, args, output_dir=None):
    adopt_init_args(self, locals())
    if output_dir is None :
      self.output_dir = os.getcwd()

  def __call__(self):
    return True

class target_with_save_result(object):
  def __init__(self, args, file_name, output_dir=None, log_file=None,
      job_title=None):
    from six import string_types
    assert isinstance(file_name, string_types)
    assert (isinstance(args, list) or isinstance(args, tuple))
    assert (output_dir is None) or isinstance(output_dir, string_types)
    assert (log_file is None) or isinstance(log_file, string_types)
    adopt_init_args(self, locals())
    if (output_dir is None):
      self.output_dir = os.getcwd()
    self._out = None

  def __call__(self):
    if (self.log_file is not None):
      log = open(self.log_file, "w")
      new_out = multi_out()
      new_out.register("log", log)
      new_out.register("stdout", sys.stdout)
      sys.stdout = new_out
      self._out = new_out
    result = self.run()
    easy_pickle.dump(self.file_name, result)
    if (self._out is not None) and (not getattr(self._out, "closed", False)):
      self._out.flush()
      # FIXME
      #self._out.close()
    return result

  def run(self):
    raise NotImplementedError()

class detached_process_driver(object):
  def __init__(self, target):
    adopt_init_args(self, locals())

  def __call__(self):
    result = self.target()
    return result

class detached_process_driver_mp(detached_process_driver):
  def __call__(self, args, kwds, child_conn):
    result = self.target()
    return result

class detached_base(object):
  def __init__(self, params):
    adopt_init_args(self, locals())
    self._accumulated_callbacks = []
    if params.prefix is None :
      params.prefix = ""
    if params.tmp_dir is not None :
      self.set_file_names(params.tmp_dir)
    elif params.output_dir is not None :
      self.set_file_names(params.output_dir)
    else :
      self.set_file_names(os.getcwd())

  def set_file_names(self, tmp_dir):
    prefix = os.path.join(tmp_dir, self.params.prefix)
    self.start_file = prefix + ".libtbx_start"
    self.stdout_file = prefix + ".libtbx_stdout"
    self.error_file =  prefix + ".libtbx_error"
    self.stop_file = prefix + ".libtbx_STOP"
    self.abort_file =  prefix + ".libtbx_abort"
    self.result_file = prefix + ".libtbx_result"
    self.info_file = prefix + ".libtbx_info"
    self.state_file =  prefix + ".libtbx_state"
    self.info_lock = self.info_file + ".LOCK"
    self.state_lock = self.state_file + ".LOCK"
    self.prefix = prefix

  def isAlive(self):
    return False

  def callback_start(self, data):
    pass

  def callback_stdout(self, data):
    pass

  def callback_error(self, error, traceback_info):
    pass

  def callback_abort(self):
    pass

  def callback_final(self, result):
    pass

  def callback_other(self, status):
    pass

class stdout_redirect(object):
  def __init__(self, handler):
    adopt_init_args(self, locals())

  def write(self, data):
    self.handler.callback_stdout(data)

  def flush(self):
    pass

  def close(self):
    pass

class detached_process_server(detached_base):
  def __init__(self, target, *args, **kwds):
    detached_base.__init__(self, *args, **kwds)
    self.target = target
    assert hasattr(self.target, "__call__")

  # XXX support for libtbx.queueing_system_utils.generic.Job
  def __call__(self, *args, **kwds):
    return self.run()

  def run(self):
    self.callback_start()
    self._stdout = multi_out()
    self._tmp_stdout = open(self.stdout_file, "w")
    self._stdout.register("Communication log", self._tmp_stdout)
    old_stdout = sys.stdout
    sys.stdout = stdout_redirect(self)
    import libtbx.callbacks
    libtbx.call_back.register_handler(self.callback_wrapper)
    try :
      return_value = self.target()
    except Abort : # FIXME why is this not working properly?
      self.callback_abort()
    except Exception as e :
      print(type(e).__name__, file=sys.stderr)
      if (type(e).__name__ == "Abort"):
        self.callback_abort()
      else :
        if e.__class__.__module__ == "Boost.Python" :
          e = RuntimeError("Boost.Python.%s: %s" % (e.__class__.__name__,
            str(e)))
        elif hasattr(e, "reset_module"):
          e.reset_module()
        traceback_str = "\n".join(traceback.format_tb(sys.exc_info()[2]))
        try:
          self.callback_error(e, traceback_str)
        except Exception as ee:
          self.callback_error(str(e), traceback_str)
    else :
      #time.sleep(1)
      self.callback_final(return_value)
    sys.stdout = old_stdout

  def callback_wrapper(self, message, data, accumulate=True, cached=True):
    if cached :
      self.callback_other(data=group_args(
        message=message,
        data=data,
        accumulate=accumulate,
        cached=cached))

  def callback_start(self, data=None):
    info = host_and_user()
    assert (info.pid is not None)
    f = open(self.start_file, "w")
    host_name = info.get_host_name()
    if (host_name is None) and (sys.platform == "darwin"):
      host_name = os.uname()[1]
    f.write("%s %d" % (host_name, info.pid))
    f.close()

  def callback_stdout(self, data):
    self._stdout.write(data)
    self._stdout.flush()
    if self.params.fsync :
      os.fsync(self._tmp_stdout.fileno())
    if os.path.isfile(self.stop_file):
      raise Abort()

  def callback_error(self, error, traceback_info):
    self.cleanup()
    easy_pickle.dump(self.error_file, (error, traceback_info))

  def callback_abort(self):
    self.cleanup()
    easy_pickle.dump(self.abort_file, True)

  def callback_final(self, result):
    self.cleanup()
    easy_pickle.dump(self.result_file, result)

  def callback_pause(self) : # TODO
    pass

  def callback_resume(self) : # TODO
    pass

  def wait_for_lock_to_vanish(self, local_wait_time = 0.1,
      max_local_wait_time = 5):
    # Catch case where jobs are very short and returns overlap...wait for lock
    # XXX looks like info_lock was never actually implemented, just set up.
    t0 = time.time()
    while os.path.isfile(self.info_lock):
        if time.time() - t0 > max_local_wait_time:
          break  # just crash
        else:
          time.sleep(local_wait_time)

  def callback_other(self, data, local_wait_time = 0.1,
     max_local_wait_time = 1):
    if not data.cached :
      return
    if data.accumulate :
      self._accumulated_callbacks.append(data)
      self.wait_for_lock_to_vanish()
      try:
        touch_file(self.info_lock)
        easy_pickle.dump(self.info_file, self._accumulated_callbacks)
        os.remove(self.info_lock)
      except Exception as e:
        pass # failed, can happen if user deletes directory
    else :
      self.wait_for_lock_to_vanish()
      try:
        touch_file(self.state_lock)
        easy_pickle.dump(self.state_file, data)
        os.remove(self.state_lock)
      except Exception as e:
        pass # failed, can happen if user deletes directory

  def cleanup(self):
    self._stdout.flush()
    self._stdout.close()

# TODO pause/resume?
class detached_process_client(detached_base):
  def __init__(self, *args, **kwds):
    detached_base.__init__(self, *args, **kwds)
    self._logfile = None
    self._info_mtime = 0.0 # time.time()
    self._state_mtime = 0.0 # time.time()
    self.running = False
    self.finished = False
    self._process_host = None
    self._process_pid = None
    self.update_progress = True

  def __del__(self):
    if self._logfile is not None:
      self._logfile.close()

  def isAlive(self):
    return (not self.finished)

  def is_started(self):
    return self.running

  def run(self):
    timeout = self.params.timeout
    while True :
      self.update()
      if self.finished :
        break
      else :
        time.sleep(timeout * 0.001)
    return True

  def update(self):
    if not self.running and os.path.exists(self.start_file):
      self.running = True
      with open(self.start_file, "r") as f:
        data = f.read()
      try :
        host, pid = data.split()
        self._process_host = host
        self._process_pid = int(pid)
      except Exception as e :
        print("Error acquiring runtime info:")
        print(e)
      self.callback_start(data)
    if self.update_progress :
      self.check_stdout()
      self.check_status()
    if os.path.exists(self.error_file):
      try :
        (error, traceback_info) = easy_pickle.load(self.error_file)
      except (EOFError, TypeError):
        pass
      else :
        self.callback_error(error, traceback_info)
    elif os.path.exists(self.abort_file):
      self.callback_abort()
    elif os.path.exists(self.result_file):
      max_retries = 5
      error_text = None
      for retry in range(max_retries):
        try:
          result = easy_pickle.load(self.result_file)
        except Exception:
          pass
        else:
          time.sleep(1)
          self.check_stdout()
          self.check_status()
          self.callback_final(result)
          break
        time.sleep(retry+1)
      if retry == max_retries-1:
        print("There was an error with loading '%s'." % self.result_file)
    else :
      self.finished = False
      return
    self.finished = True

  def check_stdout(self):
    if self._logfile is None and os.path.exists(self.stdout_file):
      self._logfile = open(self.stdout_file, "r")
      stat = os.stat(self.stdout_file)
      self._logfile.seek(stat.st_size)
    if self._logfile is not None :
      last = self._logfile.tell()
      data = self._logfile.read()
      if data == '' :
        self._logfile.seek(last)
      else :
        self.callback_stdout(data)

  def reset_logfile(self):
    if self._logfile is not None :
      self._logfile.seek(0)

  def check_status(self):
    if os.path.exists(self.info_file):
      mtime = os.path.getmtime(self.info_file)
      if mtime > self._info_mtime and not os.path.isfile(self.info_lock):
        self._info_mtime = mtime
        try :
          accumulated_status = easy_pickle.load(self.info_file)
        except KeyboardInterrupt :
          raise
        except EOFError :
          pass
        except Exception as e :
          print(e)
        else :
          n_cb = len(accumulated_status)
          n_cb_old = len(self._accumulated_callbacks)
          for i in range(n_cb_old, n_cb):
            new_cb = accumulated_status[i]
            self._accumulated_callbacks.append(new_cb)
            self.callback_other(new_cb)
    if os.path.exists(self.state_file):
      mtime = os.path.getmtime(self.state_file)
      if mtime > self._state_mtime and not os.path.isfile(self.state_lock):
        self._state_mtime = mtime
        try :
          current_status = easy_pickle.load(self.state_file)
        except KeyboardInterrupt :
          raise
        except EOFError :
          pass
        except Exception as e :
          print(e)
        else :
          self.callback_other(current_status)

  # XXX in practice the phenix GUI only sets force=True for jobs running on the
  # same machine; qdel is pretty thorough anyway.
  def abort(self, force=None):
    touch_file(self.stop_file)
    if (force) and (not None in [self._process_host, self._process_pid]):
      info = host_and_user()
      if (info.get_host_name() == self._process_host):
        os.kill(self._process_pid, signal.SIGTERM)
        self.running = False
        self.callback_abort()

  # XXX See also libtbx.thread_utils implementation
  def send_signal(self, signal_number) : # XXX experimental
    """
    Signals the process using os.kill, which despite the name, can also
    pause or resume processes on Unix.
    """
    assert (self._process_pid is not None) and (sys.platform != "win32")
    try :
      os.kill(self._process_pid, signal_number)
    except OSError as e :
      print(e)
      #self.callback_abort()
      return False
    else :
      return True

  def pause(self) : # XXX experimental, Unix only
    if (self.send_signal(signal.SIGSTOP)):
      self.callback_pause()

  def resume(self) : # XXX experimental, Unix only
    if (self.send_signal(signal.SIGCONT)):
      self.callback_resume()

  def callback_pause(self) : # TODO
    pass

  def callback_resume(self) : # TODO
    pass

  def purge_files(self):
    files = ["start","stdout","error","stop","abort","result","info","state"]
    for fn in files :
      file_name = getattr(self, "%s_file" % fn)
      if os.path.exists(file_name):
        try :
          os.remove(file_name)
        except Exception as e :
          print(e)

def touch_file(file_name):
  f = open(file_name, "w").close()

def write_params(params, file_name):
  param_phil = process_master_phil.format(python_object=params)
  f = open(file_name, "w")
  param_phil.show(out=f)
  f.close()

def write_run_script(file_name, cmds):
  f = open(file_name, "w")
  os.fchmod(f.fileno(),
    stat.S_IRUSR|stat.S_IWUSR|stat.S_IXUSR|stat.S_IRGRP|stat.S_IROTH)
  f.write("#!/bin/sh\n\n")
  use_cctbx_setpaths = True
  if "PHENIX" in os.environ and not "PHENIX_CUSTOM_ENV" in os.environ :
    env_file = os.path.join(os.environ["PHENIX"], "phenix_env.sh")
    if os.path.isfile(env_file):
      f.write("source %s\n" % env_file)
      use_cctbx_setpaths = False
  if use_cctbx_setpaths :
    f.write("source %s\n" % libtbx.env.under_build("setpaths.sh"))
  f.write("%s" % " ".join(cmds))
  f.close()

# XXX command-line launcher
def run(args):
  user_phil = []
  for arg in args :
    if os.path.isfile(arg):
      file_name = os.path.abspath(arg)
      base, ext = os.path.splitext(file_name)
      if ext in [".params", ".eff", ".def", ".phil"] :
        user_phil.append(libtbx.phil.parse(file_name=file_name))
      elif ext in [".pkl", ".pickle"] :
        input_string = "run_file = %s" % arg
        user_phil.append(libtbx.phil.parse(input_string))
    else :
      try :
        arg_phil = libtbx.phil.parse(arg)
      except RuntimeError as e :
        print(e)
      else :
        user_phil.append(arg_phil)
  working_phil = process_master_phil.fetch(sources=user_phil)
  params = working_phil.extract()
  if params.run_file is None :
    working_phil.show()
    raise Sorry("Pickled target function run_file not defined.")
  target = easy_pickle.load(params.run_file)
  server = detached_process_server(target, params=params)
  server.run()

########################################################################
# testing classes (see tst_runtime_utils.py for usage)
class simple_client(detached_process_client):
  def __init__(self, *args, **kwds):
    self.n_cb = 0
    self.out = StringIO()
    detached_process_client.__init__(self, *args, **kwds)

  def callback_error(self, error, traceback_info):
    raise error

  def callback_aborted(self):
    raise Sorry("aborted as planned.")

  def callback_stdout(self, data):
    self.out.write(data)
    #for line in data.splitlines():

  def callback_other(self, data):
    self.n_cb += 1

  def callback_final(self, result):
    self.result = result

class simple_run(object):
  def __init__(self, output_dir):
    adopt_init_args(self, locals())

  def __call__(self):
    pu_total = 0
    for run in range(0, 4):
      (x, n) = (0.1 * (run+1) , 20000)
      mu = 10.0
      pu = 0.0
      pol =[0] * 100
      r = list(range(0,100))
      libtbx.call_back("run %d" % run, None, accumulate=True)
      time.sleep(1)
      for i in range(0,n):
        for j in r:
          pol[j] = mu = (mu + 2.0) / 2.0
        su = 0.0
        for j in r:
          su = x * su + pol[j]
        pu = pu + su
      pu_total += pu
      libtbx.call_back("current_total", pu, accumulate=False)
      print("current is %f" % pu)
    return pu_total

class simple_func(object):
  def __init__(self, x):
    self.x = x
  def __call__(self):
    print(self.x)
