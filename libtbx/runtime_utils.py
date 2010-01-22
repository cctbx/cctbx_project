
import libtbx.phil
from libtbx.utils import Sorry
from libtbx import easy_pickle, thread_utils
from libtbx import adopt_init_args
import sys, os, time
import cStringIO

process_master_phil = libtbx.phil.parse("""
run_file = None
  .type = path
prefix = None
  .type = str
tmp_dir = None
  .type = path
debug = False
  .type = bool
timeout = 200
  .type = int
buffer_stdout = True
  .type = bool
""")

class detached_process_driver (object) :
  def __init__ (self, output_dir, target) :
    adopt_init_args(self, locals())

  def __call__ (self, args, kwds, child_conn) :
    os.chdir(self.output_dir)
    import libtbx.callbacks
    callback = libtbx.callbacks.piped_callback(child_conn)
    libtbx.call_back.register_handler(callback)
    result = self.target(args, kwds, child_conn)
    return result

class detached_base (object) :
  def __init__ (self, params) :
    adopt_init_args(self, locals())
    self._accumulated_callbacks = []
    if params.prefix is None :
      params.prefix = ""
    if params.tmp_dir is not None :
      self.set_file_names(params.tmp_dir)

  def set_file_names (self, tmp_dir) :
    prefix = os.path.join(tmp_dir, self.params.prefix)
    self.start_file = os.path.join(tmp_dir, prefix + ".libtbx_start")
    self.stdout_file = prefix + ".libtbx_stdout"
    self.error_file =  prefix + ".libtbx_error"
    self.abort_file =  prefix + ".libtbx_abort"
    self.result_file = prefix + ".libtbx_result"
    self.info_file = prefix + ".libtbx_info"
    self.state_file =  prefix + ".libtbx_state"
    self.info_lock = self.info_file + ".LOCK"
    self.state_lock = self.state_file + ".LOCK"
    self.prefix = prefix

  def isAlive (self) :
    return False

  def callback_stdout (self, data) :
    pass

  def callback_error (self, error, traceback_info) :
    pass

  def callback_abort (self) :
    pass

  def callback_final (self, result) :
    pass

  def callback_other (self, status) :
    pass

class detached_process_server (detached_base) :
  def __init__ (self, *args, **kwds) :
    detached_base.__init__(self, *args, **kwds)
    target = easy_pickle.load(self.params.run_file)
    assert isinstance(target, detached_process_driver)
    self.params.tmp_dir = target.output_dir
    self.set_file_names(self.params.tmp_dir)
    self.libtbx_process = thread_utils.process_with_callbacks(
      target=target,
      callback_stdout=self.callback_stdout,
      callback_final=self.callback_final,
      callback_err=self.callback_error,
      callback_abort=self.callback_abort,
      callback_other=self.callback_other,
      buffer_stdout=False) #self.params.buffer_stdout)
    f = open(self.start_file, "w", 0)
    f.write("1")
    f.close()
    self._stdout = open(self.stdout_file, "w")
    self.libtbx_process.start()

  def isAlive (self) :
    return self.libtbx_process.isAlive()

  def callback_stdout (self, data) :
    self._stdout.write(data)
    self._stdout.flush()
    os.fsync(self._stdout.fileno())

  def callback_error (self, error, traceback_info) :
    easy_pickle.dump(self.error_file, (error, traceback_info))
    self.cleanup()

  def callback_abort (self) :
    easy_pickle.dump(self.abort_file, True)
    self.cleanup()

  def callback_final (self, result) :
    easy_pickle.dump(self.result_file, result)
    self.cleanup()

  def callback_other (self, data) :
    if not data.cached :
      return
    if data.accumulate :
      self._accumulated_callbacks.append(data)
      touch_file(self.info_lock)
      easy_pickle.dump(self.info_file, self._accumulated_callbacks)
      os.remove(self.info_lock)
    else :
      touch_file(self.state_lock)
      easy_pickle.dump(self.state_file, data)
      os.remove(self.state_lock)

  def cleanup (self) :
    self._stdout.flush()
    self._stdout.close()

class detached_process_client (detached_base) :
  def __init__ (self, *args, **kwds) :
    detached_base.__init__(self, *args, **kwds)
    self._logfile = None
    self._info_mtime = time.time()
    self._state_mtime = time.time()
    self.running = False
    self.finished = False

  def isAlive (self) :
    return (not self.finished)

  def is_started (self) :
    return self.running

  def run (self) :
    timeout = self.params.timeout
    while True :
      if not self.running and os.path.exists(self.start_file) :
        self.running = True
      self.check_stdout()
      self.check_status()
      if os.path.exists(self.error_file) :
        (error, traceback_info) = easy_pickle.load(self.error_file)
        self.callback_error(error, traceback_info)
        break
      elif os.path.exists(self.abort_file) :
        self.callback_abort()
        break
      elif os.path.exists(self.result_file) :
        result = easy_pickle.load(self.result_file)
        self.check_stdout()
        self.callback_final(result)
        break
      time.sleep(timeout * 0.001)
    self.finished = True
    return True

  def check_stdout (self) :
    if self._logfile is None and os.path.exists(self.stdout_file) :
      self._logfile = open(self.stdout_file, "r", 0)
    if self._logfile is not None :
      last = self._logfile.tell()
      data = self._logfile.read()
      if data == '' :
        self._logfile.seek(last)
      else :
        self.callback_stdout(data)

  def check_status (self) :
    if os.path.exists(self.info_file) :
      mtime = os.path.getmtime(self.info_file)
      if mtime > self._info_mtime and not os.path.isfile(self.info_lock) :
        self._info_mtime = mtime
        try :
          accumulated_status = easy_pickle.load(self.info_file)
        except KeyboardInterrupt :
          raise
        except Exception, e :
          print e
        else :
          n_cb = len(accumulated_status)
          n_cb_old = len(self._accumulated_callbacks)
          for i in range(n_cb_old, n_cb) :
            new_cb = accumulated_status[i]
            self._accumulated_callbacks.append(new_cb)
            self.callback_other(new_cb)
    if os.path.exists(self.state_file) :
      mtime = os.path.getmtime(self.state_file)
      if mtime > self._state_mtime and not os.path.isfile(self.state_lock) :
        self._state_mtime = mtime
        try :
          current_status = easy_pickle.load(self.state_file)
        except KeyboardInterrupt :
          raise
        except Exception, e :
          print e
        else :
          self.callback_other(current_status)

def touch_file (file_name) :
  f = open(file_name, "w")
  f.write("1")
  f.close()

# XXX command-line launcher
def run (args) :
  user_phil = []
  for arg in args :
    if os.path.isfile(arg) :
      file_name = os.path.abspath(arg)
      user_phil.append(libtbx.phil.parse("run_file = \"%s\"" % file_name))
    else :
      try :
        arg_phil = libtbx.phil.parse(arg)
      except RuntimeError :
        pass
      else :
        user_phil.append(arg_phil)
  params = process_master_phil.fetch(sources=user_phil).extract()
  server = detached_process_server(params)

########################################################################
# testing classes (see tst_runtime_utils.py for usage)
class simple_client (detached_process_client) :
  def __init__ (self, *args, **kwds) :
    self.n_cb = 0
    self.out = cStringIO.StringIO()
    detached_process_client.__init__(self, *args, **kwds)

  def callback_error (self, error, traceback_info) :
    raise error

  def callback_aborted (self) :
    raise Sorry("aborted as planned.")

  def callback_stdout (self, data) :
    self.out.write(data)
    #for line in data.splitlines() :

  def callback_other (self, data) :
    self.n_cb += 1

  def callback_final (self, result) :
    self.result = result

class simple_run (object) :
  def __call__ (self, args, kwds, child_conn) :
    pu_total = 0
    for run in range(0, 4) :
      (x, n) = (0.1 * (run+1) , 20000)
      mu = 10.0
      pu = 0.0
      pol =[0] * 100
      r = range(0,100)
      libtbx.call_back("run %d" % run, None, accumulate=True)

      for i in range(0,n):
        for j in r:
          pol[j] = mu = (mu + 2.0) / 2.0
        su = 0.0
        for j in r:
          su = x * su + pol[j]
        pu = pu + su
      pu_total += pu
      libtbx.call_back("current_total", pu, accumulate=False)
      print "current is %f" % pu
    return pu_total
