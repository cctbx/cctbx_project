
from __future__ import absolute_import, division, print_function
from libtbx import easy_pickle
from libtbx.utils import Abort
from libtbx import object_oriented_patterns as oop
from libtbx import adopt_init_args
import libtbx.callbacks # import dependency
from six.moves import queue
import threading
import warnings
import traceback
import signal
import time
import os
import sys

class thread_with_callback_and_wait(threading.Thread):

  def __init__(self,
          run,
          callback,
          first_callback=None,
          run_args=(),
          run_kwds={}):
    self._run = run
    self._run_args = run_args
    self._run_kwds = dict(run_kwds) # copy, to avoid side effects
    self._run_kwds['callback'] = self._callback_proxy
    self._callback = callback
    self._first_callback = first_callback
    self._caller_wait = queue.Queue(0)
    self._queue = queue.Queue(0)
    threading.Thread.__init__(self)

  def run(self):
    self._run(*self._run_args, **self._run_kwds)

  def start_and_wait_for_first_callback(self):
    self.start()
    self._caller_wait.get()

  def _callback_proxy(self, *args, **kwds):
    if (self._first_callback is not None):
      result = self._first_callback(*args, **kwds)
      self._first_callback = None
      self._caller_wait.put(None)
    else:
      result = self._callback(*args, **kwds)
    if (result == False):
      return False
    last_iteration = self._queue.get()
    if (last_iteration):
      return False
    return result

  def resume(self, last_iteration=False):
    if (self.is_alive()):
      self._queue.put(last_iteration)
    return self

null_callback = oop.null()

class queue_monitor_thread(threading.Thread):
  def __init__(self, q, callback, sleep_interval=1.0):
    self.q = q
    self.cb = callback
    self._exit = False
    self._sleep_interval = sleep_interval
    threading.Thread.__init__(self)

  def run(self):
    t = self._sleep_interval
    while not self._exit :
      if (self.q.qsize() > 0):
        result = self.q.get(timeout=1)
        self.cb(result)
      else :
        time.sleep(t)

  def exit(self):
    self._exit = True

class file_monitor_thread(threading.Thread):
  def __init__(self, dir_name, file_names, callback, sleep_interval=1.0):
    assert os.path.isdir(dir_name)
    self.cb = callback
    self.dir_name = dir_name
    self.file_names = file_names
    self._exit = False
    self._sleep_interval = sleep_interval
    threading.Thread.__init__(self)

  def run(self):
    t = self._sleep_interval
    existing_files = []
    pending_files = []
    while not self._exit :
      files = os.listdir(self.dir_name)
      for file_name in files :
        if (file_name in existing_files):
          continue
        elif (file_name in pending_files):
          data = easy_pickle.load(os.path.join(self.dir_name, file_name))
          existing_files.append(file_name)
          #pending_files.remove(file_name)
          self.cb(data)
        elif (file_name in self.file_names):
          pending_files.append(file_name)
      time.sleep(t)

  def exit(self):
    self._exit = True

class simple_task_thread(threading.Thread):
  def __init__(self, thread_function, parent_window=None):
    threading.Thread.__init__(self)
    self.f = thread_function
    self.return_value = None
    self._exception = None

  def is_complete(self):
    return (self.return_value is not None)

  def exception_raised(self):
    return (self._exception is not None)

  def run(self):
    try :
      self.return_value = self.f()
    except Exception as e :
      sys.stderr.write(str(e))
      self._exception = e

  def get_error(self):
    return str(self._exception)

class child_process_message(object):
  def __init__(self, message_type, data):
    adopt_init_args(self, locals())

class child_process_pipe(object):
  def __init__(self, connection_object):
    adopt_init_args(self, locals())

  def send_confirm_abort(self, info=None):
    message = child_process_message(message_type="aborted", data=info)
    self.connection_object.send(message)

  def send(self, *args, **kwds):
    return self.connection_object.send(*args, **kwds)

  def recv(self, *args, **kwds):
    return self.connection_object.recv(*args, **kwds)

# fake filehandle, sends data up pipe to parent process
class stdout_pipe(object):
  def __init__(self, connection):
    self._c = connection
    self._data = ""

  def write(self, data):
    self._data += data
    self.flush() # this needs to be done immediately for some reason

  def flush(self):
    self._flush()

  def _flush(self):
    try :
      if self._data != "" :
        message = child_process_message(message_type="stdout", data=self._data)
        self._c.send(message)
        self._data = ""
    except Exception as e :
      sys.__stderr__.write("Exception in stdout_pipe: %s\n" % str(e))

  def close(self):
    pass

wait_before_flush = 0.2 # minimum time between send() calls
max_between_flush = 5

# this slows down the output so it won't stall a GUI
class stdout_pipe_buffered(stdout_pipe):
  def __init__(self, *args, **kwds):
    stdout_pipe.__init__(self, *args, **kwds)
    self._last_t = time.time()

  def write(self, data):
    self._data += data
    t = time.time()
    if (t >= (self._last_t + max_between_flush)):
      self.flush()

  def flush(self):
    t = time.time()
    if (t >= (self._last_t + wait_before_flush)):
      self._flush()
      self._last_t = t

  def close(self):
    pass

  def __del__(self):
    if self._data != "" :
      self._flush()

try:
  import multiprocessing
except ImportError:
  class process_with_callbacks(object):
    def __init__(self, *args, **kwds):
      raise ImportError("The multiprocessing module is not available.")
else:
  class _process_with_stdout_redirect(multiprocessing.Process):
    def __init__(self, group=None, target=None, name=None, args=(), kwargs={},
        connection=None, buffer_stdout=True):
      multiprocessing.Process.__init__(self, group, target, name, args, kwargs)
      self._c = connection
      if buffer_stdout :
        self._stdout = stdout_pipe_buffered(connection)
      else :
        self._stdout = stdout_pipe(connection)

    def run(self):
      import libtbx.callbacks # import dependency
      libtbx.call_back.add_piped_callback(self._c)
      old_stdout = sys.stdout
      sys.stdout = self._stdout
      message = None
      old_showwarning = warnings.showwarning
      warnings.showwarning = libtbx.call_back.showwarning
      try :
        try :
          if self._target :
            return_value = self._target(self._args, self._kwargs, self._c)
            message = child_process_message(message_type="return",
                                            data=return_value)
        except Abort :
          message = child_process_message(message_type="abort", data=None)
        except Exception as e :
          if e.__class__.__module__ == "Boost.Python" :
            e = RuntimeError("Boost.Python.%s: %s" % (e.__class__.__name__,
              str(e)))
          elif hasattr(e, "reset_module"):
            e.reset_module()
          traceback_str = "\n".join(traceback.format_tb(sys.exc_info()[2]))
          message = child_process_message(message_type="exception",
                                          data=(e, traceback_str))
      finally :
        if message is not None :
          self._stdout._flush()
          self._c.send(message)
        warnings.showwarning = old_showwarning
      sys.stdout = old_stdout

  #XXX: target functions must use this call signature!!!
  # I normally just write a very short wrapper function that invokes the
  # command I actually care about with args/kwds.
  def _process_target_function(args, kwds, connection):
    connection.send(True)

  # not really a process, but a wrapper for one.
  class process_with_callbacks(threading.Thread):
    def __init__(self,
        target,
        args=(),
        kwargs={},
        callback_stdout=null_callback,
        callback_final=null_callback,
        callback_err=null_callback,
        callback_abort=null_callback,
        callback_other=null_callback,
        callback_pause=null_callback,  # XXX experimental
        callback_resume=null_callback, # XXX experimental
        buffer_stdout=True,
        sleep_after_start=None):
      threading.Thread.__init__(self)
      self._target = target
      self._args = args
      self._kwargs = dict(kwargs)
      assert (hasattr(callback_stdout, "__call__") and
              hasattr(callback_final, "__call__") and
              hasattr(callback_err, "__call__") and
              hasattr(callback_abort, "__call__") and
              hasattr(callback_other, "__call__") and
              hasattr(callback_pause, "__call__") and
              hasattr(callback_resume, "__call__"))
      self._cb_stdout = callback_stdout
      self._cb_final  = callback_final
      self._cb_err    = callback_err
      self._cb_abort  = callback_abort
      self._cb_other  = callback_other
      self._cb_pause  = callback_pause
      self._cb_resume = callback_resume
      self._buffer_stdout = buffer_stdout
      self._abort = False
      self._killed = False
      self._completed = False
      self._error = False
      self._pid = None
      self._child_process = None
      assert ((sleep_after_start is None) or
              isinstance(sleep_after_start, int) or
              isinstance(sleep_after_start, float))
      self._sleep_after_start = sleep_after_start

    def abort(self, force=False):
      self._abort = True
      if (force):
        self.kill()

    # XXX experimental
    def kill(self):
      if (self._child_process is not None):
        try :
          self._child_process.terminate()
          self.join(0.1)
          self._killed = True
        except OSError as e :
          print(e)
        else :
          self._cb_abort()

    def send_signal(self, signal_number) : # XXX experimental
      """
      Signals the process using os.kill, which despite the name, can also
      pause or resume processes on Unix.
      """
      assert (self._child_process is not None)
      try :
        os.kill(self._child_process.pid, signal_number)
      except OSError as e :
        print(e)
        if (not self._child_process.is_alive()):
          self._cb_abort() # XXX not sure if this is ideal
        return False
      else :
        return True

    def pause(self):
      if sys.platform != "win32":
        if (self.send_signal(signal.SIGSTOP)):
          self._cb_pause()
      else:
        import psutil # really should use psutil globally as it is platform independent
        p = psutil.Process(self._child_process.pid)
        p.suspend()
        self._cb_pause()

    def resume(self):
      if sys.platform != "win32":
        if (self.send_signal(signal.SIGCONT)):
          self._cb_resume()
      else:
        import psutil
        p = psutil.Process(self._child_process.pid)
        p.resume()
        self._cb_resume()

    def run(self):
      if (self._sleep_after_start is not None):
        time.sleep(self._sleep_after_start)
      child_process = None
      parent_conn, child_conn = multiprocessing.Pipe()
      try :
        child_process = _process_with_stdout_redirect(
          target        = self._target,
          args          = self._args,
          kwargs        = self._kwargs,
          connection    = child_process_pipe(connection_object=child_conn),
          buffer_stdout = self._buffer_stdout)
        os.environ["OMP_NUM_THREADS"] = "1"
        self._child_process = child_process
        child_process.start()
      except Exception as e :
        sys.__stderr__.write("Error starting child process: %s\n" % str(e))
        child_process = None
      while child_process is not None :
        if (self._killed) : break
        if self._abort :
          child_process.terminate()
          self._cb_abort()
          break
        if not parent_conn.poll(1):
          continue
        pipe_output = parent_conn.recv()
        if isinstance(pipe_output, child_process_message):
          message = pipe_output
          (error, traceback_info) = (None, None)
          if message.message_type == "aborted" :
            child_process.terminate()
            self._cb_abort()
            break
          elif message.message_type == "stdout" :
            self._cb_stdout(message.data)
          elif message.message_type == "return" :
            self._cb_final(message.data)
            self._completed = True
            break
          elif message.message_type == "exception" :
            (error, traceback_info) = message.data
            self._cb_err(error, traceback_info)
            self._error = True
            child_process.terminate()
            break
        else :
          self._cb_other(pipe_output)
      if child_process is not None and child_process.is_alive():
        child_process.join()
