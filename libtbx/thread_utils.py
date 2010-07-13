import sys, traceback, time, os
import Queue
import threading
from libtbx.utils import Sorry, Abort
from libtbx import object_oriented_patterns as oop
from libtbx import group_args, adopt_init_args

class thread_with_callback_and_wait(threading.Thread):

  def __init__ (self,
          run,
          callback,
          first_callback=None,
          run_args=(),
          run_kwds={}) :
    self._run = run
    self._run_args = run_args
    self._run_kwds = dict(run_kwds) # copy, to avoid side effects
    self._run_kwds['callback'] = self._callback_proxy
    self._callback = callback
    self._first_callback = first_callback
    self._caller_wait = Queue.Queue(0)
    self._queue = Queue.Queue(0)
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
    if (self.isAlive()):
      self._queue.put(last_iteration)
    return self

null_callback = oop.null()

class child_process_message (object) :
  def __init__ (self, message_type, data) :
    adopt_init_args(self, locals())

class child_process_pipe (object) :
  def __init__ (self, connection_object) :
    adopt_init_args(self, locals())

  def send_confirm_abort (self, info=None) :
    message = child_process_message(message_type="aborted", data=info)
    self.connection_object.send(message)

  def send (self, *args, **kwds) :
    return self.connection_object.send(*args, **kwds)

  def recv (self, *args, **kwds) :
    return self.connection_object.recv(*args, **kwds)

# fake filehandle, sends data up pipe to parent process
class stdout_pipe (object) :
  def __init__ (self, connection) :
    self._c = connection
    self._data = ""

  def write (self, data) :
    self._data += data
    self.flush() # this needs to be done immediately for some reason

  def flush (self) :
    self._flush()

  def _flush (self) :
    try :
      if self._data != "" :
        message = child_process_message(message_type="stdout", data=self._data)
        self._c.send(message)
        self._data = ""
    except Exception, e :
      sys.__stderr__.write("Exception in stdout_pipe: %s\n" % str(e))

wait_before_flush = 1 # minimum time between send() calls

# this slows down the output so it won't stall a GUI
class stdout_pipe_buffered (stdout_pipe) :
  def __init__ (self, *args, **kwds) :
    stdout_pipe.__init__(self, *args, **kwds)
    self._last_t = time.time()

  def write (self, data) :
    self._data += data

  def flush (self) :
    t = time.time()
    if t >= (self._last_t + wait_before_flush) :
      self._flush()
      self._last_t = t

  def __del__ (self) :
    if self._data != "" :
      self._flush()

try:
  import multiprocessing
except ImportError:
  class process_with_callbacks (object) :
    def __init__ (self, *args, **kwds) :
      raise ImportError("The multiprocessing module is not available.")
else:
  class _process_with_stdout_redirect (multiprocessing.Process) :
    def __init__ (self, group=None, target=None, name=None, args=(), kwargs={},
        connection=None, buffer_stdout=True) :
      multiprocessing.Process.__init__(self, group, target, name, args, kwargs)
      self._c = connection
      if buffer_stdout :
        self._stdout = stdout_pipe_buffered(connection)
      else :
        self._stdout = stdout_pipe(connection)

    def run (self) :
      import libtbx.callbacks
      libtbx.call_back.add_piped_callback(self._c)
      old_stdout = sys.stdout
      sys.stdout = self._stdout
      message = None
      try :
        try :
          if self._target :
            return_value = self._target(self._args, self._kwargs, self._c)
            message = child_process_message(message_type="return",
                                            data=return_value)
        except Abort :
          message = child_process_message(message_type="abort", data=None)
        except Exception, e :
          if e.__class__.__module__ == "Boost.Python" :
            e = RuntimeError("Boost.Python.%s: %s" % (e.__class__.__name__,
              str(e)))
          elif hasattr(e.__class__, "__orig_module__") :
            e.__class__.__module__ = e.__class__.__orig_module__
          #Sorry.reset_module()
          traceback_str = "\n".join(traceback.format_tb(sys.exc_info()[2]))
          message = child_process_message(message_type="exception",
                                          data=(e, traceback_str))
      finally :
        if message is not None :
          self._stdout._flush()
          self._c.send(message)
      sys.stdout = old_stdout

  #XXX: target functions must use this call signature!!!
  # I normally just write a very short wrapper function that invokes the
  # command I actually care about with args/kwds.
  def _process_target_function (args, kwds, connection) :
    connection.send(True)

  # not really a process, but a wrapper for one.
  class process_with_callbacks (threading.Thread) :
    def __init__ (self,
        target,
        args=(),
        kwargs={},
        callback_stdout=null_callback,
        callback_final=null_callback,
        callback_err=null_callback,
        callback_abort=null_callback,
        callback_other=null_callback,
        buffer_stdout=True) :
      threading.Thread.__init__(self)
      self._target = target
      self._args = args
      self._kwargs = dict(kwargs)
      self._cb_stdout = callback_stdout
      self._cb_final  = callback_final
      self._cb_err    = callback_err
      self._cb_abort  = callback_abort
      self._cb_other  = callback_other
      self._buffer_stdout = buffer_stdout
      self._abort = False
      self._completed = False
      self._error = False

    def abort (self) :
      self._abort = True

    def run (self) :
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
        child_process.start()
      except Exception, e :
        sys.__stderr__.write("Error starting child process: %s\n" % str(e))
        child_process = None
      while child_process is not None :
        if self._abort :
          child_process.terminate()
          self._cb_abort()
          break
        if not parent_conn.poll(1) :
          continue
        pipe_output = parent_conn.recv()
        if isinstance(pipe_output, child_process_message) :
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
      if child_process is not None and child_process.is_alive() :
        child_process.join()

########################################################################
# tests
#
def exercise_threading() :

  def first_callback(i):
    collected.append(-i*10)

  def callback(i):
    collected.append(i*10)
    return (i < 5)

  def run(n=3, callback=None):
    i = 1
    while (callback is not None or i <= n):
      collected.append(i)
      if (callback is not None):
        status = callback(i)
        if (status == False):
          break
      i += 1

  collected = []
  run()
  assert collected == [1,2,3]

  for callback1,expected in [
        (None, [1,10]),
        (first_callback, [1,-10])]:
    for n_resume in xrange(7):
      collected = []
      t = thread_with_callback_and_wait(
        run=run,
        callback=callback,
        first_callback=callback1)
      if (callback1 is None):
        t.start()
      else:
        t.start_and_wait_for_first_callback()
      for i in xrange(n_resume):
        t.resume()
      t.resume(last_iteration=True).join()
      assert collected == expected
      if (n_resume < 4):
        expected.extend([n_resume+2, (n_resume+2)*10])


class _callback_handler (object) :
  def __init__ (self) :
    self._err = None
    self._result = None
    self._abort = None
    self._stdout = None
    self._other = None

  def cb_err (self, error, traceback_info) :
    self._err = error
    self._tb = traceback_info

  def cb_abort (self) :
    self._abort = True

  def cb_result (self, result) :
    self._result = result

  def cb_stdout (self, stdout) :
    if self._stdout is not None : self._stdout += str(stdout)
    else                        : self._stdout = str(stdout)

  def cb_other (self, data) :
    if self._other is not None :
      self._other.append(data)
    else :
      self._other = [data]

  def tst_error_raised (self) :
    assert self._err is not None
    assert self._result is None

  def tst_run_success (self, expected_result=None, expected_stdout=None,
      expected_other=None) :
    assert self._err is None
    assert self._abort is None
    assert self._result == expected_result
    assert self._stdout == expected_stdout
    assert self._other == expected_other

#--- test 01 : propagate args and kwds to target
def _inner_target_with_args_and_kwds (a, b, c=-1, d=-1) :
  assert a < 0 and b < 0
  assert c > 0 and d > 0
  return True

def _target_function01 (args, kwds, connection) :
  return _inner_target_with_args_and_kwds(*args, **kwds)

def tst_01 () :
  ch = _callback_handler()
  p = process_with_callbacks(
    target = _target_function01,
    args   = (-2, -2),
    kwargs = {"c": 2, "d" : 2},
    callback_stdout = ch.cb_stdout,
    callback_final  = ch.cb_result,
    callback_err    = ch.cb_err,
    callback_other  = ch.cb_other)
  p.start()
  while p.isAlive() :
    pass
  ch.tst_run_success(expected_result=True, expected_stdout=None,
    expected_other=None)

#--- test 02 : target function does not return a value
def _inner_target_no_return () :
  n = 1

def _target_function02 (args, kwds, connection) :
  return _inner_target_no_return(*args, **kwds)

def tst_02 () :
  ch = _callback_handler()
  p = process_with_callbacks(
    target = _target_function02,
    callback_stdout = ch.cb_stdout,
    callback_final  = ch.cb_result,
    callback_err    = ch.cb_err,
    callback_other  = ch.cb_other)
  p.start()
  while p.isAlive() :
    pass
  ch.tst_run_success(expected_result=None)

#--- test 03 : callbacks for stdout, runtime, result
def _target_function03 (args, kwds, connection) :
  for i in xrange(4) :
    print i
    connection.send(i)
  return 4

def tst_03 () :
  ch = _callback_handler()
  p = process_with_callbacks(
    target = _target_function03,
    callback_stdout = ch.cb_stdout,
    callback_final  = ch.cb_result,
    callback_err    = ch.cb_err,
    callback_other  = ch.cb_other)
  p.start()
  while p.isAlive() :
    pass
  tstout = \
"""0
1
2
3
"""
  ch.tst_run_success(expected_result=4, expected_stdout=tstout,
    expected_other=[0,1,2,3])

#--- test 04 : stdout consistency
def _tst_print (out=None) :
  if out is None :
    out = sys.stdout
  for i in xrange(1000) :
    out.write("%s\n" % i)
  return None

def _target_function04 (args, kwds, connection) :
  return _tst_print(*args, **kwds)

def tst_04 () :
  import cStringIO
  tmpout = cStringIO.StringIO()
  _tst_print(tmpout)

  for buffer_stdout in [True, False] :
    ch = _callback_handler()
    p = process_with_callbacks(
      target = _target_function04,
      callback_stdout = ch.cb_stdout,
      buffer_stdout = buffer_stdout)
    p.start()
    while p.isAlive() :
      pass
    assert tmpout.getvalue() == ch._stdout

#--- test 05 : propagating standard Exceptions and Sorry
def _target_function05a (args, kwds, connection) :
  raise Exception("_target_function05a")

def _target_function05b (args, kwds, connection) :
  raise Sorry("_target_function05b")

def tst_05 () :
  for f in [_target_function05a, _target_function05b] :
    ch = _callback_handler()
    p = process_with_callbacks(
      target = f,
      callback_err = ch.cb_err,
      callback_final = ch.cb_result)
    p.start()
    while p.isAlive() :
      pass
    ch.tst_error_raised()
    if f is _target_function05b :
      assert isinstance(ch._err, Sorry)

#--- test 06 : aborting
# Process.terminate() does not work on RedHat 8, but this test will still
# be successful despite taking an extra 100 seconds to finish.
def _target_function06 (args, kwds, connection) :
  for i in xrange(10) :
    print i
    time.sleep(1)

def tst_06 () :
  ch = _callback_handler()
  p = process_with_callbacks(
    target = _target_function06,
    callback_abort = ch.cb_abort)
  p.start()
  p.abort()
  while p.isAlive() :
    pass
  assert ch._abort == True

def _target_function07 (args, kwds, connection) :
  time.sleep(2)
  raise Abort()

def tst_07 () :
  ch = _callback_handler()
  p = process_with_callbacks(
    target=_target_function07,
    callback_abort=ch.cb_abort)
  p.start()
  p.abort()
  while p.isAlive() :
    pass
  assert ch._abort == True

def exercise_process () :
  #--- run all tests
  try :
    tst_01()
    tst_02()
    tst_03()
    tst_04()
    tst_05()
    tst_06()
    tst_07()
  except ImportError, e:
    print "Skipping thread_utils tests:", str(e)

if (__name__ == "__main__"):
  from libtbx import thread_utils
  thread_utils.exercise_threading()
  thread_utils.exercise_process()
  print "OK"
