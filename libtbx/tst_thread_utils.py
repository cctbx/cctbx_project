
from __future__ import absolute_import, division, print_function

from six.moves import cStringIO as StringIO

from libtbx.thread_utils import thread_with_callback_and_wait
from libtbx.thread_utils import process_with_callbacks
from libtbx.utils import Sorry, Abort
import time
import sys
from six.moves import range

########################################################################
# THREADING ONLY
#
def exercise_threading():
  collected = []

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

  run()
  assert collected == [1,2,3]

  for callback1,expected in [
        (None, [1,10]),
        (first_callback, [1,-10])]:
    for n_resume in range(7):
      collected = []
      t = thread_with_callback_and_wait(
        run=run,
        callback=callback,
        first_callback=callback1)
      if (callback1 is None):
        t.start()
      else:
        t.start_and_wait_for_first_callback()
      for i in range(n_resume):
        t.resume()
      t.resume(last_iteration=True).join()
      assert collected == expected
      if (n_resume < 4):
        expected.extend([n_resume+2, (n_resume+2)*10])

########################################################################
# THREADING + MULTIPROCESSING
#
class _callback_handler(object):
  def __init__(self):
    self._err = None
    self._result = None
    self._abort = None
    self._stdout = None
    self._other = None

  def cb_err(self, error, traceback_info):
    self._err = error
    self._tb = traceback_info

  def cb_abort(self):
    self._abort = True

  def cb_result(self, result):
    self._result = result

  def cb_stdout(self, stdout):
    if self._stdout is not None : self._stdout += str(stdout)
    else                        : self._stdout = str(stdout)

  def cb_other(self, data):
    if self._other is not None :
      self._other.append(data)
    else :
      self._other = [data]

  def tst_error_raised(self):
    assert self._err is not None
    assert self._result is None

  def tst_run_success(self, expected_result=None, expected_stdout=None,
      expected_other=None):
    assert self._err is None
    assert self._abort is None
    assert self._result == expected_result
    assert self._stdout == expected_stdout
    assert self._other == expected_other

#--- test 01 : propagate args and kwds to target
def _inner_target_with_args_and_kwds(a, b, c=-1, d=-1):
  assert a < 0 and b < 0
  assert c > 0 and d > 0
  return True

def _target_function01(args, kwds, connection):
  return _inner_target_with_args_and_kwds(*args, **kwds)

def tst_args_kwds():
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
  while p.is_alive():
    pass
  ch.tst_run_success(expected_result=True, expected_stdout=None,
    expected_other=None)

#--- test 02 : target function does not return a value
def _inner_target_no_return():
  n = 1

def _target_function02(args, kwds, connection):
  return _inner_target_no_return(*args, **kwds)

def tst_no_return_value():
  ch = _callback_handler()
  p = process_with_callbacks(
    target = _target_function02,
    callback_stdout = ch.cb_stdout,
    callback_final  = ch.cb_result,
    callback_err    = ch.cb_err,
    callback_other  = ch.cb_other)
  p.start()
  while p.is_alive():
    pass
  ch.tst_run_success(expected_result=None)

#--- test 03 : callbacks for stdout, runtime, result
def _target_function03(args, kwds, connection):
  for i in range(4):
    print(i)
    connection.send(i)
  return 4

def tst_callbacks():
  ch = _callback_handler()
  p = process_with_callbacks(
    target = _target_function03,
    callback_stdout = ch.cb_stdout,
    callback_final  = ch.cb_result,
    callback_err    = ch.cb_err,
    callback_other  = ch.cb_other)
  p.start()
  while p.is_alive():
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
def _tst_print(out=None):
  if out is None :
    out = sys.stdout
  for i in range(1000):
    out.write("%s\n" % i)
  return None

def _target_function04(args, kwds, connection):
  return _tst_print(*args, **kwds)

def tst_stdout():
  tmpout = StringIO()
  _tst_print(tmpout)

  for buffer_stdout in [True, False] :
    ch = _callback_handler()
    p = process_with_callbacks(
      target = _target_function04,
      callback_stdout = ch.cb_stdout,
      buffer_stdout = buffer_stdout)
    p.start()
    while p.is_alive():
      pass
    assert tmpout.getvalue() == ch._stdout

#--- test 05 : propagating standard Exceptions and Sorry
def _target_function05a(args, kwds, connection):
  raise Exception("_target_function05a")

def _target_function05b(args, kwds, connection):
  raise Sorry("_target_function05b")

def tst_exceptions():
  for f in [_target_function05a, _target_function05b] :
    ch = _callback_handler()
    p = process_with_callbacks(
      target = f,
      callback_err = ch.cb_err,
      callback_final = ch.cb_result)
    p.start()
    while p.is_alive():
      pass
    ch.tst_error_raised()
    if f is _target_function05b :
      assert isinstance(ch._err, Sorry)

#--- test 06 : aborting
# Process.terminate() does not work on RedHat 8, but this test will still
# be successful despite taking an extra 100 seconds to finish.
def _target_function06(args, kwds, connection):
  for i in range(10):
    print(i)
    time.sleep(1)

def tst_abort_simple():
  ch = _callback_handler()
  p = process_with_callbacks(
    target = _target_function06,
    callback_abort = ch.cb_abort)
  p.start()
  p.abort()
  while p.is_alive():
    pass
  assert ch._abort == True

def _target_function07(args, kwds, connection):
  time.sleep(2)
  raise Abort()

def tst_abort_2():
  ch = _callback_handler()
  p = process_with_callbacks(
    target=_target_function07,
    callback_abort=ch.cb_abort)
  p.start()
  p.abort()
  while p.is_alive():
    pass
  assert ch._abort == True

def _target_function08(args, kwds, connection):
  import libtbx.callbacks # import dependency

  log = StringIO()
  libtbx.call_back.set_warning_log(log)
  time.sleep(1)
  libtbx.warn("Hello, world!")
  time.sleep(1)

def tst_warn_callback():
  ch = _callback_handler()
  p = process_with_callbacks(
    target=_target_function08,
    callback_other=ch.cb_other)
  p.start()
  while p.is_alive():
    pass
  assert (ch._other[0].message == "warn")
  assert (ch._other[0].data == "Hello, world!")

#--- pause, resume, and kill
def _target_function09(args, kwds, connection):
  for i in range(10):
    print(i)
    time.sleep(1)

class _callback_handler_2(object):
  def __init__(self):
    self.lines = ""
    self.paused = False
    self.resumed = False
    self.aborted = False

  def cb_print(self, data):
    self.lines += data

  def cb_pause(self):
    self.paused = True

  def cb_resume(self):
    self.resumed = True
    self.paused = False

  def cb_abort(self):
    self.aborted = True

def tst_pause_resume_kill():
  if (sys.platform == "win32"):
    print("Skipping pause/resume test (not available on Windows)")
  else :
    ch = _callback_handler_2()
    p = process_with_callbacks(
      target=_target_function09,
      callback_stdout=ch.cb_print,
      callback_pause=ch.cb_pause,
      callback_resume=ch.cb_resume,
      callback_abort=ch.cb_abort)
    p.start()
    time.sleep(1)
    p.pause()
    assert ch.paused
    current_data = str(ch.lines)
    time.sleep(5)
    assert (ch.lines == current_data)
    p.resume()
    assert ch.resumed
    while p.is_alive():
      pass
    assert not ch.aborted
    assert (ch.lines.splitlines() == ['0','1','2','3','4','5','6','7','8','9'])
  #return
  # kill test
  ch = _callback_handler_2()
  p = process_with_callbacks(
    target=_target_function09,
    callback_stdout=ch.cb_print,
    callback_pause=ch.cb_pause,
    callback_resume=ch.cb_resume,
    callback_abort=ch.cb_abort)
  p.start()
  time.sleep(3)
  p.kill()
  assert ch.aborted

def exercise_process():
  #--- run all tests
  try :
    tst_args_kwds()
    tst_no_return_value()
    tst_callbacks()
    tst_stdout()
    tst_exceptions()
    tst_abort_simple()
    tst_abort_2()
    tst_warn_callback()
    tst_pause_resume_kill()
  except ImportError as e:
    print("Skipping thread_utils tests:", str(e))
  else :
    print("OK")

if (__name__ == "__main__"):
  exercise_threading()
  exercise_process()
