
import time
import threading
import Queue

class thread_with_callback_and_wait (threading.Thread) :
  def __init__ (self, run_function, callback_function, run_args=(),
      run_kwds={}) :
    self.f = run_function
    self.args = run_args
    self.kwds = run_kwds
    self.kwds['callback'] = self._callback
    self.cb = callback_function
    self.q = Queue.Queue(0)
    self._aborted = False
    threading.Thread.__init__(self)

  def run (self) :
    self.f(*self.args, **self.kwds)

  def _callback (self, *args, **kwds) :
    cb_status = self.cb(*args, **kwds)
    if cb_status is not None :
      return cb_status
    while True :
      go_ahead = self.q.get()
      if go_ahead :
        return True
      else :
        return False

  def resume (self) :
    self.q.put(True)

  def abort (self) :
    self.q.put(False)

def thread_test() :
  t = thread_with_callback_and_wait(tst_run, tst_cb)
  t.start()
  return t

def tst_run (callback=None) :
  for i in xrange(0, 10) :
    print "run: %d" % i
    if callback is not None :
      status = callback(i)
      if status == False :
        break
  print "done."
  return True

def tst_cb (i) :
  print "  callback: %d" % i
  return True
