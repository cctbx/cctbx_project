import Queue
import threading

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
    if (self.is_alive()):
      self._queue.put(last_iteration)
    return self

def exercise() :

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

  print "OK"

if (__name__ == "__main__"):
  exercise()
