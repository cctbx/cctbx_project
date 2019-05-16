from __future__ import absolute_import, division, print_function
from libtbx import easy_mp
import os
import random
import string
import threading
import time
from six.moves import range

class state_object():
  ''' A simple thread-safe object keeping an internal state. '''
  def __init__(self):
    self._lock = threading.Lock()
    self._state = "-init-"

  def generate_state(self):
    with self._lock:
      oldstate = self._state
      self._state = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in list(range(6)))
      print("(%s-%5s) changing state from %s to %s" % (format(id(self), "#x"), os.getpid(), oldstate, self._state))
      return self._state

  def get_state(self):
    with self._lock:
      print("(%s-%5s) object in state %s" % (format(id(self), "#x"), os.getpid(), self._state))
      return self._state

def exercise_multiprocessing(mp_nproc=1, mp_threads=1, mp_method="multiprocessing", tasks=3):
  print("Running %s test with %d processes, %d threads, %d tasks" % \
    (mp_method, mp_nproc, mp_threads, tasks))

  # Create one shared instance of the state object and extract its initial state
  master_state_object = state_object()
  initial_state = master_state_object.get_state()

  # This is a function that changes the state on the object
  def change_stored_state(task):
    time.sleep(random.random() / 2 / tasks)
    master_state_object.generate_state()

  # Call the state-changing function in parallel
  easy_mp.parallel_map(
    iterable=range(tasks),
    func=change_stored_state,
    processes=mp_nproc,
    method=mp_method)

  # Get the final state of the object
  final_state = master_state_object.get_state()

  # Did it change?
  assert initial_state != final_state

def run():
  exercise_multiprocessing(mp_nproc=1, mp_method="threading")
  print("OK")

  exercise_multiprocessing(mp_nproc=2, mp_method="threading")
  print("OK")

  exercise_multiprocessing(mp_nproc=1, mp_method="multiprocessing")
  print("OK")

  #exercise_multiprocessing(mp_nproc=2, mp_method="multiprocessing")
  #
  # This will fail. Multiprocessing creates copy-on-write copies of the
  # initial object. Later changes in the state are lost.
  # Problems can occur when an instance of state_object does not check
  # that it has been copied, and when it carries an internal reference to
  # something that cannot be copied properly, possibly a C object, eg.
  # a BZ2File()-instance or, more probably, a file handle.
  # State_objects can not directly detect that it has been copied, so
  # as a fallback it could check whether os.getpid() has changed.

if __name__ == "__main__":
  run()
