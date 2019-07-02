from __future__ import absolute_import, division, print_function
import signal
import sys

class keyboard_interrupt_handler(object):

  def __init__(self, max_n_events=3):
    self.prev_handler = signal.signal(signal.SIGINT, self.__call__)
    self.max_n_events = max_n_events
    self.n_events = 0

  def __del__(self):
    # self.disable() leads to SystemError
    if (self.prev_handler is not None):
      raise RuntimeError(
        "Internal error: interrupt not disabled before it goes out of scope.")

  def disable(self):
    if (self.prev_handler is not None):
      signal.signal(signal.SIGINT, self.prev_handler)
      self.prev_handler = None

  def __call__(self, signum, frame):
    self.n_events += 1
    sys.stdout.flush()
    sys.stderr.flush()
    print()
    print('Keyboard interrupt #%d.' % self.n_events)
    print()
    if (self.n_events == self.max_n_events):
      self.disable()
      raise KeyboardInterrupt
    print('Continuing to next checkpoint. Please wait.')
    if (self.max_n_events > 0):
      n_left = self.max_n_events-self.n_events
      print('To stop immediately send %d more keyboard interrupt%s.' % (
        n_left, ["", "s"][max(0,min(n_left-1, 1))]))
    print()
    sys.stdout.flush()
