
from __future__ import absolute_import, division, print_function
import inspect
from six.moves import range

class error(Exception):
  """ libtbx.python_frame error """


def named(name):
  """ Returns an object representing the frame of the function of the given name.
      This frame will be searched upward from the current frame

      Synopsis:

      from libtbx import pyframe

      def f(i):
        for j in range(i):
          g(j)

      def g(j):
        return k(j)

      def k(j):
        f = pyframe.named('f')
        if f.j == 3: ....

      This is particularly useful for conditional breakpoints in the WingIDE debugger.
      Let's say we have a function 'test' which tests a function 'clever_algorithm' and that
      one wants to put a breakpoint in the latter that only triggers
      for one particular test (that we know is going to fail e.g.). Let's say tests are
      identified by a local variable 'idx' in function 'test'. Then we can put a conditional
      breakpoint somewhere in function 'clever_algorithm' like so:

        pyframe.named('test').idx == 5

      The advantage of conditional tests is that they don't require modifying code for
      debugging purposes, which is always a dangerous thing to do as one may inadvertantly
      leave debugging code in.

  """
  f = inspect.currentframe()
  while f.f_code.co_name != name:
    if f.f_back is None:
      raise error("frame '%s' cannot be found among the frames calling the current frame.")
    f = f.f_back
  return frame_locals(f)


def up(n):
  """ Returns an object representing the frame of the function that is n steps up in
      the calling stack, starting from the function that called this function 'up'.
      Thus n=0 corresponds to the function F that called this, n=1 to the function
      that called F, etc.
  """
  f = inspect.currentframe()
  f = f.f_back
  caller_name = f.f_code.co_name
  for i in range(n):
    if f is None:
      raise error("asked for the frame %i steps up the function '%s'"
                  " but the calling stack is not that tall.")
    f = f.f_back
  return frame_locals(f)


class frame_locals(object):

  def __init__(self, f):
    self.__dict__.update(f.f_locals)
