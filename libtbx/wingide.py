""" Tools for the Python integrated development environment Wing IDE """
from __future__ import absolute_import, division, print_function

def log(stop=False, **kwds):
  """ Helper to ease putting watchpoints

  Synopsis:
    - set a conditional breakpoint
    - write the condition as
        some_condition and log(some_variable=some_variable, ...)
    - debug the code

  The debugger won't stop at that breakpoint but it will log name and values
  of the specified variables each time it passes over that breakpoint if
  and only if the given condition is true.

  Remark:
    Passing stop=True makes it a real breakpoint
  """
  for k,v in kwds.items():
    print("%s = %s" % (k,v))
  return stop
