from __future__ import generators

import os
if not hasattr(os.path, "devnull"):
  if os.name == "nt":
    os.path.devnull = "nul"
  else:
    os.path.devnull = "/dev/null"

__builtins__.setdefault("False", 0)
__builtins__.setdefault("True", 1)

if (__builtins__.get("bool", None) is None):
  def bool(value):
    if (value): return True
    return False
  __builtins__["bool"] = bool

if (__builtins__.get("enumerate", None) is None):
  def enumerate(iterable):
    """enumerate(iterable) -> iterator for index, value of iterable

  Return an enumerate object.  iterable must be an other object that supports
  iteration.  The enumerate object yields pairs containing a count (from
  zero) and a value yielded by the iterable argument.  enumerate is useful
  for obtaining an indexed list: (0, seq[0]), (1, seq[1]), (2, seq[2]), ..."""
    i = 0
    it = iter(iterable)
    while 1:
      yield (i, it.next())
      i += 1
  __builtins__["enumerate"] = enumerate

if (__builtins__.get("sum", None) is None):
  def sum(sequence, start=0):
    """sum(sequence, start=0) -> value

Returns the sum of a sequence of numbers (NOT strings) plus the value
of parameter 'start'.  When the sequence is empty, returns start."""
    if (len(sequence) == 0): return start
    sequence = iter(sequence)
    if (start is None):
      result = sequence.next()
    else:
      result = start
    for item in sequence:
      result += item
    return result
  __builtins__["sum"] = sum

if (__builtins__.get("sorted", None) is None):
  def sorted(iterable, cmp=None, key=None, reverse=False):
    """\
sorted(iterable, cmp=None, key=None, reverse=False) --> new sorted list"""
    assert key is None, "Not implemented."
    result = list(iterable)
    if (cmp is None): result.sort()
    else:             result.sort(cmp)
    if (reverse): result.reverse()
    return result
  __builtins__["sorted"] = sorted
