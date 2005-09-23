from __future__ import generators

__builtins__.setdefault("False", 0)
__builtins__.setdefault("True", 1)

if (__builtins__.get("bool", None) is None):
  def bool(value):
    if (value): return True
    return False
  __builtins__["bool"] = bool

if (__builtins__.get("_enumerate", None) is None):
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
  def sum(sequence):
    if (len(sequence) == 0): return 0
    sequence = iter(sequence)
    result = sequence.next()
    for item in sequence:
      result += item
    return result
  __builtins__["sum"] = sum
