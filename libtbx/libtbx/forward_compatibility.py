from __future__ import generators

__builtins__.setdefault("False", 0)
__builtins__.setdefault("True", 1)

def _enumerate(iterable):
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

__builtins__.setdefault("enumerate", _enumerate)
