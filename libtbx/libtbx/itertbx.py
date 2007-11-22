from __future__ import generators

try:
  if 0: raise
  from itertools import count, islice, izip
except:
  def count(n=0):
    while True:
      yield n
      n += 1

  def islice(iterable, *args):
    s = slice(*args)
    it = iter(xrange(s.start or 0, s.stop or sys.maxint, s.step or 1))
    nexti = it.next()
    for i, element in enumerate(iterable):
      if i == nexti:
        yield element
        nexti = it.next()

  def izip(*iterables):
    iterables = map(iter, iterables)
    while iterables:
      result = [it.next() for it in iterables]
      yield tuple(result)

def step(firstval=0, increment=1):
  val = firstval
  while True:
    yield val
    val += increment
