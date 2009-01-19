"Obsolete, do not use." # XXX backward compatibility 2009-01-18

from itertools import count, islice, izip

def step(firstval=0, increment=1):
  val = firstval
  while True:
    yield val
    val += increment
