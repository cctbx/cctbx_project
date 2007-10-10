import math

def iround(x):
  if (x < 0): return int(x-0.5)
  return int(x+.5)

def iceil(x):
  return iround(math.ceil(x))

def ifloor(x):
  return iround(math.floor(x))

def does_imply(p,q):
  """ does p => q in the sense of logical implication? """
  return not p or q

def are_equivalent(p,q):
  """ does p <=> q in the sense of logical equivalence? """
  return does_imply(p,q) and does_imply(q,p)

if (__name__ == "__main__"):
  assert iround(0) == 0
  assert iround(1.4) == 1
  assert iround(-1.4) == -1
  assert iround(1.6) == 2
  assert iround(-1.6) == -2
  assert iceil(0) == 0
  assert iceil(1.1) == 2
  assert iceil(-1.1) == -1
  assert iceil(1.9) == 2
  assert iceil(-1.9) == -1
  assert ifloor(0) == 0
  assert ifloor(1.1) == 1
  assert ifloor(-1.1) == -2
  assert ifloor(1.9) == 1
  assert ifloor(-1.9) == -2

  assert does_imply(True, True)
  assert not does_imply(True, False)
  assert does_imply(False, True)
  assert does_imply(False, False)

  assert are_equivalent(True, True)
  assert not are_equivalent(True, False)
  assert not are_equivalent(False, True)
  assert are_equivalent(False, False)

  print "OK"
