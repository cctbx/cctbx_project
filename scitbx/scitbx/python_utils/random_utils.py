import random
import os
import time
from scitbx.array_family import flex

def seed(a=None):
  if (a == None):
    try: a = long(os.getpid() * (2**16)) + long(time.time() * 2**8)
    except: pass
  random.seed(a)

def random_subset(set, n):
  assert n >= 0 and n <= len(set)
  if (n == 0): return []
  set = list(set) # deep copy
  while (len(set) > n):
    del set[random.randrange(len(set))]
  return set

def random_selection(size, n):
  assert n >= 0 and n <= size
  result = flex.bool(size, 00000)
  if (n > 0):
    set = range(size)
    while (len(set) > n):
      del set[random.randrange(len(set))]
    for i in set:
      result[i] = 0001
  return result

class weighted_choice:

  def __init__(self, weights):
    self.accumulated_weights = flex.double()
    sum = 0
    for w in weights:
      sum += w
      self.accumulated_weights.append(sum)
    self.accumulated_weights /= sum

  def next(self):
    r = random.random()
    for i,w in self.accumulated_weights.items():
      if (w >= r): return i
    return self.accumulated_weights.size()-1

if (__name__ == "__main__"):
  seed()
  print random_subset(range(5), 0)
  print random_subset(range(5), 5)
  print random_subset(range(5), 2)
  print random_subset(range(5), 3)
  print tuple(random_selection(5, 0))
  print tuple(random_selection(5, 5))
  print tuple(random_selection(5, 2))
  print tuple(random_selection(5, 3))
  for weights in ([5,5], [4,3,2,1]):
    r = weighted_choice(weights)
    hist = [0 for i in xrange(len(weights))]
    for i in xrange(10000):
      hist[r.next()] += 1
    hist = [int(round(s/1000.)) for s in hist]
    assert hist == weights, (hist, weights)
  print "OK"
