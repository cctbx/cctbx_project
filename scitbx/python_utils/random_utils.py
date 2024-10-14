from __future__ import absolute_import, division, print_function
import random
import os
import time
from scitbx.array_family import flex
from six.moves import range

def seed(a=None):
  if (a is None):
    try: a = long(os.getpid() * (2**16)) + long(time.time() * 2**8)
    except Exception: pass
  random.seed(a)

def random_subset(set, n):
  assert n >= 0 and n <= len(set)
  if (n == 0): return []
  set = list(set) # deep copy
  while (len(set) > n):
    del set[random.randrange(len(set))]
  return set

class weighted_choice(object):

  def __init__(self, weights):
    self.accumulated_weights = flex.double()
    sum = 0
    for w in weights:
      sum += w
      self.accumulated_weights.append(sum)
    self.accumulated_weights /= sum

  def next(self):
    r = random.random()
    for i,w in enumerate(self.accumulated_weights):
      if (w >= r): return i
    return self.accumulated_weights.size()-1

if (__name__ == "__main__"):
  seed()
  print(random_subset(range(5), 0))
  print(random_subset(range(5), 5))
  print(random_subset(range(5), 2))
  print(random_subset(range(5), 3))
  for weights in ([5,5], [4,3,2,1]):
    r = weighted_choice(weights)
    hist = [0 for i in range(len(weights))]
    for i in range(10000):
      hist[next(r)] += 1
    hist = [int(round(s/1000.)) for s in hist]
    assert hist == weights, (hist, weights)
  print("OK")
