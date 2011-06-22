import random
import os
import time
from scitbx.array_family import flex

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

def random_selection(n_candidates, n_keep):
  assert n_keep >= 0 and n_keep <= n_candidates
  n_discard = n_candidates - n_keep
  if (n_keep > n_discard):
    selection = flex.bool(n_candidates, True)
    if (n_discard > 0):
      _random_selection_core(selection, n_keep, False)
  else:
    selection = flex.bool(n_candidates, False)
    if (n_keep > 0):
      _random_selection_core(selection, n_discard, True)
  return selection

def _random_selection_core(selection, target_set_size, flag):
   set = range(len(selection))
   while (len(set) > target_set_size):
     i = random.randrange(len(set))
     selection[set[i]] = flag
     del set[i]

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
  print random_subset(range(5), 0)
  print random_subset(range(5), 5)
  print random_subset(range(5), 2)
  print random_subset(range(5), 3)
  for i in xrange(10):
    assert random_selection(0, 0).size() == 0
    assert random_selection(5, 0).size() == 5
    assert random_selection(3, 0).count(True) == 0
    assert random_selection(3, 3).size() == 3
    assert random_selection(5, 5).count(True) == 5
    assert random_selection(6, 2).size() == 6
    assert random_selection(6, 2).count(True) == 2
    assert random_selection(4, 3).size() == 4
    assert random_selection(4, 3).count(True) == 3
  for weights in ([5,5], [4,3,2,1]):
    r = weighted_choice(weights)
    hist = [0 for i in xrange(len(weights))]
    for i in xrange(10000):
      hist[r.next()] += 1
    hist = [int(round(s/1000.)) for s in hist]
    assert hist == weights, (hist, weights)
  print "OK"
