import random
from scitbx.array_family import flex

def random_subset(set, n):
  assert n >= 0 and n <= len(set)
  if (n == 0): return []
  set = list(set) # deep copy
  while (len(set) > n):
    del set[random.randrange(len(set))]
  return set

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
  print random_subset(range(5), 0)
  print random_subset(range(5), 5)
  print random_subset(range(5), 2)
  print random_subset(range(5), 3)
  for weights in ([5,5], [4,3,2,1]):
    r = weighted_choice(weights)
    hist = [0 for i in xrange(len(weights))]
    for i in xrange(10000):
      hist[r.next()] += 1
    hist = [int(round(s/1000.)) for s in hist]
    assert hist == weights, (hist, weights)
  print "OK"
