import random

def random_subset(set, n):
  assert n >= 0 and n <= len(set)
  if (n == 0): return []
  set = list(set) # deep copy
  while (len(set) > n):
    del set[random.randrange(len(set))]
  return set

if (__name__ == "__main__"):
  print random_subset(range(5), 0)
  print random_subset(range(5), 5)
  print random_subset(range(5), 2)
  print random_subset(range(5), 3)
