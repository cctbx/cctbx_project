class left_decomposition:

  def __init__(self, g, h):
    g = [s for s in g]
    h = [s for s in h]
    assert len(g) % len(h) == 0
    assert h[0].is_unit_mx()
    self.partition_indices = [-1] * len(g)
    self.partitions = []
    for i,gi in enumerate(g):
      if (self.partition_indices[i] != -1): continue
      self.partition_indices[i] = len(self.partitions)
      partition = [gi]
      for hj in h[1:]:
        gihj = gi.multiply(hj)
        for k in xrange(i+1,len(g)):
          if (self.partition_indices[k] != -1): continue
          gk = g[k]
          if (gk.r().num() == gihj.r().num()):
            self.partition_indices[k] = len(self.partitions)
            partition.append(gk)
            break
        else:
          raise RuntimeError("h is not a subgroup of g")
      if (len(partition) != len(h)):
        raise RuntimeError("h is not a subgroup of g")
      self.partitions.append(partition)
    if (len(self.partitions) * len(h) != len(g)):
      raise RuntimeError("h is not a subgroup of g")

def double_unique(g, h1, h2):
  g = [s for s in g]
  h1 = [s for s in h1]
  h2 = [s for s in h2]
  result = []
  done = {}
  for a in g:
    if (str(a) in done): continue
    result.append(a)
    for hi in h1:
      for hj in h2:
        b = hi.multiply(a).multiply(hj)
        done[str(b)] = None
  return result
