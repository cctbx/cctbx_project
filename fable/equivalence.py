def array_alignment(members_size, i_mbr_byte_offset_pairs):
  n = members_size
  diff_matrix = [None] * (n*(n-1))
  msg_prefix = "equivalence.array_alignment(): "
  msg_directly_conflicting = msg_prefix + "directly conflicting input"
  for p0,p1 in i_mbr_byte_offset_pairs:
    i0,a0 = p0
    i1,a1 = p1
    if (i0 == i1):
      if (a0 != a1):
        raise RuntimeError(msg_directly_conflicting)
    else:
      if (i0 < i1):
        i = i0 * n + i1
        d = a0 - a1
      else:
        i = i1 * n + i0
        d = a1 - a0
      dd = diff_matrix[i]
      if (dd is None):
        diff_matrix[i] = d
      elif (dd != d):
        raise RuntimeError(msg_directly_conflicting)
  cluster_indices = range(n)
  clusters = []
  for i in xrange(n):
    clusters.append([])
  for i0 in xrange(n-1):
    for i1 in xrange(i0+1,n):
      i = i0 * n + i1
      d = diff_matrix[i]
      if (d is not None):
        ci0 = cluster_indices[i0]
        ci1 = cluster_indices[i1]
        if (ci0 == ci1):
          continue
        if (ci0 > ci1):
          ci0, ci1 = ci1, ci0
          d *= -1
          i0, i1 = i1, i0
        if (ci1 != i1):
          continue
        c0 = clusters[ci0]
        c1 = clusters[ci1]
        if (ci0 != i0):
          for i,o in c0:
            if (i == i0):
              d += o
              break
        c0.append((ci1, d))
        cluster_indices[ci1] = ci0
        for i,o in c1:
          c0.append((i, o+d))
          cluster_indices[i] = ci0
        clusters[ci1] = []
  if (cluster_indices.count(0) != n):
    raise RuntimeError(msg_prefix + "insufficient input")
  assert len(clusters[0]) == n-1
  diffs0 = [None] * n
  diffs0[0] = 0
  for i,o in clusters[0]:
    assert i != 0
    assert diffs0[i] is None
    diffs0[i] = o
  for i0 in xrange(n-1):
    for i1 in xrange(i0+1,n):
      i = i0 * n + i1
      d = diff_matrix[i]
      if (    d is not None
          and diffs0[i1] - diffs0[i0] != d):
        raise RuntimeError(msg_prefix + "indirectly conflicting input")
  return diffs0

class cluster_unions(object):

  __slots__ = ["unions", "indices"]

  def __init__(O):
    O.unions = []
    O.indices = {}

  def add(O, key_cluster):
    curr_index = None
    for key in key_cluster:
      prev_index = O.indices.get(key)
      if (curr_index is None):
        if (prev_index is None):
          O.indices[key] = curr_index = len(O.unions)
          O.unions.append([key])
        else:
          curr_index = prev_index
      elif (prev_index is None):
        O.indices[key] = curr_index
        O.unions[curr_index].append(key)
      elif (prev_index != curr_index):
        if (prev_index < curr_index):
          curr_index, prev_index = prev_index, curr_index
        for key in O.unions[prev_index]:
          O.indices[key] = curr_index
          O.unions[curr_index].append(key)
        O.unions[prev_index] = None

  def tidy(O):
    unions = []
    indices = {}
    for union in O.unions:
      if (union is not None):
        for key in union:
          indices[key] = len(unions)
        unions.append(union)
    O.unions = unions
    O.indices = indices
