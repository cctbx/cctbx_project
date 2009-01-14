class cluster_manager(object):

  __slots__ = ["cluster_indices", "clusters", "parents"]

  def __init__(O, n_vertices):
    O.cluster_indices = range(n_vertices)
    O.clusters = []
    for i in xrange(n_vertices):
      O.clusters.append([i])
    O.parents = None

  def connect(O, i, j):
    assert O.parents is None
    ci = O.cluster_indices
    cii = ci[i]
    cij = ci[j]
    if (cii == cij): return
    ccij = O.clusters[cij]
    ccii = O.clusters[cii]
    if (len(ccij) <= len(ccii)):
      for k in ccij: ci[k] = cii
      ccii.extend(ccij)
      del ccij[:]
    else:
      for k in ccii: ci[k] = cij
      ccij.extend(ccii)
      del ccii[:]

  def refresh_indices(O):
    assert O.parents is None
    ci = O.cluster_indices
    for ic,c in enumerate(O.clusters):
      for i in c:
        ci[i] = ic

  def tidy(O):
    assert O.parents is None
    for c in O.clusters: c.sort()
    def cmp_clusters(a, b):
      if (len(a) > len(b)): return -1
      if (len(a) < len(b)): return 1
      if (len(a) != 0): return cmp(a[0], b[0])
      return 0
    O.clusters.sort(cmp_clusters)
    for ic in xrange(len(O.clusters)-1,-1,-1):
      if (len(O.clusters[ic]) != 0):
        del O.clusters[ic+1:]
        break
    O.refresh_indices()

  def cluster_edge_sets(O, edges):
    result = []
    for i in xrange(len(O.clusters)):
      result.append(set())
    ci = O.cluster_indices
    for i,j in edges:
      cii = ci[i]
      cij = ci[j]
      if (cii == cij): continue
      result[cii].add(cij)
      result[cij].add(cii)
    return result

  def merge_lones(O, edges):
    assert O.parents is None
    for cii,es in enumerate(O.cluster_edge_sets(edges=edges)):
      if (len(es) != 1): continue
      ccii = O.clusters[cii]
      if (len(ccii) != 1): continue
      cij = iter(es).next()
      O.clusters[cij].extend(ccii)
      del ccii[:]
    O.tidy()

  def construct_spanning_trees(O, edges):
    assert O.parents is None
    cluster_edge_sets = O.cluster_edge_sets(edges=edges)
    n_clusters = len(O.clusters)
    if (n_clusters == 1):
      assert len(cluster_edge_sets[0]) == 0
      O.parents = [-1]
      return
    w_max = len(O.clusters[0])
    candi = []
    for i in xrange(w_max+1):
      candi.append([])
    parents = [-1] * n_clusters
    i_new_given_i_old = [-1] * n_clusters
    i_old_given_i_new = []
    for ip in xrange(n_clusters):
      if (parents[ip] != -1): continue
      i_new_given_i_old[ip] = len(i_old_given_i_new)
      i_old_given_i_new.append(ip)
      w_max = 0
      for j in cluster_edge_sets[ip]:
        assert j != ip
        w = len(O.clusters[j])
        candi[w].append(j)
        parents[j] = ip
        if (w_max < w): w_max = w
      while True:
        kp = None
        ip = n_clusters
        cw = candi[w_max]
        for k in xrange(len(cw)):
          if (ip > cw[k]):
            kp = k
            ip = cw[k]
        if (kp is None):
          break
        del cw[kp]
        for j in cluster_edge_sets[ip]:
          assert j != ip
          if (i_new_given_i_old[j] == -1):
            w = len(O.clusters[j])
            candi[w].append(j)
            parents[j] = ip
            if (w_max < w): w_max = w
        if (i_new_given_i_old[ip] == -1):
          i_new_given_i_old[ip] = len(i_old_given_i_new)
          i_old_given_i_new.append(ip)
        for w_max in xrange(w_max,-1,-1):
          if (len(candi[w_max]) != 0):
            break
        else:
          break
    assert len(i_old_given_i_new) == len(i_new_given_i_old)
    new_clusters = []
    for i in i_old_given_i_new:
      new_clusters.append(O.clusters[i])
    del O.clusters[:]
    O.clusters.extend(new_clusters)
    O.refresh_indices()
    O.parents = [None] * len(parents)
    for i,ip in enumerate(parents):
      if (ip != -1): ip = i_new_given_i_old[ip]
      O.parents[i_new_given_i_old[i]] = ip

  def find_loop_edges(O, edges):
    assert O.parents is not None
    result = []
    ci = O.cluster_indices
    p = O.parents
    for e in edges:
      cii, cij = [ci[i] for i in e]
      if (      cii  !=   cij
          and p[cii] !=   cij
          and   cii  != p[cij]):
        result.append(e)
    return result

def construct_edge_sets(n_vertices, edges):
  result = [set() for i in xrange(n_vertices)]
  for i,j in edges:
    assert i < j
    result[i].add(j)
    result[j].add(i)
  return result

def find_loops(edge_sets, depth, loop_set, path, iv, traversing):
  path = path + [iv]
  traversing[iv] = True
  at_limit = (len(path) == depth)
  for jv in edge_sets[iv]:
    if (jv < path[0]): continue
    if (jv == path[0] and len(path) > 2):
      loop_set.update(path)
    if (at_limit): continue
    if (traversing[jv]): continue
    find_loops(edge_sets, depth, loop_set, path, jv, traversing)
  traversing[iv] = False

class construct(object):

  def __init__(O, n_vertices, edges, size_max=8):
    O.n_vertices = n_vertices
    O.edges = edges
    O.edge_sets = construct_edge_sets(n_vertices=n_vertices, edges=edges)
    O.cluster_manager = cluster_manager(n_vertices=n_vertices)
    traversing = [False] * n_vertices
    for iv in xrange(n_vertices):
      loop_set = set()
      find_loops(
        edge_sets=O.edge_sets,
        depth=size_max,
        loop_set=loop_set,
        path=[],
        iv=iv,
        traversing=traversing)
      for jv in loop_set:
        O.cluster_manager.connect(i=iv, j=jv)
    O.cluster_manager.tidy()
