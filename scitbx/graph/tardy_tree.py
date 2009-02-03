from scitbx.graph.utils import construct_edge_sets

class cluster_manager(object):

  __slots__ = [
    "cluster_indices", "clusters",
    "hinge_edges", "loop_edges",
    "loop_edge_bendings"]

  def __init__(O, n_vertices):
    O.cluster_indices = range(n_vertices)
    O.clusters = []
    for i in xrange(n_vertices):
      O.clusters.append([i])
    O.hinge_edges = None
    O.loop_edges = None
    O.loop_edge_bendings = None

  def connect(O, i, j):
    assert O.hinge_edges is None
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
    ci = O.cluster_indices
    for ic,c in enumerate(O.clusters):
      for i in c:
        ci[i] = ic

  def tidy(O):
    assert O.hinge_edges is None
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

  def merge_lones(O, edge_sets):
    assert O.hinge_edges is None
    for cii in xrange(len(O.clusters)):
      ccii = O.clusters[cii]
      if (len(ccii) != 1): continue
      i = ccii[0]
      es = edge_sets[i]
      if (len(es) != 1): continue
      j = iter(es).next()
      cij = O.cluster_indices[j]
      O.clusters[cij].extend(ccii)
      del ccii[:]
    O.tidy()

  def merge_clusters_with_multiple_connections(O, edge_sets):
    while True:
      repeat = False
      for cii in xrange(len(O.clusters)):
        while True:
          connected = set()
          multiple = set()
          for i in O.clusters[cii]:
            for j in edge_sets[i]:
              cij = O.cluster_indices[j]
              if (cij == cii): continue
              if (cij in connected): multiple.add(cij)
              else:                  connected.add(cij)
          if (len(multiple) == 0):
            break
          for cij in multiple:
            ccij = O.clusters[cij]
            for j in ccij:
              O.cluster_indices[j] = cii
            O.clusters[cii].extend(ccij)
            del ccij[:]
            repeat = True
      if (not repeat):
        break
    O.tidy()

  def construct_spanning_trees(O, edge_sets):
    assert O.hinge_edges is None
    hinge_edges = [(-1,c[0]) for c in O.clusters]
    O.loop_edges = []
    n_clusters = len(O.clusters)
    if (n_clusters < 2):
      O.hinge_edges = hinge_edges
      return
    w_max = len(O.clusters[0])
    candi = []
    for i in xrange(w_max+1):
      candi.append([])
    done = [0] * n_clusters
    new_clusters = []
    for ip in xrange(n_clusters):
      pe = hinge_edges[ip]
      if (pe[0] != -1): continue
      done[ip] = 1
      new_clusters.append(O.clusters[ip])
      w_max = 0
      for i in O.clusters[ip]:
        for j in edge_sets[i]:
          cij = O.cluster_indices[j]
          if (cij == ip): continue
          if (done[cij] != 0):
            O.loop_edges.append((i,j))
          else:
            done[cij] = -1
            w = len(O.clusters[cij])
            candi[w].append(cij)
            hinge_edges[cij] = (i,j)
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
        for i in O.clusters[ip]:
          for j in edge_sets[i]:
            cij = O.cluster_indices[j]
            if (cij == ip): continue
            if (done[cij] == 1): continue
            if (done[cij] == -1):
              O.loop_edges.append((i,j))
            else:
              done[cij] = -1
              w = len(O.clusters[cij])
              candi[w].append(cij)
              hinge_edges[cij] = (i,j)
              if (w_max < w): w_max = w
        assert done[ip] == -1
        done[ip] = 1
        new_clusters.append(O.clusters[ip])
        for w_max in xrange(w_max,-1,-1):
          if (len(candi[w_max]) != 0):
            break
        else:
          break
    assert len(new_clusters) == len(O.clusters)
    assert done.count(1) == len(done)
    del O.clusters[:]
    O.clusters.extend(new_clusters)
    O.refresh_indices()
    O.hinge_edges = [None] * n_clusters
    for pe in hinge_edges:
      cij = O.cluster_indices[pe[1]]
      assert O.hinge_edges[cij] is None
      O.hinge_edges[cij] = pe
    O.loop_edges.sort()

  def roots(O):
    assert O.hinge_edges is not None
    result = []
    for i,pe in enumerate(O.hinge_edges):
      if (pe[0] == -1):
        result.append(i)
    return result

  def tree_ids(O):
    assert O.hinge_edges is not None
    result = []
    tid = 0
    for i,pe in enumerate(O.hinge_edges):
      if (pe[0] == -1):
        result.append(tid)
        tid += 1
      else:
        result.append(result[O.cluster_indices[pe[0]]])
    return result

  def find_loop_edge_bendings(O, edge_sets):
    assert O.loop_edges is not None
    assert O.loop_edge_bendings is None
    leb = set()
    for i,j in O.loop_edges:
      for k in edge_sets[i]:
        if (k == j): continue
        assert k not in edge_sets[j]
        leb.add(tuple(sorted((j,k))))
      for k in edge_sets[j]:
        if (k == i): continue
        assert k not in edge_sets[i]
        leb.add(tuple(sorted((i,k))))
    O.loop_edge_bendings = sorted(leb)

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

  def __init__(O, n_vertices, edge_list, size_max=8):
    O.n_vertices = n_vertices
    O.edge_list = edge_list
    O.size_max = size_max
    O.edge_sets = construct_edge_sets(
      n_vertices=n_vertices, edge_list=edge_list)
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
