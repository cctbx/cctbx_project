from scitbx.graph.utils import construct_edge_sets, extract_edge_list

class cluster_manager(object):

  __slots__ = [
    "cluster_indices", "clusters",
    "merge_clusters_with_multiple_connections_passes",
    "overlapping_rigid_clusters",
    "hinge_edges", "loop_edges",
    "loop_edge_bendings"]

  def __init__(O, n_vertices):
    O.cluster_indices = range(n_vertices)
    O.clusters = []
    for i in xrange(n_vertices):
      O.clusters.append([i])
    O.merge_clusters_with_multiple_connections_passes = 0
    O.overlapping_rigid_clusters = None
    O.hinge_edges = None
    O.loop_edges = None
    O.loop_edge_bendings = None

  def connect_clusters(O, cii, cij, optimize):
    assert O.hinge_edges is None
    if (cii == cij): return
    ci = O.cluster_indices
    ccij = O.clusters[cij]
    ccii = O.clusters[cii]
    if (not optimize or len(ccij) <= len(ccii)):
      for k in ccij: ci[k] = cii
      ccii.extend(ccij)
      del ccij[:]
      return cii
    for k in ccii: ci[k] = cij
    ccij.extend(ccii)
    del ccii[:]
    return cij

  def connect_vertices(O, i, j, optimize):
    assert O.hinge_edges is None
    ci = O.cluster_indices
    return O.connect_clusters(cii=ci[i], cij=ci[j], optimize=optimize)

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

  def merge_clusters_with_multiple_connections(O, edge_sets):
    while True:
      O.merge_clusters_with_multiple_connections_passes += 1
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
            assert O.connect_clusters(cii=cii, cij=cij, optimize=False) == cii
            repeat = True
      if (not repeat):
        break
    O.tidy()

  def cluster_edge_sets(O, edge_list):
    result = []
    for i in xrange(len(O.clusters)):
      result.append(set())
    ci = O.cluster_indices
    for i,j in edge_list:
      cii = ci[i]
      cij = ci[j]
      if (cii == cij): continue
      result[cii].add(cij)
      result[cij].add(cii)
    return result

  def sort_by_overlapping_rigid_cluster_sizes(O, edge_sets):
    assert O.overlapping_rigid_clusters is None
    O.overlapping_rigid_clusters = []
    cii_orcs = []
    for cii,cluster in enumerate(O.clusters):
      c = set(cluster)
      for i in cluster:
        c.update(edge_sets[i])
      cii_orcs.append((cii, len(c)))
      O.overlapping_rigid_clusters.append(tuple(sorted(c)))
    O.overlapping_rigid_clusters.sort()
    def cmp_elems(a, b):
      if (a[1] > b[1]): return -1
      if (a[1] < b[1]): return 1
      return cmp(a[0], b[0])
    cii_orcs.sort(cmp_elems)
    new_clusters = []
    for cii,orcs in cii_orcs:
      new_clusters.append(O.clusters[cii])
    del O.clusters[:]
    O.clusters.extend(new_clusters)
    O.refresh_indices()
    return [orcs for cii,orcs in cii_orcs]

  def construct_spanning_trees(O, edge_sets):
    assert O.hinge_edges is None
    orcs = O.sort_by_overlapping_rigid_cluster_sizes(edge_sets=edge_sets)
    n_clusters = len(O.clusters)
    hinge_edges = [(-1,c[0]) for c in O.clusters]
    O.loop_edges = []
    if (n_clusters == 0): w_max = -1
    else:                 w_max = orcs[0]
    candi = []
    for i in xrange(w_max+1):
      candi.append([])
    done = [0] * n_clusters
    cluster_perm = []
    for ip in xrange(len(O.clusters)):
      he = hinge_edges[ip]
      if (he[0] != -1): continue
      done[ip] = 1
      cluster_perm.append(ip)
      def set_loop_or_hinge_edge(w_max):
        if (done[cij] != 0):
          O.loop_edges.append((i,j))
        else:
          done[cij] = -1
          w = orcs[cij]
          candi[w].append(cij)
          hinge_edges[cij] = (i,j)
          if (w_max < w): w_max = w
        return w_max
      w_max = 0
      for i in O.clusters[ip]:
        for j in edge_sets[i]:
          cij = O.cluster_indices[j]
          if (cij == ip): continue
          w_max = set_loop_or_hinge_edge(w_max=w_max)
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
            w_max = set_loop_or_hinge_edge(w_max=w_max)
        assert done[ip] == -1
        done[ip] = 1
        cluster_perm.append(ip)
        he = hinge_edges[ip]
        if (he[0] != -1):
          O.clusters[O.cluster_indices[he[0]]].append(he[1])
          O.clusters[ip].remove(he[1])
        for w_max in xrange(w_max,-1,-1):
          if (len(candi[w_max]) != 0):
            break
        else:
          break
    assert len(cluster_perm) == n_clusters
    assert done.count(1) == len(done)
    new_clusters = []
    O.hinge_edges = []
    for cii in cluster_perm:
      c = O.clusters[cii]
      if (len(c) != 0):
        new_clusters.append(sorted(c))
        O.hinge_edges.append(hinge_edges[cii])
    del O.clusters[:]
    O.clusters.extend(new_clusters)
    O.refresh_indices()
    O.loop_edges.sort()

  def roots(O):
    assert O.hinge_edges is not None
    result = []
    for i,he in enumerate(O.hinge_edges):
      if (he[0] == -1):
        result.append(i)
    return result

  def tree_ids(O):
    assert O.hinge_edges is not None
    result = []
    tid = 0
    for he in O.hinge_edges:
      if (he[0] == -1):
        result.append(tid)
        tid += 1
      else:
        result.append(result[O.cluster_indices[he[0]]])
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

def find_paths(edge_sets, iv):
  result = {}
  for jv in edge_sets[iv]:
    result[jv] = []
  result = {iv: result}
  path = []
  def depth_first_search(jv, kv):
    path.append(kv)
    def new_path_list():
      if (len(path) < 3): return []
      if (kv in path[1:-1]): return [] # XXX avoid need for this
      return [path[1:-1]]
    for lv in edge_sets[kv]:
      if (lv == jv): continue
      result_lv = result.get(lv)
      if (result_lv is None):
        result[lv] = {kv: new_path_list()}
        if (len(path) != 6):
          depth_first_search(jv=kv, kv=lv)
      else:
        result_lv_kv = result_lv.get(kv)
        if (result_lv_kv is None):
          result_lv[kv] = new_path_list()
          if (len(path) != 6):
            depth_first_search(jv=kv, kv=lv)
        elif (    len(path) > 2
              and kv not in path[1:-1]): # XXX avoid need for this
          result_lv_kv.append(path[1:-1])
    path.pop()
  depth_first_search(jv=-1, kv=iv)
  return result

class construct(object):

  def __init__(O, n_vertices, edge_list):
    O.n_vertices = n_vertices
    O.edge_list = edge_list
    O.edge_sets = construct_edge_sets(
      n_vertices=n_vertices, edge_list=edge_list)
    O.cluster_manager = cluster_manager(n_vertices=n_vertices)
    O._find_paths()
    O.cluster_manager.tidy()
    O.find_cluster_loop_repeats = None

  def _find_paths(O):
    for iv in xrange(O.n_vertices):
      jv_kv_paths = find_paths(edge_sets=O.edge_sets, iv=iv)
      for jv,kv_paths in jv_kv_paths.items():
        if (jv == iv):
          for kv,paths in kv_paths.items():
            if (len(paths) == 0): continue
            O.cluster_manager.connect_vertices(i=iv, j=kv, optimize=True)
            for path in paths:
              assert len(path) != 0
              for lv in path:
                O.cluster_manager.connect_vertices(i=iv, j=lv, optimize=True)
        elif (len(kv_paths) > 2):
          paths_by_length = [[], [], [], [], []]
          for kv,paths in kv_paths.items():
            if (len(paths) == 0): # XXX adjust find_paths to include []?
              if (kv != iv):
                paths_by_length[0].append((kv, []))
            else:
              for path in paths:
                assert len(path) < 5
                paths_by_length[len(path)].append((kv, path))
          n_l1_l2_l3_lt_7 = 0
          for l1 in xrange(5):
            for l2 in xrange(min(5,7-l1)):
              for l3 in xrange(min(5,7-l1-l2)):
                assert l1+l2+l3 < 7
                n_l1_l2_l3_lt_7 += 1
                for kv1,path1 in paths_by_length[l1]:
                  if (iv in path1): continue
                  if (jv in path1): continue
                  for kv2,path2 in paths_by_length[l2]:
                    if (kv2 == kv1): continue
                    if (iv in path2): continue
                    if (jv in path2): continue
                    for kv3,path3 in paths_by_length[l3]:
                      if (kv3 == kv1): continue
                      if (kv3 == kv2): continue
                      if (iv in path3): continue
                      if (jv in path3): continue
                      sp1 = set(path1)
                      sp2 = set(path2)
                      sp3 = set(path3)
                      if (not sp1.isdisjoint(sp2)): continue
                      if (not sp1.isdisjoint(sp3)): continue
                      if (not sp2.isdisjoint(sp3)): continue
                      # XXX next six conditions needed?
                      if (kv1 in sp2): continue
                      if (kv1 in sp3): continue
                      if (kv2 in sp1): continue
                      if (kv2 in sp3): continue
                      if (kv3 in sp1): continue
                      if (kv3 in sp2): continue
                      for lvs in [[jv, kv1, kv2, kv3], path1, path2, path3]:
                        for lv in lvs:
                          O.cluster_manager.connect_vertices(
                            i=iv, j=lv, optimize=True)
          assert n_l1_l2_l3_lt_7 == 72

  def find_cluster_loops(O):
    assert O.find_cluster_loop_repeats is None
    O.find_cluster_loop_repeats = -1
    cm = O.cluster_manager
    while True:
      O.find_cluster_loop_repeats += 1
      cm.merge_clusters_with_multiple_connections(edge_sets=O.edge_sets)
      ces = cm.cluster_edge_sets(edge_list=O.edge_list)
      cel = extract_edge_list(edge_sets=ces)
      ctt = construct(n_vertices=len(cm.clusters), edge_list=cel)
      ccm = ctt.cluster_manager
      ccm.merge_clusters_with_multiple_connections(edge_sets=ctt.edge_sets)
      if (len(ccm.clusters) == len(cm.clusters)):
        break
      for cc in ccm.clusters:
        cii = cc[0]
        for cij in cc[1:]:
          cii = cm.connect_clusters(cii=cii, cij=cij, optimize=True)
      cm.tidy()
    return O

  def finalize(O):
    O.find_cluster_loops()
    cm = O.cluster_manager
    cm.construct_spanning_trees(edge_sets=O.edge_sets)
    cm.find_loop_edge_bendings(edge_sets=O.edge_sets)
    return O
