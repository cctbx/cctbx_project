from __future__ import division
from scitbx.graph.utils import \
  construct_edge_sets, extract_edge_list, sub_edge_list, tree_marking
from libtbx import slots_getstate_setstate
import math

class cluster_manager(slots_getstate_setstate):

  __slots__ = [
    "fixed_vertex_lists",
    "cluster_indices",
    "clusters",
    "merge_clusters_with_multiple_connections_passes",
    "hinge_edges",
    "loop_edges",
    "loop_edge_bendings",
    "fixed_hinges"]

  def __init__(O,
        n_vertices,
        all_in_one_rigid_body=False,
        fixed_vertex_lists=()):
    O.fixed_vertex_lists = fixed_vertex_lists
    if (len(O.fixed_vertex_lists) == 0):
      if (all_in_one_rigid_body):
        O.cluster_indices = [0] * n_vertices
        O.clusters = [range(n_vertices)]
      else:
        O.cluster_indices = range(n_vertices)
        O.clusters = []
        for i in xrange(n_vertices):
          O.clusters.append([i])
    else:
      assert not all_in_one_rigid_body # not implemented
      O.cluster_indices = [-1] * n_vertices
      O.clusters = []
      for fixed_vertices in O.fixed_vertex_lists:
        assert len(fixed_vertices) != 0
        for i in fixed_vertices:
          assert O.cluster_indices[i] == -1
          O.cluster_indices[i] = len(O.clusters)
        O.clusters.append(list(fixed_vertices))
      for i in xrange(n_vertices):
        if (O.cluster_indices[i] != -1): continue
        O.cluster_indices[i] = len(O.clusters)
        O.clusters.append([i])
    O.merge_clusters_with_multiple_connections_passes = 0
    O.hinge_edges = None
    O.loop_edges = None
    O.loop_edge_bendings = None
    O.fixed_hinges = None

  def all_in_one_rigid_body(O):
    return len(O.clusters) == 1

  def show_summary(O, out=None, prefix=""):
    from libtbx.utils import xlen, plural_s
    import sys
    if (out is None): out = sys.stdout
    print >> out, prefix+"number of fixed vertex lists:", \
      len(O.fixed_vertex_lists)
    print >> out, prefix+"number of fixed vertices:", \
      sum([len(fixed_vertices) for fixed_vertices in O.fixed_vertex_lists])
    print >> out, prefix+"number of clusters:", len(O.clusters)
    print >> out, prefix+"merge clusters with multiple connections: %d pass%s"\
      % plural_s(O.merge_clusters_with_multiple_connections_passes, "es")
    print >> out, prefix+"number of hinge edges:", xlen(O.hinge_edges)
    print >> out, prefix+"number of loop edges:", xlen(O.loop_edges)
    print >> out, prefix+"number of loop edge bendings:", \
      xlen(O.loop_edge_bendings)
    print >> out, prefix+"number of fixed hinges:", xlen(O.fixed_hinges)
    return O

  def fixed_vertices_given_cluster_index_dict(O):
    result = {}
    ci = O.cluster_indices
    for fixed_vertices in O.fixed_vertex_lists:
      i = ci[fixed_vertices[0]]
      assert i not in result
      result[i] = fixed_vertices
    return result

  def connect_clusters(O, cii, cij, optimize):
    assert O.hinge_edges is None
    if (cii == cij): return None
    lfvl = len(O.fixed_vertex_lists)
    if (cii < lfvl and cij < lfvl):
      raise RuntimeError(
        "connect_clusters():"
        " fixed vertex lists in same connected tree.")
    ci = O.cluster_indices
    ccij = O.clusters[cij]
    ccii = O.clusters[cii]
    if ((not optimize
         or len(ccij) <= len(ccii)
         or cii < lfvl)
        and (cij >= lfvl or lfvl == 0)):
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
      if (len(a) != 0 and len(b) != 0):
        fa = O.cluster_indices[a[0]] < len(O.fixed_vertex_lists)
        fb = O.cluster_indices[b[0]] < len(O.fixed_vertex_lists)
        if (fa):
          if (not fb): return -1
        else:
          if (fb): return 1
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

  def overlapping_rigid_clusters(O, edge_sets):
    assert O.hinge_edges is None
    result = []
    for cluster in O.clusters:
      def add_cluster_and_connected_vertices():
        c = set(cluster)
        for i in cluster:
          c.update(edge_sets[i])
        result.append(tuple(sorted(c)))
      if (len(cluster) != 1):
        add_cluster_and_connected_vertices()
      else:
        i = cluster[0]
        esi = edge_sets[i]
        if (len(esi) != 1):
          add_cluster_and_connected_vertices()
        else:
          j = list(esi)[0]
          if (len(edge_sets[j]) == 1 and j > i):
            result.append((i,j))
    return result

  def determine_weighted_order_for_construct_spanning_tree(O, edge_sets):
    fixed_vertex_info = [0] * len(O.clusters)
    for fixed_vertices in O.fixed_vertex_lists:
      lvf = len(fixed_vertices)
      assert lvf != 0
      i = fixed_vertices[0]
      cii = O.cluster_indices[i]
      fixed_vertex_info[cii] = lvf
      if (lvf == 1):
        for j in edge_sets[i]:
          cij = O.cluster_indices[j]
          if (cij == cii): continue
          if (fixed_vertex_info[cij] != 0):
            raise RuntimeError(
              "determine_weighted_order_for_construct_spanning_tree():"
              " fixed vertex lists in same connected tree.")
          fixed_vertex_info[cij] = -1
    cii_orcs = []
    for cii,cluster in enumerate(O.clusters):
      c = set(cluster)
      for i in cluster:
        c.update(edge_sets[i])
      cii_orcs.append((cii, len(c)))
    def cmp_elems(a, b):
      fa = abs(fixed_vertex_info[a[0]])
      fb = abs(fixed_vertex_info[b[0]])
      if (fa > fb): return -1
      if (fa < fb): return 1
      if (a[1] > b[1]): return -1
      if (a[1] < b[1]): return 1
      return cmp(a[0], b[0])
    cii_orcs.sort(cmp_elems)
    return cii_orcs, fixed_vertex_info

  def construct_spanning_trees(O, edge_sets):
    assert O.hinge_edges is None
    if (edge_sets is None):
      assert O.all_in_one_rigid_body()
      O.hinge_edges = [(-1,0)]
      O.loop_edges = []
      return
    cii_orcs, fixed_vertex_info = \
      O.determine_weighted_order_for_construct_spanning_tree(
        edge_sets=edge_sets)
    n_clusters = len(O.clusters)
    hinge_edges = [(-1,c[0]) for c in O.clusters]
    O.loop_edges = []
    if (n_clusters == 0): w_max = -1
    else:                 w_max = max([orcs for cii,orcs in cii_orcs])
    candi = []
    for i in xrange(w_max+1):
      candi.append([])
    cii_wo_given_cii = [n_clusters] * n_clusters
    for ip_wo in xrange(n_clusters):
      ip = cii_orcs[ip_wo][0]
      cii_wo_given_cii[ip] = ip_wo
    done = [0] * n_clusters
    cluster_perm = []
    for ip_wo in xrange(n_clusters):
      ip = cii_orcs[ip_wo][0]
      he = hinge_edges[ip]
      if (he[0] != -1): continue
      have_fixed = [(fixed_vertex_info[ip] > 0)]
      done[ip] = 1
      cluster_perm.append(ip)
      def set_loop_or_hinge_edge(w_max):
        if (fixed_vertex_info[cij] > 0):
          if (have_fixed[0]):
            raise RuntimeError(
              "construct_spanning_trees():"
              " fixed vertex lists in same connected tree.")
          have_fixed[0] = True
        if (done[cij] != 0):
          O.loop_edges.append((i,j))
        else:
          done[cij] = -1
          cij_wo = cii_wo_given_cii[cij]
          w = cii_orcs[cij_wo][1]
          candi[w].append(cij_wo)
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
        ip_wo = n_clusters
        cw = candi[w_max]
        for k in xrange(len(cw)):
          if (ip_wo > cw[k]):
            kp = k
            ip_wo = cw[k]
        if (kp is None):
          break
        del cw[kp]
        ip = cii_orcs[ip_wo][0]
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
    if (edge_sets is None):
      assert O.all_in_one_rigid_body()
      O.loop_edge_bendings = []
      return
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

  def fix_near_singular_hinges(O, sites, angular_tolerance_deg):
    assert O.loop_edge_bendings is not None
    assert O.fixed_hinges is None
    O.fixed_hinges = []
    if (O.all_in_one_rigid_body()): return
    if (sites is None): return
    if (hasattr(sites, "accessor")):
      from scitbx import matrix
      sites = matrix.col_list(sites)
    abs_cos_limit = abs(math.cos(math.radians(angular_tolerance_deg)))
    for jc in xrange(len(O.clusters)-1,-1,-1):
      hi,hj = O.hinge_edges[jc]
      if (hi == -1):
        continue
      pivot = sites[hi]
      axis = sites[hj] - pivot
      for i in O.clusters[jc]:
        abs_cos = abs(axis.cos_angle(sites[i] - pivot, value_if_undefined=1))
        if (abs_cos < abs_cos_limit):
          break
      else:
        O.fixed_hinges.append(O.hinge_edges[jc])
        del O.hinge_edges[jc]
        ic = O.cluster_indices[hj]
        O.clusters[ic].extend(O.clusters[jc])
        O.clusters[ic].sort()
        ci = O.cluster_indices
        for i in O.clusters[jc]:
          ci[i] = ic
        del O.clusters[jc]
        for ic in xrange(jc,len(O.clusters)):
          for i in O.clusters[ic]:
            ci[i] = ic
    O.fixed_hinges.sort()

  def edge_classifier(O):
    return edge_classifier(cluster_manager=O)

class edge_classifier(object):

  def __init__(O, cluster_manager):
    O.cluster_manager = cluster_manager
    O.hinge_edge_set = set()
    for e in cluster_manager.hinge_edges:
      if (e[0] == -1): continue
      O.hinge_edge_set.add(tuple(sorted(e)))
    O.loop_edge_set = set([tuple(sorted(e))
      for e in cluster_manager.loop_edges])
    assert len(O.hinge_edge_set.intersection(O.loop_edge_set)) == 0
    O.fixed_hinge_set = set(cluster_manager.fixed_hinges)
    assert len(O.hinge_edge_set.intersection(O.fixed_hinge_set)) == 0
    assert len(O.loop_edge_set.intersection(O.fixed_hinge_set)) == 0

  def __call__(O, edge):
    edge = tuple(sorted(edge))
    if (edge in O.fixed_hinge_set): return "fixed"
    if (edge in O.hinge_edge_set): return "hinge"
    if (edge in O.loop_edge_set): return "loop"
    cm = O.cluster_manager
    cii, cij = [cm.cluster_indices[i] for i in edge]
    if (cii == cij and cm.hinge_edges[cii][0] == -1): return "base"
    return "intra"

class find_paths(object):

  def __init__(O, edge_sets):
    O.edge_sets = edge_sets
    O.in_path = [False] * len(O.edge_sets)

  def search_from(O, iv):
    edge_sets = O.edge_sets
    in_path = O.in_path
    loops = {}
    dendrites = {}
    path = []
    def depth_first_search(jv, kv):
      path.append(kv)
      in_path[kv] = True
      closing = False
      for lv in edge_sets[kv]:
        if (lv == jv): continue
        if (lv == iv):
          loops.setdefault(path[0], []).append(path[1:])
          closing = True
        elif (in_path[lv]):
          closing = True
      if (not closing and len(path) != 6):
        for lv in edge_sets[kv]:
          if (lv == jv): continue
          dendrites.setdefault(lv, []).append(set(path))
          depth_first_search(jv=kv, kv=lv)
      path.pop()
      in_path[kv] = False
    for jv in edge_sets[iv]:
      depth_first_search(jv=iv, kv=jv)
    return loops, dendrites

class construct(slots_getstate_setstate):

  __slots__ = [
    "n_vertices",
    "edge_list",
    "edge_sets",
    "cluster_manager",
    "external_clusters_connect_count",
    "find_cluster_loop_repeats"]

  def __init__(O,
        n_vertices=None,
        sites=None,
        edge_list=None,
        external_clusters=None,
        fixed_vertices=None,
        fixed_vertex_lists=None,
        near_singular_hinges_angular_tolerance_deg=5):
    assert [n_vertices, sites].count(None) == 1
    assert edge_list is not None
    assert [fixed_vertices, fixed_vertex_lists].count(None) != 0
    if (sites is not None):
      n_vertices = len(sites)
    O.n_vertices = n_vertices
    all_in_one_rigid_body = (edge_list == "all_in_one_rigid_body")
    if (all_in_one_rigid_body):
      assert external_clusters is None
      O.edge_list = None
      O.edge_sets = None
    else:
      O.edge_list = edge_list
      O.edge_sets = construct_edge_sets(
        n_vertices=n_vertices, edge_list=edge_list)
    if (fixed_vertex_lists is None):
      if (fixed_vertices is None or len(fixed_vertices) == 0):
        fixed_vertex_lists = ()
      else:
        assert O.edge_sets is not None # not implemented
        fixed_vertex_lists = tree_marking(edge_sets=O.edge_sets) \
          .partitions_of(vertex_indices=fixed_vertices)
    O.cluster_manager = cluster_manager(
      n_vertices=n_vertices,
      all_in_one_rigid_body=all_in_one_rigid_body,
      fixed_vertex_lists=fixed_vertex_lists)
    if (not all_in_one_rigid_body):
      O._find_paths()
      O._process_external_clusters(clusters=external_clusters)
    O.cluster_manager.tidy()
    O.find_cluster_loop_repeats = None
    if (sites is not None):
      O.build_tree()
      O.fix_near_singular_hinges(
        sites=sites,
        angular_tolerance_deg=near_singular_hinges_angular_tolerance_deg)

  def show_summary(O, vertex_labels, out=None, prefix=""):
    from libtbx.utils import xlen, plural_s
    import sys
    if (out is None): out = sys.stdout
    if (vertex_labels is None):
      fmt = "%%0%dd" % len(str(max(0, O.n_vertices-1)))
      vertex_labels = [fmt % i for i in xrange(O.n_vertices)]
    else:
      assert len(vertex_labels) == O.n_vertices
    print >> out, prefix+"number of vertices:", O.n_vertices
    print >> out, prefix+"number of edges:", xlen(O.edge_list)
    if (O.find_cluster_loop_repeats is None):
      print >> out, prefix+"find cluster loops: None"
    else:
      print >> out, prefix+"find cluster loops: %d repeat%s" % \
        plural_s(O.find_cluster_loop_repeats)
    cm = O.cluster_manager
    cm.show_summary(out=out, prefix=prefix)
    if (cm.fixed_hinges is not None):
      for i,j in cm.fixed_hinges:
        print >> out, prefix+"tardy fixed hinge:", vertex_labels[i]
        print >> out, prefix+"                  ", vertex_labels[j]
    return O

  def extract_edge_list(O):
    return extract_edge_list(edge_sets=O.edge_sets)

  def _find_paths(O):
    fp = find_paths(edge_sets=O.edge_sets)
    for iv in xrange(O.n_vertices):
      loops, dendrites = fp.search_from(iv=iv)
      #
      for jv,loops_through_jv in loops.items():
        have_small = False
        l5s = []
        for loop_through_jv in loops_through_jv:
          if (len(loop_through_jv) < 5):
            for kv in loop_through_jv:
              O.cluster_manager.connect_vertices(i=iv, j=kv, optimize=True)
            have_small = True
          else:
            l5s.append(loop_through_jv)
        if (have_small):
          for loop_through_jv in l5s:
            for kv in loop_through_jv:
              O.cluster_manager.connect_vertices(i=iv, j=kv, optimize=True)
      #
      for jv,sps_to_jv in dendrites.items():
        if (len(sps_to_jv) < 3): continue
        sps_by_length = [None, [], [], [], [], []]
        for sp_to_jv in sps_to_jv:
          sps_by_length[len(sp_to_jv)].append(sp_to_jv)
        n_l1_l2_l3_lt_10 = 0
        for l1 in xrange(1,6):
          for l2 in xrange(1,min(6,10-l1)):
            for l3 in xrange(1,min(6,10-l1-l2)):
              n_l1_l2_l3_lt_10 += 1
              for sp1 in sps_by_length[l1]:
                for sp2 in sps_by_length[l2]:
                  sp12 = sp1.union(sp2)
                  if (len(sp12) != len(sp1) + len(sp2)): continue
                  for sp3 in sps_by_length[l3]:
                    sp123 = sp12.union(sp3)
                    if (len(sp123) != len(sp12) + len(sp3)): continue
                    O.cluster_manager.connect_vertices(
                      i=iv, j=jv, optimize=True)
                    for kv in sp123:
                      O.cluster_manager.connect_vertices(
                        i=iv, j=kv, optimize=True)
        assert n_l1_l2_l3_lt_10 == 72

  def _process_external_clusters(O, clusters):
    O.external_clusters_connect_count = 0
    if (clusters is None): return
    cv = O.cluster_manager.connect_vertices
    for cluster in clusters:
      sub = sub_edge_list(edge_sets=O.edge_sets, vertex_indices=cluster)
      sub_edge_sets = sub.edge_sets()
      for i_sub,j_sub in sub.edge_list:
        if (len(sub_edge_sets[i_sub]) == 1): continue
        if (len(sub_edge_sets[j_sub]) == 1): continue
        if (cv(i=cluster[i_sub], j=cluster[j_sub], optimize=True) is not None):
          O.external_clusters_connect_count += 1

  def find_cluster_loops(O):
    assert O.find_cluster_loop_repeats is None
    if (O.edge_sets is None):
      O.find_cluster_loop_repeats = 0
      return
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

  def build_tree(O):
    O.find_cluster_loops()
    cm = O.cluster_manager
    cm.construct_spanning_trees(edge_sets=O.edge_sets)
    cm.find_loop_edge_bendings(edge_sets=O.edge_sets)
    return O

  def fix_near_singular_hinges(O, sites, angular_tolerance_deg=5):
    O.cluster_manager.fix_near_singular_hinges(
      sites=sites, angular_tolerance_deg=angular_tolerance_deg)
    return O

  def rmsd_calculator(O):
    return O.rmsd_calculation

  def rmsd_calculation(O, sites_cart_1, sites_cart_2):
    return sites_cart_1.rms_difference(sites_cart_2.select(
      O.rmsd_permutation(
        sites_cart_1=sites_cart_1,
        sites_cart_2=sites_cart_2)))

  def rmsd_permutation(O, sites_cart_1, sites_cart_2):
    "simple, limited handling of flipped sites"
    assert sites_cart_1.size() == len(O.edge_sets)
    assert sites_cart_2.size() == len(O.edge_sets)
    from scitbx.array_family import flex
    result = flex.size_t_range(len(O.edge_sets))
    for i,esi in enumerate(O.edge_sets):
      if (len(esi) not in[2, 3]): continue
      n1 = flex.size_t()
      for j in esi:
        if (len(O.edge_sets[j]) == 1):
          n1.append(j)
      if (len(n1) != 2): continue
      n1_rev = flex.size_t(reversed(n1))
      pair_1 = sites_cart_1.select(n1)
      rmsd_1 = pair_1.rms_difference(sites_cart_2.select(n1))
      rmsd_2 = pair_1.rms_difference(sites_cart_2.select(n1_rev))
      if (rmsd_2 < rmsd_1*(1-1e-6)):
        result.set_selected(n1, n1_rev)
    return result

  def viewer_lines_with_colors_legend(O, include_loop_edge_bendings):
    result = [
      "Edge colors:",
      "  turquoise: intra-base-cluster, six degrees of freedom",
      "  green:     rotatable bond, one degree of freedom",
      "  blue:      intra-cluster",
      "  red:       loop edge (restrained only)",
      "  pink:      near singular rotatable bond (fixed)"]
    if (include_loop_edge_bendings): result.append(
      "  purple:    loop bending edge (restrained only)")
    return result

  def viewer_lines_with_colors(O, include_loop_edge_bendings):
    result = []
    colors = {
      "base":  (0,1,1),
      "hinge": (0,1,0),
      "intra": (0,0,1),
      "loop":  (1,0,0),
      "fixed": (1,182/255,193/255)}
    ec = O.cluster_manager.edge_classifier()
    for line in O.edge_list:
      result.append((line, colors[ec(edge=line)]))
    if (include_loop_edge_bendings):
      for line in O.cluster_manager.loop_edge_bendings:
        result.append((line, (0.5,0,0.5)))
    return result
