from scitbx.graph.tardy_tree import \
  cluster_manager, construct_edge_sets, construct
from StringIO import StringIO
import sys

def exercise_cluster_manager():
  cm = cluster_manager(n_vertices=0)
  assert cm.cluster_indices == []
  assert cm.clusters == []
  cm = cluster_manager(n_vertices=5)
  for p in xrange(2):
    assert cm.cluster_indices == [0,1,2,3,4]
    assert cm.clusters == [[0],[1],[2],[3],[4]]
    cm.connect_vertices(i=0, j=0, optimize=True)
  for p in xrange(2):
    cm.connect_vertices(i=1, j=3, optimize=True)
    assert cm.cluster_indices == [0,1,2,1,4]
    assert cm.clusters == [[0],[1,3],[2],[],[4]]
  for p in xrange(2):
    cm.connect_vertices(i=0, j=3, optimize=True)
    assert cm.cluster_indices == [1,1,2,1,4]
    for q in xrange(2):
      assert cm.clusters == [[],[1,3,0],[2],[],[4]]
      cm.refresh_indices()
  cm.connect_vertices(i=2, j=4, optimize=True)
  assert cm.clusters == [[],[1,3,0],[2,4],[],[]]
  assert cm.cluster_indices == [1,1,2,1,2]
  cm.connect_vertices(i=2, j=3, optimize=True)
  assert cm.clusters == [[],[1,3,0,2,4],[],[],[]]
  assert cm.cluster_indices == [1,1,1,1,1]
  cm.tidy()
  assert cm.clusters == [[0,1,2,3,4]]
  assert cm.cluster_indices == [0,0,0,0,0]
  #
  cm = cluster_manager(n_vertices=6)
  cm.connect_vertices(i=3, j=0, optimize=True)
  cm.connect_vertices(i=2, j=4, optimize=True)
  cm.connect_vertices(i=1, j=2, optimize=True)
  cm.tidy()
  assert cm.clusters == [[1,2,4],[0,3],[5]]
  assert cm.cluster_indices == [1,0,0,1,0,2]
  edges = [(0,1), (0,2), (3,4), (4,5)]
  cm.connect_vertices(i=4, j=5, optimize=True)
  cm.tidy()
  assert cm.clusters == [[1,2,4,5],[0,3]]
  assert cm.cluster_indices == [1,0,0,1,0,0]
  es = construct_edge_sets(n_vertices=6, edge_list=edges)
  cm.construct_spanning_trees(edge_sets=es)
  assert cm.clusters == [[0,1,2,4,5],[3]]
  assert cm.cluster_indices == [0,0,0,1,0,0]
  assert cm.hinge_edges == [(-1,1), (1,0)]
  assert cm.loop_edges == [(2,0), (4,3)]
  assert cm.roots() == [0]
  assert cm.tree_ids() == [0,0]
  cm.find_loop_edge_bendings(edge_sets=es)
  assert cm.loop_edge_bendings == [(1,2), (3,5)]
  #
  cm.cluster_indices = [1,0,0,1,0,0]
  cm.clusters = [[1,2,4,5],[0,3]]
  cm.hinge_edges = None
  cm.loop_edges = None
  cm.loop_edge_bendings = None
  cm.merge_clusters_with_multiple_connections(edge_sets=es)
  assert cm.clusters == [[0,1,2,3,4,5]]
  assert cm.cluster_indices == [0,0,0,0,0,0]

class test_case_data(object):

  def __init__(O,
        art,
        n_vertices,
        edge_list,
        clusters1,
        hinge_edges1,
        roots1,
        tree_ids1,
        clusters1_5=None,
        hinge_edges1_5=None,
        roots1_5=None,
        tree_ids1_5=None,
        loop_edges1_5=[],
        loop_edge_bendings1_5=[],
        clusters2_5=None,
        hinge_edges2_5=None,
        loop_edges2_5=None,
        loop_edge_bendings2_5=None):
    assert art[0] == "\n"
    assert (hinge_edges1_5 is None) == (clusters1_5 is None)
    assert (roots1_5 is None) == (hinge_edges1_5 is None)
    assert (tree_ids1_5 is None) == (hinge_edges1_5 is None)
    assert (len(loop_edges1_5) == 0) == (hinge_edges1_5 is None)
    assert (len(loop_edge_bendings1_5) == 0) == (len(loop_edges1_5) == 0)
    if (clusters2_5 is None):
      assert hinge_edges2_5 is None
      clusters2_5 = clusters1_5
      hinge_edges2_5 = hinge_edges1_5
    else:
      assert hinge_edges2_5 is not None
      assert clusters2_5 != clusters1_5
      assert hinge_edges2_5 != hinge_edges1_5
    assert (loop_edges2_5 is None) == (hinge_edges2_5 is None)
    assert (loop_edge_bendings2_5 is None) == (loop_edges2_5 is None)
    if (loop_edges2_5 is None): loop_edges2_5 = []
    if (loop_edge_bendings2_5 is None): loop_edge_bendings2_5 = []
    O.art = art[1:]
    O.n_vertices = n_vertices
    O.edge_list = edge_list
    O.clusters1 = clusters1
    O.hinge_edges1 = hinge_edges1
    O.roots1 = roots1
    O.tree_ids1 = tree_ids1
    O.clusters1_5 = clusters1_5
    O.hinge_edges1_5 = hinge_edges1_5
    O.roots1_5 = roots1_5
    O.tree_ids1_5 = tree_ids1_5
    O.loop_edges1_5 = loop_edges1_5
    O.loop_edge_bendings1_5 = loop_edge_bendings1_5
    O.clusters2_5 = clusters2_5
    O.hinge_edges2_5 = hinge_edges2_5
    O.loop_edges2_5 = loop_edges2_5
    O.loop_edge_bendings2_5 = loop_edge_bendings2_5

test_cases = [
  test_case_data(
    art=r"""
0
""",
    n_vertices=1,
    edge_list=[],
    clusters1=[[0]],
    hinge_edges1=[(-1,0)],
    roots1=[0],
    tree_ids1=[0]),
  test_case_data(
    art=r"""
0 - 1
""",
    n_vertices=2,
    edge_list=[(0,1)],
    clusters1=[[0], [1]],
    hinge_edges1=[(-1,0)],
    roots1=[0],
    tree_ids1=[0]),
  test_case_data(
    art=r"""
6-membered loop
""",
    n_vertices=6,
    edge_list=[(0,1), (0,5), (1,2), (2,3), (3,4), (4,5)],
    clusters1=[[0, 1, 2, 3, 4, 5]],
    hinge_edges1=[(-1,0)],
    roots1=[0],
    tree_ids1=[0],
    clusters1_5=[[0], [1], [2], [3], [4], [5]],
    hinge_edges1_5=[(-1, 0), (0,1), (1, 2), (2, 3)],
    roots1_5=[0],
    tree_ids1_5=[0, 0, 0, 0],
    loop_edges1_5=[(4,5)],
    loop_edge_bendings1_5=[(0,4), (3,5)],
    loop_edges2_5=[(4,5)],
    loop_edge_bendings2_5=[(0,4), (3,5)]),
  test_case_data(
    art=r"""
8-membered loop
""",
    n_vertices=8,
    edge_list=[(0,1), (0,7), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7)],
    clusters1=[[0, 1, 2, 3, 4, 5, 6, 7]],
    hinge_edges1=[(-1,0)],
    roots1=[0],
    tree_ids1=[0],
    clusters1_5=[[0], [1], [2], [3], [4], [5], [6], [7]],
    hinge_edges1_5=[(-1,0), (0,1), (1,2), (2,3), (3,4), (4,5)],
    roots1_5=[0],
    tree_ids1_5=[0, 0, 0, 0, 0, 0],
    loop_edges1_5=[(6,7)],
    loop_edge_bendings1_5=[(0,6), (5,7)],
    loop_edges2_5=[(6,7)],
    loop_edge_bendings2_5=[(0,6), (5,7)]),
  test_case_data(
    art=r"""
           7
           |
  5        1
   \     /   \
    4 - 0     2 - 8
   /         /
  6        3
           |
           9
""",
    n_vertices=10,
    edge_list=[
      (0,1), (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6)],
    clusters1=[[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]],
    hinge_edges1=[(-1,1), (1,2), (1,0), (0,4), (2,3)],
    roots1=[0],
    tree_ids1=[0] * 5),
  test_case_data(
    art=r"""
           7
           |
  5        1
   \     /   \
    4 - 0     2 - 8
   /         //
  6        3 /
           |/
           9
""",
    n_vertices=10,
    edge_list=[
      (0,1), (0,4),
      (1,2), (1,7),
      (2,3), (2,8), (2,9),
      (3,9),
      (4,5), (4,6)],
    clusters1=[[2,3,9], [0], [1], [4], [5], [6], [7], [8]],
    hinge_edges1=[(-1,2), (2,1), (1,0), (0,4)],
    roots1=[0],
    tree_ids1=[0] * 4),
  test_case_data(
    art=r"""
           7
           |
  5        1
   \     /   \
    4 - 0     2 - 8
   /     \   /
  6        3
           |
           9
""",
    n_vertices=10,
    edge_list=[
      (0,1), (0,3), (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6)],
    clusters1=[[0,1,2,3], [4], [5], [6], [7], [8], [9]],
    hinge_edges1=[(-1,0), (0,4)],
    roots1=[0],
    tree_ids1=[0] * 2),
  test_case_data(
    art=r"""
           7
           |
  5        1--\
   \     /   \ \-\
    4 - 0     2 - 8
   /     \   /
  6        3
           |
           9
""",
    n_vertices=10,
    edge_list=[
      (0,1), (0,3), (0,4),
      (1,7), (1,8),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6)],
    clusters1=[[0, 1, 2, 3, 8], [4], [5], [6], [7], [9]],
    hinge_edges1=[(-1,0), (0,4)],
    roots1=[0],
    tree_ids1=[0] * 2),
  test_case_data(
    art=r"""
           7
           | \
  5        1   \
   \     /       \
    4 - 0     2 - 8
   /     \   /
  6        3
           |
           9
""",
    n_vertices=10,
    edge_list=[
      (0,1), (0,3), (0,4),
      (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6),
      (7,8)],
    clusters1=[[0, 1, 2, 3, 7, 8], [4], [5], [6], [9]],
    hinge_edges1=[(-1,0), (0,4)],
    roots1=[0],
    tree_ids1=[0] * 2,
    clusters1_5=[[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]],
    hinge_edges1_5=[(-1, 0), (0, 3), (0, 4), (0, 1), (3, 2)],
    roots1_5=[0],
    tree_ids1_5=[0] * 5,
    loop_edges1_5=[(7,8)],
    loop_edge_bendings1_5=[(1,8), (2,7)],
    loop_edges2_5=[(7,8)],
    loop_edge_bendings2_5=[(1,8), (2,7)]),
  test_case_data(
    art=r"""
           7
           | \
  5        1   \
   \     /   \   \
    4 - 0     2 - 8
   /     \   /
  6        3
           |
           9
""",
    n_vertices=10,
    edge_list=[
      (0,1), (0,3), (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6),
      (7,8)],
    clusters1=[[0, 1, 2, 3, 7, 8], [4], [5], [6], [9]],
    hinge_edges1=[(-1,0), (0,4)],
    roots1=[0],
    tree_ids1=[0] * 2),
  test_case_data(
    art=r"""
           7
           |
  5        1
   \     /   \
    4 - 0     2 - 8
   /     \   /   /
  6        3   /
           | /
           9
""",
    n_vertices=10,
    edge_list=[
      (0,1), (0,3), (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6),
      (8,9)],
    clusters1=[[0, 1, 2, 3, 8, 9], [4], [5], [6], [7]],
    hinge_edges1=[(-1,0), (0,4)],
    roots1=[0],
    tree_ids1=[0] * 2),
  test_case_data(
    art=r"""
           7
           | \
  5        1   \
   \     /   \   \
    4 - 0     2 - 8
   /     \   /   /
  6        3   /
           | /
           9
""",
    n_vertices=10,
    edge_list=[
      (0,1), (0,3), (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6),
      (7,8),
      (8,9)],
    clusters1=[[0, 1, 2, 3, 7, 8, 9], [4], [5], [6]],
    hinge_edges1=[(-1,0), (0,4)],
    roots1=[0],
    tree_ids1=[0] * 2),
  test_case_data(
    art=r"""
           7
           |
  5        1
  |\     /   \
  | 4 - 0     2 - 8
  |/     \   /
  6        3
           |
           9
""",
    n_vertices=10,
    edge_list=[
      (0,1), (0,3), (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6),
      (5,6)],
    clusters1=[[0, 1, 2, 3], [4, 5, 6], [7], [8], [9]],
    hinge_edges1=[(-1,0), (0,4)],
    roots1=[0],
    tree_ids1=[0] * 2),
  test_case_data(
    art=r"""
0   1
""",
    n_vertices=2,
    edge_list=[],
    clusters1=[[0], [1]],
    hinge_edges1=[(-1,0), (-1,1)],
    roots1=[0, 1],
    tree_ids1=[0, 1]),
  test_case_data(
    art=r"""
0 - 2   1
""",
    n_vertices=3,
    edge_list=[(0,2)],
    clusters1=[[0], [1], [2]],
    hinge_edges1=[(-1,0), (-1,1)],
    roots1=[0, 2],
    tree_ids1=[0, 1]),
  test_case_data(
    art=r"""
           7
           | \
  5        1   \
  |\         \   \
  | 4 - 0     2 - 8
  |/         /
  6        3
           |
           9
""",
    n_vertices=10,
    edge_list=[
      (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6),
      (5,6),
      (7,8)],
    clusters1=[[1, 2, 7, 8], [4, 5, 6], [0], [3], [9]],
    hinge_edges1=[(-1,1), (2,3), (-1,4)],
    roots1=[0, 3],
    tree_ids1=[0, 0, 1]),
  test_case_data(
    art=r"""
           7
           | \
  5        1   \
             \   \
    4 - 0     2 - 8
   /     \   /
  6 ------ 3
           |
           9
""",
    n_vertices=10,
    edge_list=[
      (0,3), (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,6), (3,9),
      (4,6),
      (7,8)],
    clusters1=[[0, 3, 4, 6], [1, 2, 7, 8], [5], [9]],
    hinge_edges1=[(-1,0), (3,2), (-1,5)],
    roots1=[0, 3],
    tree_ids1=[0, 0, 1]),
  test_case_data(
    art=r"""
    3 - 6 ----- 4 - 0
   /     \     /     \
  8       1   5       9
    \   /       \   /
      2 --------  7
""",
    n_vertices=10,
    edge_list=[
      (0,4), (0,9),
      (1,2), (1,6),
      (2,7), (2,8),
      (3,6), (3,8),
      (4,5), (4,6),
      (5,7),
      (7,9)],
    clusters1=[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]],
    hinge_edges1=[(-1,0)],
    roots1=[0],
    tree_ids1=[0],
    clusters1_5=[[0, 4, 5, 7, 9], [1, 2, 3, 6, 8]],
    hinge_edges1_5=[(-1,0), (4,6)],
    roots1_5=[0],
    tree_ids1_5=[0, 0],
    loop_edges1_5=[(7,2)],
    loop_edge_bendings1_5=[(1,7), (2,5), (2,9), (7,8)],
    clusters2_5=[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]],
    hinge_edges2_5=[(-1,0)],
    loop_edges2_5=[],
    loop_edge_bendings2_5=[]),
  ]

def run(args):
  assert args in [[], ["--verbose"]]
  verbose = "--verbose" in args
  if (verbose):
    out = sys.stdout
  else:
    out = StringIO()
  exercise_cluster_manager()
  for tc in test_cases:
    print >> out, tc.art
    for rigid_loop_size_max in [8, 5]:
      print >> out, "rigid_loop_size_max:", rigid_loop_size_max
      if (rigid_loop_size_max == 5 and tc.clusters1_5 is not None):
       tc_c1, tc_he1, tc_r1, tc_tid1, tc_le1, tc_leb1, \
       tc_c2, tc_he2, tc_le2, tc_leb2 = \
         tc.clusters1_5, tc.hinge_edges1_5, \
         tc.roots1_5, tc.tree_ids1_5, \
         tc.loop_edges1_5, tc.loop_edge_bendings1_5, \
         tc.clusters2_5, tc.hinge_edges2_5, \
         tc.loop_edges2_5, tc.loop_edge_bendings2_5
      else:
       tc_c1, tc_he1, tc_r1, tc_tid1, tc_le1, tc_leb1, \
       tc_c2, tc_he2, tc_le2, tc_leb2 = \
         tc.clusters1, tc.hinge_edges1, \
         tc.roots1, tc.tree_ids1, [], [], \
         tc.clusters1, tc.hinge_edges1, [], []
      def assert_same(label, have, expected):
        print >> out, label, have
        if (expected is not None):
          if (have != expected):
            print >> out, "expected:", expected
          assert have == expected, "Note: --verbose for details"
      #
      tt = construct(
        n_vertices=tc.n_vertices,
        edge_list=tc.edge_list,
        rigid_loop_size_max=rigid_loop_size_max)
      cm = tt.cluster_manager
      assert_same("c1:", cm.clusters, tc_c1)
      cm.construct_spanning_trees(edge_sets=tt.edge_sets)
      print >> out, "c1t:", cm.clusters
      assert_same("he1:", cm.hinge_edges, tc_he1)
      assert_same("le1:", cm.loop_edges, tc_le1)
      tid = cm.tree_ids()
      assert_same("tid1:", tid, tc_tid1)
      cm.find_loop_edge_bendings(edge_sets=tt.edge_sets)
      assert_same("leb1:", cm.loop_edge_bendings, tc_leb1)
      #
      tt = construct(
        n_vertices=tc.n_vertices,
        edge_list=tc.edge_list,
        rigid_loop_size_max=rigid_loop_size_max)
      cm = tt.cluster_manager
      cm.merge_clusters_with_multiple_connections(edge_sets=tt.edge_sets)
      assert_same("c2:", cm.clusters, tc_c2)
      cm.construct_spanning_trees(edge_sets=tt.edge_sets)
      print >> out, "c2t:", cm.clusters
      assert_same("he2:", cm.hinge_edges, tc_he2)
      assert_same("le2:", cm.loop_edges, tc_le2)
      cm.find_loop_edge_bendings(edge_sets=tt.edge_sets)
      assert_same("leb2:", cm.loop_edge_bendings, tc_leb2)
      #
      print >> out
  #
  from scitbx.graph import tst_tardy_pdb
  for tc in tst_tardy_pdb.test_cases:
    tt = tc.tardy_tree_construct()
    tx = construct(
      n_vertices=tt.n_vertices,
      edge_list=tt.edge_list,
      rigid_loop_size_max=6)
    tx.find_cluster_loops()
    tx.finalize()
  #
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
