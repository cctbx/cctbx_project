from scitbx.graph.tardy_tree import cluster_manager, construct
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
    cm.connect(i=0, j=0)
  for p in xrange(2):
    cm.connect(i=1, j=3)
    assert cm.cluster_indices == [0,1,2,1,4]
    assert cm.clusters == [[0],[1,3],[2],[],[4]]
  for p in xrange(2):
    cm.connect(i=0, j=3)
    assert cm.cluster_indices == [1,1,2,1,4]
    for q in xrange(2):
      assert cm.clusters == [[],[1,3,0],[2],[],[4]]
      cm.refresh_indices()
  cm.connect(i=2, j=4)
  assert cm.clusters == [[],[1,3,0],[2,4],[],[]]
  assert cm.cluster_indices == [1,1,2,1,2]
  cm.connect(i=2, j=3)
  assert cm.clusters == [[],[1,3,0,2,4],[],[],[]]
  assert cm.cluster_indices == [1,1,1,1,1]
  cm.tidy()
  assert cm.clusters == [[0,1,2,3,4]]
  assert cm.cluster_indices == [0,0,0,0,0]
  #
  cm = cluster_manager(n_vertices=6)
  cm.connect(i=3, j=0)
  cm.connect(i=2, j=4)
  cm.connect(i=1, j=2)
  cm.tidy()
  assert cm.clusters == [[1,2,4],[0,3],[5]]
  assert cm.cluster_indices == [1,0,0,1,0,2]
  edges = [(0,1), (0,2), (3,4), (4,5)]
  ces = cm.cluster_edge_sets(edges=edges)
  assert [sorted(e) for e in ces] == [[1, 2], [0], [0]]
  cm.merge_lones(edges=edges)
  assert cm.clusters == [[1,2,4,5],[0,3]]
  assert cm.cluster_indices == [1,0,0,1,0,0]
  ces = cm.cluster_edge_sets(edges=edges)
  assert [sorted(e) for e in ces] == [[1], [0]]
  cm.construct_spanning_trees(edges=edges)
  assert cm.parents == [-1,0]
  assert cm.roots() == [0]
  assert cm.tree_ids() == [0,0]
  assert cm.find_loop_edges(edges=edges) == []

class test_case_data(object):

  def __init__(O,
        art,
        n_vertices,
        edges,
        clusters1,
        parents1,
        roots1,
        tree_ids1,
        clusters2,
        parents2,
        clusters1_5=None,
        parents1_5=None,
        roots1_5=None,
        tree_ids1_5=None,
        loop_edges1_5=[],
        clusters2_5=None,
        parents2_5=None,
        loop_edges2_5=[]):
    assert art[0] == "\n"
    assert (tree_ids1_5 is None) == (parents1_5 is None)
    assert (roots1_5 is None) == (parents1_5 is None)
    O.art = art[1:]
    O.n_vertices = n_vertices
    O.edges = edges
    O.clusters1 = clusters1
    O.parents1 = parents1
    O.roots1 = roots1
    O.tree_ids1 = tree_ids1
    O.clusters2 = clusters2
    O.parents2 = parents2
    O.clusters1_5 = clusters1_5
    O.parents1_5 = parents1_5
    O.roots1_5 = roots1_5
    O.tree_ids1_5 = tree_ids1_5
    O.loop_edges1_5 = loop_edges1_5
    O.clusters2_5 = clusters2_5
    O.parents2_5 = parents2_5
    O.loop_edges2_5 = loop_edges2_5

test_cases = [
  test_case_data(
    art=r"""
0
""",
    n_vertices=1,
    edges=[],
    clusters1=[[0]],
    parents1=[-1],
    roots1=[0],
    tree_ids1=[0],
    clusters2=[[0]],
    parents2=[-1]),
  test_case_data(
    art=r"""
0 - 1
""",
    n_vertices=2,
    edges=[(0,1)],
    clusters1=[[0], [1]],
    parents1=[-1, 0],
    roots1=[0],
    tree_ids1=[0, 0],
    clusters2=[[0, 1]],
    parents2=[-1]),
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
    edges=[
      (0,1), (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6)],
    clusters1=[[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]],
    parents1=[-1, 0, 1, 2, 0, 4, 4, 1, 2, 3],
    roots1=[0],
    tree_ids1=[0] * 10,
    clusters2=[[4, 5, 6], [1, 7], [2, 8], [3, 9], [0]],
    parents2=[-1, 0, 1, 2, 3]),
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
    edges=[
      (0,1), (0,4),
      (1,2), (1,7),
      (2,3), (2,8), (2,9),
      (3,9),
      (4,5), (4,6)],
    clusters1=[[2,3,9], [0], [1], [4], [5], [6], [7], [8]],
    parents1=[-1, 0, 1, 2, 3, 3, 1, 0],
    roots1=[0],
    tree_ids1=[0] * 8,
    clusters2=[[2, 3, 8, 9], [4, 5, 6], [1, 7], [0]],
    parents2=[-1, 0, 1, 2]),
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
    edges=[
      (0,1), (0,3), (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6)],
    clusters1=[[0,1,2,3], [4], [5], [6], [7], [8], [9]],
    parents1=[-1, 0, 1, 1, 0, 0, 0],
    roots1=[0],
    tree_ids1=[0] * 7,
    clusters2=[[0, 1, 2, 3, 7, 8, 9], [4, 5, 6]],
    parents2=[-1,0]),
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
    edges=[
      (0,1), (0,3), (0,4),
      (1,7), (1,8),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6)],
    clusters1=[[0, 1, 2, 3, 8], [4], [5], [6], [7], [9]],
    parents1=[-1, 0, 1, 1, 0, 0],
    roots1=[0],
    tree_ids1=[0] * 6,
    clusters2=[[0, 1, 2, 3, 7, 8, 9], [4, 5, 6]],
    parents2=[-1,0]),
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
    edges=[
      (0,1), (0,3), (0,4),
      (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6),
      (7,8)],
    clusters1=[[0, 1, 2, 3, 7, 8], [4], [5], [6], [9]],
    parents1=[-1, 0, 1, 1, 0],
    roots1=[0],
    tree_ids1=[0] * 5,
    clusters2=[[0, 1, 2, 3, 7, 8, 9], [4, 5, 6]],
    parents2=[-1,0],
    clusters1_5=[[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]],
    parents1_5=[-1, 0, 0, 2, 0, 4, 4, 1, 7, 2],
    roots1_5=[0],
    tree_ids1_5=[0] * 10,
    loop_edges1_5=[(2,8)],
    clusters2_5=[[4, 5, 6], [3, 9], [0], [1], [2], [7], [8]],
    parents2_5=[-1, 0, 1, 1, 2, 3, 5],
    loop_edges2_5=[(2,8)]),
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
    edges=[
      (0,1), (0,3), (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6),
      (7,8)],
    clusters1=[[0, 1, 2, 3, 7, 8], [4], [5], [6], [9]],
    parents1=[-1, 0, 1, 1, 0],
    roots1=[0],
    tree_ids1=[0] * 5,
    clusters2=[[0, 1, 2, 3, 7, 8, 9], [4, 5, 6]],
    parents2=[-1,0]),
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
    edges=[
      (0,1), (0,3), (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6),
      (8,9)],
    clusters1=[[0, 1, 2, 3, 8, 9], [4], [5], [6], [7]],
    parents1=[-1, 0, 1, 1, 0],
    roots1=[0],
    tree_ids1=[0] * 5,
    clusters2=[[0, 1, 2, 3, 7, 8, 9], [4, 5, 6]],
    parents2=[-1,0]),
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
    edges=[
      (0,1), (0,3), (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6),
      (7,8),
      (8,9)],
    clusters1=[[0, 1, 2, 3, 7, 8, 9], [4], [5], [6]],
    parents1=[-1, 0, 1, 1],
    roots1=[0],
    tree_ids1=[0] * 4,
    clusters2=[[0, 1, 2, 3, 7, 8, 9], [4, 5, 6]],
    parents2=[-1,0]),
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
    edges=[
      (0,1), (0,3), (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6),
      (5,6)],
    clusters1=[[0, 1, 2, 3], [4, 5, 6], [7], [8], [9]],
    parents1=[-1, 0, 0, 0, 0],
    roots1=[0],
    tree_ids1=[0] * 5,
    clusters2=[[0, 1, 2, 3, 7, 8, 9], [4, 5, 6]],
    parents2=[-1,0]),
  test_case_data(
    art=r"""
0   1
""",
    n_vertices=2,
    edges=[],
    clusters1=[[0], [1]],
    parents1=[-1, -1],
    roots1=[0, 1],
    tree_ids1=[0, 1],
    clusters2=[[0], [1]],
    parents2=[-1, -1]),
  test_case_data(
    art=r"""
0 - 2   1
""",
    n_vertices=3,
    edges=[(0,2)],
    clusters1=[[0], [1], [2]],
    parents1=[-1, 0, -1],
    roots1=[0, 2],
    tree_ids1=[0, 0, 1],
    clusters2=[[0, 2], [1]],
    parents2=[-1, -1]),
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
    edges=[
      (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,9),
      (4,5), (4,6),
      (5,6),
      (7,8)],
    clusters1=[[1, 2, 7, 8], [4, 5, 6], [0], [3], [9]],
    parents1=[-1, 0, 1, -1, 3],
    roots1=[0, 3],
    tree_ids1=[0, 0, 0, 1, 1],
    clusters2=[[0, 4, 5, 6], [1, 2, 7, 8], [3, 9]],
    parents2=[-1, -1, 1]),
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
    edges=[
      (0,3), (0,4),
      (1,2), (1,7),
      (2,3), (2,8),
      (3,6), (3,9),
      (4,6),
      (7,8)],
    clusters1=[[0, 3, 4, 6], [1, 2, 7, 8], [5], [9]],
    parents1=[-1, 0, 0, -1],
    roots1=[0, 3],
    tree_ids1=[0, 0, 0, 1],
    clusters2=[[0, 3, 4, 6, 9], [1, 2, 7, 8], [5]],
    parents2=[-1, 0, -1]),
  ]

def exercise_tyr_with_h():
  # matching PDB file:
  #  scitbx/rigid_body/proto/tst_molecules.py (svn rev. 8373):
  #    simulation_tyr_with_h()
  tt = construct(
    n_vertices=21,
    edges=[
      (0, 1), (0, 2), (0, 12), (1, 3), (1, 6), (2, 4), (2, 7), (3, 5),
      (3, 8), (4, 5), (4, 9), (5, 10), (10, 11), (12, 13), (12, 14),
      (12, 16), (15, 16), (15, 20), (16, 17), (16, 18), (17, 19)])
  tt.cluster_manager.merge_lones(edges=tt.edges)
  tt.cluster_manager.construct_spanning_trees(edges=tt.edges)
  assert tt.cluster_manager.clusters == [
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    [12, 13, 14], [10, 11], [16, 18], [15, 20], [17, 19]]
  assert tt.cluster_manager.parents == [-1, 0, 0, 1, 3, 3]

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
    for loop_size_max in [8, 5]:
      print >> out, "loop_size_max:", loop_size_max
      if (loop_size_max == 5 and tc.clusters1_5 is not None):
       tc_c1, tc_p1, tc_r1, tc_tid1, tc_le1, tc_c2, tc_p2, tc_le2 = \
         tc.clusters1_5, tc.parents1_5, \
         tc.roots1_5, tc.tree_ids1_5, tc.loop_edges1_5, \
         tc.clusters2_5, tc.parents2_5, tc.loop_edges2_5
      else:
       tc_c1, tc_p1, tc_r1, tc_tid1, tc_le1, tc_c2, tc_p2, tc_le2 = \
         tc.clusters1, tc.parents1, \
         tc.roots1, tc.tree_ids1, [], \
         tc.clusters2, tc.parents2, []
      def assert_same(label, have, expected):
        print >> out, label, have
        if (expected is not None):
          if (have != expected):
            print >> out, "expected:", expected
          assert have == expected, "Note: --verbose for details"
      tt = construct(
        n_vertices=tc.n_vertices, edges=tc.edges, size_max=loop_size_max)
      assert_same("c1:", tt.cluster_manager.clusters, tc_c1)
      tt.cluster_manager.construct_spanning_trees(edges=tc.edges)
      print >> out, "c1p:", tt.cluster_manager.clusters
      assert_same("p1:", tt.cluster_manager.parents, tc_p1)
      tid = tt.cluster_manager.tree_ids()
      assert_same("tid1:", tid, tc_tid1)
      le = tt.cluster_manager.find_loop_edges(edges=tc.edges)
      assert_same("le1:", le, tc_le1)
      #
      tt = construct(
        n_vertices=tc.n_vertices, edges=tc.edges, size_max=loop_size_max)
      tt.cluster_manager.merge_lones(edges=tc.edges)
      assert_same("c2:", tt.cluster_manager.clusters, tc_c2)
      tt.cluster_manager.construct_spanning_trees(edges=tc.edges)
      print >> out, "c2p:", tt.cluster_manager.clusters
      assert_same("p2:", tt.cluster_manager.parents, tc_p2)
      le = tt.cluster_manager.find_loop_edges(edges=tc.edges)
      assert_same("le2:", le, tc_le2)
  exercise_tyr_with_h()
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
