"""
ext.girvan_newman(n_vertices, edges, threshold): the edges left by boost's
betweenness_centrality_clustering (bc_clustering_threshold(threshold, g,
normalize=False)) on an undirected setS/vecS adjacency list built from the
edges in order: duplicates, reversed pairs and self-pairs ignored; while the
largest (halved, undirected) edge betweenness is >= threshold, the first
edge in insertion order with that value is removed.
"""
import os, random
import libtbx.load_env
from scitbx.array_family import flex
from mmtbx.pair_interaction import pair_interaction
from mmtbx.pair_interaction import ext
from boost_adaptbx import graph
from boost_adaptbx.graph import clustering_algorithm

def reference(n_vertices, edges, threshold):
  """Remaining edges (sorted pairs) from the boost algorithm."""
  g = graph.adjacency_list(graph_type="undirected", vertex_type="vector",
    edge_type="set")
  vertices = [g.add_vertex() for i in range(n_vertices)]
  for i, j in edges:
    g.add_edge(vertex1=vertices[i], vertex2=vertices[j], weight=1)
  clustering_algorithm.betweenness_centrality_clustering(graph=g,
    threshold=threshold)
  index = dict((v, k) for k, v in enumerate(vertices))
  result = []
  for e in g.edges():
    a, b = index[g.source(edge=e)], index[g.target(edge=e)]
    result.append((min(a, b), max(a, b)))
  return sorted(result)

def as_flex(edges):
  result = flex.int()
  for i, j in edges:
    result.append(i); result.append(j)
  return result

def from_flex(a):
  return [(a[2*k], a[2*k+1]) for k in range(a.size()//2)]

def check(n_vertices, edges, threshold):
  expected = reference(n_vertices, edges, threshold)
  got = from_flex(ext.girvan_newman(n_vertices=n_vertices,
    edges=as_flex(edges), threshold=threshold))
  assert sorted((min(i,j), max(i,j)) for i, j in got) == expected, (
    n_vertices, edges, threshold, got, expected)
  # input order kept (the first occurrence of a pair, in either orientation)
  order = {}
  for k, (i, j) in enumerate(edges):
    order.setdefault((min(i, j), max(i, j)), k)
  keys = [order[(min(i, j), max(i, j))] for i, j in got]
  assert keys == sorted(keys), keys
  return got

def exercise_random():
  random.seed(7)
  for trial in range(30):
    n = random.randint(8, 60)
    edges = []
    for v in range(n):
      for k in range(random.randint(1, 3)):
        w = random.randint(0, n-1)
        if(w != v): edges.append((v, w))
    for threshold in (2, 4, 9):
      check(n, edges, threshold)

def exercise_degenerate():
  """Duplicates, reversed pairs, self-pairs, isolated vertices, no edges."""
  edges = [(0,1), (1,0), (0,1), (2,2), (1,2), (2,3), (3,4), (4,1), (5,6)]
  for threshold in (1, 2, 4):
    check(9, edges, threshold)
  assert ext.girvan_newman(n_vertices=3, edges=flex.int(), threshold=2).size() == 0
  got = from_flex(ext.girvan_newman(n_vertices=4, edges=as_flex([(1,1),(2,2)]),
    threshold=1))
  assert got == [(1,1), (2,2)], got # self-pairs stay (boost keeps them)

def exercise_chain():
  """A path of 20 vertices at thresholds 9..4, chained as clustering.py does."""
  edges = [(i, i+1) for i in range(19)]
  remaining = edges
  for threshold in (9, 8, 7, 6, 5, 4):
    remaining = check(20, remaining, threshold)

def exercise_1yjp():
  """The interaction graph of 1yjp (qrefine data), thresholds 9..4 chained."""
  d = libtbx.env.find_in_repositories("qrefine/tests/unit/data_files")
  if(d is None): return
  import iotbx.pdb
  ph = iotbx.pdb.input(file_name=os.path.join(d, "1yjp.pdb")).construct_hierarchy()
  n = len(list(ph.residue_groups()))
  il = pair_interaction.run(ph)
  edges = [(i-1, j-1) for i, j in il]
  remaining = edges
  for threshold in (9, 8, 7, 6, 5, 4):
    remaining = check(n, remaining, threshold)

def run():
  exercise_random()
  exercise_degenerate()
  exercise_chain()
  exercise_1yjp()

if(__name__=="__main__"):
  run()
  print("OK")
