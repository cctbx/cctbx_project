def exercise_specific():
  from libtbx.topological_sort import stable
  from libtbx.topological_sort import strongly_connected_components as scc
  connections = [
    ("a", ["b", "d"]),
    ("b", []),
    ("c", []),
    ("d", ["c"])]
  node_list = stable(connections=connections)
  assert node_list == ["b", "c", "d", "a"]
  assert scc(dict(connections)) == []
  assert scc(dict(connections), omit_single_node_components=False) \
      == [("b",), ("c",), ("d",), ("a",)]
  connections = [
    ("a", ["d", "b"]),
    ("b", []),
    ("c", []),
    ("d", ["c"])]
  node_list = stable(connections=connections)
  assert node_list == ["b", "c", "d", "a"]
  assert scc(dict(connections)) == []
  connections = [
    (0, [1]),
    (1, [2]),
    (2, [1,3]),
    (3, [3])]
  node_list = stable(connections=connections)
  assert node_list == [3,2,1,0]
  assert scc(dict(connections)) == [(1, 2)]
  connections = [
    ("a", ["d", "b", "e"]),
    ("b", []),
    ("c", []),
    ("d", ["c"])]
  node_list = stable(connections=connections)
  assert node_list == ["b", "c", "d", "e", "a"]
  assert scc(dict(connections)) == []
  #
  assert len(scc(
    successors_by_node={
      "node1": ["successor1", "successor2"],
      "node2": ["successor1", "successor3"]},
    omit_single_node_components=False)) == 5

def exercise_random(rng, n_nodes):
  # meant to discover endless loops or similar bugs
  connections = []
  for i_node in xrange(n_nodes):
    if (rng.randrange(10) > 7):
      continue
    n_del = max(int(n_nodes*0.6), rng.randrange(n_nodes))
    deps = range(n_nodes)
    for i_del in xrange(n_del):
      i = rng.randrange(len(deps))
      del deps[i]
    connections.append((i_node, deps))
  from libtbx.topological_sort import stable
  stable(connections=connections)
  #
  from libtbx.topological_sort import strongly_connected_components as scc
  from libtbx.topological_sort import find_path
  sbn = dict(connections)
  components = scc(successors_by_node=sbn)
  for component in components:
    for a in component:
      for b in component:
        path = find_path(successors_by_node=sbn, from_node=a, to_node=b)
        assert path is not None

def run(args):
  assert len(args) == 0
  exercise_specific()
  import random
  random.seed(0)
  for i_trial in xrange(10):
    exercise_random(rng=random, n_nodes=10)
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
