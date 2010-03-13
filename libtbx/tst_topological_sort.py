def exercise_specific():
  from libtbx.topological_sort import stable
  connections = [
    ("a", ["b", "d"]),
    ("b", []),
    ("c", []),
    ("d", ["c"])]
  node_list = stable(connections=connections)
  assert node_list == ["b", "c", "d", "a"]
  connections = [
    ("a", ["d", "b"]),
    ("b", []),
    ("c", []),
    ("d", ["c"])]
  node_list = stable(connections=connections)
  assert node_list == ["b", "c", "d", "a"]
  connections = [
    (0, [1]),
    (1, [2]),
    (2, [1,3]),
    (3, [3])]
  node_list = stable(connections=connections)
  assert node_list == [3,2,1,0]
  connections = [
    ("a", ["d", "b", "e"]),
    ("b", []),
    ("c", []),
    ("d", ["c"])]
  node_list = stable(connections=connections)
  assert node_list == ["b", "c", "d", "e", "a"]

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
