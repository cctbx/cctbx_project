def stable(connections):
  ranks = {}
  for node, deps in connections:
    assert node not in ranks
    ranks[node] = len(ranks)
  deps_by_node = {}
  for node, deps in connections:
    deps_by_node[node] = deps
    for d in deps:
      if (d not in ranks):
        ranks[d] = len(ranks)
  lower_bounds = {}
  node_list = []
  def process(dependent_node, node):
    if (node in lower_bounds):
      return
    if (dependent_node is None):
      lower_bounds[node] = len(node_list)
      node_list.append(node)
    else:
      n = len(node_list)
      i = lower_bounds[dependent_node]
      while (i < n):
        if (node_list[i] == dependent_node):
          break
        i += 1
      else:
        raise AssertionError
      lower_bounds[node] = i
      node_list.insert(i, node)
    deps = deps_by_node.get(node)
    if (deps is not None):
      del deps_by_node[node]
      for rank,dependency in sorted([(ranks[d],d) for d in deps]):
        process(dependent_node=node, node=dependency)
  for node, deps in connections:
    process(dependent_node=None, node=node)
  return node_list
