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

def strongly_connected_components(
      successors_by_node,
      omit_single_node_components=True,
      low_infinite=2**30):
  """
successors_by_node = {
  "node1": ["successor1", "successor2"],
  "node2": ["successor1", "successor3"]
}

http://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
http://www.logarithmic.net/pfh-files/blog/01208083168/sort.py

Original implementation (by Paul Harrison), modified to accommodate
successors that do not appear as a key in successors_by_node.
  """
  result = []
  stack = []
  low = {}
  def visit(node):
    if (node in low):
      return
    num = len(low)
    low[node] = num
    stack_pos = len(stack)
    stack.append(node)
    for successor in successors_by_node.get(node, []):
      visit(successor)
      low[node] = min(low[node], low[successor])
    if (num == low[node]):
      component = tuple(stack[stack_pos:])
      del stack[stack_pos:]
      if (not omit_single_node_components or len(component) != 1):
        result.append(component)
      for item in component:
        low[item] = low_infinite
  for node in successors_by_node:
    visit(node)
  return result

def find_path(successors_by_node, from_node, to_node):
  visited = set()
  path = []
  def depth_first_search(node):
    visited.add(node)
    for successor in successors_by_node.get(node, []):
      if (successor == to_node):
        return True
      if (successor not in visited):
        path.append(successor)
        if (depth_first_search(successor)):
          return True
        path.pop()
    return False
  if (depth_first_search(from_node)):
    return path
  return None
