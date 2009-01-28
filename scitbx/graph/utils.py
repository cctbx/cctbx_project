def construct_edge_sets(n_vertices, edge_list):
  result = [set() for i in xrange(n_vertices)]
  for i,j in edge_list:
    assert i < j
    result[i].add(j)
    result[j].add(i)
  return result

def bond_bending_edge_sets(edge_sets):
  result = []
  for edge_set in edge_sets:
    result.append(set(edge_set))
  for i,edge_set in enumerate(edge_sets):
    for j in edge_set:
      if (j < i): continue
      for k in edge_sets[j]:
        if (k == i): continue
        result[i].add(k)
      for k in edge_sets[i]:
        if (k == j): continue
        result[j].add(k)
  return result

def extract_edge_list(edge_sets):
  result = []
  for i,edge_set in enumerate(edge_sets):
    for j in sorted(edge_set):
      if (j < i): continue
      result.append((i,j))
  return result
