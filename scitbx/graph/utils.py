def construct_edge_sets(n_vertices, edge_list):
  result = [set() for i in xrange(n_vertices)]
  for i,j in edge_list:
    assert i < j
    result[i].add(j)
    result[j].add(i)
  return result

def extract_edge_list(edge_sets):
  result = []
  for i,edge_set in enumerate(edge_sets):
    for j in sorted(edge_set):
      if (j < i): continue
      result.append((i,j))
  return result

def bond_bending_edge_sets(edge_sets):
  result = [set(edge_set) for edge_set in edge_sets]
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

def potential_implied_one_way_edge_sets(edge_sets, bond_bending_edge_sets):
  result = [set() for i in xrange(len(edge_sets))]
  for i,edge_set in enumerate(edge_sets):
    for j in edge_set:
      if (j == i): continue
      for k in edge_sets[j]:
        if (k == j): continue
        if (k == i): continue
        for l in edge_sets[k]:
          if (l == i): continue
          if (l in bond_bending_edge_sets[i]): continue
          if (i < l): result[i].add(l)
          else:       result[l].add(i)
  return result

def potential_implied_edge_list(edge_sets, bond_bending_edge_sets):
  return extract_edge_list(
    edge_sets=potential_implied_one_way_edge_sets(
      edge_sets=edge_sets, bond_bending_edge_sets=bond_bending_edge_sets))
