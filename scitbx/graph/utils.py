from libtbx import slots_getstate_setstate

def construct_edge_sets(n_vertices, edge_list, assert_i_lt_j=True):
  result = [set() for i in xrange(n_vertices)]
  for i,j in edge_list:
    if (assert_i_lt_j):
      assert i < j
    else:
      assert i != j
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

def bond_bending_edge_sets(edge_sets, omit_bonds=False):
  if (omit_bonds):
    result = [set() for i in xrange(len(edge_sets))]
  else:
    result = [set(edge_set) for edge_set in edge_sets]
  for i,edge_set in enumerate(edge_sets):
    for j in edge_set:
      if (j < i): continue
      for k in edge_sets[j]:
        if (k == i): continue
        if (omit_bonds and k in edge_sets[i]): continue
        result[i].add(k)
      for k in edge_sets[i]:
        if (k == j): continue
        if (omit_bonds and k in edge_sets[j]): continue
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

class sub_edge_list(object):

  def __init__(O, edge_sets, vertex_indices):
    O.vertex_indices = vertex_indices
    O.reindexing_dict = {}
    for i_sub,i in enumerate(vertex_indices):
      O.reindexing_dict[i] = i_sub
    assert len(O.reindexing_dict) == len(vertex_indices)
    O.edge_list = []
    ea = O.edge_list.append
    for i_sub,i in enumerate(vertex_indices):
      for j in edge_sets[i]:
        if (i > j): continue
        assert i != j
        j_sub = O.reindexing_dict.get(j)
        if (j_sub is None): continue
        if (i_sub < j_sub): ea((i_sub, j_sub))
        else:               ea((j_sub, i_sub))

  def edge_sets(O):
    return construct_edge_sets(
      n_vertices=len(O.vertex_indices),
      edge_list=O.edge_list)

class tree_marking(slots_getstate_setstate):

  __slots__ = ["edge_sets", "n_trees", "indices"]

  def __init__(O, edge_sets):
    n_vertices = len(edge_sets)
    indices = [n_vertices] * n_vertices
    i_tree = 0
    for i_root in xrange(n_vertices):
      if (indices[i_root] == n_vertices):
        follow = [i_root]
        while (len(follow) != 0):
          i = follow.pop()
          if (indices[i] == n_vertices):
            indices[i] = i_tree
            follow.extend(edge_sets[i])
        i_tree += 1
    O.edge_sets = edge_sets
    O.n_trees = i_tree
    O.indices = indices

  def partitions_of(O, vertex_indices):
    result = []
    ti = O.indices
    nt = O.n_trees
    ri = [nt] * nt
    for i in vertex_indices:
      tii = ti[i]
      rii = ri[tii]
      if (rii == nt):
        ri[tii] = len(result)
        result.append([i])
      else:
        result[rii].append(i)
    return result
