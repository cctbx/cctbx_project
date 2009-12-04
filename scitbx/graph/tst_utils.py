from scitbx.graph import utils
import sys

def exercise_tree_marking():
  tree_marking = utils.tree_marking
  el = []
  for n_vertices in xrange(5):
    es = utils.construct_edge_sets(n_vertices=n_vertices, edge_list=el)
    ti = tree_marking(edge_sets=es).indices
    assert ti == range(n_vertices)
  el = [(0,1)]
  for n_vertices in xrange(2,5):
    es = utils.construct_edge_sets(n_vertices=n_vertices, edge_list=el)
    ti = tree_marking(edge_sets=es).indices
    assert ti == [0] + range(n_vertices-1)
  el = [(2,5)]
  es = utils.construct_edge_sets(n_vertices=7, edge_list=el)
  ti = tree_marking(edge_sets=es).indices
  assert ti == [0,1,2,3,4,2,5]
  for n_vertices in xrange(7):
    el = [(i,i+1) for i in xrange(n_vertices-1)]
    es = utils.construct_edge_sets(n_vertices=n_vertices, edge_list=el)
    ti = tree_marking(edge_sets=es).indices
    assert ti == [0] * n_vertices
  el = [(0,1),(2,3)]
  es = utils.construct_edge_sets(n_vertices=4, edge_list=el)
  ti = tree_marking(edge_sets=es).indices
  assert ti == [0,0,1,1]
  el = [(1,2),(0,3)]
  es = utils.construct_edge_sets(n_vertices=4, edge_list=el)
  ti = tree_marking(edge_sets=es).indices
  assert ti == [0,1,1,0]
  el = [
    (8,9),(7,9),(3,7),(8,11),(8,12),(12,14),(13,15),(7,13),
    (1,6),(4,6),
    (0,5)]
  es = utils.construct_edge_sets(n_vertices=16, edge_list=el)
  tm = tree_marking(edge_sets=es)
  assert tm.indices == [0,1,2,3,1,0,1,3,3,3,4,3,3,3,3,3]
  assert tm.partitions_of(vertex_indices=()) == []
  for i in xrange(16):
    assert tm.partitions_of(vertex_indices=[i]) == [[i]]
  assert tm.partitions_of(vertex_indices=[14,3]) == [[14,3]]
  assert tm.partitions_of(vertex_indices=[5,10]) == [[5],[10]]
  assert tm.partitions_of(vertex_indices=range(16)) == [
    [0,5], [1,4,6], [2], [3,7,8,9,11,12,13,14,15], [10]]

def run(args):
  assert len(args) == 0
  #
  el = []
  n_vertices = 0
  es = utils.construct_edge_sets(n_vertices=n_vertices, edge_list=el)
  assert len(es) == 0
  bbes = utils.bond_bending_edge_sets(edge_sets=es)
  assert len(bbes) == 0
  bbel = utils.extract_edge_list(edge_sets=bbes)
  assert len(bbel) == 0
  piel = utils.potential_implied_edge_list(
    edge_sets=es, bond_bending_edge_sets=bbes)
  assert len(piel) == 0
  bbesob = utils.bond_bending_edge_sets(edge_sets=es, omit_bonds=True)
  assert len(bbesob) == 0
  #
  el = [(0,1), (1,2), (2,3)]
  n_vertices = 4
  es = utils.construct_edge_sets(n_vertices=n_vertices, edge_list=el)
  assert len(es) == n_vertices
  bbes = utils.bond_bending_edge_sets(edge_sets=es)
  assert len(bbes) == n_vertices
  bbel = utils.extract_edge_list(edge_sets=bbes)
  assert bbel == [(0,1), (0,2), (1,2), (1,3), (2,3)]
  assert set(el).issubset(set(bbel))
  piel = utils.potential_implied_edge_list(
    edge_sets=es, bond_bending_edge_sets=bbes)
  assert piel == [(0,3)]
  assert set(bbel).isdisjoint(set(piel))
  bbesob = utils.bond_bending_edge_sets(edge_sets=es, omit_bonds=True)
  assert len(bbesob) == n_vertices
  bbelob = utils.extract_edge_list(edge_sets=bbesob)
  assert bbelob == [(0,2), (1,3)]
  #
  el = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7)]
  n_vertices = 8
  es = utils.construct_edge_sets(n_vertices=n_vertices, edge_list=el)
  bbes = utils.bond_bending_edge_sets(edge_sets=es)
  bbel = utils.extract_edge_list(edge_sets=bbes)
  assert bbel == [
    (0, 1), (0, 2),
    (1, 2), (1, 3),
    (2, 3), (2, 4),
    (3, 4), (3, 5),
    (4, 5), (4, 6),
    (5, 6), (5, 7),
    (6, 7)]
  assert set(el).issubset(set(bbel))
  piel = utils.potential_implied_edge_list(
    edge_sets=es, bond_bending_edge_sets=bbes)
  assert piel == [(0,3), (1,4), (2,5), (3,6), (4,7)]
  assert set(bbel).isdisjoint(set(piel))
  bbesob = utils.bond_bending_edge_sets(edge_sets=es, omit_bonds=True)
  bbelob = utils.extract_edge_list(edge_sets=bbesob)
  assert bbelob == [
    (0, 2),
    (1, 3),
    (2, 4),
    (3, 5),
    (4, 6),
    (5, 7)]
  #
  el = [(0,1), (1,2), (2,3), (0,3)]
  n_vertices = 4
  es = utils.construct_edge_sets(n_vertices=n_vertices, edge_list=el)
  bbes = utils.bond_bending_edge_sets(edge_sets=es)
  bbel = utils.extract_edge_list(edge_sets=bbes)
  assert bbel == [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
  assert set(el).issubset(set(bbel))
  piel = utils.potential_implied_edge_list(
    edge_sets=es, bond_bending_edge_sets=bbes)
  assert piel == []
  bbesob = utils.bond_bending_edge_sets(edge_sets=es, omit_bonds=True)
  bbelob = utils.extract_edge_list(edge_sets=bbesob)
  assert bbelob == [(0, 2), (1, 3)]
  #
  el = [(0,1), (1,2), (2,3), (3,4), (4,5), (0,5)]
  n_vertices = 6
  es = utils.construct_edge_sets(n_vertices=n_vertices, edge_list=el)
  bbes = utils.bond_bending_edge_sets(edge_sets=es)
  bbel = utils.extract_edge_list(edge_sets=bbes)
  assert bbel == [
    (0,1), (0,2), (0,4), (0,5),
    (1,2), (1,3), (1,5),
    (2,3), (2,4),
    (3,4), (3,5),
    (4,5)]
  assert set(el).issubset(set(bbel))
  piel = utils.potential_implied_edge_list(
    edge_sets=es, bond_bending_edge_sets=bbes)
  assert piel == [(0,3), (1,4), (2,5)]
  assert set(bbel).isdisjoint(set(piel))
  bbesob = utils.bond_bending_edge_sets(edge_sets=es, omit_bonds=True)
  bbelob = utils.extract_edge_list(edge_sets=bbesob)
  assert bbelob == [
    (0,2), (0,4),
    (1,3), (1,5),
    (2,4),
    (3,5)]
  #
  el = [(0,1), (1,2), (2,3), (2,6), (3,4), (3,7), (4,5), (0,5)]
  n_vertices = 8
  es = utils.construct_edge_sets(n_vertices=n_vertices, edge_list=el)
  bbes = utils.bond_bending_edge_sets(edge_sets=es)
  bbel = utils.extract_edge_list(edge_sets=bbes)
  assert bbel == [
    (0,1), (0,2), (0,4), (0,5),
    (1,2), (1,3), (1,5), (1,6),
    (2,3), (2,4), (2,6), (2,7),
    (3,4), (3,5), (3,6), (3,7),
    (4,5), (4,7)]
  assert set(el).issubset(set(bbel))
  piel = utils.potential_implied_edge_list(
    edge_sets=es, bond_bending_edge_sets=bbes)
  assert piel == [(0,3), (0,6), (1,4), (1,7), (2,5), (4,6), (5,7), (6,7)]
  assert set(bbel).isdisjoint(set(piel))
  bbesob = utils.bond_bending_edge_sets(edge_sets=es, omit_bonds=True)
  bbelob = utils.extract_edge_list(edge_sets=bbesob)
  assert bbelob == [
    (0,2), (0,4),
    (1,3), (1,5), (1,6),
    (2,4), (2,7),
    (3,5), (3,6),
    (4,7)]
  #
  el = [(0,1), (1,2), (2,3), (2,6), (3,4), (3,7), (4,5), (0,5)]
  n_vertices = 8
  es = utils.construct_edge_sets(n_vertices=n_vertices, edge_list=el)
  sub = utils.sub_edge_list(edge_sets=es, vertex_indices=[])
  assert len(sub.edge_list) == 0
  assert len(sub.edge_sets()) == 0
  sub = utils.sub_edge_list(edge_sets=es, vertex_indices=[1])
  assert len(sub.edge_list) == 0
  assert len(sub.edge_sets()) == 1
  sub = utils.sub_edge_list(edge_sets=es, vertex_indices=[1,3])
  assert len(sub.edge_list) == 0
  assert len(sub.edge_sets()) == 2
  sub = utils.sub_edge_list(edge_sets=es, vertex_indices=[1,0])
  assert sub.edge_list == [(0,1)]
  assert len(sub.edge_sets()) == 2
  assert sub.reindexing_dict == {0: 1, 1: 0}
  sub = utils.sub_edge_list(edge_sets=es, vertex_indices=[6,7,3,2])
  assert sub.edge_list == [(1,2), (2,3), (0,3)]
  assert len(sub.edge_sets()) == 4
  assert sub.reindexing_dict == {2: 3, 3: 2, 6: 0, 7: 1}
  #
  exercise_tree_marking()
  #
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
