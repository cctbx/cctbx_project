from scitbx.graph import utils
import sys

def run(args):
  assert len(args) == 0
  #
  el = []
  es = utils.construct_edge_sets(n_vertices=0, edge_list=el)
  assert len(es) == 0
  bbes = utils.bond_bending_edge_sets(edge_sets=es)
  assert len(bbes) == 0
  eel = utils.extract_edge_list(edge_sets=bbes)
  assert len(eel) == 0
  #
  el = [(0,1), (1,2), (2,3)]
  es = utils.construct_edge_sets(n_vertices=4, edge_list=el)
  assert len(es) == 4
  bbes = utils.bond_bending_edge_sets(edge_sets=es)
  assert len(bbes) == 4
  eel = utils.extract_edge_list(edge_sets=bbes)
  assert eel == [(0,1), (0,2), (1,2), (1,3), (2,3)]
  #
  el = [(0,1), (1,2), (2,3), (0,3)]
  es = utils.construct_edge_sets(n_vertices=4, edge_list=el)
  bbes = utils.bond_bending_edge_sets(edge_sets=es)
  eel = utils.extract_edge_list(edge_sets=bbes)
  assert eel == [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
  #
  el = [(0,1), (1,2), (2,3), (3,4), (4,5), (0,5)]
  es = utils.construct_edge_sets(n_vertices=6, edge_list=el)
  bbes = utils.bond_bending_edge_sets(edge_sets=es)
  eel = utils.extract_edge_list(edge_sets=bbes)
  assert eel == [
    (0,1), (0,2), (0,4), (0,5),
    (1,2), (1,3), (1,5),
    (2,3), (2,4),
    (3,4), (3,5),
    (4,5)]
  #
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
