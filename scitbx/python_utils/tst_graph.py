from __future__ import absolute_import, division, print_function
from scitbx.python_utils import graph_tools as gt

def tst_graph():
  ll = gt.graph()
  ll.insert_node( name = 'a',
                  node_object = 'ploep_a',
                  edge_object = { 'b': 'a_to_b', 'c': 'a_to_c' } )

  ll.insert_node( name = 'b',
                  node_object = 'ploep_b',
                  edge_object = { 'c': 'b_to_c' } )

  ll.insert_node( name = 'c',
                  node_object = 'ploep_c',
                  edge_object = {  } )




  path_found = ll.find_all_paths('a', 'c')
  path_possible =  [['a', 'c'], ['a', 'b', 'c']]
  for path in path_possible:
    assert (path in path_found)

  shortest = ll.find_shortest_path('a', 'c')
  assert (shortest == ['a','c'])

  assert ( (ll.o == { 'a':['b','c'], 'b':['c'],'c':[] }) or
           (ll.o == { 'a':['c','b'], 'b':['c'],'c':[] }) )

  kk = gt.graph()
  kk.insert_node( name = 'a',
                  node_object = 'ploep_a',
                  edge_object = { 'b': 'a_to_b', 'c': 'a_to_c' } )

  kk.insert_node( name = 'b',
                  node_object = 'ploep_b',
                  edge_object = { 'c': 'b_to_c' } )

  kk.insert_node( name = 'c',
                  node_object = 'ploep_c',
                  edge_object = {  } )

  assert ll.is_equivalent_to( kk )

  kk.insert_node( name = 'b',
                  node_object=None,
                  edge_object={ 'a': 'b_to_a'} )

  assert ll.is_contained_in( kk )
  assert not kk.is_contained_in( ll )

  ll.assert_is_clean()
  kk.assert_is_clean()

  kk.remove_node( 'b' )
  kk.assert_is_clean()

  kk.insert_node( name = 'a',
                  node_object = 'ploep_a',
                  edge_object = { 'b': 'a_to_b', 'c': 'a_to_c' } )

  kk.insert_node( name = 'b',
                  node_object = 'ploep_b',
                  edge_object = { 'c': 'b_to_c' } )

  kk.insert_node( name = 'd',
                  node_object = 'ploep_a',
                  edge_object = { 'b': 'd_to_b', 'c': 'd_to_c' } )

  kk.insert_node( name = 'e',
                  node_object = 'ploep_b',
                  edge_object = { 'd': 'e_to_d' } )

  kk.insert_node( name = 'c',
                  node_object = 'ploep_b',
                  edge_object = { 'd': 'c_to_d' } )

  a_paths = [['a', 'c', 'd'], ['a', 'b', 'c']]
  c_paths = [['c', 'd', 'b']]
  b_paths = [['b', 'c', 'd']]
  e_paths = [['e', 'd', 'c'], ['e', 'd', 'b']]
  d_paths = [['d', 'b', 'c']]
  for ps, solution in [('a', a_paths), ('b', b_paths), ('c', c_paths),
                       ('d', d_paths), ('e', e_paths)]:
    for p in kk.find_paths_of_length(ps, 3):
      assert p in solution



def run():
  tst_graph()


  print('OK')

if (__name__ == "__main__"):
  run()
