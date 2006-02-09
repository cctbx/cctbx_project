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


def run():
  tst_graph()


  print 'OK'

if (__name__ == "__main__"):
  run()
