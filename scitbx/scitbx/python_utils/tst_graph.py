from scitbx.python_utils import graph_tools as gt

def tst_graph():
  ref_in_graph  = { 'a':[], 'b':['a','c'],'c':['a'] }
  ref_out_graph = { 'a':['b','c'], 'b':[],'c':['b'] }

  g = gt.graph()

  g.insert_node( name = 'a',
                 in_names = [],
                 out_names = [ 'b', 'c' ]
                 )

  g.insert_node( name = 'c',
                 in_names = [],
                 out_names = [ 'b' ]
                )

  assert ( g.i == ref_in_graph )
  assert ( g.o == ref_out_graph )

  k = gt.graph()

  k.insert_node( name = 'a',
                 in_names = [],
                 out_names = [ 'b', 'c' ]
                 )
  k.insert_node( name = 'b',
                 in_names = ['c'],
                 out_names = [  ]
                )


  assert ( k.i == ref_in_graph )
  assert ( k.o == ref_out_graph )

  g.assert_is_clean()
  k.assert_is_clean()

  assert k.is_contained_in( g )
  assert k.is_equivalent_to( g )
  assert g.is_equivalent_to( k )

  k.insert_node( name = 'b',
                 in_names = [],
                 out_names = ['a']
                )

  assert ( not k.is_contained_in(g) )
  assert ( g.is_contained_in(k) )
  assert ( not g.is_equivalent_to(k) )
  assert ( not k.is_equivalent_to(g) )

  path_found = g.find_all_paths('a', 'b')
  path_ref = [['a', 'b'], ['a', 'c', 'b']]
  for path in path_ref:
    assert (path in path_found)

  shortest = g.find_shortest_path('a', 'b')
  assert (shortest==[ 'a', 'b' ])

def run():
 tst_graph()


 print 'OK'

if (__name__ == "__main__"):
  run()
