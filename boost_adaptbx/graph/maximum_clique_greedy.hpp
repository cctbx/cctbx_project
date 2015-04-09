#ifndef BOOST_ADAPTBX_GRAPH_MAXIMUM_CLIQUE_GREEDY_H
#define BOOST_ADAPTBX_GRAPH_MAXIMUM_CLIQUE_GREEDY_H

#include <boost_adaptbx/graph/vertex_map.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/tuple/tuple.hpp>

#include <vector>
#include <set>
#include <queue>
#include <algorithm>
#include <iterator>
#include <utility>

namespace boost_adaptbx
{
namespace graph
{
namespace greedy
{

class maximum_clique_exception
{};

class bad_vertex_exception : public maximum_clique_exception
{};

template< typename AdjacencyGraph >
class partition
{
public:
  typedef AdjacencyGraph graph_type;
  typedef boost::graph_traits< AdjacencyGraph > graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::adjacency_iterator adjacency_iterator;

  typedef std::set< vertex_descriptor > vertex_set_type;
  typedef typename vertex_set_type::iterator vertex_set_iterator;
  typedef typename vertex_set_type::const_iterator vertex_set_const_iterator;

private:
  vertex_set_type clique_;
  vertex_set_type neighbours_;

public:
  partition(graph_type const& graph)
    : neighbours_( boost::vertices( graph ).first, boost::vertices( graph ).second )
  {}

  vertex_set_type const& clique() const
  {
    return clique_;
  }

  vertex_set_type const& neighbours() const
  {
    return neighbours_;
  }

  void add(graph_type const& graph, vertex_descriptor vertex)
  {
    if ( neighbours_.find( vertex ) == neighbours_.end() )
    {
      throw bad_vertex_exception();
    }

    clique_.insert( vertex );

    adjacency_iterator v, vend;
    boost::tie( v, vend ) = boost::adjacent_vertices( vertex, graph );
    vertex_set_type adjacents( v, vend );

    vertex_set_iterator it1 = neighbours_.begin();
    vertex_set_iterator it2 = adjacents.begin();

    while ( ( it1 != neighbours_.end() ) && ( it2 != adjacents.end() ) )
    {
      if ( *it1 < *it2 )
      {
        neighbours_.erase( it1++ );
      }
      else if ( *it2 < *it1 )
      {
        ++it2;
      }
      else
      {
        ++it1;
        ++it2;
      }
    }

    neighbours_.erase( it1, neighbours_.end() );
  }
};

template< typename Partition >
struct degree_scorer
{
  typedef Partition partition_type;
  typedef typename partition_type::vertex_set_type vertex_set_type;
  typedef typename vertex_set_type::size_type set_size_type;
  typedef set_size_type result_type;

public:
  degree_scorer()
  {}

  result_type operator ()(partition_type const& partition) const
  {
    return result_type( partition.neighbours().size() );
  }
};

template< typename AdjacencyGraph, typename Partition >
struct exdegree_scorer
{
  typedef AdjacencyGraph graph_type;
  typedef boost::graph_traits< AdjacencyGraph > graph_traits;
  typedef typename graph_traits::adjacency_iterator adjacency_iterator;
  typedef typename graph_traits::vertices_size_type vertices_size_type;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::vertex_iterator vertex_iterator;

  typedef Partition partition_type;
  typedef typename partition_type::vertex_set_type vertex_set_type;
  typedef typename vertex_set_type::size_type set_size_type;
  typedef typename vertex_set_type::const_iterator vertices_set_const_iterator;
  typedef std::pair< set_size_type, vertices_size_type > result_type;

private:
  graph_type const& graph_;

public:
  exdegree_scorer(graph_type const& graph)
    : graph_( graph )
  {}

  result_type operator ()(partition_type const& partition) const
  {
    vertices_size_type vertices_count( 0 );
    adjacency_iterator av, avend;
    vertex_set_type const& neighbours = partition.neighbours();

    for (
      vertices_set_const_iterator it = neighbours.begin();
       it != neighbours.end();
       ++it
       )
    {
      boost::tie( av, avend ) = boost::adjacent_vertices( *it, graph_ );
      vertices_count += std::distance( av, avend );
    }

    return result_type( neighbours.size(), vertices_count );
  }
};

template< typename AdjacencyGraph, typename Partition >
struct exdegree_scorer_cached
{
  typedef AdjacencyGraph graph_type;
  typedef boost::graph_traits< AdjacencyGraph > graph_traits;
  typedef typename graph_traits::adjacency_iterator adjacency_iterator;
  typedef typename graph_traits::vertices_size_type vertices_size_type;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::vertex_iterator vertex_iterator;

  typedef vertex_map::generic_vertex_map< graph_type, vertices_size_type >
    vertex_exdegree_mapping_type;
  typedef typename vertex_exdegree_mapping_type::property_map_type
    vertex_exdegree_property_map_type;

  typedef Partition partition_type;
  typedef typename partition_type::vertex_set_type vertex_set_type;
  typedef typename vertex_set_type::size_type set_size_type;
  typedef typename vertex_set_type::const_iterator vertices_set_const_iterator;
  typedef std::pair< set_size_type, vertices_size_type > result_type;

private:
  vertex_exdegree_mapping_type exdegree_for_;

public:
  exdegree_scorer_cached(graph_type const& graph) : exdegree_for_( graph )
  {
    vertex_iterator v, vend;
    adjacency_iterator av, avend;
    vertex_exdegree_property_map_type ve_propmap( exdegree_for_.get() );

    for ( boost::tie( v, vend ) = boost::vertices( graph ); v != vend; ++v )
    {
      boost::tie( av, avend ) = boost::adjacent_vertices( *v, graph );
      boost::put( ve_propmap, *v, std::distance( av, avend ) );
    }
  }

  result_type operator ()(partition_type const& partition) const
  {
    vertices_size_type vertices_count( 0 );
    vertex_set_type const& neighbours = partition.neighbours();
    vertex_exdegree_property_map_type ve_propmap( exdegree_for_.get() );

    for (
      vertices_set_const_iterator it = neighbours.begin();
      it != neighbours.end();
      ++it
      )
    {
      vertices_count += boost::get( ve_propmap, *it );
    }

    result_type( neighbours.size(), vertices_count );
  }
};

template< typename Pair >
struct greater_on_second_type
{
  typedef typename Pair::second_type score_type;

  bool operator ()(Pair const& lhs, Pair const& rhs) const
  {
    return lhs.second > rhs.second;
  }
};

template< typename AdjacencyGraph, typename Scorer >
std::vector< partition< AdjacencyGraph > >
maximum_clique(
  AdjacencyGraph const& graph,
  Scorer const& scorer = Scorer(),
  std::size_t maxsol = 0
  )
{
  typedef boost::graph_traits< AdjacencyGraph > graph_traits;
  typedef typename graph_traits::vertex_iterator vertex_iterator;

  typedef partition< AdjacencyGraph > partition_type;
  typedef typename partition_type::vertex_set_type vertex_set_type;
  typedef typename partition_type::vertex_set_const_iterator vertex_set_const_iterator;

  typedef std::vector< partition_type > partition_list;
  typedef typename partition_list::const_iterator partition_const_iterator;

  typedef Scorer scorer_type;
  typedef typename scorer_type::result_type score_type;
  typedef std::pair< partition_type, score_type > partition_score_data;
  typedef greater_on_second_type< partition_score_data > ranker_type;
  typedef std::priority_queue<
    partition_score_data,
    std::vector< partition_score_data >,
    ranker_type >
    partition_data_priority_queue;

  typedef std::set< vertex_set_type > clique_set_type;

  if ( maxsol == 0 )
  {
    vertex_iterator v, vend;
    boost::tie( v, vend ) = boost::vertices( graph );
    maxsol = std::distance( v, vend );
  }
  assert ( 0 < maxsol );

  partition_list solutions;
  solutions.push_back( partition_type( graph ) );
  partition_data_priority_queue solution_queue;
  clique_set_type known_cliques;

  while ( true )
  {
    assert ( solution_queue.empty() );

    for (
      partition_const_iterator pit = solutions.begin();
      pit != solutions.end();
      ++pit
      )
    {
      vertex_set_type const& neighbours = pit->neighbours();

      for (
        vertex_set_const_iterator nit = neighbours.begin();
        nit != neighbours.end();
        ++nit
        )
      {
        partition_type nextpart( *pit );
        nextpart.add( graph, *nit );
        score_type nextscore = scorer( nextpart );

        if ( solution_queue.size() < maxsol )
        {
          if ( known_cliques.insert( nextpart.clique() ).second )
          {
            solution_queue.push( partition_score_data( nextpart, nextscore ) );
          }
        }

        else if ( solution_queue.top().second < nextscore )
        {
          if ( known_cliques.insert( nextpart.clique() ).second )
          {
            known_cliques.erase( solution_queue.top().first.clique() );
            solution_queue.pop();
            solution_queue.push( partition_score_data( nextpart, nextscore ) );
          }
        }
      }
    }

    assert ( solution_queue.size() <= maxsol );

    if ( solution_queue.size() == 0 )
    {
      break;
    }

    solutions.clear();

    while ( !solution_queue.empty() )
    {
      solutions.push_back( solution_queue.top().first );
      solution_queue.pop();
    }

    known_cliques.clear();
  }

  return solutions;
}

} // namespace greedy
} // namespace graph
} // namespace boost_adaptbx

#endif // BOOST_ADAPTBX_GRAPH_MAXIMUM_CLIQUE_GREEDY_H
