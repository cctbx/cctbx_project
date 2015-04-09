#ifndef BOOST_ADAPTBX_GRAPH_MAXIMUM_CLIQUE_RASCAL_H
#define BOOST_ADAPTBX_GRAPH_MAXIMUM_CLIQUE_RASCAL_H

#include <boost_adaptbx/graph/vertex_map.hpp>

#include <boost/graph/sequential_vertex_coloring.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>

#include <vector>
#include <set>
#include <stack>
#include <algorithm>
#include <iterator>

namespace boost_adaptbx
{
namespace graph
{

template< typename VertexDescriptor, typename SizeType >
class rascal_state
{
public:
  typedef VertexDescriptor vertex_descriptor_type;
  typedef SizeType size_type;
  typedef std::set< vertex_descriptor_type > vertex_group_type;
  typedef std::vector< vertex_group_type > vertex_partition_type;

  typedef typename vertex_partition_type::const_iterator const_partition_iterator;

private:
  vertex_partition_type m_partition;
  size_type m_upper_bound_kwp;

public:
  rascal_state()
    : m_upper_bound_kwp( 0 )
  {};

  rascal_state(size_type const& upper_bound_kwp)
    : m_partition( upper_bound_kwp ), m_upper_bound_kwp( upper_bound_kwp )
  {};

  rascal_state(vertex_partition_type const& partition, size_type const& upper_bound_kwp)
    : m_partition( partition ), m_upper_bound_kwp( upper_bound_kwp )
  {};

  size_type& upper_bound_kwp()
  {
    return m_upper_bound_kwp;
  }

  size_type const& upper_bound_kwp() const
  {
    return m_upper_bound_kwp;
  }

  size_type upper_bound_ka() const
  {
    return m_partition.size();
  }

  size_type upper_bound() const
  {
    return std::min( upper_bound_kwp(), upper_bound_ka() );
  }

  vertex_partition_type& partition()
  {
    return m_partition;
  }

  vertex_partition_type const& partition() const
  {
    return m_partition;
  }
};


struct size_sort_predicate
{
  template< typename Container >
  bool operator ()(Container const& left, Container const& right) const
  {
    return left.size() > right.size();
  }
};

struct initial_partition_by_vertex_coloring
{
  template< typename VertexListGraph >
  rascal_state<
    typename boost::graph_traits< VertexListGraph >::vertex_descriptor,
    typename boost::graph_traits< VertexListGraph >::vertices_size_type
    >
  operator ()(VertexListGraph const& g) const
  {
    typedef boost::graph_traits< VertexListGraph > graph_traits;
    typedef typename graph_traits::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits::vertex_iterator vertex_iterator;
    typedef typename graph_traits::vertices_size_type vertices_size_type;

    typedef vertex_map::generic_vertex_map< VertexListGraph, vertices_size_type >
      color_map_type;
    typedef typename color_map_type::property_map_type color_property_map_type;

    typedef rascal_state< vertex_descriptor, vertices_size_type > state_type;
    typedef typename state_type::vertex_group_type vertex_group_type;
    typedef typename state_type::vertex_partition_type vertex_partition_type;

    color_map_type color_map( g );
    color_property_map_type color_property_map( color_map.get() );
    vertices_size_type color_count = boost::sequential_vertex_coloring(
      g,
      color_property_map
      );

    state_type result( color_count );

    vertex_iterator v, vend;
    vertex_partition_type& partition = result.partition();

    for( boost::tie( v, vend ) = boost::vertices( g ); v != vend; ++v )
    {
      vertex_descriptor const& vertex = *v;
      partition[ boost::get( color_property_map, vertex ) ].insert( vertex );
    }

    std::stable_sort( partition.begin(), partition.end(), size_sort_predicate() );
    return result;
  }
};

template< typename Descriptor >
class partial_graph_selection_predicate
{
public:
  typedef std::set< Descriptor > descriptor_set_type;
  typedef boost::shared_ptr< descriptor_set_type > descriptor_set_ptr_type;

private:
  descriptor_set_ptr_type m_descriptor_set_ptr;

public:
  partial_graph_selection_predicate()
    : m_descriptor_set_ptr( new descriptor_set_type() )
  {};

  partial_graph_selection_predicate(descriptor_set_ptr_type const& ptr)
    : m_descriptor_set_ptr( ptr )
  {};

  template< typename InputIterator >
  partial_graph_selection_predicate(InputIterator begin, InputIterator end)
    : m_descriptor_set_ptr( new descriptor_set_type( begin, end ) )
  {}

  bool operator ()(Descriptor const& desc) const
  {
    return m_descriptor_set_ptr->find( desc ) != m_descriptor_set_ptr->end();
  }

  descriptor_set_type& selection() const
  {
    return *m_descriptor_set_ptr;
  }
};

template< typename Graph1, typename Graph2, typename InputIterator >
void
selected_subgraph(Graph1 const& graph, Graph2& subgraph, InputIterator begin, InputIterator end)
{
  typedef Graph1 graph_type;
  typedef boost::graph_traits< graph_type > graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::vertex_iterator vertex_iterator;
  typedef typename graph_traits::edge_descriptor edge_descriptor;
  typedef typename graph_traits::out_edge_iterator out_edge_iterator;
  typedef typename graph_traits::vertices_size_type vertices_size_type;

  typedef partial_graph_selection_predicate< vertex_descriptor > vertex_predicate_type;
  typedef typename vertex_predicate_type::descriptor_set_type vertex_descriptor_set_type;
  typedef partial_graph_selection_predicate< edge_descriptor > edge_predicate_type;
  typedef typename edge_predicate_type::descriptor_set_type edge_descriptor_set_type;

  typedef boost::filtered_graph< graph_type, edge_predicate_type, vertex_predicate_type >
    filtered_graph_type;
  typedef vertex_map::index_map< filtered_graph_type > index_map_type;

  vertex_predicate_type keep_vertices_pred( begin, end );
  vertex_descriptor_set_type const& keep_vertices_set = keep_vertices_pred.selection();

  edge_predicate_type keep_edges_pred;
  edge_descriptor_set_type& keep_edges_set = keep_edges_pred.selection();

  for (
    typename vertex_descriptor_set_type::const_iterator it = keep_vertices_set.begin();
    it != keep_vertices_set.end();
    ++it
    )
  {
    out_edge_iterator e, eend;

    for( boost::tie( e, eend ) = boost::out_edges( *it, graph ); e != eend; ++e )
    {
      if ( keep_vertices_set.find( boost::target( *e, graph ) ) != keep_vertices_set.end() )
      {
        keep_edges_set.insert( *e );
      }
    }
  }

  filtered_graph_type filtgraph(
    graph,
    keep_edges_pred,
    keep_vertices_pred
    );

  index_map_type index_map( filtgraph );

  boost::copy_graph(
    filtgraph,
    subgraph,
    boost::vertex_index_map( index_map.get() )
    );
}

struct upper_bound_by_chromatic_number
{
  template< typename VertexListGraph, typename Partition >
  typename boost::graph_traits< VertexListGraph >::vertices_size_type
  operator ()(VertexListGraph const& graph, Partition const& partition) const
  {
    typedef VertexListGraph graph_type;
    typedef boost::graph_traits< graph_type > graph_traits;
    typedef typename graph_traits::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits::vertices_size_type vertices_size_type;

    typedef vertex_map::generic_vertex_map< VertexListGraph, vertices_size_type >
      color_map_type;

    std::vector< vertex_descriptor > vertices;

    for (
      typename Partition::const_iterator it = partition.begin();
      it != partition.end();
      ++it
      )
    {
      std::copy( it->begin(), it->end(), std::back_inserter( vertices ) );
    }

    graph_type subgraph;
    selected_subgraph( graph, subgraph, vertices.begin(), vertices.end() );

    color_map_type color_map( graph );
    return boost::sequential_vertex_coloring( subgraph, color_map.get() );
  }

};

struct empty_size_predicate
{
  template< typename Container >
  bool operator ()(Container const& con) const
  {
    return con.empty();
  }
};

template< typename VertexListGraph, typename Partition, typename SizeType >
void
move_vertices_from_excess_partitions_if_possible(
  VertexListGraph const& graph,
  Partition& partition,
  SizeType lower_bound
  )
{
  typedef VertexListGraph graph_type;
  typedef boost::graph_traits< graph_type > graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::vertex_iterator vertex_iterator;
  typedef typename graph_traits::edge_descriptor edge_descriptor;
  typedef typename graph_traits::adjacency_iterator adjacency_iterator;
  typedef typename graph_traits::vertices_size_type vertices_size_type;
  typedef typename Partition::value_type vertex_group_type;
  typedef typename Partition::iterator partition_iterator;
  typedef typename vertex_group_type::iterator vertex_group_iterator;

  typedef std::set< vertex_descriptor > vertex_descriptor_set;

  partition_iterator divisor = partition.begin() + std::min( lower_bound, partition.size() );

  for ( partition_iterator pit_e = divisor; pit_e != partition.end(); ++pit_e )
  {
    vertex_group_iterator next = pit_e->begin();

    for ( vertex_group_iterator vgit = next; vgit != pit_e->end(); vgit = next )
    {
      ++next;
      adjacency_iterator v, vend;
      boost::tie( v, vend ) = boost::adjacent_vertices( *vgit, graph );
      vertex_descriptor_set neighbours( v, vend );

      for ( partition_iterator pit_m = partition.begin(); pit_m != divisor; ++pit_m )
      {
        vertex_descriptor_set intersect;
        std::set_intersection(
          neighbours.begin(),
          neighbours.end(),
          pit_m->begin(),
          pit_m->end(),
          std::inserter( intersect, intersect.end() )
          );

        if ( intersect.empty() )
        {
          pit_m->insert( *vgit );
          pit_e->erase( vgit );
          break;
        }
      }
    }
  }

  partition.erase(
    std::remove_if( divisor, partition.end(), empty_size_predicate() ),
    partition.end()
    );
  std::stable_sort( partition.begin(), partition.end(), size_sort_predicate() );
}

template<
  typename AdjacencyGraph,
  typename InitialPartition,
  typename UpperBound,
  typename UserCallback
  >
void
maximum_clique_rascal(
  AdjacencyGraph const& graph,
  InitialPartition const& initial_partition,
  UpperBound const& upper_bound,
  UserCallback callback,
  typename boost::graph_traits< AdjacencyGraph >::vertices_size_type lower_bound = 1
  )
{
  typedef boost::graph_traits< AdjacencyGraph > graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::vertices_size_type vertices_size_type;
  typedef typename graph_traits::adjacency_iterator adjacency_iterator;
  typedef rascal_state< vertex_descriptor, vertices_size_type > state_type;
  typedef typename state_type::vertex_group_type vertex_group_type;

  vertices_size_type max_clique_size = lower_bound + 1;
  std::vector< vertex_descriptor > branch_points;
  std::stack< state_type > states;
  states.push( initial_partition( graph ) );
  branch_points.push_back( graph_traits::null_vertex() );

  while ( true )
  {
    state_type& state = states.top();

    if (
      !state.partition().empty()
      && ( max_clique_size <= branch_points.size() + state.upper_bound() )
      )
    {
      typename vertex_group_type::iterator selit = state.partition().back().begin();
      vertex_descriptor selected = *selit;
      state.partition().back().erase( selit );
      branch_points.push_back( selected );
      max_clique_size = std::max( max_clique_size, branch_points.size() );

      adjacency_iterator v, vend;
      boost::tie( v, vend ) = boost::adjacent_vertices( selected, graph );
      std::set< vertex_descriptor > adjacents( v, vend );
      states.push( state_type() );
      state_type& next = states.top();

      for (
        typename state_type::const_partition_iterator pit = state.partition().begin();
        pit != state.partition().end() - 1;
        ++pit
        )
      {
        vertex_group_type igroup;
        std::set_intersection(
          adjacents.begin(),
          adjacents.end(),
          pit->begin(),
          pit->end(),
          std::inserter( igroup, igroup.end() )
          );

        if ( ! igroup.empty() )
        {
          next.partition().push_back( igroup );
        }
      }

      if ( state.partition().back().empty() )
      {
        state.partition().pop_back();
      }

      std::stable_sort(
        next.partition().begin(),
        next.partition().end(),
        size_sort_predicate()
        );
      move_vertices_from_excess_partitions_if_possible(
        graph,
        next.partition(),
        max_clique_size - branch_points.size()
        );
      next.upper_bound_kwp() = upper_bound( graph, next.partition() );
    }

    else
    {
      if ( max_clique_size <= branch_points.size() )
      {
        callback( branch_points.begin() + 1, branch_points.end() );
      }

      do
      {
        states.pop();
        branch_points.pop_back();
      }
      while (
        !states.empty()
        && !( max_clique_size <= branch_points.size() + states.top().upper_bound() )
        );

      if ( !states.empty() )
      {
        state_type& top = states.top();
        top.upper_bound_kwp() =  upper_bound( graph, top.partition() );
      }

      else
      {
        break;
      }
    }
  }
}


} // namespace graph
} // namespace boost_adaptbx

#endif // BOOST_ADAPTBX_GRAPH_MAXIMUM_CLIQUE_RASCAL_H
