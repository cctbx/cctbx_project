#ifndef BOOST_ADAPTBX_GRAPH_VERTEX_MAP_H
#define BOOST_ADAPTBX_GRAPH_VERTEX_MAP_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

#include <map>
#include <vector>

namespace boost_adaptbx
{
namespace vertex_map
{

namespace detail
{

template< typename Graph, typename Value >
struct associative_vertex_map_impl
{
  typedef boost::graph_traits< Graph > graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor_type;
  typedef std::map< vertex_descriptor_type, Value > storage_type;
  typedef boost::associative_property_map< storage_type > property_map_type;
  typedef Value map_value_type;

  storage_type data_for_;
  property_map_type vertex_map_;

  associative_vertex_map_impl(Graph const& graph) : vertex_map_( data_for_ )
  {}

  property_map_type get()
  {
    return vertex_map_;
  }
};

template< typename Graph, typename Value >
struct serial_vertex_map_impl
{
  typedef std::vector< Value > storage_type;
  typedef typename storage_type::pointer property_map_type;
  typedef Value map_value_type;

  storage_type vertex_map_;

  serial_vertex_map_impl(Graph const& graph)
    : vertex_map_( boost::num_vertices( graph ) )
  {}

  property_map_type get()
  {
    return &vertex_map_[0];
  }
};

template< typename Graph, typename PropertyTag >
struct internal_vertex_map_impl
{
  typedef typename boost::mpl::if_<
    boost::is_const< Graph >,
    typename boost::property_map< Graph, PropertyTag >::const_type,
    typename boost::property_map< Graph, PropertyTag >::type
    >::type property_map_type;
  typedef typename boost::property_traits< property_map_type >::value_type
    map_value_type;

  property_map_type vertex_map_;

  internal_vertex_map_impl(Graph const& graph)
    : vertex_map_( boost::get( PropertyTag(), graph ) )
  {}

  property_map_type get()
  {
    return vertex_map_;
  }
};

template< typename Graph >
struct associative_index_map_impl
{
  typedef boost::graph_traits< Graph > graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor_type;
  typedef typename graph_traits::vertex_iterator vertex_iterator_type;
  typedef typename graph_traits::vertices_size_type vertex_index_type;
  typedef associative_vertex_map_impl< Graph, vertex_index_type > vertex_map_type;
  typedef typename vertex_map_type::property_map_type property_map_type;

  vertex_map_type vertex_map_;

  associative_index_map_impl(Graph const& graph) :vertex_map_( graph )
  {
    property_map_type vim( vertex_map_.get() );
    vertex_index_type index = 0;
    vertex_iterator_type di, dj;

    for ( boost::tie( di, dj ) = boost::vertices( graph ) ; di != dj; ++di )
    {
      boost::put( vim, *di, index++ );
    }
  }

  property_map_type get()
  {
    return vertex_map_.get();
  }
};

template< typename Graph >
struct property_index_map_impl
{
  typedef internal_vertex_map_impl< Graph, boost::vertex_index_t > vertex_map_type;
  typedef typename vertex_map_type::map_value_type vertex_index_type;
  typedef typename vertex_map_type::property_map_type property_map_type;

  vertex_map_type vertex_map_;

  property_index_map_impl(Graph const& graph) : vertex_map_( graph )
  {}

  property_map_type get()
  {
    return vertex_map_.get();
  }
};

} // namespace detail

template< typename Graph, typename Value >
struct generic_vertex_map
{
  typedef Graph graph_type;
  typedef Value value_type;
  typedef detail::associative_vertex_map_impl< graph_type, value_type > vertex_map_impl_type;
  typedef typename vertex_map_impl_type::property_map_type property_map_type;

  vertex_map_impl_type vertex_map_impl_;

  generic_vertex_map(Graph const& graph) : vertex_map_impl_( graph )
  {}

  property_map_type get()
  {
    return vertex_map_impl_.get();
  }
};

template<
  typename OutEdgeList,
  typename VertexList,
  typename Directed,
  typename VertexProperties,
  typename EdgeProperties,
  typename GraphProperties,
  typename EdgeList,
  typename Value
  >
struct generic_vertex_map<
  boost::adjacency_list<
    OutEdgeList,
    VertexList,
    Directed,
    VertexProperties,
    EdgeProperties,
    GraphProperties,
    EdgeList
    >,
  Value
  >
{
  typedef boost::adjacency_list<
    OutEdgeList,
    VertexList,
    Directed,
    VertexProperties,
    EdgeProperties,
    GraphProperties,
    EdgeList
    > graph_type;
  typedef Value value_type;
  typedef typename boost::mpl::if_<
    boost::is_same< typename graph_type::vertex_list_selector, boost::vecS >,
    detail::serial_vertex_map_impl< graph_type, value_type >,
    detail::associative_vertex_map_impl< graph_type, value_type >
    >::type vertex_map_impl_type;
  typedef typename vertex_map_impl_type::property_map_type property_map_type;

  vertex_map_impl_type vertex_map_impl_;

  generic_vertex_map(graph_type const& graph) : vertex_map_impl_( graph )
  {}

  property_map_type get()
  {
    return vertex_map_impl_.get();
  }
};

template< typename Graph >
struct index_map
{
  typedef Graph graph_type;
  typedef detail::associative_index_map_impl< graph_type > index_map_impl_type;
  typedef typename index_map_impl_type::property_map_type property_map_type;

  index_map_impl_type index_map_impl_;

  index_map(Graph const& graph) : index_map_impl_( graph )
  {}

  property_map_type get()
  {
    return index_map_impl_.get();
  }
};

template<
  typename OutEdgeList,
  typename VertexList,
  typename Directed,
  typename VertexProperties,
  typename EdgeProperties,
  typename GraphProperties,
  typename EdgeList
  >
struct index_map<
  boost::adjacency_list<
    OutEdgeList,
    VertexList,
    Directed,
    VertexProperties,
    EdgeProperties,
    GraphProperties,
    EdgeList
    >
  >
{
  typedef boost::adjacency_list<
    OutEdgeList,
    VertexList,
    Directed,
    VertexProperties,
    EdgeProperties,
    GraphProperties,
    EdgeList
    > graph_type;
  typedef typename boost::mpl::if_<
    boost::is_same< typename graph_type::vertex_list_selector, boost::vecS >,
    detail::property_index_map_impl< graph_type >,
    detail::associative_index_map_impl< graph_type >
    >::type index_map_impl_type;

  typedef typename index_map_impl_type::property_map_type property_map_type;

  index_map_impl_type index_map_impl_;

  index_map(graph_type const& graph) : index_map_impl_( graph )
  {}

  property_map_type get()
  {
    return index_map_impl_.get();
  }
};

} // namespace vertex_map

namespace edge_map
{

template< typename Graph, typename Value >
struct generic_edge_map
{
  typedef Graph graph_type;
  typedef Value value_type;
  typedef boost::graph_traits< Graph > graph_traits;
  typedef typename graph_traits::edge_descriptor edge_descriptor_type;
  typedef std::map< edge_descriptor_type, Value > storage_type;
  typedef boost::associative_property_map< storage_type > property_map_type;

  storage_type data_for_;
  property_map_type edge_map_;

  generic_edge_map(Graph const& graph) : edge_map_( data_for_ )
  {}

  property_map_type get()
  {
    return edge_map_;
  }
};

} // namespace edge_map
} // namespace boost_adaptbx

#endif // BOOST_ADAPTBX_GRAPH_VERTEX_MAP_H
