#include <boost_adaptbx/graph/graph_type.hpp>
#include <boost_adaptbx/graph/graph_export_adaptor.hpp>

#include <boost_adaptbx/exporting.hpp>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/iterator.hpp>

#include <boost/graph/graph_traits.hpp>

#include <boost/mpl/transform.hpp>
#include <boost/mpl/string.hpp>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#include <boost/functional/hash.hpp>

#include <string>

namespace boost_adaptbx
{
namespace
{

template< typename Graph >
struct graph_export
{
  typedef boost::graph_traits< Graph > graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::edge_descriptor edge_descriptor;
  typedef typename graph_traits::vertex_iterator vertex_iterator;
  typedef typename graph_traits::adjacency_iterator adjacency_iterator;
  typedef typename graph_traits::edge_iterator edge_iterator;
  typedef typename graph_traits::out_edge_iterator out_edge_iterator;
  typedef graph_export_adaptor::vertex_descriptor_converter< vertex_descriptor > converter;
  typedef boost::transform_iterator< converter, vertex_iterator > transformed_vertex_iterator;
  typedef boost::transform_iterator< converter, adjacency_iterator > transformed_adjacency_iterator;

  static
  typename converter::type
  add_vertex(Graph& graph, boost::python::object const& prop)
  {
    return converter::forward(
      boost::add_vertex( graph_type::vertex_property( prop ), graph )
      );
  }

  static
  boost::python::tuple
  add_edge(
    Graph& graph,
    typename converter::type u,
    typename converter::type v,
    boost::python::object const& prop
    )
  {
    std::pair< edge_descriptor, bool > res = boost::add_edge(
      converter::backward( u ),
      converter::backward( v ),
      graph_type::edge_property( prop ),
      graph
      );
    return boost::python::make_tuple( res.first, res.second );
  }

  static
  typename converter::type
  source_vertex(Graph const& graph, edge_descriptor edge)
  {
    return converter::forward( boost::source( edge, graph ) );
  }

  static
  typename converter::type
  target_vertex(Graph const& graph, edge_descriptor edge)
  {
    return converter::forward( boost::target( edge, graph ) );
  }

  static
  transformed_vertex_iterator
  vertex_iterator_begin(Graph const& graph)
  {
    return transformed_vertex_iterator(
      boost::vertices( graph ).first,
      converter()
      );
  }

  static
  transformed_vertex_iterator
  vertex_iterator_end(Graph const& graph)
  {
    return transformed_vertex_iterator(
      boost::vertices( graph ).second,
      converter()
      );
  }

  static
  transformed_adjacency_iterator
  adjacent_vertex_iterator_begin(Graph const& graph, typename converter::type vertex)
  {
    return transformed_adjacency_iterator(
      boost::adjacent_vertices( converter::backward( vertex ), graph ).first,
      converter()
      );
  }

  static
  transformed_adjacency_iterator
  adjacent_vertex_iterator_end(Graph const& graph, typename converter::type vertex)
  {
    return transformed_adjacency_iterator(
      boost::adjacent_vertices( converter::backward( vertex ), graph ).second,
      converter()
      );
  }

  static
  boost::python::object
  adjacent_vertices(Graph const& graph, typename converter::type vertex)
  {
    return boost::python::range< boost::python::default_call_policies, Graph const>(
      boost::bind( &graph_export< Graph >::adjacent_vertex_iterator_begin, _1, vertex ),
      boost::bind( &graph_export< Graph >::adjacent_vertex_iterator_end, _1, vertex )
      )( boost::cref( graph ) );
  }

  static
  boost::python::object
  vertex_label(Graph const& graph, typename converter::type vertex)
  {
    return boost::get( boost::vertex_name_t(), graph, converter::backward( vertex ) );
  }

  static
  void
  set_vertex_label(
    Graph& graph,
    typename converter::type vertex,
    boost::python::object label
    )
  {
    boost::put( boost::vertex_name_t(), graph, converter::backward( vertex ), label );
  }

  static
  boost::python::object
  edge_weight(Graph const& graph, edge_descriptor edge)
  {
    return boost::get( boost::edge_weight_t(), graph, edge );
  }

  static
  void
  set_edge_weight(
    Graph& graph,
    edge_descriptor edge,
    boost::python::object weight
    )
  {
    boost::put( boost::edge_weight_t(), graph, edge, weight );
  }


  static
  edge_iterator
  edge_iterator_begin(Graph const& graph)
  {
    return boost::edges( graph ).first;
  }

  static
  edge_iterator
  edge_iterator_end(Graph const& graph)
  {
    return boost::edges( graph ).second;
  }

  static
  out_edge_iterator
  out_edge_iterator_begin(Graph const& graph, typename converter::type vertex)
  {
    return boost::out_edges( converter::backward( vertex ), graph ).first;
  }

  static
  out_edge_iterator
  out_edge_iterator_end(Graph const& graph, typename converter::type vertex)
  {
    return boost::out_edges( converter::backward( vertex ), graph ).second;
  }

  static
  boost::python::object
  out_edges(Graph const& graph, typename converter::type vertex)
  {
    return boost::python::range< boost::python::default_call_policies, Graph const>(
      boost::bind( &graph_export< Graph >::out_edge_iterator_begin, _1, vertex ),
      boost::bind( &graph_export< Graph >::out_edge_iterator_end, _1, vertex )
      )( boost::cref( graph ) );
  }

  static
  void
  remove_vertex(Graph& graph, typename converter::type vertex)
  {
    vertex_descriptor vd = converter::backward( vertex );
    boost::clear_vertex( vd, graph );
    boost::remove_vertex( vd, graph );
  }

  static
  std::size_t
  num_vertices(Graph const& graph)
  {
    return boost::num_vertices( graph );
  }

  static
  std::size_t
  num_edges(Graph const& graph)
  {
    return boost::num_edges( graph );
  }

  static
  void
  remove_edge(Graph& graph, edge_descriptor edge)
  {
    boost::remove_edge( edge, graph );
  }

  static
  void
  process(std::string const& name)
  {
    using namespace boost::python;

    class_< Graph >( ( "graph_" + name ).c_str(), no_init )
      .def( init<>() )
      .def( "vertices", boost::python::range( vertex_iterator_begin, vertex_iterator_end ) )
      .def( "source", source_vertex, arg( "edge" ) )
      .def( "target", target_vertex, arg( "edge" ) )
      .def( "adjacent_vertices", adjacent_vertices, arg( "vertex" ) )
      .def( "edges", boost::python::range( edge_iterator_begin, edge_iterator_end ) )
      .def( "out_edges", out_edges, arg( "vertex" ) )
      .def( "vertex_label", vertex_label, arg( "vertex" ) )
      .def( "set_vertex_label", set_vertex_label, ( arg( "vertex" ), arg( "label" ) ) )
      .def( "edge_weight", edge_weight, arg( "edge" ) )
      .def( "set_edge_weight", set_edge_weight, ( arg( "edge" ), arg( "weight" ) ) )
      .def( "add_vertex", add_vertex, arg( "label" ) = object() )
      .def(
        "add_edge",
        &add_edge,
        ( arg( "vertex1" ), arg( "vertex2" ), arg( "weight" ) = object() )
        )
      .def( "remove_vertex", remove_vertex, arg( "vertex" ) )
      .def( "remove_edge", remove_edge, arg( "edge" ) )
      .def( "num_vertices", num_vertices )
      .def( "num_edges", num_edges )
      ;
  }
};

struct graph_exporter
{
  template< typename Export >
  void operator ()(boost::mpl::identity< Export > myexport) const
  {
    typedef typename Export::first graph_type;
    typedef typename Export::second name_type;

    graph_export< graph_type >::process(
      std::string( boost::mpl::c_str< name_type >::value )
      );
  }
};

struct get_graph_vertex_descriptor_type
{
  template< class Export >
  struct apply
  {
    typedef typename Export::first graph_type;
    typedef typename boost::graph_traits< graph_type >::vertex_descriptor type;
  };
};

struct get_graph_edge_descriptor_type
{
  template< class Export >
  struct apply
  {
    typedef typename Export::first graph_type;
    typedef typename boost::graph_traits< graph_type >::edge_descriptor type;
  };
};

} // namespace <anonymous>

namespace exporting
{

template< typename DirectedCatg, typename VertexDesc >
struct python_type_export_traits<
  boost::detail::edge_desc_impl< DirectedCatg, VertexDesc >
  >
{
  typedef boost::detail::edge_desc_impl< DirectedCatg, VertexDesc > edge_descriptor;

  static
  std::size_t edge_descriptor_hash_value(edge_descriptor const& edge)
  {
    return boost::hash_value( edge.get_property() );
  }

  static
  void process(boost::python::class_< edge_descriptor >& myclass)
  {
    myclass.def( "__hash__", edge_descriptor_hash_value );
  }
};

} // namespace exporting

} // namespace boost_adaptbx

BOOST_PYTHON_MODULE(boost_adaptbx_graph_ext)
{
  boost_adaptbx::exporting::class_list< boost_adaptbx::graph_type::exports >::process(
    boost_adaptbx::graph_exporter()
    );

  //
  // Dependent types
  //
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wmultichar"

  // Vertex descriptors
  typedef boost::mpl::transform<
    boost_adaptbx::graph_type::exports,
    boost_adaptbx::get_graph_vertex_descriptor_type
    >::type graph_vertex_descriptor_types;
  typedef boost_adaptbx::exporting::unique_type_set<
    graph_vertex_descriptor_types
    >::type unique_graph_vertex_descriptor_types;
  typedef boost_adaptbx::exporting::novel_python_type<
     unique_graph_vertex_descriptor_types
     >::type novel_unique_graph_vertex_descriptor_types;

  boost_adaptbx::exporting::class_list< novel_unique_graph_vertex_descriptor_types >::process(
    boost_adaptbx::exporting::python_type_export< boost::mpl::string< 'vert', 'ex_' > >()
    );

  // Edge descriptors
  typedef boost::mpl::transform<
    boost_adaptbx::graph_type::exports,
    boost_adaptbx::get_graph_edge_descriptor_type
    >::type graph_edge_descriptor_types;
  typedef boost_adaptbx::exporting::unique_type_set<
    graph_edge_descriptor_types
    >::type unique_graph_edge_descriptor_types;

  boost_adaptbx::exporting::class_list< unique_graph_edge_descriptor_types >::process(
    boost_adaptbx::exporting::python_type_export<
      boost::mpl::string< 'edge', '_' >
      >().enable_equality_operators()
    );

  #pragma clang diagnostic pop
}
