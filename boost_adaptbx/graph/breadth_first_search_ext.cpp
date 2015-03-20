#include <boost_adaptbx/graph/graph_type.hpp>
#include <boost_adaptbx/graph/graph_export_adaptor.hpp>
#include <boost_adaptbx/graph/vertex_map.hpp>

#include <boost_adaptbx/exporting.hpp>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include <boost/ref.hpp>

#include <map>
#include <vector>

namespace boost_adaptbx
{

template< typename Graph >
class bfs_visitor_adaptor
{
public:
  typedef Graph graph_type;
  typedef boost::graph_traits< graph_type > graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::edge_descriptor edge_descriptor;
  typedef graph_export_adaptor::vertex_descriptor_converter< vertex_descriptor > converter;

private:
  boost::python::object m_visitor;

public:
  bfs_visitor_adaptor(boost::python::object visitor) : m_visitor( visitor ) {};
  ~bfs_visitor_adaptor() {};

  void call_python_method(const char* name, edge_descriptor e, Graph const& g)
  {
    boost::python::object attr = m_visitor.attr( name );
    attr(e, boost::cref( g ) );
  }

  void call_python_method(const char* name, vertex_descriptor v, Graph const& g)
  {
    boost::python::object attr = m_visitor.attr( name );
    attr(converter::forward( v ), boost::cref( g ) );
  }

  void initialize_vertex(vertex_descriptor v, Graph const& g)
  {
    call_python_method( "initialize_vertex", v, g );
  }

  void discover_vertex(vertex_descriptor v, Graph const& g)
  {
    call_python_method( "discover_vertex", v, g );
  }

  void examine_vertex(vertex_descriptor v, Graph const& g)
  {
    call_python_method( "examine_vertex", v, g );
  }

  void finish_vertex(vertex_descriptor v, Graph const& g)
  {
    call_python_method( "finish_vertex", v, g );
  }

  void examine_edge(edge_descriptor v, Graph const& g)
  {
    call_python_method( "examine_edge", v, g );
  }

  void tree_edge(edge_descriptor v, Graph const& g)
  {
    call_python_method( "tree_edge", v, g );
  }

  void non_tree_edge(edge_descriptor v, Graph const& g)
  {
    call_python_method( "non_tree_edge", v, g );
  }

  void gray_target(edge_descriptor v, Graph const& g)
  {
    call_python_method( "gray_target", v, g );
  }

  void black_target(edge_descriptor v, Graph const& g)
  {
    call_python_method( "black_target", v, g );
  }
};

template< typename Graph >
struct breadth_first_search_export
{
  typedef Graph graph_type;
  typedef boost::graph_traits< graph_type > graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef graph_export_adaptor::vertex_descriptor_converter< vertex_descriptor > converter;
  typedef vertex_map::index_map< Graph > index_map_type;

  static void breadth_first_search(
    graph_type const& graph,
    typename converter::type vertex,
    boost::python::object vis
    )
  {
    using namespace boost;
    index_map_type index_map( graph );
    boost::breadth_first_search(
      graph,
      converter::backward( vertex ),
      visitor( bfs_visitor_adaptor< graph_type >( vis ) ).
      vertex_index_map( index_map.get() )
      );
  }

  static void process()
  {
    using namespace boost::python;

    def(
      "breadth_first_search",
      breadth_first_search,
      ( arg( "graph" ), arg( "vertex" ), arg( "visitor" ) )
      );
  }
};

struct bfs_exporter
{
  template< typename Export >
  void operator ()(boost::mpl::identity< Export > myexport) const
  {
    typedef typename Export::first graph_type;
    typedef typename Export::second name_type;

    breadth_first_search_export< graph_type >::process();
  }
};

} // namespace boost_adaptbx

BOOST_PYTHON_MODULE(boost_adaptbx_graph_breadth_first_search_ext)
{
  boost_adaptbx::exporting::class_list< boost_adaptbx::graph_type::exports >::process(
    boost_adaptbx::bfs_exporter()
    );
}
