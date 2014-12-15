#include <boost_adaptbx/graph/graph_type.hpp>
#include <boost_adaptbx/graph/graph_export_adaptor.hpp>

#include <boost_adaptbx/exporting.hpp>

#include <boost/python/module.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/list.hpp>
#include <boost/python/def.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/property_map/property_map.hpp>

#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

#include <map>
#include <vector>

namespace boost_adaptbx
{
namespace
{

template< typename Graph >
struct connected_components_export
{
  typedef boost::graph_traits< Graph > graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::vertices_size_type vertices_size_type;
  typedef std::map< vertex_descriptor, vertices_size_type > component_map_type;
  typedef boost::associative_property_map< component_map_type >
    component_property_map_type;
  typedef std::map< vertex_descriptor, boost::default_color_type > color_map_type;
  typedef boost::associative_property_map< color_map_type >
    color_property_map_type;
  typedef graph_export_adaptor::vertex_descriptor_converter< vertex_descriptor > converter;

  static boost::python::list connected_components(Graph const& graph)
  {
    using namespace boost;
    component_map_type component_map;
    color_map_type color_map;
    int num = boost::connected_components(
      graph,
      component_property_map_type( component_map ),
      boost::color_map( color_property_map_type( color_map ) )
      );

    boost::python::list result;

    for(
      typename component_map_type::const_iterator it = component_map.begin();
      it != component_map.end();
      ++it
      )
    {
      result.append( boost::python::make_tuple( converter::forward( it->first ), it->second ) );
    }

    return result;
  }

  static void process()
  {
    using namespace boost::python;

    def( "connected_components", connected_components, arg( "graph" ) );
  }
};

template< typename EdgeList, typename VertexProperty, typename EdgeProperty >
struct connected_components_export<
  boost::adjacency_list< EdgeList, boost::vecS, boost::undirectedS, VertexProperty, EdgeProperty >
  >
{
  typedef boost::adjacency_list< EdgeList, boost::vecS, boost::undirectedS, VertexProperty, EdgeProperty > graph_type;
  typedef boost::graph_traits< graph_type > graph_traits;
  typedef typename graph_traits::vertices_size_type vertices_size_type;
  typedef std::vector< vertices_size_type > component_map_type;

  static boost::python::list connected_components(graph_type const& graph)
  {
    using namespace boost;
    component_map_type component_map( boost::num_vertices( graph ) );
    int num = boost::connected_components(
      graph,
      &component_map[0]
      );

    boost::python::list result;

    for(
      typename component_map_type::size_type i = 0;
      i != component_map.size();
      ++i
      )
    {
      result.append( boost::python::make_tuple( i, component_map[i] ) );
    }

    return result;
  }

  static void process()
  {
    using namespace boost::python;

    def( "connected_components", connected_components, arg( "graph" ) );
  }
};

struct connected_component_algorithm_exporter
{
  template< typename Export >
  void operator ()(boost::mpl::identity< Export > myexport) const
  {
    typedef typename Export::first graph_type;
    typedef typename Export::second name_type;

    typedef boost::graph_traits< graph_type > graph_traits;
    typedef typename boost::mpl::if_<
      boost::is_same<
        typename  graph_traits::directed_category,
        boost::undirected_tag
        >,
      connected_components_export< graph_type >,
      graph_export_adaptor::no_export< graph_type >
      >::type exporter_type;
    exporter_type::process();
  }
};

} // namespace <anonymous>
} // namespace boost_adaptbx

BOOST_PYTHON_MODULE(boost_adaptbx_graph_connected_component_algorithm_ext)
{
  boost_adaptbx::exporting::class_list< boost_adaptbx::graph_type::exports >::process(
    boost_adaptbx::connected_component_algorithm_exporter()
    );
}
