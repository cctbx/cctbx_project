#include <boost_adaptbx/graph/graph_type.hpp>
#include <boost_adaptbx/graph/graph_export_adaptor.hpp>
#include <boost_adaptbx/graph/vertex_map.hpp>

#include <boost_adaptbx/exporting.hpp>

#include <boost/python/module.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/list.hpp>
#include <boost/python/def.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/connected_components.hpp>

#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

namespace boost_adaptbx
{
namespace
{

template< typename Graph >
struct connected_components_export
{
  typedef boost::graph_traits< Graph > graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::vertex_iterator vertex_iterator;
  typedef typename graph_traits::vertices_size_type vertex_index_type;
  typedef graph_export_adaptor::vertex_descriptor_converter< vertex_descriptor > converter;
  typedef vertex_map::generic_vertex_map< Graph, vertex_index_type > component_map_type;
  typedef typename component_map_type::property_map_type component_property_map_type;
  typedef vertex_map::index_map< Graph > index_map_type;

  static boost::python::list connected_components(Graph const& graph)
  {
    using namespace boost;
    index_map_type index_map( graph );
    component_map_type component_map( graph );
    component_property_map_type component_property_map( component_map.get() );

    int num = boost::connected_components(
      graph,
      component_property_map,
      vertex_index_map( index_map.get()  )
      );

    vertex_iterator di, dj;
    boost::python::list result;

    for( boost::tie( di, dj ) = boost::vertices( graph ); di != dj; ++di )
    {
      result.append(
        boost::python::make_tuple(
          converter::forward( *di ),
          boost::get( component_property_map, *di )
          )
        );
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
