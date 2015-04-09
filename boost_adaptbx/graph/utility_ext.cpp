#include <boost_adaptbx/graph/graph_type.hpp>
#include <boost_adaptbx/graph/graph_export_adaptor.hpp>
#include <boost_adaptbx/graph/vertex_map.hpp>

#include <boost_adaptbx/exporting.hpp>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/dict.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/copy.hpp>

#include <vector>

namespace boost_adaptbx
{
namespace
{

template< typename Graph >
struct basic_operation_export
{
  typedef vertex_map::index_map< Graph > index_map_type;
  typedef typename index_map_type::property_map_type index_property_map_type;

  typedef boost::graph_traits< Graph > graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef graph_export_adaptor::vertex_descriptor_converter< vertex_descriptor > converter;

  static
  void
  copy_graph(Graph const& source, Graph& target)
  {
    index_map_type index_map( source );

    boost::copy_graph(
      source,
      target,
      boost::vertex_index_map( index_map.get() )
      );
  }

  static
  boost::python::dict
  copy_graph_and_map_vertices(Graph const& source, Graph& target)
  {
    index_map_type index_map( source );
    index_property_map_type index_propmap( index_map.get() );

    typedef typename graph_traits::vertex_iterator vertex_iterator;
    typedef boost::iterator_property_map<
      typename std::vector< vertex_descriptor >::iterator,
      index_property_map_type,
      vertex_descriptor,
      vertex_descriptor&
      > iso_map_type;

    std::vector< vertex_descriptor > iso_values( boost::num_vertices( source ) );
    iso_map_type mapV( iso_values.begin(), index_propmap );

    boost::copy_graph(
      source,
      target,
      boost::vertex_index_map( index_propmap ).
      orig_to_copy( mapV )
      );

    boost::python::dict result;

    vertex_iterator vi, vj;

    for( boost::tie( vi, vj ) = boost::vertices( source ); vi != vj; ++vi )
    {
      result[ converter::forward( *vi ) ] = converter::forward( boost::get( mapV, *vi ) );
    }

    return result;
  }

  static void process()
  {
    using namespace boost::python;

    def( "copy_graph", copy_graph, ( arg( "source" ), arg( "target" ) ) );
    def(
      "copy_graph_and_map_vertices",
      copy_graph_and_map_vertices,
      ( arg( "source" ), arg( "target" ) )
      );
  }
};

struct basic_operation_exporter
{
  template< typename Export >
  void operator ()(boost::mpl::identity< Export > myexport) const
  {
    typedef typename Export::first graph_type;
    typedef typename Export::second name_type;
    typedef basic_operation_export< graph_type > exporter_type;
    exporter_type::process();
  }
};

} // namespace <anonymous>
} // namespace boost_adaptbx

BOOST_PYTHON_MODULE(boost_adaptbx_graph_utility_ext)
{
  boost_adaptbx::exporting::class_list< boost_adaptbx::graph_type::exports >::process(
    boost_adaptbx::basic_operation_exporter()
    );
}
