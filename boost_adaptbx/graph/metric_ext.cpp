#include <boost_adaptbx/graph/graph_type.hpp>
#include <boost_adaptbx/graph/graph_export_adaptor.hpp>
#include <boost_adaptbx/graph/vertex_map.hpp>

#include <boost_adaptbx/exporting.hpp>

#include <boost/python/module.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/def.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/betweenness_centrality.hpp>

namespace boost_adaptbx
{
namespace
{

template< typename Graph >
struct brandes_betweenness_centrality_export
{
  typedef boost::graph_traits< Graph > graph_traits;
  typedef typename graph_traits::vertex_iterator vertex_iterator;
  typedef typename graph_traits::edge_iterator edge_iterator;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef graph_export_adaptor::vertex_descriptor_converter< vertex_descriptor > converter;
  typedef vertex_map::index_map< Graph > index_map_type;
  typedef vertex_map::generic_vertex_map< Graph, double > vertex_centrality_map_type;
  typedef edge_map::generic_edge_map< Graph, double > edge_centrality_map_type;
  typedef edge_centrality_map_type weight_map_type;
  typedef typename weight_map_type::property_map_type weight_property_map_type;

  static boost::python::tuple brandes_betweenness_centrality(Graph const& graph)
  {
    using namespace boost;
    index_map_type index_map( graph );
    vertex_centrality_map_type vertex_centrality_map( graph );
    edge_centrality_map_type edge_centrality_map( graph );
    weight_map_type weight_map( graph );
    weight_property_map_type wpropmap( weight_map.get() );
    edge_iterator ei, ej;
    vertex_iterator vi, vj;

    double w;

    for( boost::tie( ei, ej ) = boost::edges( graph ); ei != ej; ++ei )
    {
      w = boost::python::extract< double >( boost::get( boost::edge_weight, graph, *ei ) );
      boost::put( wpropmap, *ei, w);
    }

    boost::brandes_betweenness_centrality(
      graph,
      centrality_map( vertex_centrality_map.get() ).
      edge_centrality_map( edge_centrality_map.get() ).
      vertex_index_map( index_map.get()  ).
      weight_map( wpropmap )
      );

    boost::python::dict vertex_centrality_dict;
    boost::python::dict edge_centrality_dict;

    for( boost::tie( vi, vj ) = boost::vertices( graph ); vi != vj; ++vi )
    {
      vertex_centrality_dict[ converter::forward( *vi ) ] = boost::get(
        vertex_centrality_map.get(),
        *vi
        );
    }

    for( boost::tie( ei, ej ) = boost::edges( graph ); ei != ej; ++ei )
    {
      edge_centrality_dict[ *ei ] = boost::get( edge_centrality_map.get(), *ei );
    }

    return boost::python::make_tuple( vertex_centrality_dict, edge_centrality_dict );
  }

  static void process()
  {
    using namespace boost::python;

    def( "brandes_betweenness_centrality", brandes_betweenness_centrality, arg( "graph" ) );
  }
};

struct metric_exporter
{
  template< typename Export >
  void operator ()(boost::mpl::identity< Export > myexport) const
  {
    typedef typename Export::first graph_type;
    typedef typename Export::second name_type;

    brandes_betweenness_centrality_export< graph_type >::process();
  }
};

} // namespace <anonymous>
} // namespace boost_adaptbx

BOOST_PYTHON_MODULE(boost_adaptbx_graph_metric_ext)
{
  boost_adaptbx::exporting::class_list< boost_adaptbx::graph_type::exports >::process(
    boost_adaptbx::metric_exporter()
    );
}
