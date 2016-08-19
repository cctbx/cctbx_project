#include <boost_adaptbx/graph/graph_type.hpp>
#include <boost_adaptbx/graph/graph_export_adaptor.hpp>
#include <boost_adaptbx/graph/vertex_map.hpp>

#include <boost_adaptbx/exporting.hpp>

#include <boost/python/module.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/def.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/bc_clustering.hpp>


namespace boost_adaptbx
{
namespace
{

template< typename Graph >
struct betweenness_centrality_clustering_export
{
  typedef boost::graph_traits< Graph > graph_traits;
  typedef typename graph_traits::edge_iterator edge_iterator;
  typedef edge_map::generic_edge_map< Graph, double > edge_centrality_map_type;
  typedef vertex_map::index_map< Graph > index_map_type;

  static boost::python::dict betweenness_centrality_clustering(
    Graph& graph,
    double const& threshold
    )
  {
    using namespace boost;
    index_map_type index_map( graph );
    edge_centrality_map_type edge_centrality_map( graph );

    boost::betweenness_centrality_clustering(
      graph,
      boost::bc_clustering_threshold< double >( threshold, graph, false ),
      edge_centrality_map.get(),
      index_map.get()
      );

    boost::python::dict edge_centrality_dict;
    edge_iterator ei, ej;

    for( boost::tie( ei, ej ) = boost::edges( graph ); ei != ej; ++ei )
    {
      edge_centrality_dict[ *ei ] = boost::get( edge_centrality_map.get(), *ei );
    }

    return edge_centrality_dict;
  }

  static void process()
  {
    using namespace boost::python;

    def(
      "betweenness_centrality_clustering",
      betweenness_centrality_clustering,
      ( arg( "graph" ), arg( "threshold" ) )
      );
  }
};

struct clustering_algorithm_exporter
{
  template< typename Export >
  void operator ()(boost::mpl::identity< Export > myexport) const
  {
    typedef typename Export::first graph_type;
    typedef typename Export::second name_type;

    betweenness_centrality_clustering_export< graph_type >::process();
  }
};

} // namespace <anonymous>
} // namespace boost_adaptbx

BOOST_PYTHON_MODULE(boost_adaptbx_graph_clustering_algorithm_ext)
{
  boost_adaptbx::exporting::class_list< boost_adaptbx::graph_type::exports >::process(
    boost_adaptbx::clustering_algorithm_exporter()
    );
}
