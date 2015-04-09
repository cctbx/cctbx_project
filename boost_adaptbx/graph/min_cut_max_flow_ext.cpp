#include <boost_adaptbx/graph/graph_type.hpp>
#include <boost_adaptbx/graph/graph_export_adaptor.hpp>
#include <boost_adaptbx/graph/vertex_map.hpp>

#include <boost_adaptbx/exporting.hpp>

#include <boost/python/module.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/def.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/stl_iterator.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/stoer_wagner_min_cut.hpp>
#include <boost/graph/one_bit_color_map.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/properties.hpp>

#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

#include <iostream>

namespace boost_adaptbx
{

template< typename Graph >
struct minimum_cut_export
{
  typedef boost::graph_traits< Graph > graph_traits;
  typedef typename graph_traits::vertex_iterator vertex_iterator;
  typedef typename graph_traits::edge_iterator edge_iterator;
  typedef vertex_map::index_map< Graph > index_map_type;
  typedef typename index_map_type::property_map_type index_property_map_type;
  typedef boost::one_bit_color_map< index_property_map_type > parity_map_type;
  typedef edge_map::generic_edge_map< Graph, double > edge_map_type;
  typedef typename edge_map_type::property_map_type edge_property_map_type;

  static boost::python::tuple stoer_wagner_minimum_cut(Graph const& graph)
  {
    index_map_type index_map( graph );
    parity_map_type parities = boost::make_one_bit_color_map(
      boost::num_vertices( graph ),
      index_map.get()
      );

    edge_map_type edge_map( graph );
    edge_property_map_type epropmap( edge_map.get() );
    edge_iterator ei, ej;
    double w;

    for( boost::tie( ei, ej ) = boost::edges( graph ); ei != ej; ++ei )
    {
      w = boost::python::extract< double >( boost::get( boost::edge_weight, graph, *ei ) );
      boost::put( epropmap, *ei, w);
    }

    w = boost::stoer_wagner_min_cut(
      graph,
      epropmap,
      boost::parity_map( parities ).
      vertex_index_map( index_map.get() )
      );

    vertex_iterator di, dj;
    boost::python::list result;

    for( boost::tie( di, dj ) = boost::vertices( graph ); di != dj; ++di )
    {
      result.append( bool( boost::get( parities, *di ) ) );
    }

    return boost::python::make_tuple( w, result );
  }

  static void process()
  {
    using namespace boost::python;

    def( "stoer_wagner_min_cut", stoer_wagner_minimum_cut, arg( "graph" ) );
  }
};

struct min_cut_exporter
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
      minimum_cut_export< graph_type >,
      graph_export_adaptor::no_export< graph_type >
      >::type exporter_type;
    exporter_type::process();
  }
};

template< typename Graph >
struct maximum_flow_export
{
  typedef boost::graph_traits< Graph > graph_traits;
  typedef typename graph_traits::vertex_iterator vertex_iterator;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::edge_iterator edge_iterator;
  typedef typename graph_traits::edge_descriptor edge_descriptor;

  typedef vertex_map::index_map< Graph > index_map_type;
  typedef typename index_map_type::property_map_type index_property_map_type;

  typedef edge_map::generic_edge_map< Graph, double > capacity_map_type;
  typedef typename capacity_map_type::property_map_type capacity_property_map_type;

  typedef edge_map::generic_edge_map< Graph, edge_descriptor > reverse_edge_map_type;
  typedef typename reverse_edge_map_type::property_map_type reverse_edge_property_map_type;

  typedef vertex_map::generic_vertex_map< Graph, boost::default_color_type > color_map_type;
  typedef typename color_map_type::property_map_type color_property_map_type;

  typedef graph_export_adaptor::vertex_descriptor_converter< vertex_descriptor > converter;

  static boost::python::tuple boykov_kolmogorov_max_flow(
    Graph const& graph,
    boost::python::dict reverse_edge_map,
    typename converter::type source,
    typename converter::type sink
    )
  {
    using namespace boost;
    using namespace boost::python;

    index_map_type index_map( graph );

    capacity_map_type capmap( graph );
    capacity_property_map_type cappropmap( capmap.get() );
    edge_iterator ei, ej;

    for( tie( ei, ej ) = edges( graph ); ei != ej; ++ei )
    {
      boost::put( cappropmap, *ei, extract< double >( get( edge_weight, graph, *ei ) ) );
    }

    reverse_edge_map_type revmap( graph );
    reverse_edge_property_map_type revpropmap( revmap.get() );

    stl_input_iterator< object > end;

    for( stl_input_iterator< object > it( reverse_edge_map.iteritems() ); it!=end; ++it)
    {
      edge_descriptor ed = extract< edge_descriptor >( (*it)[0] );
      edge_descriptor red = extract< edge_descriptor >( (*it)[1] );
      put( revpropmap, ed, red );
    }

    color_map_type colormap( graph );
    color_property_map_type colorpropmap( colormap.get() );
    capacity_map_type resicapmap( graph );

    double w = boost::boykov_kolmogorov_max_flow(
        graph,
        cappropmap,
        resicapmap.get(),
        revpropmap,
        colorpropmap,
        index_map.get(),
        converter::backward( source ),
        converter::backward( sink )
        );

    vertex_iterator di, dj;
    boost::python::list result;
    default_color_type const black( color_traits< default_color_type >::black() );
    default_color_type current;

    for( boost::tie( di, dj ) = boost::vertices( graph ); di != dj; ++di )
    {
      current = get( colorpropmap, *di );
      result.append( current == black );
    }

    return boost::python::make_tuple( w, result );
  }

  static void process()
  {
    using namespace boost::python;

    def(
      "boykov_kolmogorov_max_flow",
      boykov_kolmogorov_max_flow,
      ( arg( "graph" ), arg( "reverse_edge_map" ), arg( "source" ), arg( "sink" ) )
      );
  }
};

struct max_flow_exporter
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
        boost::directed_tag
        >,
      maximum_flow_export< graph_type >,
      graph_export_adaptor::no_export< graph_type >
      >::type exporter_type;
    exporter_type::process();
  }
};

} // namespace boost_adaptbx

BOOST_PYTHON_MODULE(boost_adaptbx_graph_min_cut_max_flow_ext)
{
  boost_adaptbx::exporting::class_list< boost_adaptbx::graph_type::exports >::process(
    boost_adaptbx::min_cut_exporter()
    );
  boost_adaptbx::exporting::class_list< boost_adaptbx::graph_type::exports >::process(
    boost_adaptbx::max_flow_exporter()
    );
}
