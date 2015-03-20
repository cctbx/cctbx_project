#include <boost_adaptbx/graph/graph_type.hpp>
#include <boost_adaptbx/graph/graph_export_adaptor.hpp>
#include <boost_adaptbx/graph/vertex_map.hpp>

#include <boost_adaptbx/exporting.hpp>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/mcgregor_common_subgraphs.hpp>
#include <boost/property_map/property_map.hpp>

#include <map>
#include <iterator>

namespace boost_adaptbx
{
namespace
{

template< typename Graph >
class mcg_python_callback_adaptor
{
public:
  typedef boost::graph_traits< Graph > graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::vertices_size_type vertices_size_type;

  typedef graph_export_adaptor::vertex_descriptor_converter< vertex_descriptor > converter;

private:
  const Graph& m_graph1;
  const Graph& m_graph2;
  boost::python::object m_callable;

public:

mcg_python_callback_adaptor(
  Graph const& graph1,
  Graph const& graph2,
  boost::python::object callable
  )
  : m_graph1( graph1 ), m_graph2( graph2 ), m_callable( callable )
{}

~mcg_python_callback_adaptor()
{}

template< typename VertexMap >
bool operator()(
  VertexMap cmap_1_to_2,
  VertexMap cmap_2_to_1,
  vertices_size_type subgraph_size
  )
{
  boost::python::list pairs;

  BGL_FORALL_VERTICES_T(vertex1, m_graph1, Graph)
  {
    vertex_descriptor vertex2 = boost::get( cmap_1_to_2, vertex1 );

    if ( vertex2 != graph_traits::null_vertex() )
    {
      pairs.append(
        boost::python::make_tuple( converter::forward( vertex1 ), converter::forward( vertex2 ) )
        );
    }
  }

  return m_callable( pairs );
}

};

template< typename PropertyMap >
class python_property_equivalent
{
private:
  const PropertyMap m_property_map1;
  const PropertyMap m_property_map2;
  boost::python::object m_callable;

public:

python_property_equivalent(
  PropertyMap const property_map1,
  PropertyMap const property_map2,
  boost::python::object callable
  )
  : m_property_map1( property_map1 ), m_property_map2( property_map2 ),
    m_callable( callable )
{}

template< typename Item >
bool operator()(Item const item1, Item const item2)
{
  typename PropertyMap::value_type property1 = get(m_property_map1, item1);
  typename PropertyMap::value_type property2 = get(m_property_map2, item2);

  return m_callable( property1, property2 );
}

};

template< typename PropertyMap >
python_property_equivalent< PropertyMap >
make_python_property_equivalent(
  PropertyMap const property_map1,
  PropertyMap const property_map2,
  boost::python::object callable
  )
{
  return python_property_equivalent< PropertyMap >(
    property_map1,
    property_map2,
    callable
    );
}

template< typename Graph >
struct mcgregor_common_subgraphs_unique_export
{
  typedef Graph graph_type;
  typedef boost::graph_traits< Graph > graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::vertices_size_type vertices_size_type;
  typedef vertex_map::index_map< graph_type > index_map_type;

  static
  void
  mcgregor_common_subgraphs_unique(
    graph_type const& graph1,
    graph_type  const& graph2,
    boost::python::object vertex_equality,
    boost::python::object edge_equality,
    boost::python::object callback
    )
  {
    // Vertex equality
    typedef typename boost::property_map< graph_type, boost::vertex_name_t>::const_type vertex_name_map_t;
    typedef python_property_equivalent< vertex_name_map_t > vertex_comp_t;
    vertex_comp_t vertex_comp = make_python_property_equivalent(
      boost::get( boost::vertex_name_t(), graph1 ),
      boost::get( boost::vertex_name_t(), graph2 ),
      vertex_equality
      );

    // Edge equality
    typedef typename boost::property_map< graph_type, boost::edge_weight_t >::const_type edge_weight_map_t;
    typedef python_property_equivalent< edge_weight_map_t > edge_comp_t;
    edge_comp_t edge_comp = make_python_property_equivalent(
      boost::get( boost::edge_weight_t(), graph1 ),
      boost::get( boost::edge_weight_t(), graph2 ),
      edge_equality
      );

    // Callback
    mcg_python_callback_adaptor< graph_type > mcg_callback( graph1, graph2, callback );

    // Call
    index_map_type index_map1( graph1 ), index_map2( graph2 );
    boost::mcgregor_common_subgraphs_unique(
      graph1,
      graph2,
      index_map1.get(),
      index_map2.get(),
      edge_comp,
      vertex_comp,
      true,
      mcg_callback
      );
  }

  static void process()
  {
    using namespace boost::python;

    def(
      "mcgregor_common_subgraphs_unique",
      mcgregor_common_subgraphs_unique,
      ( arg( "graph1" ), arg( "graph2" ), arg( "vertex_equality" ),
        arg( "edge_equality" ), arg( "callback" ) )
      );
  }
};

struct graph_structure_comparison_exporter
{
  template< typename Export >
  void operator ()(boost::mpl::identity< Export > myexport) const
  {
    typedef typename Export::first graph_type;
    typedef typename Export::second name_type;

    mcgregor_common_subgraphs_unique_export< graph_type >::process();
  }
};

} // namespace <anonymnous>
} // namespace boost_adaptbx

BOOST_PYTHON_MODULE(boost_adaptbx_graph_graph_structure_comparison_ext)
{
  boost_adaptbx::exporting::class_list< boost_adaptbx::graph_type::exports >::process(
    boost_adaptbx::graph_structure_comparison_exporter()
    );
}
