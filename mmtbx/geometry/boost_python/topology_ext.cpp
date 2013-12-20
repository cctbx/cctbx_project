#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/mcgregor_common_subgraphs.hpp>

#include <string>
#include <utility>

namespace mmtbx
{
namespace geometry
{
namespace graph
{
namespace
{

template< typename GraphFirst, typename GraphSecond >
class mcg_python_callback_adaptor
{
private:
  const GraphFirst& m_graph1;
  const GraphSecond& m_graph2;
  boost::python::object m_callable;

public:

mcg_python_callback_adaptor(
  GraphFirst const& graph1,
  GraphSecond const& graph2,
  boost::python::object callable
  )
  : m_graph1( graph1 ), m_graph2( graph2 ), m_callable( callable )
{}

~mcg_python_callback_adaptor()
{}

template< typename CMapFirstToSecond, typename CMapSecondToFirst >
bool operator()(
  CMapFirstToSecond cmap_1_to_2,
  CMapSecondToFirst cmap_2_to_1,
  typename boost::graph_traits< GraphFirst >::vertices_size_type subgraph_size
  )
{
  typedef typename boost::graph_traits< GraphSecond >::vertex_descriptor vertex_descriptor;
  boost::python::list pairs;

  BGL_FORALL_VERTICES_T(vertex1, m_graph1, GraphFirst)
  {
    vertex_descriptor vertex2 = boost::get( cmap_1_to_2, vertex1 );

    if ( vertex2 != boost::graph_traits< GraphSecond >::null_vertex() )
    {
      pairs.append( boost::python::make_tuple( vertex1, vertex2 ) );
    }
  }

  return m_callable( pairs );
}

};

template< typename PropertyMapFirst, typename PropertyMapSecond >
class python_property_equivalent
{
private:
  const PropertyMapFirst m_property_map1;
  const PropertyMapSecond m_property_map2;
  boost::python::object m_callable;

public:

python_property_equivalent(
  PropertyMapFirst const property_map1,
  PropertyMapSecond const property_map2,
  boost::python::object callable
  )
  : m_property_map1( property_map1 ), m_property_map2( property_map2 ),
    m_callable( callable )
{}

template< typename ItemFirst, typename ItemSecond >
bool operator()(ItemFirst const item1, ItemSecond const item2)
{
  typename PropertyMapFirst::value_type property1 = get(m_property_map1, item1);
  typename PropertyMapSecond::value_type property2 = get(m_property_map2, item2);

  return m_callable( property1, property2 );
}

};

template< typename PropertyMapFirst, typename PropertyMapSecond >
python_property_equivalent< PropertyMapFirst, PropertyMapSecond >
make_python_property_equivalent(
  PropertyMapFirst const property_map1,
  PropertyMapSecond const property_map2,
  boost::python::object callable
  )
{
  return python_property_equivalent< PropertyMapFirst, PropertyMapSecond >(
    property_map1,
    property_map2,
    callable
    );
}

template< typename GraphType >
void
mcgregor_subgraphs(
  GraphType const& graph1,
  GraphType const& graph2,
  boost::python::object vertex_equality,
  boost::python::object edge_equality,
  boost::python::object callback
  )
{
  // Vertex equality
  typedef typename boost::property_map< GraphType, boost::vertex_name_t>::const_type vertex_name_map_t;
  typedef python_property_equivalent< vertex_name_map_t, vertex_name_map_t > vertex_comp_t;
  vertex_comp_t vertex_comp = make_python_property_equivalent(
    boost::get( boost::vertex_name_t(), graph1 ),
    boost::get( boost::vertex_name_t(), graph2 ),
    vertex_equality
    );

  // Edge equality
  typedef typename boost::property_map< GraphType, boost::edge_weight_t >::const_type edge_weight_map_t;
  typedef python_property_equivalent< edge_weight_map_t, edge_weight_map_t > edge_comp_t;
  edge_comp_t edge_comp = make_python_property_equivalent(
    boost::get( boost::edge_weight_t(), graph1 ),
    boost::get( boost::edge_weight_t(), graph2 ),
    edge_equality
    );

  // Callback
  mcg_python_callback_adaptor< GraphType, GraphType > mcg_callback(
    graph1,
    graph2,
    callback
    );

  boost::mcgregor_common_subgraphs_unique(
    graph1,
    graph2,
    boost::get( boost::vertex_index_t(), graph1 ),
    boost::get( boost::vertex_index_t(), graph2 ),
    edge_comp,
    vertex_comp,
    true,
    mcg_callback
    );
}

template< typename LabelType, typename WeightType >
struct python_exports
{

typedef boost::property< boost::edge_weight_t, WeightType > edge_property;
typedef boost::property< boost::vertex_name_t, LabelType > vertex_property;
typedef boost::adjacency_list<
  boost::setS,
  boost::vecS,
  boost::undirectedS,
  vertex_property,
  edge_property
  >
  graph_type;
typedef boost::graph_traits< graph_type > graph_traits_type;
typedef typename graph_traits_type::vertex_descriptor vertex_descriptor;
typedef typename graph_traits_type::edge_descriptor edge_descriptor;

static
vertex_descriptor
add_vertex(graph_type& graph, LabelType const& prop)
{
  return boost::add_vertex( vertex_property( prop ), graph );
}

static
boost::python::tuple
add_edge(
  graph_type& graph,
  vertex_descriptor u,
  vertex_descriptor v,
  WeightType const& prop
  )
{
  std::pair< edge_descriptor, bool > res =
    boost::add_edge( u, v, edge_property( prop ), graph );
  return boost::python::make_tuple( res.first, res.second );
}


static void wrap()
{
  typedef boost::property< boost::edge_weight_t, WeightType > edge_property;
  typedef boost::property< boost::vertex_name_t, LabelType > vertex_property;
  typedef boost::adjacency_list<
    boost::setS,
    boost::vecS,
    boost::undirectedS,
    vertex_property,
    edge_property
    >
    graph_type;

  using namespace boost::python;

  class_< edge_descriptor >( "edge_descriptor", no_init )
    ;

  class_< graph_type >( "graph", no_init )
    .def( init<>() )
    .def( "add_vertex", &add_vertex, arg( "label" ) )
    .def( "add_edge", &add_edge, ( arg( "vertex1" ), arg( "vertex2" ), arg( "weight" ) ) )
    ;

  void (*mcsg_unique_ptr)(
    graph_type const&,
    graph_type const&,
    boost::python::object,
    boost::python::object,
    boost::python::object
    ) =
    &mcgregor_subgraphs< graph_type >;

  def(
    "mcgregor_common_subgraphs_unique",
    mcsg_unique_ptr,
    ( arg( "graph1" ), arg( "graph2" ),
      arg( "vertex_equality" ), arg( "edge_equality" ), arg( "callback" ) )
    );
}

};

} // namespace <anonymnous>
} // namespace graph
} // namespace graph
} // namespace mmtbx

BOOST_PYTHON_MODULE(mmtbx_geometry_topology_ext)
{
  mmtbx::geometry::graph::python_exports< boost::python::object, boost::python::object >::wrap();
}
