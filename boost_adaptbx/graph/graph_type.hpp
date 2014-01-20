#ifndef BOOST_ADAPTBX_GRAPH_TYPE_H
#define BOOST_ADAPTBX_GRAPH_TYPE_H

#include <boost/python/object.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/string.hpp>

namespace boost_adaptbx
{
namespace graph_type
{

typedef boost::property< boost::edge_weight_t, boost::python::object > edge_property;
typedef boost::property< boost::vertex_name_t, boost::python::object > vertex_property;

typedef boost::adjacency_list<
  boost::setS,
  boost::listS,
  boost::undirectedS,
  vertex_property,
  edge_property
  >
  adjacency_list_undirected_listS_setS_type;

typedef boost::adjacency_list<
  boost::setS,
  boost::vecS,
  boost::undirectedS,
  vertex_property,
  edge_property
  >
  adjacency_list_undirected_vecS_setS_type;

typedef boost::adjacency_list<
  boost::vecS,
  boost::listS,
  boost::undirectedS,
  vertex_property,
  edge_property
  >
  adjacency_list_undirected_listS_vecS_type;

typedef boost::adjacency_list<
  boost::vecS,
  boost::vecS,
  boost::undirectedS,
  vertex_property,
  edge_property
  >
  adjacency_list_undirected_vecS_vecS_type;

typedef boost::adjacency_list<
  boost::vecS,
  boost::vecS,
  boost::directedS,
  vertex_property,
  edge_property
  >
  adjacency_list_directed_vecS_vecS_type;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmultichar"
typedef boost::mpl::vector<
  boost::mpl::pair<
    adjacency_list_undirected_listS_setS_type,
    boost::mpl::string< 'adj', 'list', '_', 'und', 'ir_', 'list', '_', 'set' >
    >,
  boost::mpl::pair<
    adjacency_list_undirected_vecS_setS_type,
    boost::mpl::string< 'adj', 'list', '_', 'und', 'ir_', 'vect', '_', 'set' >
    >,
  boost::mpl::pair<
    adjacency_list_undirected_listS_vecS_type,
    boost::mpl::string< 'adj', 'list', '_', 'und', 'ir_', 'list', '_', 'vect' >
    >,
  boost::mpl::pair<
    adjacency_list_undirected_vecS_vecS_type,
    boost::mpl::string< 'adj', 'list', '_', 'und', 'ir_', 'vect', '_', 'vect' >
    >,
  boost::mpl::pair<
    adjacency_list_directed_vecS_vecS_type,
    boost::mpl::string< 'adj', 'list', '_', 'dir_', 'vect', '_', 'vect' >
    >
  > exports;
#pragma clang diagnostic pop

} // namespace graph_type
} // namespace boost_adaptbx

#endif // BOOST_ADAPTBX_GRAPH_TYPE_H
