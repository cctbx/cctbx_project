#ifndef BOOST_ADAPTBX_GRAPH_EXPORT_ADAPTOR_H
#define BOOST_ADAPTBX_GRAPH_EXPORT_ADAPTOR_H

#include <utility>

namespace boost_adaptbx
{
namespace graph_export_adaptor
{

template< typename Graph >
struct no_export
{
  static void process()
  {}
};

template< typename VertexDescriptor >
struct vertex_descriptor_converter
{
  typedef VertexDescriptor type;
  typedef type result_type;

  static type forward(VertexDescriptor vd)
  {
    return vd;
  }

  static VertexDescriptor backward(type vd)
  {
    return vd;
  }

  inline type operator ()(VertexDescriptor const& vd) const
  {
    return vd;
  }
};

template<>
struct vertex_descriptor_converter< void* >
{
  typedef size_t type;
  typedef type result_type;

  static type forward(void* vd)
  {
    return reinterpret_cast< type >( vd );
  }

  static void* backward(type vd)
  {
    return reinterpret_cast< void* >( vd );
  }

  inline type operator ()(void* vd) const
  {
    return forward( vd );
  }
};

template< typename VertexRange, typename OutputIterator >
void
fill_vertex_index_map(VertexRange range, OutputIterator out)
{
  std::size_t index = 0;

  for (;range.first != range.second; ++range.first)
  {
    *(out++) = std::make_pair( *range.first, index++ );
  }
}

} // namespace graph_export_adaptor
} // namespace boost_adaptbx

#endif // BOOST_ADAPTBX_GRAPH_EXPORT_ADAPTOR_H
