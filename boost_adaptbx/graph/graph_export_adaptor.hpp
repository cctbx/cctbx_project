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
  typedef VertexDescriptor argument_type;
  typedef VertexDescriptor type;
  typedef type result_type;

  static result_type forward(argument_type vd)
  {
    return vd;
  }

  static argument_type backward(result_type vd)
  {
    return vd;
  }

  inline result_type operator ()(argument_type const& vd) const
  {
    return vd;
  }
};

template<>
struct vertex_descriptor_converter< void* >
{
  typedef void* argument_type;
  typedef size_t type;
  typedef type result_type;

  static result_type forward(argument_type vd)
  {
    return reinterpret_cast< result_type >( vd );
  }

  static argument_type backward(result_type vd)
  {
    return reinterpret_cast< void* >( vd );
  }

  inline result_type operator ()(argument_type vd) const
  {
    return forward( vd );
  }
};

template< typename VertexDescriptor >
struct vertex_descriptor_backconverter
{
  typedef vertex_descriptor_converter< VertexDescriptor > converter;
  typedef typename converter::result_type argument_type;
  typedef typename converter::argument_type result_type;

  inline result_type operator ()(argument_type vd) const
  {
    return converter::backward( vd );
  }
};

} // namespace graph_export_adaptor
} // namespace boost_adaptbx

#endif // BOOST_ADAPTBX_GRAPH_EXPORT_ADAPTOR_H
