#ifndef SUFFIXTREE_PYTHON_EDGE_EXPORTS_HPP_
#define SUFFIXTREE_PYTHON_EDGE_EXPORTS_HPP_

#include <boost/python/class.hpp>
#include <boost/python/object.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/list.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/def.hpp>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/functional.hpp>

#include <scitbx/suffixtree/boost_python/iterator_wrapper.hpp>
#include <scitbx/suffixtree/edge.hpp>
#include <scitbx/suffixtree/iterator.hpp>

#include <iostream>

namespace scitbx
{
namespace suffixtree
{
namespace python
{

template<
  typename Glyph,
  typename Index,
  typename WordLength,
  typename SuffixLabel,
  template< typename, typename > class NodeAdapter
  >
struct edge_exports
{
  typedef edge::Edge< Glyph, Index, WordLength, SuffixLabel, NodeAdapter > edge_type;
  typedef typename edge_type::ptr_type ptr_type;
  typedef typename edge_type::weak_ptr_type weak_ptr_type;
  typedef typename edge_type::index_type index_type;
  typedef typename edge_type::suffix_label_type suffix_label_type;
  typedef typename edge_type::node_type node_type;
  typedef typename edge_type::word_length_type word_length_type;
  typedef typename edge_type::glyph_type glyph_type;

  typedef iterator::PreOrder< edge_type > preorder_iterator;
  typedef iterator::PostOrder< edge_type > postorder_iterator;

  typedef python_iterator< preorder_iterator > preorder_python_iterator;
  typedef python_iterator< postorder_iterator > postorder_python_iterator;

  static void throw_python_key_error()
  {
    PyErr_SetString(PyExc_KeyError, "Key not found");
    boost::python::throw_error_already_set();
  }

  static index_type get_start(ptr_type const& edge_ptr)
  {
    return edge_ptr->get_start();
  }

  static void set_start(ptr_type const& edge_ptr, index_type const& start)
  {
    edge_ptr->set_start( start );
  }

  static index_type get_stop(ptr_type const& edge_ptr)
  {
    return edge_ptr->get_stop();
  }

  static suffix_label_type get_suffix_label(ptr_type const& edge_ptr)
  {
    return edge_ptr->label();
  }

  static bool node_contains(ptr_type const& edge_ptr, glyph_type const& key)
  {
    try
    {
      return edge_ptr->find( key ) != edge_ptr->end();
    }
    catch ( bad_edge_type& e )
    {
      return false;
    }
  }

  static bool node_empty(ptr_type const& edge_ptr)
  {
    return edge_ptr->empty();
  }

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreturn-type"

  static ptr_type node_get_item(ptr_type const& edge_ptr, glyph_type const& key)
  {
    try
    {
      return edge_ptr->get_child_with_label( key );
    }

    catch ( nonexistent& e )
    {
      throw_python_key_error();
    }

    catch ( bad_edge_type& e )
    {
      throw_python_key_error();
    }
  }

#pragma clang diagnostic pop

  static void node_set_item(
    ptr_type const& edge_ptr,
    glyph_type const& key,
    ptr_type const& value
    )
  {
    edge_ptr->attach_child( value, key );
  }

  static boost::python::list node_keys(ptr_type const& edge_ptr)
  {
    boost::python::list result;

    try
    {
      for(
        typename edge_type::const_iterator it = edge_ptr->begin();
        it != edge_ptr->end();
        ++it
        )
      {
        const glyph_type& k = it->first;
        result.append( k );
      }
    }
    catch ( bad_edge_type& e )
    {}

    return result;
  }

  static boost::python::list node_values(ptr_type const& edge_ptr)
  {
    boost::python::list result;

    try
    {
      for(
        typename edge_type::const_iterator it = edge_ptr->begin();
        it != edge_ptr->end();
        ++it
        )
      {
        ptr_type const& v = it->second;
        result.append( v );
      }
    }
    catch ( bad_edge_type& e )
    {}

    return result;
  }

  static boost::python::object const get_parent(ptr_type const& edge_ptr)
  {
    try
    {
      return boost::python::object( edge_ptr->get_parent() );
    }
    catch ( unavailable& e )
    {
      return boost::python::object();
    }
  }

  static void set_parent(ptr_type const& edge_ptr, ptr_type const& parent)
  {
    edge_ptr->parent() = weak_ptr_type( parent );
  }

  static boost::python::object const get_suffix(ptr_type const& edge_ptr)
  {
    try
    {
      return boost::python::object( edge_ptr->get_suffix() );
    }
    catch ( unavailable& e )
    {
      return boost::python::object();
    }
  }

  static void set_suffix(ptr_type const& edge_ptr, ptr_type const& suffix)
  {
    edge_ptr->suffix() = suffix;
  }

  static bool is_root(ptr_type const& edge_ptr)
  {
    return edge_ptr->is_root();
  }

  static bool is_leaf(ptr_type const& edge_ptr)
  {
    return edge_ptr->is_leaf();
  }

  static std::size_t calculate_hash(ptr_type const& edge_ptr)
  {
    return boost::hash_value( edge_ptr );
  }

  static preorder_python_iterator get_preorder_range(ptr_type const& root)
  {
    return preorder_python_iterator(
      preorder_iterator::begin( root ),
      preorder_iterator::end( root )
      );
  }


  static postorder_python_iterator get_postorder_range(ptr_type const& root)
  {
    return postorder_python_iterator(
      postorder_iterator::begin( root ),
      postorder_iterator::end( root )
      );
  }
/*
  static bool operator_eq_nonconst(ptr_type const& lhs, ptr_type const& rhs)
  {
    return lhs == rhs;
  }

  static bool operator_eq_mixed(ptr_type const& lhs, const_ptr_type const& rhs)
  {
    return lhs == rhs;
  }

  static bool operator_eq_const(const_ptr_type const& lhs, const_ptr_type const& rhs)
  {
    return lhs == rhs;
  }
*/
  static void wrap()
  {
    using namespace boost::python;

    preorder_python_iterator::wrap( "preorder_iteration_range" );
    postorder_python_iterator::wrap( "postorder_iteration_range" );

    class_< ptr_type >( "edge", no_init )
      .def(
        "root",
        edge_type::root
        )
      .staticmethod( "root" )
      .def(
        "branch",
        edge_type::branch,
        ( arg( "start" ), arg( "stop" ) )
        )
      .staticmethod( "branch" )
      .def(
        "leaf",
        edge_type::leaf,
        ( arg( "start" ), arg( "length" ), arg( "label" ) )
        )
      .staticmethod( "leaf" )
      .add_property( "start",get_start, set_start )
      .add_property( "stop", get_stop )
      .add_property( "label", get_suffix_label )
      .add_property( "parent", get_parent, set_parent )
      .add_property( "suffix", get_suffix, set_suffix )
      .def( "is_root", is_root )
      .def( "is_leaf", is_leaf )
      .def( "is_empty", node_empty )
      .def( "__contains__", node_contains, arg( "key" ) )
      .def( "__getitem__", node_get_item, arg( "key" ) )
      .def( "__setitem__", node_set_item, ( arg( "key" ), arg( "value" ) ) )
      .def( "keys", node_keys )
      .def( "values", node_values )
      .def( "preorder_iteration", get_preorder_range )
      .def( "postorder_iteration", get_postorder_range )
      .def( self == self )
      .def( self != self )
      .def( "__hash__", calculate_hash )
      ;
  }
};

} // namespace python
} // namespace suffixtree
} // namespace scitbx

#endif // SUFFIXTREE_PYTHON_EDGE_EXPORTS_HPP_
