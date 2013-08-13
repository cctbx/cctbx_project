#ifndef SUFFIXTREE_PYTHON_TREE_EXPORTS_HPP_
#define SUFFIXTREE_PYTHON_TREE_EXPORTS_HPP_

#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>

#include <scitbx/suffixtree/tree.hpp>
#include <scitbx/suffixtree/builder.hpp>

namespace scitbx
{
namespace suffixtree
{
namespace python
{

template<
  typename Word,
  typename SuffixLabel,
  template< typename, typename > class NodeAdaptor
  >
struct tree_exports
{
  typedef Tree< Word, SuffixLabel, NodeAdaptor > tree_type;

  static void wrap()
  {
    using namespace boost::python;

    class_< tree_type >( "tree", no_init )
      .def( init<>() )
      .add_property( "root", &tree_type::root )
      .add_property(
        "word",
        make_function(
          &tree_type::word,
          return_value_policy< copy_const_reference >()
          )
        )
      .add_property( "in_construction", &tree_type::in_construction )
      ;
  }
};

template< typename Tree >
struct ukkonen_builder_exports
{
  typedef Tree tree_type;
  typedef builder::Ukkonen< Tree > builder_type;

  static void wrap()
  {
    using namespace boost::python;

    class_< builder_type >( "ukkonen", no_init )
      .def( init< tree_type const& >( arg( "tree" ) ) )
      .add_property( "is_attached", &builder_type::is_attached )
      .add_property( "is_valid", &builder_type::is_valid )
      .def( "append", &builder_type::push_back, arg( "glyph" ) )
      .def( "detach", &builder_type::detach )
      ;
  }
};

} // namespace python
} // namespace suffixtree
} // namespace scitbx

#endif // SUFFIXTREE_PYTHON_TREE_EXPORTS_HPP_
