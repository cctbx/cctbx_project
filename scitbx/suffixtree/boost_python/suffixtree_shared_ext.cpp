#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace scitbx
{
namespace suffixtree
{
namespace
{

template< typename ValueType >
const ValueType&
dereference(const boost::shared_ptr< ValueType >& ptr)
{
  return *ptr;
}

struct shared_python_exports
{
  static void wrap()
  {
    using namespace boost::python;

    typedef boost::shared_ptr< const std::size_t > length_descriptor_type;

    class_< length_descriptor_type >( "length_descriptor", no_init )
      .def(
        "__call__",
        &dereference< typename length_descriptor_type::element_type >,
        return_value_policy< copy_const_reference >()
        )
      ;
  }
};

} // namespace anonymous
} // namespace suffixtree
} // namespace scitbx

BOOST_PYTHON_MODULE(scitbx_suffixtree_shared_ext)
{
  scitbx::suffixtree::shared_python_exports::wrap();
}
