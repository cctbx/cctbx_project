#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_value_policy.hpp>

#include <scitbx/array_family/shared.h>
#include <scitbx/optional_copy/conversion.h>

namespace scitbx { namespace boost_python {

namespace optional_copy_conversions {

  struct test {
    optional_copy<af::shared<unsigned> > a, b;

    test(af::shared<unsigned> const &array)
      : a(array)
    {}

    static void wrap() {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      class_<test>("test", no_init)
        .def(init<af::shared<unsigned> const&>())
        .add_property("a", make_getter(&test::a, rbv))
        .add_property("b", make_getter(&test::b, rbv))
        ;
    }
  };

}}} // scitbx::boost_python::optional_copy_conversions


BOOST_PYTHON_MODULE(scitbx_optional_copy_ext)
{
  using namespace scitbx;
  boost_python::optional_copy_conversions::to_python<af::shared<unsigned> >();
  boost_python::optional_copy_conversions::test::wrap();
}
