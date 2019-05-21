#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/reference_existing_object.hpp>
#include <boost/python/return_value_policy.hpp>
#include <cctbx/eltbx/chemical_elements.h>
#include <scitbx/boost_python/array_as_list.h>

namespace cctbx { namespace eltbx { namespace chemical_elements {
namespace boost_python {

namespace {

  boost::python::object
  proper_caps_list()
  {
    boost::python::object result = scitbx::boost_python::array_as_list(
      proper_caps, sizeof(proper_caps) / sizeof(char*) - 1);
    return result;
  }

  boost::python::object
  proper_upper_list()
  {
    boost::python::object result = scitbx::boost_python::array_as_list(
      proper_upper, sizeof(proper_upper) / sizeof(char*) - 1);
    return result;
  }

  void init_module()
  {
    using namespace boost::python;
    typedef return_value_policy<reference_existing_object> reo;
    def("proper_caps_list", proper_caps_list);
    def("proper_upper_list", proper_upper_list);
    def("proper_caps_set", proper_caps_set, reo());
    def("proper_upper_set", proper_upper_set, reo());
    def("proper_and_isotopes_upper_set", proper_and_isotopes_upper_set, reo());
  }

} // namespace <anonymous>
}}}} // namespace cctbx::eltbx::chemical_elements::boost_python

BOOST_PYTHON_MODULE(cctbx_eltbx_chemical_elements_ext)
{
  cctbx::eltbx::chemical_elements::boost_python::init_module();
}
