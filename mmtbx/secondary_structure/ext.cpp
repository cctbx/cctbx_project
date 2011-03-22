#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <mmtbx/secondary_structure/identify.hpp>

namespace mmtbx { namespace secondary_structure {
namespace {
  void wrap_identification ()
  {
    using namespace boost::python;
    def("delta_distance_squared", delta_distance_squared, (
      arg("base_1"),
      arg("base_2"),
      arg("name_1"),
      arg("name_2"),
      arg("atom_names"),
      arg("sites_cart"),
      arg("distance_ideal")));
  }
}

namespace boost_python {
  void wrap_sec_str ()
  {
    wrap_identification();
  }
}
}}

BOOST_PYTHON_MODULE(mmtbx_secondary_structure_ext)
{
  mmtbx::secondary_structure::boost_python::wrap_sec_str();
}
