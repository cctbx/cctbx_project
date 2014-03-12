#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <cctbx/french_wilson.h>

namespace cctbx {
namespace {

  void init_module()
  {
    using namespace boost::python;

    def("expectEFWacen",
      (double(*)
        (double,
         double)) expectEFWacen, (
      arg("eosq"),
      arg("sigesq")));

  }

} // namespace <anonymous>
} // namespace cctbx::maptbx::boost_python

BOOST_PYTHON_MODULE(cctbx_ext)
{
  cctbx::init_module();
}
