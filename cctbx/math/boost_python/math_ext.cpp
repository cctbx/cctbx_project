#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <cctbx/math/cos_sin_table.h>
#include <boost/optional.hpp>

namespace cctbx { namespace math { namespace boost_python {
namespace {

  struct cos_sin_table_wrappers
  {
    typedef cos_sin_table<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("cos_sin_table", no_init)
        .def(init<int, optional<bool> >())
        .def("n_points", &w_t::n_points)
        .def("get", &w_t::get)
      ;
    }
  };

  void init_module()
  {
    cos_sin_table_wrappers::wrap();
  }

} // namespace <anonymous>
}}} // namespace cctbx::math::boost_python

BOOST_PYTHON_MODULE(cctbx_math_ext)
{
  cctbx::math::boost_python::init_module();
}
