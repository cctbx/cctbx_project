#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/math/erf.h>
#include <scitbx/math/eigensystem.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>

namespace scitbx { namespace math {
namespace boost_python {

  void wrap_gaussian();

namespace {

  struct eigensystem_real_symmetric_wrappers
  {
    typedef eigensystem::real_symmetric<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("eigensystem_real_symmetric", no_init)
        .def(init<af::const_ref<double, af::c_grid<2> > const&,
                  optional<double> >())
        .def(init<scitbx::sym_mat3<double> const&,
                  optional<double> >())
        .def("vectors", &w_t::vectors)
        .def("values", &w_t::values)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;
    def("erf", (double(*)(double const&)) erf);
    def("erfc", (double(*)(double const&)) erfc);
    def("erfcx", (double(*)(double const&)) erfcx);

    eigensystem_real_symmetric_wrappers::wrap();

    wrap_gaussian();
  }

}}}} // namespace scitbx::math::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_math_ext)
{
  scitbx::math::boost_python::init_module();
}
