#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/array_family/linear_interpolation.h>
#include <boost/python/def.hpp>
#include <boost/python/overloads.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {
  BOOST_PYTHON_FUNCTION_OVERLOADS(
    linear_interpolation_overloads, linear_interpolation, 3, 4)
}

  void wrap_flex_linear_interpolation()
  {
    using namespace boost::python;

    def("linear_interpolation",
      (float(*)(
        const_ref<float> const&,
        const_ref<float> const&,
        float const&,
        float const&))
      linear_interpolation, linear_interpolation_overloads());
    def("linear_interpolation",
      (double(*)(
        const_ref<double> const&,
        const_ref<double> const&,
        double const&,
        double const&))
      linear_interpolation, linear_interpolation_overloads());
    def("linear_interpolation",
      (shared<float>(*)(
        const_ref<float> const&,
        const_ref<float> const&,
        const_ref<float> const&,
        float const&))
      linear_interpolation, linear_interpolation_overloads());
    def("linear_interpolation",
      (shared<double>(*)(
        const_ref<double> const&,
        const_ref<double> const&,
        const_ref<double> const&,
        double const&))
      linear_interpolation, linear_interpolation_overloads());
  }

}}} // namespace scitbx::af::boost_python
