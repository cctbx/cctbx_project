#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <scitbx/array_family/boost_python/flex_wrapper_complex.h>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  versa<std::complex<double>, flex_grid<> >
  mul_ac_ar(
    versa<std::complex<double>, flex_grid<> > const& a1,
    versa<             double , flex_grid<> > const& a2)
  {
    return a1 * a2;
  }

} // namespace <anonymous>

  void wrap_flex_complex_double()
  {
    flex_wrapper<std::complex<double> >::numeric_common(
      "complex_double", boost::python::scope())
      .def_pickle(flex_pickle_single_buffered<std::complex<double> >())
      .def("__mul__", mul_ac_ar)
      .def("__rmul__", mul_ac_ar)
    ;
    flex_wrapper_complex_functions<double>::wrap(boost::python::scope());
  }

}}} // namespace scitbx::af::boost_python
