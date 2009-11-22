#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/array_family/linear_interpolation.h>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  template <typename ElementType>
  void
  linear_interpolation_wrapper()
  {
    using namespace boost::python;
    def("linear_interpolation",
      (ElementType(*)(
        const_ref<ElementType> const&,
        const_ref<ElementType> const&,
        ElementType const&,
        ElementType const&)) linear_interpolation, (
          arg("table_x"),
          arg("table_y"),
          arg("x"),
          arg("tolerance")=1e-6));
    def("linear_interpolation",
      (shared<ElementType>(*)(
        const_ref<ElementType> const&,
        const_ref<ElementType> const&,
        const_ref<ElementType> const&,
        ElementType const&)) linear_interpolation, (
          arg("table_x"),
          arg("table_y"),
          arg("x"),
          arg("tolerance")=1e-6));
  }

} // namespace <anonymous>

  void wrap_flex_linear_interpolation()
  {
    linear_interpolation_wrapper<float>();
    linear_interpolation_wrapper<double>();
  }

}}} // namespace scitbx::af::boost_python
