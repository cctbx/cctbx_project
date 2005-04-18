#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/matrix/tensor_rank_2.h>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

namespace scitbx { namespace matrix { namespace boost_python {

namespace {

  template <typename NumTypeG>
  struct gradient_average_cache_wrappers
  {
    typedef tensor_rank_2::gradient_average_cache<NumTypeG> w_t;

    static void
    wrap(const char* python_name)
    {
      using namespace boost::python;
      class_<w_t>(python_name, no_init)
        .def(init<>())
        .def("accumulate", &w_t::accumulate, (arg_("a")))
        .def("average",
          (sym_mat3<double>(w_t::*)(
            sym_mat3<double> const&,
            const double&) const)
              &w_t::average, (arg_("g"), arg_("denominator")))
      ;
    }
  };

  void wrap_tensor_rank_2()
  {
    using namespace boost::python;

    def("tensor_rank_2_gradient_transform",
      (sym_mat3<double>(*)(
        mat3<double> const& a,
        sym_mat3<double> const& g))
      tensor_rank_2::gradient_transform, (arg_("a"), arg_("g")));

    gradient_average_cache_wrappers<int>::wrap(
      "tensor_rank_2_gradient_average_cache_int");
  }

}}} // namespace matrix::boost_python::<anonymous>

namespace math { namespace boost_python {

  void wrap_tensor_rank_2()
  {
    matrix::boost_python::wrap_tensor_rank_2();
  }

}}} // namespace scitbx::math::boost_python
