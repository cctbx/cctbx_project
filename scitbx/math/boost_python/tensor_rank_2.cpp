#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/matrix/tensor_rank_2.h>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

namespace scitbx { namespace matrix { namespace boost_python {

namespace {

  void wrap_tensor_rank_2()
  {
    using namespace boost::python;

    def("tensor_rank_2_gradient_transform",
      (sym_mat3<double>(*)(
        mat3<double> const& a,
        sym_mat3<double> const& g))
      tensor_rank_2::gradient_transform, (arg_("a"), arg_("g")));

    def("tensor_rank_2_gradient_transform_matrix",
      (af::versa<double, af::c_grid<2> >(*)(
        mat3<double> const& a))
      tensor_rank_2::gradient_transform_matrix, (arg_("a")));
  }

}}} // namespace matrix::boost_python::<anonymous>

namespace math { namespace boost_python {

  void wrap_tensor_rank_2()
  {
    matrix::boost_python::wrap_tensor_rank_2();
  }

}}} // namespace scitbx::math::boost_python
