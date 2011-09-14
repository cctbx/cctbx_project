#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/serialization/single_buffered.h>
#include <scitbx/matrix/transpose_multiply.h>
#include <scitbx/math/utils.h>
#include <scitbx/vec2.h>
#include <scitbx/mat2.h>
#include <boost/python/make_constructor.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_arg.hpp>
#include "flex_helpers.h"

namespace scitbx { namespace serialization { namespace single_buffered {

  inline
  char* to_string(char* start, vec2<double> const& value)
  {
    return
      to_string(to_string(start, value[0]), value[1]);
  }

  template <>
  struct from_string<vec2<double> >
  {
    from_string(const char* start)
    {
      end = start;
      for(std::size_t i=0;i<2;i++) {
        from_string<double> proxy(end);
        value[i] = proxy.value;
        end = proxy.end;
      }
    }

    vec2<double> value;
    const char* end;
  };

}}} // namespace scitbx::serialization::single_buffered

#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>

namespace scitbx { namespace af {
namespace {

  flex<vec2<double> >::type*
  join(
    af::const_ref<double> const& x,
    af::const_ref<double> const& y)
  {
    SCITBX_ASSERT(y.size() == x.size());
    af::shared<vec2<double> > result((af::reserve(x.size())));
    for(std::size_t i=0;i<x.size();i++) {
      result.push_back(vec2<double>(x[i],y[i]));
    }
    return new flex<vec2<double> >::type(result, result.size());
  }

  flex<vec2<double> >::type*
  from_double(
    af::const_ref<double> const& x)
  {
    SCITBX_ASSERT(x.size() % 2 == 0);
    std::size_t result_size = x.size() / 2;
    af::shared<vec2<double> > result((af::reserve(result_size)));
    const double* d = x.begin();
    for(std::size_t i=0;i<result_size;i++) {
      result.push_back(vec2<double>(d));
      d += 2;
    }
    return new flex<vec2<double> >::type(result, result.size());
  }

  flex_double
  as_double(flex<vec2<double> >::type const& a)
  {
    SCITBX_ASSERT(a.accessor().is_trivial_1d());
    flex_double result(a.size()*2, init_functor_null<double>());
    double* r = result.begin();
    const_ref<vec2<double> > a_ref = a.const_ref().as_1d();
    for(std::size_t i=0;i<a_ref.size();i++) {
      for(std::size_t j=0;j<2;j++) {
        *r++ = a_ref[i][j];
      }
    }
    return result;
  }

  vec2<double>
  vec2_min(flex<vec2<double> >::type const& a)
  {
    SCITBX_ASSERT(!a.accessor().is_padded());
    vec2<double> result(0,0);
    af::const_ref<vec2<double>, af::flex_grid<> > a_ref = a.const_ref();
    if (a_ref.size() > 0) {
      result = a_ref[0];
      for(std::size_t i=1;i<a_ref.size();i++) {
        result.each_update_min(a_ref[i]);
      }
    }
    return result;
  }

  vec2<double>
  vec2_max(flex<vec2<double> >::type const& a)
  {
    SCITBX_ASSERT(!a.accessor().is_padded());
    vec2<double> result(0,0);
    af::const_ref<vec2<double>, af::flex_grid<> > a_ref = a.const_ref();
    if (a_ref.size() > 0) {
      result = a_ref[0];
      for(std::size_t i=1;i<a_ref.size();i++) {
        result.each_update_max(a_ref[i]);
      }
    }
    return result;
  }

  vec2<double>
  mean_weighted_a_a(
    af::const_ref<vec2<double> > const& self,
    af::const_ref<double> const& weights)
  {
    return af::mean_weighted(self, weights);
  }

  af::shared<vec2<double> >
  mul_a_scalar(
    af::const_ref<vec2<double> > const& a,
    double f)
  {
    af::shared<vec2<double> > result((af::reserve(a.size())));
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(a[i] * f);
    }
    return result;
  }

  void
  imul_a_scalar(
    af::ref<vec2<double> > const& a,
    double f)
  {
    for(std::size_t i=0;i<a.size();i++) a[i] *= f;
  }

  af::shared<vec2<double> >
  div_a_as(
    af::ref<vec2<double> > const& lhs,
    af::ref<double> const& rhs)
  {
    SCITBX_ASSERT(lhs.size() == rhs.size());
    af::shared<vec2<double> > result((af::reserve(lhs.size())));
    for(std::size_t i=0;i<lhs.size();i++) {
      SCITBX_ASSERT(rhs[i] != 0);
      result.push_back(lhs[i] / rhs[i]);
    }
    return result;
  }

  af::shared<vec2<double> >
  mul_a_mat2(
    af::const_ref<vec2<double> > const& a,
    mat2<double> const& m)
  {
    af::shared<vec2<double> > result((af::reserve(a.size())));
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(a[i] * m);
    }
    return result;
  }

  af::shared<vec2<double> >
  rmul_a_mat2(
    af::const_ref<vec2<double> > const& a,
    mat2<double> const& m)
  {
    mat2<double> m_transposed = m.transpose();
    af::shared<vec2<double> > result((af::reserve(a.size())));
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(a[i] * m_transposed);
    }
    return result;
  }

  af::shared<double>
  dot_a_a(
    af::const_ref<vec2<double> > const& lhs,
    af::const_ref<vec2<double> > const& rhs)
  {
    SCITBX_ASSERT(lhs.size() == rhs.size());
    af::shared<double> result((af::reserve(lhs.size())));
    for(std::size_t i=0;i<lhs.size();i++) {
      result.push_back(lhs[i] * rhs[i]);
    }
    return result;
  }

  af::shared<double>
  dot_a(
    af::const_ref<vec2<double> > const& lhs)
  {
    af::shared<double> result((af::reserve(lhs.size())));
    for(std::size_t i=0;i<lhs.size();i++) {
      result.push_back(lhs[i] * lhs[i]);
    }
    return result;
  }

  double
  sum_sq_(
    af::const_ref<vec2<double> > const& self)
  {
    double result = 0;
    for(std::size_t i=0;i<self.size();i++) {
      result += self[i] * self[i];
    }
    return result;
  }

  double
  norm_(
    af::const_ref<vec2<double> > const& self)
  {
    return std::sqrt(sum_sq_(self));
  }

  double
  max_distance(
    af::const_ref<vec2<double> > const& lhs,
    af::const_ref<vec2<double> > const& rhs)
  {
    SCITBX_ASSERT(lhs.size() == rhs.size());
    if (lhs.size() == 0) return 0;
    double max_length_sq = 0;
    for(std::size_t i=0;i<lhs.size();i++) {
      math::update_max(max_length_sq, (lhs[i]-rhs[i]).length_sq());
    }
    return std::sqrt(max_length_sq);
  }

  double
  rms_difference(
    af::const_ref<vec2<double> > const& lhs,
    af::const_ref<vec2<double> > const& rhs)
  {
    SCITBX_ASSERT(lhs.size() == rhs.size());
    if (lhs.size() == 0) return 0;
    double sum_length_sq = 0;
    for(std::size_t i=0;i<lhs.size();i++) {
      sum_length_sq += (lhs[i]-rhs[i]).length_sq();
    }
    return std::sqrt(sum_length_sq / lhs.size());
  }

  double
  rms_length(
    af::const_ref<vec2<double> > const& lhs)
  {
    if (lhs.size() == 0) return 0;
    double sum_length_sq = 0;
    for(std::size_t i=0;i<lhs.size();i++) {
      sum_length_sq += lhs[i].length_sq();
    }
    return std::sqrt(sum_length_sq / lhs.size());
  }

} // namespace <anonymous>

namespace boost_python {

  template <>
  struct flex_default_element<vec2<double> >
  {
    static vec2<double>
    get() { return vec2<double>(0,0); }
  };

  void wrap_flex_vec2_double()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef flex_wrapper<vec2<double> > f_w;
    f_w::plain("vec2_double")
      .def_pickle(flex_pickle_single_buffered<vec2<double>,
        2*pickle_size_per_element<double>::value>())
      .def("__init__", make_constructor(join))
      .def("__init__", make_constructor(from_double))
      .def("as_double", as_double)
      .def("add_selected",
        (object(*)(
          object const&,
          af::const_ref<std::size_t> const&,
          af::const_ref<vec2<double> > const&)) add_selected_unsigned_a,
        (arg("self"), arg("indices"), arg("values")))
      .def("min", vec2_min)
      .def("max", vec2_max)
      .def("sum", f_w::sum_a)
      .def("mean", f_w::mean_a)
      .def("mean_weighted", mean_weighted_a_a, (arg("self"), arg("weights")))
      .def("__add__", f_w::add_a_s)
      .def("__add__", f_w::add_a_a)
      .def("__iadd__", f_w::iadd_a_s)
      .def("__iadd__", f_w::iadd_a_a)
      .def("__sub__", f_w::sub_a_s)
      .def("__sub__", f_w::sub_a_a)
      .def("__isub__", f_w::isub_a_s)
      .def("__mul__", mul_a_scalar)
      .def("__rmul__", mul_a_scalar)
      .def("__imul__", imul_a_scalar, return_self<>())
      .def("__div__", div_a_as)
      .def("__truediv__", div_a_as)
      .def("__mul__", mul_a_mat2)
      .def("__rmul__", rmul_a_mat2)
      .def("dot", dot_a_a)
      .def("dot", dot_a)
      .def("transpose_multiply",
        (mat2<double>(*)(
          af::const_ref<vec2<double> > const&,
          af::const_ref<vec2<double> > const&)) matrix::transpose_multiply)
      .def("sum_sq", sum_sq_)
      .def("norm", norm_)
      .def("max_distance", max_distance)
      .def("rms_difference", rms_difference)
      .def("rms_length", rms_length)
    ;
  }

}}} // namespace scitbx::af::boost_python
