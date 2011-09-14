#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/serialization/single_buffered.h>
#include <scitbx/matrix/transpose_multiply.h>
#include <scitbx/math/utils.h>
#include <boost/python/make_constructor.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_arg.hpp>
#include <boost/format.hpp>
#include "flex_helpers.h"

namespace scitbx { namespace serialization { namespace single_buffered {

  inline
  char* to_string(char* start, vec3<double> const& value)
  {
    return
      to_string(to_string(to_string(start, value[0]), value[1]), value[2]);
  }

  template <>
  struct from_string<vec3<double> >
  {
    from_string(const char* start)
    {
      end = start;
      for(std::size_t i=0;i<3;i++) {
        from_string<double> proxy(end);
        value[i] = proxy.value;
        end = proxy.end;
      }
    }

    vec3<double> value;
    const char* end;
  };

}}} // namespace scitbx::serialization::single_buffered

#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>

namespace scitbx { namespace af {
namespace {

  flex<vec3<double> >::type*
  join(
    af::const_ref<double> const& x,
    af::const_ref<double> const& y,
    af::const_ref<double> const& z)
  {
    SCITBX_ASSERT(y.size() == x.size());
    SCITBX_ASSERT(z.size() == x.size());
    af::shared<vec3<double> > result((af::reserve(x.size())));
    for(std::size_t i=0;i<x.size();i++) {
      result.push_back(vec3<double>(x[i],y[i],z[i]));
    }
    return new flex<vec3<double> >::type(result, result.size());
  }

  flex<vec3<double> >::type*
  from_double(
    af::const_ref<double> const& x)
  {
    SCITBX_ASSERT(x.size() % 3 == 0);
    std::size_t result_size = x.size() / 3;
    af::shared<vec3<double> > result((af::reserve(result_size)));
    const double* d = x.begin();
    for(std::size_t i=0;i<result_size;i++) {
      result.push_back(vec3<double>(d));
      d += 3;
    }
    return new flex<vec3<double> >::type(result, result.size());
  }

  boost::python::tuple
  part_names()
  {
    return boost::python::make_tuple("x", "y", "z");
  }

  boost::python::tuple
  parts(
    versa<vec3<double>, flex_grid<> > const& O)
  {
    tiny<versa<double, flex_grid<> >, 3> result;
    std::size_t n = O.size();
    for(std::size_t i=0;i<3;i++) {
      result[i].resize(O.accessor());
      for(std::size_t j=0;j<n;j++) {
        result[i][j] = O[j][i];
      }
    }
    return boost::python::make_tuple(result[0], result[1], result[2]);
  }

  flex_double
  as_double(flex<vec3<double> >::type const& a)
  {
    SCITBX_ASSERT(a.accessor().is_trivial_1d());
    flex_double result(a.size()*3, init_functor_null<double>());
    double* r = result.begin();
    const_ref<vec3<double> > a_ref = a.const_ref().as_1d();
    for(std::size_t i=0;i<a_ref.size();i++) {
      for(std::size_t j=0;j<3;j++) {
        *r++ = a_ref[i][j];
      }
    }
    return result;
  }

  vec3<double>
  vec3_min(flex<vec3<double> >::type const& a)
  {
    SCITBX_ASSERT(!a.accessor().is_padded());
    vec3<double> result(0,0,0);
    af::const_ref<vec3<double>, af::flex_grid<> > a_ref = a.const_ref();
    if (a_ref.size() > 0) {
      result = a_ref[0];
      for(std::size_t i=1;i<a_ref.size();i++) {
        result.each_update_min(a_ref[i]);
      }
    }
    return result;
  }

  vec3<double>
  vec3_max(flex<vec3<double> >::type const& a)
  {
    SCITBX_ASSERT(!a.accessor().is_padded());
    vec3<double> result(0,0,0);
    af::const_ref<vec3<double>, af::flex_grid<> > a_ref = a.const_ref();
    if (a_ref.size() > 0) {
      result = a_ref[0];
      for(std::size_t i=1;i<a_ref.size();i++) {
        result.each_update_max(a_ref[i]);
      }
    }
    return result;
  }

  vec3<double>
  mean_weighted_a_a(
    af::const_ref<vec3<double> > const& self,
    af::const_ref<double> const& weights)
  {
    return af::mean_weighted(self, weights);
  }

  af::shared<vec3<double> >
  mul_a_scalar(
    af::const_ref<vec3<double> > const& a,
    double f)
  {
    af::shared<vec3<double> > result((af::reserve(a.size())));
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(a[i] * f);
    }
    return result;
  }

  af::shared<vec3<double> >
  mul_a_a_scalar(
    af::const_ref<vec3<double> > const& lhs,
    af::const_ref<double> const& rhs)
  {
    SCITBX_ASSERT(lhs.size() == rhs.size());
    af::shared<vec3<double> > result((af::reserve(lhs.size())));
    for(std::size_t i=0;i<lhs.size();i++) {
      result.push_back(lhs[i] * rhs[i]);
    }
    return result;
  }

  void
  imul_a_scalar(
    af::ref<vec3<double> > const& a,
    double f)
  {
    for(std::size_t i=0;i<a.size();i++) a[i] *= f;
  }

  af::shared<vec3<double> >
  div_a_as(
    af::ref<vec3<double> > const& lhs,
    af::ref<double> const& rhs)
  {
    SCITBX_ASSERT(lhs.size() == rhs.size());
    af::shared<vec3<double> > result((af::reserve(lhs.size())));
    for(std::size_t i=0;i<lhs.size();i++) {
      SCITBX_ASSERT(rhs[i] != 0);
      result.push_back(lhs[i] / rhs[i]);
    }
    return result;
  }

  af::shared<vec3<double> >
  mul_a_mat3(
    af::const_ref<vec3<double> > const& a,
    mat3<double> const& m)
  {
    af::shared<vec3<double> > result((af::reserve(a.size())));
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(a[i] * m);
    }
    return result;
  }

  af::shared<vec3<double> >
  rmul_a_mat3(
    af::const_ref<vec3<double> > const& a,
    mat3<double> const& m)
  {
    mat3<double> m_transposed = m.transpose();
    af::shared<vec3<double> > result((af::reserve(a.size())));
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(a[i] * m_transposed);
    }
    return result;
  }

  af::shared<double>
  dot_a_a(
    af::const_ref<vec3<double> > const& lhs,
    af::const_ref<vec3<double> > const& rhs)
  {
    SCITBX_ASSERT(lhs.size() == rhs.size());
    af::shared<double> result((af::reserve(lhs.size())));
    for(std::size_t i=0;i<lhs.size();i++) {
      result.push_back(lhs[i] * rhs[i]);
    }
    return result;
  }

  af::shared<double>
  dot_a_s(
    af::const_ref<vec3<double> > const& lhs,
    vec3<double> rhs)
  {
    af::shared<double> result((af::reserve(lhs.size())));
    for(std::size_t i=0;i<lhs.size();i++) {
      result.push_back(lhs[i] * rhs);
    }
    return result;
  }

  af::shared<double>
  dot_a(
    af::const_ref<vec3<double> > const& lhs)
  {
    af::shared<double> result((af::reserve(lhs.size())));
    for(std::size_t i=0;i<lhs.size();i++) {
      result.push_back(lhs[i] * lhs[i]);
    }
    return result;
  }

  af::shared<vec3<double> >
  cross_a_a(
    af::const_ref<vec3<double> > const& lhs,
    af::const_ref<vec3<double> > const& rhs)
  {
    SCITBX_ASSERT(lhs.size() == rhs.size());
    af::shared<vec3<double> > result((af::reserve(lhs.size())));
    for(std::size_t i=0;i<lhs.size();i++) {
      result.push_back(lhs[i].cross(rhs[i]));
    }
    return result;
  }

  af::shared<double>
  norms_(
    af::const_ref<vec3<double> > const& lhs)
  {
    af::shared<double> result((af::reserve(lhs.size())));
    for(std::size_t i=0;i<lhs.size();i++) {
      result.push_back(std::sqrt(lhs[i] * lhs[i]));
    }
    return result;
  }

  double
  sum_sq_(
    af::const_ref<vec3<double> > const& self)
  {
    double result = 0;
    for(std::size_t i=0;i<self.size();i++) {
      result += self[i] * self[i];
    }
    return result;
  }

  double
  norm_(
    af::const_ref<vec3<double> > const& self)
  {
    return std::sqrt(sum_sq_(self));
  }

  af::shared<vec3<double> >
  each_normalize(
    af::const_ref<vec3<double> > const& a,
    bool raise_if_length_zero=true)
  {
    af::shared<vec3<double> > result(a.begin(), a.end());
    vec3<double>* r = result.begin();
    std::size_t n_zero = 0;
    for(std::size_t i=0;i<a.size();i++) {
      double length = r[i].length();
      if (length == 0) n_zero++;
      else r[i] *= (1 / length);
    }
    if (n_zero != 0 && raise_if_length_zero) {
      throw std::runtime_error((boost::format(
        "flex.vec3_double.each_normalize():"
        " number of vectors with length zero: %lu of %lu")
          % n_zero % a.size()).str());
    }
    return result;
  }

  double
  max_distance(
    af::const_ref<vec3<double> > const& lhs,
    af::const_ref<vec3<double> > const& rhs)
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
    af::const_ref<vec3<double> > const& lhs,
    af::const_ref<vec3<double> > const& rhs)
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
    af::const_ref<vec3<double> > const& lhs)
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
  struct flex_default_element<vec3<double> >
  {
    static vec3<double>
    get() { return vec3<double>(0,0,0); }
  };

  void wrap_flex_vec3_double()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef flex_wrapper<vec3<double> > f_w;
    f_w::plain("vec3_double")
      .def_pickle(flex_pickle_single_buffered<vec3<double>,
        3*pickle_size_per_element<double>::value>())
      .def("__init__", make_constructor(join))
      .def("__init__", make_constructor(from_double))
      .def("part_names", part_names)
      .staticmethod("part_names")
      .def("parts", parts)
      .def("as_double", as_double)
      .def("add_selected",
        (object(*)(
          object const&,
          af::const_ref<std::size_t> const&,
          af::const_ref<vec3<double> > const&)) add_selected_unsigned_a, (
            arg("indices"), arg("values")))
      .def("min", vec3_min)
      .def("max", vec3_max)
      .def("sum", f_w::sum_a)
      .def("mean", f_w::mean_a)
      .def("mean_weighted", mean_weighted_a_a, (arg("weights")))
      .def("__add__", f_w::add_a_s)
      .def("__add__", f_w::add_a_a)
      .def("__iadd__", f_w::iadd_a_s)
      .def("__iadd__", f_w::iadd_a_a)
      .def("__sub__", f_w::sub_a_s)
      .def("__sub__", f_w::sub_a_a)
      .def("__isub__", f_w::isub_a_s)
      .def("__mul__", mul_a_scalar)
      .def("__rmul__", mul_a_scalar)
      .def("__mul__", mul_a_a_scalar)
      .def("__rmul__", mul_a_a_scalar)
      .def("__imul__", imul_a_scalar, return_self<>())
      .def("__div__", div_a_as)
      .def("__truediv__", div_a_as)
      .def("__mul__", mul_a_mat3)
      .def("__rmul__", rmul_a_mat3)
      .def("dot", dot_a_s)
      .def("dot", dot_a_a)
      .def("dot", dot_a)
      .def("cross", cross_a_a)
      .def("norms", norms_)
      .def("transpose_multiply",
        (mat3<double>(*)(
          af::const_ref<vec3<double> > const&,
          af::const_ref<vec3<double> > const&)) matrix::transpose_multiply)
      .def("sum_sq", sum_sq_)
      .def("norm", norm_)
      .def("each_normalize", each_normalize, (
        arg("raise_if_length_zero")=true))
      .def("max_distance", max_distance)
      .def("rms_difference", rms_difference)
      .def("rms_length", rms_length)
    ;
  }

}}} // namespace scitbx::af::boost_python
