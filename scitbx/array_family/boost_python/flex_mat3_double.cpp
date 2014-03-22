#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/serialization/single_buffered.h>
#include <scitbx/matrix/transpose_multiply.h>
#include <scitbx/math/utils.h>
#include <scitbx/mat3.h>
#include <boost/python/make_constructor.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_arg.hpp>
#include <boost/format.hpp>
#include <scitbx/array_family/boost_python/flex_helpers.h>

namespace scitbx { namespace serialization { namespace single_buffered {

  inline
  char* to_string(char* start, scitbx::mat3<double> const& value)
  {
    return
      to_string(to_string(to_string(
        to_string(to_string(to_string(
          to_string(to_string(to_string(
            start,
              value[0]), value[1]), value[2]),
              value[3]), value[4]), value[5]),
              value[6]), value[7]), value[8]);
  }

  template <>
  struct from_string<scitbx::mat3<double> >
  {
    from_string(const char* start)
    {
      end = start;
      for(std::size_t i=0;i<9;i++) {
        from_string<double> proxy(end);
        value[i] = proxy.value;
        end = proxy.end;
      }
    }

    scitbx::mat3<double> value;
    const char* end;
  };

}}} // namespace scitbx::serialization::single_buffered

#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>

namespace scitbx { namespace af { namespace boost_python {

  template <>
  struct flex_default_element<scitbx::mat3<double> >
  {
    static scitbx::mat3<double>
    get() { return scitbx::mat3<double>(0,0,0,0,0,0,0,0,0); }
  };

}}}

namespace scitbx { namespace af {
namespace {

  scitbx::af::flex<scitbx::mat3<double> >::type*
  from_double(
    scitbx::af::const_ref<double> const& x)
  {
    SCITBX_ASSERT(x.size() % 9 == 0);
    std::size_t result_size = x.size() / 9;
    scitbx::af::shared<scitbx::mat3<double> > result(result_size,
      scitbx::af::init_functor_null< scitbx::mat3<double> >());
    const double* d = x.begin();
    for(std::size_t i=0;i<result_size;i++) {
      for (std::size_t j=0;j<9;++j) {
        result[i][j] = *d++;
      }
    }
    return new scitbx::af::flex<scitbx::mat3<double> >::type(result, result.size());
  }

  scitbx::af::flex_double
  as_double(scitbx::af::flex<scitbx::mat3<double> >::type const& a)
  {
    SCITBX_ASSERT(a.accessor().is_trivial_1d());
    scitbx::af::flex_double result(a.size()*9, scitbx::af::init_functor_null<double>());
    double* r = result.begin();
    scitbx::af::const_ref<scitbx::mat3<double> > a_ref = a.const_ref().as_1d();
    for(std::size_t i=0;i<a_ref.size();i++) {
      for(std::size_t j=0;j<9;j++) {
        *r++ = a_ref[i][j];
      }
    }
    return result;
  }

  scitbx::af::shared<mat3<double> >
  mul_a_scalar(
    scitbx::af::const_ref<mat3<double> > const& a,
    double f)
  {
    scitbx::af::shared<mat3<double> > result((af::reserve(a.size())));
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(a[i] * f);
    }
    return result;
  }

  scitbx::af::shared<mat3<double> >
  mul_a_a_scalar(
    scitbx::af::const_ref<mat3<double> > const& lhs,
    scitbx::af::const_ref<double> const& rhs)
  {
    SCITBX_ASSERT(lhs.size() == rhs.size());
    scitbx::af::shared<mat3<double> > result((af::reserve(lhs.size())));
    for(std::size_t i=0;i<lhs.size();i++) {
      result.push_back(lhs[i] * rhs[i]);
    }
    return result;
  }

  scitbx::af::shared<vec3<double> >
  mul_a_vec3(
    scitbx::af::const_ref<mat3<double> > const& a,
    vec3<double> const& m)
  {
    af::shared<vec3<double> > result((af::reserve(a.size())));
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(a[i] * m);
    }
    return result;
  }

  scitbx::af::shared<vec3<double> >
  mul_a_a_vec3(
    scitbx::af::const_ref<mat3<double> > const& lhs,
    scitbx::af::const_ref<vec3<double> > const& rhs)
  {
    SCITBX_ASSERT(lhs.size() == rhs.size());
    af::shared<vec3<double> > result((af::reserve(lhs.size())));
    for(std::size_t i=0;i<lhs.size();i++) {
      result.push_back(lhs[i] * rhs[i]);
    }
    return result;
  }

  scitbx::af::shared<mat3<double> >
  mul_a_mat3(
    scitbx::af::const_ref<mat3<double> > const& a,
    mat3<double> const& m)
  {
    af::shared<mat3<double> > result((af::reserve(a.size())));
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(a[i] * m);
    }
    return result;
  }

  scitbx::af::shared<mat3<double> >
  mul_a_a_mat3(
    scitbx::af::const_ref<mat3<double> > const& lhs,
    scitbx::af::const_ref<mat3<double> > const& rhs)
  {
    SCITBX_ASSERT(lhs.size() == rhs.size());
    af::shared<mat3<double> > result((af::reserve(lhs.size())));
    for(std::size_t i=0;i<lhs.size();i++) {
      result.push_back(lhs[i] * rhs[i]);
    }
    return result;
  }

} // namespace <anonymous>

namespace boost_python {

  void wrap_flex_mat3_double()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef scitbx::af::boost_python::flex_wrapper<scitbx::mat3<double> > f_w;
    f_w::plain("mat3_double")
      .def_pickle(scitbx::af::boost_python::flex_pickle_single_buffered<
        scitbx::mat3<double>,
        9*scitbx::af::boost_python::pickle_size_per_element<double>::value>())
      .def("__init__", make_constructor(from_double))
      .def("__mul__", mul_a_scalar)
      .def("__mul__", mul_a_a_scalar)
      .def("__mul__", mul_a_vec3)
      .def("__mul__", mul_a_a_vec3)
      .def("__mul__", mul_a_mat3)
      .def("__mul__", mul_a_a_mat3)
      .def("as_double", as_double);
    ;
  }

}}} // namespace scitbx::af::boost_python
