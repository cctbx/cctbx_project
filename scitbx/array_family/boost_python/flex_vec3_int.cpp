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
  char* to_string(char* start, vec3<int> const& value)
  {
    return
      to_string(to_string(to_string(start, value[0]), value[1]), value[2]);
  }

  template <>
  struct from_string<vec3<int> >
  {
    from_string(const char* start)
    {
      end = start;
      for(std::size_t i=0;i<3;i++) {
        from_string<int> proxy(end);
        value[i] = proxy.value;
        end = proxy.end;
      }
    }

    vec3<int> value;
    const char* end;
  };

}}} // namespace scitbx::serialization::single_buffered

#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>

namespace scitbx { namespace af {
namespace {

  flex<vec3<int> >::type*
  join(
    af::const_ref<int> const& x,
    af::const_ref<int> const& y,
    af::const_ref<int> const& z)
  {
    SCITBX_ASSERT(y.size() == x.size());
    SCITBX_ASSERT(z.size() == x.size());
    af::shared<vec3<int> > result((af::reserve(x.size())));
    for(std::size_t i=0;i<x.size();i++) {
      result.push_back(vec3<int>(x[i],y[i],z[i]));
    }
    return new flex<vec3<int> >::type(result, result.size());
  }

  flex<vec3<int> >::type*
  from_int(
    af::const_ref<int> const& x)
  {
    SCITBX_ASSERT(x.size() % 3 == 0);
    std::size_t result_size = x.size() / 3;
    af::shared<vec3<int> > result((af::reserve(result_size)));
    const int* d = x.begin();
    for(std::size_t i=0;i<result_size;i++) {
      result.push_back(vec3<int>(d));
      d += 3;
    }
    return new flex<vec3<int> >::type(result, result.size());
  }

  flex_int
  as_int(flex<vec3<int> >::type const& a)
  {
    SCITBX_ASSERT(a.accessor().is_trivial_1d());
    flex_int result(a.size()*3, init_functor_null<int>());
    int* r = result.begin();
    const_ref<vec3<int> > a_ref = a.const_ref().as_1d();
    for(std::size_t i=0;i<a_ref.size();i++) {
      for(std::size_t j=0;j<3;j++) {
        *r++ = a_ref[i][j];
      }
    }
    return result;
  }

  af::shared<vec3<double> >
  as_vec3_double(flex<vec3<int> >::type const& a)
  {
    af::shared<vec3<double> > result((af::reserve(a.size())));
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(vec3<double>(a[i]));
    }
    return result;
  }

} // namespace <anonymous>

namespace boost_python {

  template <>
  struct flex_default_element<vec3<int> >
  {
    static vec3<int>
    get() { return vec3<int>(0,0,0); }
  };

  void wrap_flex_vec3_int()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef flex_wrapper<vec3<int> > f_w;
    f_w::plain("vec3_int")
      .def_pickle(flex_pickle_single_buffered<vec3<int>,
        3*pickle_size_per_element<int>::value>())
      .def("__init__", make_constructor(join))
      .def("__init__", make_constructor(from_int))
      .def("as_int", as_int)
      .def("as_vec3_double", as_vec3_double)
      .def("add_selected",
        (object(*)(
          object const&,
          af::const_ref<std::size_t> const&,
          af::const_ref<vec3<int> > const&)) add_selected_unsigned_a, (
            arg("indices"), arg("values")))
    ;
  }

}}} // namespace scitbx::af::boost_python
