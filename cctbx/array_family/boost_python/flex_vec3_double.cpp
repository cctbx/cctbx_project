/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (R.W. Grosse-Kunstleve)
 */

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/boost_python/pickle_single_buffered.h>
#include <scitbx/vec3.h>
#include <cctbx/error.h>

namespace scitbx { namespace boost_python { namespace pickle_single_buffered {

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

}}} // namespace scitbx::boost_python::pickle_single_buffered

#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>

namespace scitbx { namespace af { namespace boost_python {

  namespace {

    flex_double
    as_double(flex<vec3<double> >::type const& a)
    {
      CCTBX_ASSERT(a.accessor().is_trivial_1d());
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

    flex<vec3<double> >::type
    from_double(flex<vec3<double> >::type& vec, flex_double const& dbl)
    {
      CCTBX_ASSERT(vec.size() == 0);
      CCTBX_ASSERT(dbl.size() % 3 == 0);
      std::size_t vec_size = dbl.size() / 3;
      shared<vec3<double> > v = vec.as_base_array();
      v.reserve(vec_size);
      const double* d = dbl.begin();
      for(std::size_t i=0;i<vec_size;i++) {
        v.push_back(vec3<double>(d));
        d += 3;
      }
      vec.resize(flex_grid<>(vec_size));
      return vec;
    }

  } // namespace <anonymous>

  void wrap_flex_vec3_double()
  {
    flex_wrapper<vec3<double> >::plain("vec3_double")
      .def_pickle(flex_pickle_single_buffered<vec3<double>,
        3*pickle_size_per_element<double>::value>())
      .def("as_double", as_double)
      .def("from_double", from_double)
    ;
  }

}}} // namespace scitbx::af::boost_python
