#include <cctbx/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/boost_python/pickle_single_buffered.h>

namespace scitbx { namespace boost_python { namespace pickle_single_buffered {

  inline
  char* to_string(char* start, sym_mat3<double> const& value)
  {
    for(std::size_t i=0;i<6;i++) start = to_string(start, value[i]);
    return start;
  }

  template <>
  struct from_string<sym_mat3<double> >
  {
    from_string(const char* start)
    {
      end = start;
      for(std::size_t i=0;i<6;i++) {
        from_string<double> proxy(end);
        value[i] = proxy.value;
        end = proxy.end;
      }
    }

    sym_mat3<double> value;
    const char* end;
  };

}}} // namespace scitbx::boost_python::pickle_single_buffered

#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>

namespace scitbx { namespace af { namespace boost_python {

  namespace {

    flex_double
    as_double(const_ref<sym_mat3<double> > const& a)
    {
      flex_double result(a.size()*6, init_functor_null<double>());
      double* r = result.begin();
      for(std::size_t i=0;i<a.size();i++) {
        for(std::size_t j=0;j<6;j++) {
          *r++ = a[i][j];
        }
      }
      return result;
    }

    flex<sym_mat3<double> >::type
    from_double(flex<sym_mat3<double> >::type& vec, flex_double const& dbl)
    {
      CCTBX_ASSERT(vec.size() == 0);
      CCTBX_ASSERT(dbl.size() % 6 == 0);
      std::size_t vec_size = dbl.size() / 6;
      shared<sym_mat3<double> > v = vec.as_base_array();
      v.reserve(vec_size);
      const double* d = dbl.begin();
      for(std::size_t i=0;i<vec_size;i++) {
        v.push_back(sym_mat3<double>(d));
        d += 6;
      }
      vec.resize(flex_grid<>(vec_size));
      return vec;
    }

  } // namespace <anonymous>

  void wrap_flex_sym_mat3_double()
  {
    flex_wrapper<sym_mat3<double> >::plain("sym_mat3_double")
      .def_pickle(flex_pickle_single_buffered<sym_mat3<double>,
        6*pickle_size_per_element<double>::value>())
      .def("as_double", as_double)
      .def("from_double", from_double)
    ;
  }

}}} // namespace scitbx::af::boost_python
