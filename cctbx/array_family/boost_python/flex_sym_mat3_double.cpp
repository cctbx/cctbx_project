#include <cctbx/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/boost_python/pickle_single_buffered.h>
#include <boost/python/make_constructor.hpp>

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

    flex<sym_mat3<double> >::type*
    from_double(
      af::const_ref<double> const& x)
    {
      SCITBX_ASSERT(x.size() % 6 == 0);
      std::size_t result_size = x.size() / 6;
      af::shared<sym_mat3<double> > result((af::reserve(result_size)));
      const double* d = x.begin();
      for(std::size_t i=0;i<result_size;i++) {
        result.push_back(sym_mat3<double>(d));
        d += 6;
      }
      return new flex<sym_mat3<double> >::type(result, result.size());
    }

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

  } // namespace <anonymous>

  void wrap_flex_sym_mat3_double()
  {
    flex_wrapper<sym_mat3<double> >::plain("sym_mat3_double")
      .def_pickle(flex_pickle_single_buffered<sym_mat3<double>,
        6*pickle_size_per_element<double>::value>())
      .def("__init__", boost::python::make_constructor(from_double))
      .def("as_double", as_double)
    ;
  }

}}} // namespace scitbx::af::boost_python
