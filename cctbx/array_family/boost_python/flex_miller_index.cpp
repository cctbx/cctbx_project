#include <cctbx/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/boost_python/pickle_single_buffered.h>

namespace scitbx { namespace boost_python { namespace pickle_single_buffered {

  inline
  char* to_string(char* start, cctbx::miller::index<> const& value)
  {
    return
      to_string(to_string(to_string(start, value[0]), value[1]), value[2]);
  }

  template <>
  struct from_string<cctbx::miller::index<> >
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

    cctbx::miller::index<> value;
    const char* end;
  };

}}} // namespace scitbx::boost_python::pickle_single_buffered

#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  af::shared<vec3<double> >
  as_vec3_double(af::const_ref<cctbx::miller::index<> > const& a)
  {
    af::shared<vec3<double> > result((af::reserve(a.size())));
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(vec3<double>(a[i]));
    }
    return result;
  }

} // namespace <anonymous>

  void wrap_flex_miller_index(boost::python::object const& flex_root_scope)
  {
    using namespace cctbx;

    flex_wrapper<miller::index<> >::ordered("miller_index", flex_root_scope)
      .def("__neg__", flex_wrapper<miller::index<> >::neg_a)
      .def_pickle(flex_pickle_single_buffered<miller::index<>,
        3*pickle_size_per_element<miller::index<>::value_type>::value>())
      .def("as_vec3_double", as_vec3_double)
    ;
  }

}}} // namespace scitbx::af::boost_python
