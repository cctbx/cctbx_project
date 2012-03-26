#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/serialization/single_buffered.h>
#include <boost/python/make_constructor.hpp>
#include <tbxx/error_utils.hpp>

namespace scitbx { namespace serialization { namespace single_buffered {

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

}}} // namespace scitbx::serialization::single_buffered

#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  flex<cctbx::miller::index<> >::type*
  join(
    af::const_ref<int> const& h,
    af::const_ref<int> const& k,
    af::const_ref<int> const& l)
  {
    SCITBX_ASSERT(k.size() == h.size());
    SCITBX_ASSERT(l.size() == h.size());
    af::shared<cctbx::miller::index<> > result((af::reserve(h.size())));
    for(std::size_t i=0;i<h.size();i++) {
      result.push_back(cctbx::miller::index<>(h[i],k[i],l[i]));
    }
    return new flex<cctbx::miller::index<> >::type(result, result.size());
  }

  af::shared<vec3<double> >
  as_vec3_double(af::const_ref<cctbx::miller::index<> > const& a)
  {
    af::shared<vec3<double> > result((af::reserve(a.size())));
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(vec3<double>(a[i]));
    }
    return result;
  }

  double
  fourier_transform_real_part_at_x(
    af::const_ref<cctbx::miller::index<> > const& miller_indices,
    af::const_ref<std::complex<double> > const& fourier_coeffs,
    af::tiny<double, 3> const& x)
  {
    TBXX_ASSERT(fourier_coeffs.size() == miller_indices.size());
    double result = 0;
    double mtp = -constants::two_pi;
    for(std::size_t i=0;i<miller_indices.size();i++) {
      cctbx::miller::index<> const& h = miller_indices[i];
      double mtphx = mtp * (h[0]*x[0] + h[1]*x[1] + h[2]*x[2]);
      double a = std::cos(mtphx);
      double b = std::sin(mtphx);
      std::complex<double> const& fc = fourier_coeffs[i];
      result += fc.real() * a - fc.imag() * b;
    }
    return result;
  }

} // namespace <anonymous>

  void wrap_flex_miller_index(boost::python::object const& flex_root_scope)
  {
    namespace bp = boost::python;
    namespace miller = cctbx::miller;

    flex_wrapper<miller::index<> >::ordered("miller_index", flex_root_scope)
      .def("__init__", bp::make_constructor(join))
      .def("__neg__", flex_wrapper<miller::index<> >::neg_a)
      .def_pickle(flex_pickle_single_buffered<miller::index<>,
        3*pickle_size_per_element<miller::index<>::value_type>::value>())
      .def("as_vec3_double", as_vec3_double)
      .def("fourier_transform_real_part_at_x",
        fourier_transform_real_part_at_x, (
          bp::arg("fourier_coeffs"), bp::arg("x")))
    ;
  }

}}} // namespace scitbx::af::boost_python
