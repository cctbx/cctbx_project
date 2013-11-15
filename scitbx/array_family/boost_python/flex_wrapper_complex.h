#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_WRAPPER_COMPLEX_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_WRAPPER_COMPLEX_H

#include <scitbx/array_family/versa_algebra.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/constants.h>
#include <scitbx/array_family/boost_python/utils.h>
#include <boost/python/def.hpp>
#include <boost/python/scope.hpp>

namespace scitbx { namespace af { namespace boost_python {

  template <typename FloatType>
  struct flex_wrapper_complex_functions
  {
    static versa<FloatType, flex_grid<> >
    real_complex(versa<std::complex<FloatType>, flex_grid<> > const& a)
    {
      return real(a);
    }

    static versa<FloatType, flex_grid<> >
    imag_complex(versa<std::complex<FloatType>, flex_grid<> > const& a)
    {
      return imag(a);
    }

    static versa<std::complex<FloatType>, flex_grid<> >
    conj_complex(versa<std::complex<FloatType>, flex_grid<> > const& a)
    {
      return conj(a);
    }

    static versa<FloatType, flex_grid<> >
    abs_complex(versa<std::complex<FloatType>, flex_grid<> > const& a)
    {
      shared_plain<FloatType> result(a.size(), init_functor_null<FloatType>());
      for(std::size_t i=0;i<a.size();i++) result[i] = std::abs(a[i]);
      return versa<FloatType, flex_grid<> >(result, a.accessor());
    }

    static versa<FloatType, flex_grid<> >
    arg_complex_2(versa<std::complex<FloatType>, flex_grid<> > const& a,
                  bool deg)
    {
      shared_plain<FloatType> result(a.size(), init_functor_null<FloatType>());
      for(std::size_t i=0;i<a.size();i++) {
        result[i] = std::arg(a[i]);
        if (deg) result[i] /= constants::pi_180;
      }
      return versa<FloatType, flex_grid<> >(result, a.accessor());
    }

    static versa<FloatType, flex_grid<> >
    arg_complex_1(versa<std::complex<FloatType>, flex_grid<> > const& a)
    {
      return arg_complex_2(a, false);
    }

    static versa<FloatType, flex_grid<> >
    norm_complex(versa<std::complex<FloatType>, flex_grid<> > const& a)
    {
      return norm(a);
    }

    static versa<std::complex<FloatType>, flex_grid<> >
    polar_complex_r_r_3(
      versa<FloatType, flex_grid<> > const& rho,
      versa<FloatType, flex_grid<> > const& theta,
      bool deg)
    {
      if (rho.accessor() != theta.accessor()) {
        raise_incompatible_arrays();
      }
      shared_plain<std::complex<FloatType> > result(
        rho.size(), init_functor_null<std::complex<FloatType> >());
      if (deg) {
        for(std::size_t i=0;i<rho.size();i++) {
          SCITBX_ASSERT(rho[i] >= 0)(rho[i]);
          result[i] = std::polar(rho[i], theta[i] * constants::pi_180);
        }
      }
      else {
        for(std::size_t i=0;i<rho.size();i++) {
          SCITBX_ASSERT(rho[i] >= 0)(rho[i]);
          result[i] = std::polar(rho[i], theta[i]);
        }
      }
      return versa<std::complex<FloatType>, flex_grid<> >(
        result, rho.accessor());
    }

    static versa<std::complex<FloatType>, flex_grid<> >
    polar_complex_r_r_2(
      versa<FloatType, flex_grid<> > const& rho,
      versa<FloatType, flex_grid<> > const& theta)
    {
      return polar_complex_r_r_3(rho, theta, false);
    }

    static versa<std::complex<FloatType>, flex_grid<> >
    polar_complex_c_r_3(
      versa<std::complex<FloatType>, flex_grid<> > const& rho,
      versa<FloatType, flex_grid<> > const& theta,
      bool deg)
    {
      if (rho.accessor() != theta.accessor()) {
        raise_incompatible_arrays();
      }
      shared_plain<std::complex<FloatType> > result(
        rho.size(), init_functor_null<std::complex<FloatType> >());
      if (deg) {
        for(std::size_t i=0;i<rho.size();i++) {
          result[i] = std::polar(std::abs(rho[i]), theta[i]*constants::pi_180);
        }
      }
      else {
        for(std::size_t i=0;i<rho.size();i++) {
          result[i] = std::polar(std::abs(rho[i]), theta[i]);
        }
      }
      return versa<std::complex<FloatType>, flex_grid<> >(
        result, rho.accessor());
    }

    static versa<std::complex<FloatType>, flex_grid<> >
    polar_complex_c_r_2(
      versa<std::complex<FloatType>, flex_grid<> > const& rho,
      versa<FloatType, flex_grid<> > const& theta)
    {
      return polar_complex_c_r_3(rho, theta, false);
    }

    static versa<std::complex<FloatType>, flex_grid<> >
    polar_complex_r_c(
      versa<FloatType, flex_grid<> > const& rho,
      versa<std::complex<FloatType>, flex_grid<> > const& theta)
    {
      if (rho.accessor() != theta.accessor()) {
        raise_incompatible_arrays();
      }
      shared_plain<std::complex<FloatType> > result(
      rho.size(), init_functor_null<std::complex<FloatType> >());
      for(std::size_t i=0;i<rho.size();i++) {
        SCITBX_ASSERT(rho[i] >= 0)(rho[i]);
        result[i] = std::polar(rho[i], std::arg(theta[i]));
      }
      return versa<std::complex<FloatType>, flex_grid<> >(
        result, rho.accessor());
    }

    static versa<std::complex<FloatType>, flex_grid<> >
    polar_complex_c_c(
      versa<std::complex<FloatType>, flex_grid<> > const& rho,
      versa<std::complex<FloatType>, flex_grid<> > const& theta)
    {
      if (rho.accessor() != theta.accessor()) {
        raise_incompatible_arrays();
      }
      shared_plain<std::complex<FloatType> > result(
      rho.size(), init_functor_null<std::complex<FloatType> >());
      for(std::size_t i=0;i<rho.size();i++) {
        result[i] = std::polar(std::abs(rho[i]), std::arg(theta[i]));
      }
      return versa<std::complex<FloatType>, flex_grid<> >(
        result, rho.accessor());
    }

    static versa<std::complex<FloatType>, flex_grid<> >
    polar_complex_rs_r(
      FloatType const& rho,
      versa<FloatType, flex_grid<> > const& theta)
    {
      shared_plain<std::complex<FloatType> > result(
        theta.size(), init_functor_null<std::complex<FloatType> >());
      for(std::size_t i=0;i<theta.size();i++) {
        SCITBX_ASSERT(rho >= 0)(rho);
        result[i] = std::polar(rho, theta[i]);
      }
      return versa<std::complex<FloatType>, flex_grid<> >(
        result, theta.accessor());
    }

    static versa<std::complex<FloatType>, flex_grid<> >
    polar_complex_r_rs(
      versa<FloatType, flex_grid<> > const& rho,
      FloatType const& theta)
    {
      shared_plain<std::complex<FloatType> > result(
        rho.size(), init_functor_null<std::complex<FloatType> >());
      for(std::size_t i=0;i<rho.size();i++) {
        SCITBX_ASSERT(rho[i] >= 0)(rho[i]);
        result[i] = std::polar(rho[i], theta);
      }
      return versa<std::complex<FloatType>, flex_grid<> >(
        result, rho.accessor());
    }

    static void
    wrap(boost::python::object const& flex_root_scope)
    {
      boost::python::scope local_scope(flex_root_scope);
      boost::python::def("real", real_complex);
      boost::python::def("imag", imag_complex);
      boost::python::def("conj", conj_complex);
      boost::python::def("abs", abs_complex);
      boost::python::def("arg", arg_complex_2);
      boost::python::def("arg", arg_complex_1);
      boost::python::def("norm", norm_complex);
      boost::python::def("polar", polar_complex_r_r_3);
      boost::python::def("polar", polar_complex_r_r_2);
      boost::python::def("polar", polar_complex_c_r_3);
      boost::python::def("polar", polar_complex_c_r_2);
      boost::python::def("polar", polar_complex_r_c);
      boost::python::def("polar", polar_complex_c_c);
      boost::python::def("polar", polar_complex_rs_r);
      boost::python::def("polar", polar_complex_r_rs);
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_WRAPPER_COMPLEX_H
