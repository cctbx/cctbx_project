// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Mar 2002: modified copy of shared_storagemodule.cpp (rwgk)
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/array_family/tiny_bpl.h>
#include <cctbx/array_family/shared_bpl.h>
#include <cctbx/math/linear_regression.h>
#include <cctbx/math/array_utils.h>

#include <cctbx/miller_bpl.h>
#include <cctbx/hendrickson_lattman_bpl.h>

#include <cctbx/sgtbx/matrix.h>
#include <cctbx/sftbx/xray_scatterer.h>

namespace {

  cctbx::af::shared<double>
  py_reinterpret_complex_as_real(
    cctbx::af::shared<std::complex<double> > a)
  {
    return cctbx::af::shared<double>(a.handle());
  }

  cctbx::af::shared<std::complex<double> >
  py_reinterpret_real_as_complex(
    cctbx::af::shared<double> a)
  {
    return cctbx::af::shared<std::complex<double> >(a.handle());
  }

  void
  py_square(cctbx::af::shared<double> a)
  {
    for(std::size_t i=0;i<a.size();i++) a[i] *= a[i];
  }

  void
  py_set_if_less_than(
    cctbx::af::shared<double> a,
    double threshold_value,
    double imposed_value)
  {
    for(std::size_t i=0;i<a.size();i++) {
      if (a[i] < threshold_value) a[i] = imposed_value;
    }
  }

  // to preserve VC6 compatibility we are not using shared_algebra.h
  cctbx::af::shared<double>
  py_abs_complex(cctbx::af::shared<std::complex<double> > a)
  {
    cctbx::af::shared<double> result;
    result.reserve(a.size());
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(std::abs(a[i]));
    }
    return result;
  }
  cctbx::af::shared<double>
  py_arg_complex(
    cctbx::af::shared<std::complex<double> > a,
    bool deg = false)
  {
    cctbx::af::shared<double> result;
    result.reserve(a.size());
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(std::arg(a[i]));
      if (deg) result[i] /= cctbx::constants::pi_180;
    }
    return result;
  }
  cctbx::af::shared<double>
  py_norm_complex(cctbx::af::shared<std::complex<double> > a)
  {
    cctbx::af::shared<double> result;
    result.reserve(a.size());
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(std::norm(a[i]));
    }
    return result;
  }
  cctbx::af::shared<std::complex<double> >
  py_polar_complex(
    cctbx::af::shared<double> rho,
    cctbx::af::shared<double> theta,
    bool deg = false)
  {
    cctbx_assert(rho.size() == theta.size());
    cctbx::af::shared<std::complex<double> > result;
    result.reserve(rho.size());
    if (deg) {
      for(std::size_t i=0;i<rho.size();i++) {
        result.push_back(std::polar(rho[i], theta[i]));
      }
    }
    else {
      for(std::size_t i=0;i<rho.size();i++) {
        result.push_back(
          std::polar(rho[i], theta[i] * cctbx::constants::pi_180));
      }
    }
    return result;
  }

  template <typename FloatType>
  struct ex_linear_regression : cctbx::math::linear_regression<FloatType>
  {
    typedef cctbx::math::linear_regression<FloatType> base_type;

    ex_linear_regression() {}
    ex_linear_regression(const cctbx::af::shared<FloatType>& x,
                         const cctbx::af::shared<FloatType>& y)
      : base_type(x.const_ref(), y.const_ref())
    {}
    ex_linear_regression(const cctbx::af::shared<FloatType>& x,
                         const cctbx::af::shared<FloatType>& y,
                         const FloatType& epsilon)
      : base_type(x.const_ref(), y.const_ref(), epsilon)
    {}
    bool is_well_defined() const { return base_type::is_well_defined(); }
    const FloatType& b() const { return base_type::b(); }
    const FloatType& m() const { return base_type::m(); }
    const FloatType& cc() const { return base_type::cc(); }
  };

  template <typename FloatType>
  struct ex_statistics : cctbx::math::array_statistics<FloatType>
  {
    typedef cctbx::math::array_statistics<FloatType> base_type;

    ex_statistics() {}
    ex_statistics(const cctbx::af::shared<FloatType>& x)
      : base_type(x.const_ref())
    {}
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1300) // VC++ 7.0
    template <typename OtherFloatType>
    ex_statistics(const cctbx::af::shared<OtherFloatType>& x)
      : base_type(cctbx::af::shared<FloatType>(x.begin(), x.end()).const_ref())
    {}
#endif
    const FloatType& min() const { return base_type::min(); }
    const FloatType& max() const { return base_type::max(); }
    const FloatType& mean() const { return base_type::mean(); }
    const FloatType& mean2() const { return base_type::mean2(); }
    const FloatType& sigma() const { return base_type::sigma(); }
  };

# include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    python::import_converters<cctbx::sgtbx::RTMx>
    py_RTMx("cctbx_boost.sgtbx", "RTMx");

    typedef cctbx::sftbx::XrayScatterer<
      double, cctbx::eltbx::CAASF_WK1995> XrayScatterer;
    python::import_converters<XrayScatterer>
    py_XrayScatterer("cctbx_boost.sftbx", "XrayScatterer");

#define WRAP_TYPE(python_name, element_type) \
    cctbx::af::wrap_shared<element_type >::run(this_module, python_name)

    WRAP_TYPE("bool", bool);
    WRAP_TYPE("int", int);
    WRAP_TYPE("long", long);
    WRAP_TYPE("float", float);
    WRAP_TYPE("double", double);
    WRAP_TYPE("complex_double", std::complex<double>);
    WRAP_TYPE("std_string", std::string);

    WRAP_TYPE("Miller_Index", cctbx::Miller::Index);
    WRAP_TYPE("hendrickson_lattman", cctbx::hendrickson_lattman<double>);
    WRAP_TYPE("RTMx", cctbx::sgtbx::RTMx);
    WRAP_TYPE("XrayScatterer", XrayScatterer);

    WRAP_TYPE("double3", cctbx::af::double3);

    typedef std::size_t size_t;
    WRAP_TYPE("size_t", size_t);
    typedef cctbx::af::tiny<size_t, 2> tiny_size_t_2;
    WRAP_TYPE("tiny_size_t_2", tiny_size_t_2);

    this_module.def(py_reinterpret_complex_as_real,
      "reinterpret_complex_as_real");
    this_module.def(py_reinterpret_real_as_complex,
      "reinterpret_real_as_complex");

    this_module.def(py_square, "square");
    this_module.def(py_set_if_less_than, "set_if_less_than");
    this_module.def(py_abs_complex, "abs");
    this_module.def(py_arg_complex, "arg");
    this_module.def(
      (cctbx::af::shared<double>
       (*)(cctbx::af::shared<std::complex<double> >))
      py_arg_complex, "arg");
    this_module.def(py_norm_complex, "norm");
    this_module.def(py_polar_complex, "polar");
    this_module.def(
      (cctbx::af::shared<std::complex<double> >
       (*)(cctbx::af::shared<double>, cctbx::af::shared<double>))
      py_polar_complex, "polar");

    class_builder<ex_linear_regression<double> >
    py_linear_regression(this_module, "linear_regression");
    class_builder<ex_statistics<double> >
    py_statistics(this_module, "statistics");

    class_builder<cctbx::Miller::join_sets>
    py_join_sets(this_module, "join_sets");

    py_linear_regression.def(constructor<>());
    py_linear_regression.def(constructor<
      const cctbx::af::shared<double>&,
      const cctbx::af::shared<double>&>());
    py_linear_regression.def(constructor<
      const cctbx::af::shared<double>&,
      const cctbx::af::shared<double>&,
      const double&>());
    py_linear_regression.def(&ex_linear_regression<double>::is_well_defined,
      "is_well_defined");
    py_linear_regression.def(&ex_linear_regression<double>::b, "b");
    py_linear_regression.def(&ex_linear_regression<double>::m, "m");
    py_linear_regression.def(&ex_linear_regression<double>::cc, "cc");

    py_statistics.def(constructor<>());
    py_statistics.def(constructor<const cctbx::af::shared<double>&>());
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1300) // VC++ 7.0
    py_statistics.def(constructor<const cctbx::af::shared<float>&>());
#endif
    py_statistics.def(&ex_statistics<double>::min, "min");
    py_statistics.def(&ex_statistics<double>::max, "max");
    py_statistics.def(&ex_statistics<double>::mean, "mean");
    py_statistics.def(&ex_statistics<double>::mean2, "mean2");
    py_statistics.def(&ex_statistics<double>::sigma, "sigma");

    py_join_sets.def(constructor<>());
    py_join_sets.def(constructor<
      cctbx::af::shared<cctbx::Miller::Index>,
      cctbx::af::shared<cctbx::Miller::Index> >());
    py_join_sets.def(&cctbx::Miller::join_sets::pairs, "pairs");
    py_join_sets.def(&cctbx::Miller::join_sets::singles, "singles");
    py_join_sets.def(&cctbx::Miller::join_sets::have_singles, "have_singles");
  }

}

BOOST_PYTHON_MODULE_INIT(shared)
{
  boost::python::module_builder this_module("shared");
  init_module(this_module);
}
