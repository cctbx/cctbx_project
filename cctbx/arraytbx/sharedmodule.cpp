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
#include <cctbx/math/utils.h>

#include <cctbx/miller_bpl.h>

#include <cctbx/sgtbx/matrix.h>
#include <cctbx/xray_scatterer.h>

namespace {

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
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
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

    WRAP_TYPE("int", int);
    WRAP_TYPE("long", long);
    WRAP_TYPE("float", float);
    WRAP_TYPE("double", double);
    WRAP_TYPE("complex_double", std::complex<double>);
    WRAP_TYPE("std_string", std::string);

    WRAP_TYPE("Miller_Index", cctbx::Miller::Index);
    WRAP_TYPE("RTMx", cctbx::sgtbx::RTMx);
    WRAP_TYPE("XrayScatterer", XrayScatterer);

    WRAP_TYPE("double3", cctbx::af::double3);

    class_builder<ex_linear_regression<double> >
    py_linear_regression(this_module, "linear_regression");
    class_builder<ex_statistics<double> >
    py_statistics(this_module, "statistics");

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
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
    py_statistics.def(constructor<const cctbx::af::shared<float>&>());
#endif
    py_statistics.def(&ex_statistics<double>::min, "min");
    py_statistics.def(&ex_statistics<double>::max, "max");
    py_statistics.def(&ex_statistics<double>::mean, "mean");
    py_statistics.def(&ex_statistics<double>::mean2, "mean2");
    py_statistics.def(&ex_statistics<double>::sigma, "sigma");
  }

}

BOOST_PYTHON_MODULE_INIT(shared)
{
  boost::python::module_builder this_module("shared");
  init_module(this_module);
}
