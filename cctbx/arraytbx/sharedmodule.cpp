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

namespace cctbx { namespace af {

  boost::python::ref shared_miller_index_getstate(
    shared<miller::Index> const& a);
  void shared_miller_index_setstate(
    shared<miller::Index>& a,
    boost::python::ref state);
  boost::python::ref shared_hendrickson_lattman_double_getstate(
    shared<hendrickson_lattman<double> > const& a);
  void shared_hendrickson_lattman_double_setstate(
    shared<hendrickson_lattman<double> >& a,
    boost::python::ref state);
  boost::python::ref shared_xray_scatterer_double_wk1995_getstate(
    shared<sftbx::XrayScatterer<double, eltbx::CAASF_WK1995> > const& a);
  void shared_xray_scatterer_double_wk1995_setstate(
    shared<sftbx::XrayScatterer<double, eltbx::CAASF_WK1995> >& a,
    boost::python::ref state);

  template <>
  struct shared_pickle<miller::Index>
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(shared_miller_index_getstate, "__getstate__");
      class_bldr.def(shared_miller_index_setstate, "__setstate__");
    }
  };

  template <>
  struct shared_pickle<hendrickson_lattman<double> >
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(
        shared_hendrickson_lattman_double_getstate, "__getstate__");
      class_bldr.def(
        shared_hendrickson_lattman_double_setstate, "__setstate__");
    }
  };

  template <>
  struct shared_pickle<sftbx::XrayScatterer<double, eltbx::CAASF_WK1995> >
  {
    template <typename ClassBuilderType>
    static void def(ClassBuilderType& class_bldr)
    {
      class_bldr.def(
        shared_xray_scatterer_double_wk1995_getstate, "__getstate__");
      class_bldr.def(
        shared_xray_scatterer_double_wk1995_setstate, "__setstate__");
    }
  };

}} // namespace cctbx::af

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
  py_in_place_pow(cctbx::af::shared<double> a, double exponent)
  {
    for(std::size_t i=0;i<a.size();i++) a[i] = std::pow(a[i], exponent);
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

#define WRAP_PLAIN(python_name, element_type) \
    cctbx::af::shared_wrapper<element_type >::plain( \
      this_module, python_name)
#define WRAP_NUMERIC(python_name, element_type) \
    cctbx::af::shared_wrapper<element_type >::numeric( \
      this_module, python_name)
#define WRAP_INTEGER(python_name, element_type) \
    cctbx::af::shared_wrapper<element_type >::integer( \
      this_module, python_name)
#define WRAP_COMPLEX(python_name, element_type) \
    cctbx::af::shared_wrapper<element_type >::complex( \
      this_module, python_name)

    // bool is wrapped here to enable boolean operators for numeric types
    cctbx::af::shared_wrapper<bool>::logical(this_module, "bool");

    // double is wrapped here to enable .as_double() for the other types
    WRAP_NUMERIC("double", double);

    WRAP_INTEGER("int", int);
    WRAP_NUMERIC("float", float);
    WRAP_COMPLEX("complex_double", std::complex<double>);

#ifndef FAST_COMPILE
    WRAP_INTEGER("long", long);
    WRAP_PLAIN("std_string", std::string);

    WRAP_PLAIN("miller_Index", cctbx::miller::Index);
    WRAP_PLAIN("hendrickson_lattman", cctbx::hendrickson_lattman<double>);
    WRAP_PLAIN("RTMx", cctbx::sgtbx::RTMx);
    WRAP_PLAIN("XrayScatterer", XrayScatterer);
    WRAP_PLAIN("double3", cctbx::af::double3);

    typedef std::size_t size_t;
    WRAP_PLAIN("size_t", size_t);
    typedef cctbx::af::tiny<size_t, 2> tiny_size_t_2;
    WRAP_PLAIN("tiny_size_t_2", tiny_size_t_2);
#endif

    this_module.def(py_reinterpret_complex_as_real,
      "reinterpret_complex_as_real");
    this_module.def(py_reinterpret_real_as_complex,
      "reinterpret_real_as_complex");

    this_module.def(py_in_place_pow, "in_place_pow");
    this_module.def(py_set_if_less_than, "set_if_less_than");

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
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1300) // VC++ 7.0
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
