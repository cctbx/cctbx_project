// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created, fragments from sharedmodule.cpp (R.W. Grosse-Kunstleve)
 */

#include <iostream>

#include <boost/python/cross_module.hpp>
#include <cctbx/array_family/tiny_bpl.h>
#include <cctbx/array_family/tiny_types.h>
#include <cctbx/math/linear_regression.h>
#include <cctbx/array_family/loops.h>
#include <cctbx/math/array_utils.h>
#include <cctbx/maps/accessors.h>

#include <cctbx/array_family/flex_bpl.h>

namespace cctbx { namespace af { namespace bpl { namespace { namespace {

  void import_flex()
  {
    CCTBX_ARRAY_FAMILY_FLEX_IMPORT(float, "float")
    CCTBX_ARRAY_FAMILY_FLEX_IMPORT(double, "double")
  }

}}}}} // namespace cctbx::af::bpl<anonymous><anonymous>

CCTBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(float)
CCTBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(double)

namespace {

  typedef cctbx::af::versa<float, cctbx::af::flex_grid<> > flex_float;
  typedef cctbx::af::versa<double, cctbx::af::flex_grid<> > flex_double;

  void
  py_in_place_pow(flex_double a, double exponent)
  {
    for(std::size_t i=0;i<a.size();i++) a[i] = std::pow(a[i], exponent);
  }

  void
  py_set_if_less_than(
    flex_double a,
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
    ex_linear_regression(cctbx::af::shared<FloatType> x,
                         cctbx::af::shared<FloatType> y)
      : base_type(x.const_ref(), y.const_ref())
    {}
    ex_linear_regression(cctbx::af::shared<FloatType> x,
                         cctbx::af::shared<FloatType> y,
                         FloatType const& epsilon)
      : base_type(x.const_ref(), y.const_ref(), epsilon)
    {}
    bool is_well_defined() const { return base_type::is_well_defined(); }
    FloatType const& b() const { return base_type::b(); }
    FloatType const& m() const { return base_type::m(); }
    FloatType const& cc() const { return base_type::cc(); }
  };

  template <typename FloatType>
  struct ex_statistics : cctbx::math::array_statistics<FloatType>
  {
    typedef cctbx::math::array_statistics<FloatType> base_type;

    ex_statistics() {}
    ex_statistics(
      cctbx::af::versa<FloatType, cctbx::af::flex_grid<> > x)
      : base_type(x.as_base_array().const_ref())
    {
      cctbx_assert(!x.accessor().is_padded());
    }
    ex_statistics(
      cctbx::af::versa<float, cctbx::af::flex_grid<> > x)
      : base_type(cctbx::af::shared_plain<FloatType>(
          cctbx::af::adapt(x)).const_ref())
    {
      cctbx_assert(!x.accessor().is_padded());
    }
    FloatType const& min() const { return base_type::min(); }
    FloatType const& max() const { return base_type::max(); }
    FloatType const& mean() const { return base_type::mean(); }
    FloatType const& mean2() const { return base_type::mean2(); }
    FloatType const& sigma() const { return base_type::sigma(); }
  };

  template <typename FloatType>
  struct map_utils
  {
    static
    void
    inplace_unpad(
      cctbx::af::versa<FloatType, cctbx::af::flex_grid<> >& map)
    {
      using namespace cctbx;
      cctbx_assert(map.accessor().nd() == 3);
      cctbx_assert(map.accessor().is_0_based());
      if (map.accessor().layout().size() == 0) return;
      af::long3 m_real(af::adapt(map.accessor().grid()));
      af::long3 n_real(af::adapt(map.accessor().layout()));
      cctbx_assert(n_real[0] == m_real[0]);
      cctbx_assert(n_real[1] == m_real[1]);
      cctbx_assert(n_real[2] <= m_real[2]);
      if (n_real[2] == m_real[2]) return;
      af::ref<FloatType, af::grid<3> > map_n(map.begin(), n_real);
      af::ref<FloatType, af::grid<3> > map_m(map.begin(), m_real);
      af::nested_loop<af::long3> loop(n_real);
      for(af::long3 const& point = loop(); !loop.over(); loop.incr()) {
        map_n(point) = map_m(point);
      }
      map.resize(af::flex_grid<>(af::adapt(n_real)));
    }
  };

  // PyMol support, based on code by N.K. Sauter
  template <typename InpFloatType,
            typename OutFloatType>
  struct as_CObjectZYX
  {
    static std::size_t out_size(
      cctbx::af::long3 const& first, cctbx::af::long3 const& last)
    {
      std::size_t result = 1;
      for(std::size_t i=0;i<3;i++) {
        cctbx_assert(last[i] >= first[i]);
        result *= (last[i] - first[i] + 1);
      }
      return result;
    }

    static
    boost::python::ref
    convert(
      cctbx::af::versa<InpFloatType, cctbx::af::flex_grid<> > a,
      cctbx::af::long3 const& first,
      cctbx::af::long3 const& last,
      bool apply_sigma_scaling)
    {
      using namespace cctbx;
      cctbx_assert(a.accessor().nd() == 3);
      cctbx_assert(a.accessor().is_0_based());
      cctbx_assert(!a.accessor().is_padded());
      math::array_statistics<InpFloatType> map_statistics;
      InpFloatType mean = 0;
      InpFloatType sigma = 0;
      if (apply_sigma_scaling) {
        math::array_statistics<InpFloatType>
        map_statistics(a.as_base_array().const_ref());
        mean = map_statistics.mean();
        sigma = map_statistics.sigma();
        if (sigma == 0) sigma = 1;
      }
      OutFloatType* out_mem = reinterpret_cast<OutFloatType*>(
        malloc(out_size(first, last) * sizeof(OutFloatType)));
      OutFloatType* out_ptr = out_mem;
      af::ref<InpFloatType, maps::grid_p1<3> > a3d(
        a.begin(), af::adapt(a.accessor().grid()));
      af::long3 out_pt;
      for (out_pt[2] = first[2]; out_pt[2] <= last[2]; out_pt[2]++) {
      for (out_pt[1] = first[1]; out_pt[1] <= last[1]; out_pt[1]++) {
      for (out_pt[0] = first[0]; out_pt[0] <= last[0]; out_pt[0]++) {
        InpFloatType val = a3d(out_pt);
        if (apply_sigma_scaling) val = (val - mean) / sigma;
        *out_ptr++ = static_cast<OutFloatType>(val);
      }}}
      return boost::python::ref(PyCObject_FromVoidPtr(out_mem, free));
    }
  };

# include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    cctbx::af::bpl::import_flex();

    class_builder<ex_linear_regression<double> >
    py_linear_regression(this_module, "linear_regression");
    class_builder<ex_statistics<double> >
    py_statistics(this_module, "statistics");

    this_module.def(py_in_place_pow, "in_place_pow");
    this_module.def(py_set_if_less_than, "set_if_less_than");

    py_linear_regression.def(constructor<>());
    py_linear_regression.def(constructor<
      cctbx::af::shared<double>,
      cctbx::af::shared<double> >());
    py_linear_regression.def(constructor<
      cctbx::af::shared<double>,
      cctbx::af::shared<double>,
      double const&>());
    py_linear_regression.def(&ex_linear_regression<double>::is_well_defined,
      "is_well_defined");
    py_linear_regression.def(&ex_linear_regression<double>::b, "b");
    py_linear_regression.def(&ex_linear_regression<double>::m, "m");
    py_linear_regression.def(&ex_linear_regression<double>::cc, "cc");

    py_statistics.def(constructor<>());
    py_statistics.def(constructor<flex_double>());
    py_statistics.def(constructor<flex_float>());
    py_statistics.def(&ex_statistics<double>::min, "min");
    py_statistics.def(&ex_statistics<double>::max, "max");
    py_statistics.def(&ex_statistics<double>::mean, "mean");
    py_statistics.def(&ex_statistics<double>::mean2, "mean2");
    py_statistics.def(&ex_statistics<double>::sigma, "sigma");

    this_module.def(
      map_utils<float>::inplace_unpad, "inplace_unpad");
    this_module.def(
      map_utils<double>::inplace_unpad, "inplace_unpad");
    this_module.def(
      as_CObjectZYX<float, float>::convert, "as_CObjectZYXfloat");
    this_module.def(
      as_CObjectZYX<double, float>::convert, "as_CObjectZYXfloat");
  }

} // namespace

BOOST_PYTHON_MODULE_INIT(flex_utils)
{
  boost::python::module_builder this_module("flex_utils");
  init_module(this_module);
}
