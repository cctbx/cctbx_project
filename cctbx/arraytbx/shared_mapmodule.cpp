// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/array_family/tiny_bpl.h>
#include <cctbx/error.h>
#include <cctbx/array_family/tiny_types.h>
#include <cctbx/array_family/shared.h>
#include <cctbx/array_family/reductions.h>
#include <cctbx/array_family/loops.h>
#include <cctbx/math/array_utils.h>
#include <cctbx/maps/accessors.h>

namespace {

  using namespace cctbx;

  template <typename FloatType>
  struct map_utils
  {
    static
    void
    inplace_unpad(
      af::shared<FloatType> map,
      af::long3 const& n_real,
      af::long3 const& m_real)
    {
      cctbx_assert(af::product(m_real.const_ref()) == map.size());
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
      map.resize(af::product(n_real.const_ref()));
    }
  };

  // PyMol support, based on code by N.K. Sauter
  template <typename InpFloatType,
            typename OutFloatType>
  struct as_CObjectZYX
  {
    static std::size_t out_size(af::long3 const& first, af::long3 const& last)
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
      af::shared<InpFloatType> a,
      af::long3 const& gridding,
      af::long3 const& first,
      af::long3 const& last,
      bool apply_sigma_scaling)
    {
      cctbx_assert(af::product(gridding.const_ref()) == a.size());
      math::array_statistics<InpFloatType> map_statistics;
      InpFloatType mean = 0;
      InpFloatType sigma = 0;
      if (apply_sigma_scaling) {
        math::array_statistics<InpFloatType>
        map_statistics(a.const_ref());
        mean = map_statistics.mean();
        sigma = map_statistics.sigma();
        if (sigma == 0) sigma = 1;
      }
      OutFloatType* out_mem = reinterpret_cast<OutFloatType*>(
        malloc(out_size(first, last) * sizeof(OutFloatType)));
      OutFloatType* out_ptr = out_mem;
      af::ref<InpFloatType, maps::grid_p1<3> > a3d(a.begin(), gridding);
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

    python::import_converters<af::shared<float> >
    py_sh_float("cctbx_boost.arraytbx.shared", "float");

    python::import_converters<af::shared<double> >
    py_sh_double("cctbx_boost.arraytbx.shared", "double");

    this_module.def(
      map_utils<float>::inplace_unpad, "inplace_unpad");
    this_module.def(
      map_utils<double>::inplace_unpad, "inplace_unpad");
    this_module.def(
      as_CObjectZYX<float, float>::convert, "as_CObjectZYXfloat");
    this_module.def(
      as_CObjectZYX<double, float>::convert, "as_CObjectZYXfloat");
  }

}

BOOST_PYTHON_MODULE_INIT(shared_map)
{
  boost::python::module_builder this_module("shared_map");
  init_module(this_module);
}
