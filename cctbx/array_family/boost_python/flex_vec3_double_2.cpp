/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Apr: Created (R.W. Grosse-Kunstleve)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/mat3.h>
#include <cctbx/math/utils.h>
#include <cctbx/error.h>

namespace scitbx { namespace af { namespace boost_python {

  namespace {

    vec3<double>
    vec3_min(flex<vec3<double> >::type const& a)
    {
      CCTBX_ASSERT(!a.accessor().is_padded());
      vec3<double> result(0,0,0);
      af::const_ref<vec3<double>, af::flex_grid<> > a_ref = a.const_ref();
      if (a_ref.size() > 0) {
        result = a_ref[0];
        for(std::size_t i=0;i<a_ref.size();i++) {
          for(std::size_t j=0;j<3;j++) {
            cctbx::math::update_min(result[j], a_ref[i][j]);
          }
        }
      }
      return result;
    }

    vec3<double>
    vec3_max(flex<vec3<double> >::type const& a)
    {
      CCTBX_ASSERT(!a.accessor().is_padded());
      vec3<double> result(0,0,0);
      af::const_ref<vec3<double>, af::flex_grid<> > a_ref = a.const_ref();
      if (a_ref.size() > 0) {
        result = a_ref[0];
        for(std::size_t i=0;i<a_ref.size();i++) {
          for(std::size_t j=0;j<3;j++) {
            cctbx::math::update_max(result[j], a_ref[i][j]);
          }
        }
      }
      return result;
    }

    af::shared<vec3<double> >
    mul_a_mat3(
      af::const_ref<vec3<double> > const& a,
      mat3<double> const& m)
    {
      af::shared<vec3<double> > result;
      for(std::size_t i=0;i<a.size();i++) {
        result.push_back(a[i] * m);
      }
      return result;
    }

    af::shared<vec3<double> >
    rmul_a_mat3(
      af::const_ref<vec3<double> > const& a,
      mat3<double> const& m)
    {
      mat3<double> m_transposed = m.transpose();
      af::shared<vec3<double> > result;
      for(std::size_t i=0;i<a.size();i++) {
        result.push_back(a[i] * m_transposed);
      }
      return result;
    }

  } // namespace <anonymous>

  void wrap_vec3_double_2(flex_wrapper<vec3<double> >::class_f_t class_object)
  {
    class_object
      .def("min", vec3_min)
      .def("max", vec3_max)
      .def("__add__", flex_wrapper<vec3<double> >::add_a_s)
      .def("__iadd__", flex_wrapper<vec3<double> >::iadd_a_s)
      .def("__mul__", mul_a_mat3)
      .def("__rmul__", rmul_a_mat3)
    ;
  }

}}} // namespace scitbx::af::boost_python
