// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Mar 2002: modified copy of cctbx/maps/utils.h (rwgk)
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MATH_ARRAY_UTILS_H
#define CCTBX_MATH_ARRAY_UTILS_H

#include <cctbx/fixes/cmath>
#include <cctbx/array_family/reductions.h>
#include <cctbx/array_family/shared.h>

namespace cctbx { namespace math {

  template <typename FloatType>
  class array_statistics
  {
    public:
      array_statistics() {}
      array_statistics(const af::const_ref<FloatType>& a)
      {
        m_min = af::min(a);
        m_max = af::max(a);
        m_mean = af::mean(a);
        af::shared<FloatType> a2(a.begin(), a.end());
        af::ref<FloatType> a2r = a2.ref();
        for(std::size_t i=0;i<a2r.size();i++) a2r[i] *= a2r[i];
        m_mean2 = af::mean(a2r);
        m_sigma = m_mean2 - m_mean * m_mean;
        if (m_sigma < FloatType(0)) m_sigma = 0;
        m_sigma = std::sqrt(m_sigma);
      }
      const FloatType& min() const { return m_min; }
      const FloatType& max() const { return m_max; }
      const FloatType& mean() const { return m_mean; }
      const FloatType& mean2() const { return m_mean2; }
      const FloatType& sigma() const { return m_sigma; }
    protected:
      FloatType m_min;
      FloatType m_max;
      FloatType m_mean;
      FloatType m_mean2;
      FloatType m_sigma;
  };

}} // namespace cctbx::math

#endif // CCTBX_MATH_ARRAY_UTILS_H
