// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MAPS_UTILS_H
#define CCTBX_MAPS_UTILS_H

#include <cctbx/fixes/cmath>
#include <cctbx/vector/standard.h>

namespace cctbx { namespace maps {

  template <typename ValueType>
  class statistics
  {
    public:
      statistics() {}
      template <typename VectorType>
      statistics(const VectorType& vec)
      {
        using namespace cctbx::vector;
        m_min = vector::min(vec);
        m_max = vector::max(vec);
        m_mean = vector::mean(vec);
        m_mean2 = vector::mean(vec * vec);
        m_sigma = m_mean2 - m_mean * m_mean;
        if (m_sigma < ValueType(0)) m_sigma = 0;
        m_sigma = std::sqrt(m_sigma);
      }
      const ValueType& min() const { return m_min; }
      const ValueType& max() const { return m_max; }
      const ValueType& mean() const { return m_mean; }
      const ValueType& mean2() const { return m_mean2; }
      const ValueType& sigma() const { return m_sigma; }
    protected:
      ValueType m_min;
      ValueType m_max;
      ValueType m_mean;
      ValueType m_mean2;
      ValueType m_sigma;
  };

}} // namespace cctbx::maps

#endif // CCTBX_MAPS_UTILS_H
