// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Dec 2001: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_VECTOR_REDUCTIONS_H
#define CCTBX_VECTOR_REDUCTIONS_H

namespace cctbx { namespace vector {

  template <typename VectorType>
  typename VectorType::value_type
  sum(const VectorType& v)
  {
    typename VectorType::value_type result = 0;
    for(std::size_t i=0;i<v.size();i++) result += v[i];
    return result;
  }

  template <typename VectorType>
  typename VectorType::value_type
  mean(const VectorType& v)
  {
    return sum(v) / v.size();
  }

  template <typename VectorTypeW, typename VectorTypeV>
  typename VectorTypeV::value_type
  weighted_mean(const VectorTypeW& weights, const VectorTypeV& values)
  {
    return sum(weights * values) / sum(weights);
  }

}} // namespace cctbx::vector

#endif // CCTBX_VECTOR_REDUCTIONS_H
