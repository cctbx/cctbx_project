/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Dec 2001: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_VECTOR_UNARY_H
#define CCTBX_VECTOR_UNARY_H

#include <complex>

namespace cctbx { namespace vector {

  template <typename T> struct return_type {};

  template <typename ComplexVectorType, typename ResultVectorType>
  ResultVectorType complex_abs(const ComplexVectorType& v,
                               return_type<ResultVectorType>)
  {
    ResultVectorType result(v.size());
    for(std::size_t i=0;i<v.size();i++)
      result[i] = v[i].real() * v[i].real() + v[i].imag() * v[i].imag();
    return result;
  }

}} // namespace cctbx::vector

#endif // CCTBX_VECTOR_UNARY_H
